
#include "convdiff_core.h"

const int NUM_INTEG_POINTS_EACH_AXIS = 3;
const double MAT_EPSILON             = 1.0e-10;
const double MAT_MAX_ITER            = 100000;


double manusol_get_sol(
		double x,
		double y,
		double z,
		double t)
{
	double val = sin( 0.25*x ) * sin( 0.5*y ) * sin( 1.0*z );
	
	return val;
}


void manusol_get_conv_vel(
		double a[3],
		double x[3])
{
	a[0] = 0.0;
	a[1] = 0.0;
	a[2] = 0.0;
}


double manusol_get_diff_coef(
		double x[3])
{
	double val = 1.0;

	return val;
}


double manusol_get_source(
		double x[3],
		double a[3],
		double k)
{
	double dk_dx[3];
	dk_dx[0] = 0.0;
	dk_dx[1] = 0.0;
	dk_dx[2] = 0.0;

	double val = 
		a[0]     * ( 0.25*cos( 0.25*x[0] ) *     sin( 0.5*x[1] ) * sin( 1.0*x[2] ) ) + 
		a[1]     * (      sin( 0.25*x[0] ) * 0.5*cos( 0.5*x[1] ) * sin( 1.0*x[2] ) ) + 
		a[2]     * (      sin( 0.25*x[0] ) *     sin( 0.5*x[1] ) * cos( 1.0*x[2] ) ) - 
		dk_dx[0] * ( 0.25*cos( 0.25*x[0] ) *     sin( 0.5*x[1] ) * sin( 1.0*x[2] ) ) - 
		dk_dx[1] * (      sin( 0.25*x[0] ) * 0.5*cos( 0.5*x[1] ) * sin( 1.0*x[2] ) ) - 
		dk_dx[2] * (      sin( 0.25*x[0] ) *     sin( 0.5*x[1] ) * cos( 1.0*x[2] ) ) - 
		k * (-(0.25*0.25+0.5*0.5+1.0*1.0)*sin( 0.25*x[0] ) * sin( 0.5*x[1] ) * sin( 1.0*x[2] ));
	
	return val;
}


void manusol_set_theo_sol(
		FE_DATA* fe,
		double*  theo_sol,
		double   t)
{
	for(int i=0; i<(fe->total_num_nodes); i++) {
		theo_sol[i] = manusol_get_sol(fe->x[i][0], fe->x[i][1], fe->x[i][2], t);
	}
}


void manusol_set_source(
		FE_DATA* fe,
		double*  source)
{
	for(int i=0; i<(fe->total_num_nodes); i++) {
		double a[3];  double k;
		manusol_get_conv_vel(a, fe->x[i]);
		k = manusol_get_diff_coef(fe->x[i]);
		source[i] = manusol_get_source(fe->x[i], a, k);
	}
}


void output_result_file_vtk(
		FE_DATA*       fe,
		NODAL_VALUES*  vals,
		const char*    filename,
		const char*    directory)
{
	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, filename, directory);

	switch( fe->local_num_nodes ) {
		case 4:
			BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_TETRA);
			break;

		case 8:
			BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_HEXAHEDRON);
			break;
	}

	fprintf(fp, "POINT_DATA %d\n", fe->total_num_nodes);
	BB_vtk_write_point_vals_scalar(fp, vals->T, fe->total_num_nodes, "temperature");

	// for manufactured solution
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, vals->theo_sol, vals->T);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "abs_error");
	BB_vtk_write_point_vals_scalar(fp, vals->theo_sol, fe->total_num_nodes, "theoretical");
	
	double* source;
	source = BB_std_calloc_1d_double(source, fe->total_num_nodes);
	manusol_set_source(fe, source);
	BB_vtk_write_point_vals_scalar(fp, source, fe->total_num_nodes, "source");
	BB_std_free_1d_double(source, fe->total_num_nodes);

	fclose(fp);

}


void output_files(
		FE_SYSTEM* sys)
{
	output_result_file_vtk(
			&(sys->fe),
			&(sys->vals),
			OUTPUT_FILENAME_VTK,
			sys->cond.directory);

	BBFE_write_ascii_nodal_vals_scalar(
			&(sys->fe),
			sys->vals.T,
			OUTPUT_FILENAME_ASCII_TEMP,
			sys->cond.directory);

	/**** for manufactured solution ****/
	double* source;
	source = BB_std_calloc_1d_double(source, sys->fe.total_num_nodes);
	manusol_set_source(&(sys->fe), source);

	BBFE_write_ascii_nodal_vals_scalar(
			&(sys->fe),
			source,
			OUTPUT_FILENAME_ASCII_RHS,
			sys->cond.directory);
	double L2_error = BBFE_elemmat_equivval_relative_L2_error_scalar(
			&(sys->fe),
			&(sys->basis),
			0.0,
			sys->vals.T,
			manusol_get_sol);
	printf("L2 error: %e\n", L2_error);
	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, "l2_error.txt", sys->cond.directory);
	fprintf(fp, "%e\n", L2_error);
	fclose(fp);

	BB_std_free_1d_double(source, sys->fe.total_num_nodes);
	/***********************************/
}


void set_element_mat_vec(
		MONOLIS*     monolis,
		FE_DATA*     fe,
		FE_3D_BASIS* basis)
{
	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip      , basis->num_integ_points);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip , basis->num_integ_points);

	double** local_x;
	local_x   = BB_std_calloc_2d_double(local_x  , fe->local_num_nodes, 3);
	
	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(
				Jacobian_ip,
				basis->num_integ_points,
				e,
				fe);

		double vol = BBFE_std_integ_calc_volume(
				basis->num_integ_points,
				basis->integ_weight,
				Jacobian_ip);

		double h_e = cbrt(vol);
		
		for(int i=0; i<(fe->local_num_nodes); i++) {
			local_x[i][0] = fe->x[ fe->conn[e][i] ][0];
			local_x[i][1] = fe->x[ fe->conn[e][i] ][1];
			local_x[i][2] = fe->x[ fe->conn[e][i] ][2];
		}

		for(int i=0; i<(fe->local_num_nodes); i++) {
			for(int j=0; j<(fe->local_num_nodes); j++) {

				for(int p=0; p<(basis->num_integ_points); p++) {
					val_ip[p] = 0.0;
					
					double x_ip[3];
					BBFE_std_mapping_vector_value_integ_point_3d(
							x_ip,
							fe->local_num_nodes,
							local_x,
							basis->N[p]);

					double a_ip[3];  double k_ip;
					manusol_get_conv_vel(a_ip, x_ip);
					k_ip = manusol_get_diff_coef(x_ip);

					val_ip[p] += BBFE_elemmat_convdiff_mat_conv(
							basis->N[p][i],
							fe->geo[e][p].grad_N[j],
							a_ip);

					val_ip[p] -= BBFE_elemmat_convdiff_mat_diff(
							fe->geo[e][p].grad_N[i],
							fe->geo[e][p].grad_N[j],
							k_ip);

					double tau = BBFE_elemmat_convdiff_stab_coef(
							k_ip, a_ip, h_e);
					val_ip[p] += BBFE_elemmat_convdiff_mat_stab_conv(
							fe->geo[e][p].grad_N[i],
							fe->geo[e][p].grad_N[j],
							a_ip, tau);

				}

				double integ_val = BBFE_std_integ_calc(
						basis->num_integ_points,
						val_ip,
						basis->integ_weight,
						Jacobian_ip);

				monolis_add_scalar_to_sparse_matrix(
						monolis,
						integ_val,
						fe->conn[e][i], fe->conn[e][j], 0, 0);
			}
		}
	}

	// for manufactured solution
	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(
				Jacobian_ip,
				basis->num_integ_points,
				e,
				fe);

		double vol = BBFE_std_integ_calc_volume(
				basis->num_integ_points,
				basis->integ_weight,
				Jacobian_ip);

		double h_e = cbrt(vol);

		for(int i=0; i<(fe->local_num_nodes); i++) {
			local_x[i][0] = fe->x[ fe->conn[e][i] ][0];
			local_x[i][1] = fe->x[ fe->conn[e][i] ][1];
			local_x[i][2] = fe->x[ fe->conn[e][i] ][2];
		}

		for(int i=0; i<(fe->local_num_nodes); i++) {
			for(int p=0; p<(basis->num_integ_points); p++) {
				val_ip[p] = 0.0;

				double x_ip[3];
				BBFE_std_mapping_vector_value_integ_point_3d(
						x_ip,
						fe->local_num_nodes,
						local_x,
						basis->N[p]);

				double a_ip[3];  double k_ip;
				manusol_get_conv_vel(a_ip, x_ip);
				k_ip = manusol_get_diff_coef(x_ip);
				double source_ip = manusol_get_source(x_ip, a_ip, k_ip);
				val_ip[p] += BBFE_elemmat_convdiff_vec_source(
						basis->N[p][i],
						source_ip);

				double tau = BBFE_elemmat_convdiff_stab_coef(
						k_ip, a_ip, h_e);
				val_ip[p] += BBFE_elemmat_convdiff_vec_stab_source(
						fe->geo[e][p].grad_N[i],
						a_ip,
						tau,
						source_ip);

			}

			double integ_val = BBFE_std_integ_calc(
					basis->num_integ_points,
					val_ip,
					basis->integ_weight,
					Jacobian_ip);

			monolis->mat.B[ fe->conn[e][i] ] += integ_val;
		}
	}

	BB_std_free_1d_double(val_ip,     basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x,   fe->local_num_nodes, 3);
}


int main (
		int argc,
		char* argv[])
{
	printf("\n");

	FE_SYSTEM sys;

	monolis_global_initialize();
	double t1 = monolis_get_time();

	BBFE_convdiff_pre(&sys, argc, argv, 
			NUM_INTEG_POINTS_EACH_AXIS, true);

	manusol_set_theo_sol(&(sys.fe), sys.vals.theo_sol, 0.0);

	/****************** solver ********************/
	BBFE_elemmat_set_Jacobi_mat(
			&(sys.fe),
			&(sys.basis));
	BBFE_elemmat_set_shapefunc_derivative(
			&(sys.fe),
			&(sys.basis));
	set_element_mat_vec(
			&(sys.monolis),
			&(sys.fe),
			&(sys.basis));

	// for manufactured solution
	BBFE_manusol_set_bc_scalar(
			&(sys.fe),
			&(sys.bc),
			sys.vals.theo_sol,
			0.0);

	BBFE_sys_monowrap_set_Dirichlet_bc(
			&(sys.monolis),
			sys.fe.total_num_nodes,
			BLOCK_SIZE,
			&(sys.bc),
			sys.monolis.mat.B);

	BBFE_sys_monowrap_solve(
			&(sys.monolis),
			sys.vals.T,
			monolis_iter_BiCGSTAB_noprec,
			monolis_prec_DIAG,
			MAT_MAX_ITER,
			MAT_EPSILON);
	/**********************************************/

	output_files(&sys);

	monolis_finalize(&(sys.monolis));
	monolis_global_finalize();

	double t2 = monolis_get_time();
	printf("** Total time: %f\n", t2 - t1);

	printf("\n");

	return 0;
}
