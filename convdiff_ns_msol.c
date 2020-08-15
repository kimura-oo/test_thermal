
#include "convdiff_core.h"

const int NUM_INTEG_POINTS_EACH_AXIS = 3;
const double MAT_EPSILON             = 1.0e-10;
const double MAT_MAX_ITER            = 10000;
const int BUFFER_SIZE                = 10000;

const double DELTA = 1.0E-06;

static const char* OUTPUT_FILENAME_VTK        = "result_%06d.vtk";
static const char* OUTPUT_FILENAME_ASCII_TEMP = "temparature_%06d.dat";
static const char* OUTPUT_FILENAME_ASCII_SOURCE = "source_%06d.dat";

const double DT = 0.01;
const double FIN_T = 10.0;
const int OUTPUT_INTERVAL = 10;


typedef struct
{
	double* T;
	double* error;
	double* theo_sol;

} VALUES;


typedef struct
{
	const char* directory;

} CONDITIONS;


typedef struct
{
	BBFE_BASIS   basis;
	BBFE_DATA    fe;
	BBFE_BC      bc;
	MONOLIS      monolis;
	
	CONDITIONS   cond;
	VALUES       vals;

	MONOLIS      monolis0; // for nonsteady analysis

} FE_SYSTEM;


double manusol_get_sol(
		double x,
		double y,
		double z,
		double t)
{
	double val = sin( 0.25*x ) * sin( 0.5*y ) * sin( 1.0*z ) * (0.5*t + sin( 1.0*t ));

	return val;
}


void manusol_get_conv_vel(
		double v[3],
		double x[3])
{
	double val = 1.0/sqrt(3) * 0.0;
	v[0] = 1.0 + x[0]*x[0];
	v[1] = 1.0 + x[1]*x[1];
	v[2] = 1.0 + x[2]*x[2];
}


double manusol_get_diff_coef(
		double x[3])
{
	double val = (2.0 + sin(1.0*x[0]) * sin(0.5*x[1]) * sin(0.25*x[2]));
	//double val = 1.0;

	return val;
}


double manusol_get_source(
		double x[3],
		double t,
		double a,
		double v[3],
		double k)
{
	double invD  = 1.0/DELTA;
	double invD2 = invD*invD;

	double dT_dt = ( manusol_get_sol(x[0], x[1], x[2], t+DELTA) -  
		manusol_get_sol(x[0], x[1], x[2], t-DELTA) ) * (invD/2.0);

	double dk_dx[3];
	double x_p[3];  double x_m[3];
	x_p[0] = x[0]+DELTA;  x_p[1] = x[1];  x_p[2] = x[2];
	x_m[0] = x[0]-DELTA;  x_m[1] = x[1];  x_m[2] = x[2];
	dk_dx[0] = ( manusol_get_diff_coef(x_p) - manusol_get_diff_coef(x_m) ) * (invD/2.0);
	x_p[0] = x[0];  x_p[1] = x[1]+DELTA;  x_p[2] = x[2];
	x_m[0] = x[0];  x_m[1] = x[1]-DELTA;  x_m[2] = x[2];
	dk_dx[1] = ( manusol_get_diff_coef(x_p) - manusol_get_diff_coef(x_m) ) * (invD/2.0);	
	x_p[0] = x[0];  x_p[1] = x[1];  x_p[2] = x[2]+DELTA;
	x_m[0] = x[0];  x_m[1] = x[1];  x_m[2] = x[2]-DELTA;
	dk_dx[2] = ( manusol_get_diff_coef(x_p) - manusol_get_diff_coef(x_m) )* (invD/2.0);

	double dT_dx[3];
	dT_dx[0] =  ( manusol_get_sol(x[0]+DELTA, x[1], x[2], t) -  
		manusol_get_sol(x[0]-DELTA, x[1], x[2], t) ) * (invD/2.0);
	dT_dx[1] =  ( manusol_get_sol(x[0], x[1]+DELTA, x[2], t) -  
		manusol_get_sol(x[0], x[1]-DELTA, x[2], t) ) * (invD/2.0);
	dT_dx[2] =  ( manusol_get_sol(x[0], x[1], x[2]+DELTA, t) -  
		manusol_get_sol(x[0], x[1], x[2]-DELTA, t) ) * (invD/2.0);

	double d2T_dx2[3];
	d2T_dx2[0] = ( 
			    manusol_get_sol(x[0]+DELTA, x[1], x[2], t) + 
			    manusol_get_sol(x[0]-DELTA, x[1], x[2], t) -
			2.0*manusol_get_sol(x[0]      , x[1], x[2], t) ) * invD2;
	d2T_dx2[1] = ( 
			    manusol_get_sol(x[0], x[1]+DELTA, x[2], t) + 
			    manusol_get_sol(x[0], x[1]-DELTA, x[2], t) -
			2.0*manusol_get_sol(x[0], x[1]      , x[2], t) ) * invD2;
	d2T_dx2[2] = ( 
			    manusol_get_sol(x[0], x[1], x[2]+DELTA, t) + 
			    manusol_get_sol(x[0], x[1], x[2]-DELTA, t) -
			2.0*manusol_get_sol(x[0], x[1], x[2]      , t) ) * invD2;

	double val = 
		a * (dT_dt + v[0]*dT_dx[0] + v[1]*dT_dx[1] + v[2]*dT_dx[2]) - 
		(dk_dx[0]*dT_dx[0] + dk_dx[1]*dT_dx[1] + dk_dx[2]*dT_dx[2]) - 
		k * (d2T_dx2[0] + d2T_dx2[1] + d2T_dx2[2]);

	return val;
}


void manusol_set_theo_sol(
		BBFE_DATA* fe,
		double*  theo_sol,
		double   t)
{
	for(int i=0; i<(fe->total_num_nodes); i++) {
		theo_sol[i] = manusol_get_sol(fe->x[i][0], fe->x[i][1], fe->x[i][2], t);
	}
}


void manusol_set_source(
		BBFE_DATA* fe,
		double*  source,
		double   t)
{
	for(int i=0; i<(fe->total_num_nodes); i++) {
		double v[3];  double k;
		manusol_get_conv_vel(v, fe->x[i]);
		k = manusol_get_diff_coef(fe->x[i]);
		source[i] = manusol_get_source(fe->x[i], t, 1.0, v, k);
	}
}


void memory_allocation_nodal_values(
		VALUES*         vals,
		const int       total_num_nodes)
{
	vals->T        = BB_std_calloc_1d_double(vals->T,     total_num_nodes);
	vals->error    = BB_std_calloc_1d_double(vals->error, total_num_nodes);
	vals->theo_sol = BB_std_calloc_1d_double(vals->error, total_num_nodes);
}


void output_result_file_vtk(
		BBFE_DATA*       fe,
		VALUES*        vals,
		const char*    filename,
		const char*    directory,
		double         t)
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
	manusol_set_source(fe, source, t);
	BB_vtk_write_point_vals_scalar(fp, source, fe->total_num_nodes, "source");
	BB_std_free_1d_double(source, fe->total_num_nodes);

	fclose(fp);

}


void output_files(
		FE_SYSTEM* sys,
		int file_num,
		double t)
{
	char fname_vtk[BUFFER_SIZE];
	char fname_tem[BUFFER_SIZE];
	char fname_sou[BUFFER_SIZE];
	snprintf(fname_vtk, BUFFER_SIZE, OUTPUT_FILENAME_VTK, file_num);
	snprintf(fname_tem, BUFFER_SIZE, OUTPUT_FILENAME_ASCII_TEMP, file_num);
	snprintf(fname_sou, BUFFER_SIZE, OUTPUT_FILENAME_ASCII_SOURCE, file_num);

	output_result_file_vtk(
			&(sys->fe),
			&(sys->vals),
			fname_vtk,
			sys->cond.directory,
			t);

	BBFE_write_ascii_nodal_vals_scalar(
			&(sys->fe),
			sys->vals.T,
			fname_tem,
			sys->cond.directory);

	/**** for manufactured solution ****/
	double* source;
	source = BB_std_calloc_1d_double(source, sys->fe.total_num_nodes);
	manusol_set_source(&(sys->fe), source, t);

	BBFE_write_ascii_nodal_vals_scalar(
			&(sys->fe),
			source,
			fname_sou,
			sys->cond.directory);
	double L2_error = BBFE_elemmat_equivval_relative_L2_error_scalar(
			&(sys->fe),
			&(sys->basis),
			t,
			sys->vals.T,
			manusol_get_sol);
	printf("%s L2 error: %e\n", CODENAME, L2_error);
	FILE* fp;
	fp = BBFE_sys_write_add_fopen(fp, "l2_error.txt", sys->cond.directory);
	fprintf(fp, "%e %e\n", t, L2_error);
	fclose(fp);

	BB_std_free_1d_double(source, sys->fe.total_num_nodes);
	/***********************************/
}


void set_element_mat(
		MONOLIS*     monolis,
		BBFE_DATA*     fe,
		BBFE_BASIS* basis,
		VALUES*      vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	local_x   = BB_std_calloc_2d_double(local_x, nl, 3);

	double** x_ip;  double** v_ip;  double* k_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	k_ip = BB_std_calloc_1d_double(k_ip, np);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);
		
		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector_value_integ_point_3d(x_ip[p], nl, local_x, basis->N[p]);
			manusol_get_conv_vel(v_ip[p], x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);

		}

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;

					double tau = BBFE_elemmat_convdiff_stab_coef_ns(
							k_ip[p], v_ip[p], 1.0, h_e, DT);

					val_ip[p] += BBFE_elemmat_convdiff_mat_conv(
							basis->N[p][i], fe->geo[e][p].grad_N[j], v_ip[p]);

					val_ip[p] -= BBFE_elemmat_convdiff_mat_diff(
							fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], k_ip[p]);

					val_ip[p] += BBFE_elemmat_convdiff_mat_stab_conv(
							fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], v_ip[p], tau);

					val_ip[p] += 1.0/DT * BBFE_elemmat_convdiff_mat_mass(
								basis->N[p][i], basis->N[p][j], 1.0);

					val_ip[p] += 1.0/DT *  BBFE_elemmat_convdiff_mat_stab_mass(
								fe->geo[e][p].grad_N[i], basis->N[p][j], 1.0, v_ip[p], tau);

				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				monolis_add_scalar_to_sparse_matrix(
						monolis, integ_val,
						fe->conn[e][i], fe->conn[e][j], 0, 0);
			}
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x,   fe->local_num_nodes, 3);
	
	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
}


void set_element_vec(
		MONOLIS*     monolis,
		BBFE_DATA*     fe,
		BBFE_BASIS* basis,
		VALUES*      vals,
		double       t)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;  double* local_T;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);
	local_T = BB_std_calloc_1d_double(local_T, nl);

	double** x_ip;  double** v_ip;  double* k_ip;  double* T_ip;  double* f_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	k_ip = BB_std_calloc_1d_double(k_ip, np);
	T_ip = BB_std_calloc_1d_double(T_ip, np);
	f_ip = BB_std_calloc_1d_double(f_ip, np);

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

		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x,   e, 3);
		BBFE_elemmat_set_local_array_scalar(local_T, fe, vals->T, e);
		
		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector_value_integ_point_3d(x_ip[p], nl, local_x, basis->N[p]);
			manusol_get_conv_vel(v_ip[p], x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);
			T_ip[p] = BBFE_std_mapping_scalar_value_integ_point_3d(nl, local_T, basis->N[p]);
			f_ip[p] = manusol_get_source(x_ip[p], t, 1.0, v_ip[p], k_ip[p]);
		}

		for(int i=0; i<nl; i++) {
			for(int p=0; p<np; p++) {
				val_ip[p] = 0.0;

				double tau = BBFE_elemmat_convdiff_stab_coef_ns(
						k_ip[p], v_ip[p], 1.0, h_e, DT);

				val_ip[p] += BBFE_elemmat_convdiff_vec_source(
						basis->N[p][i], f_ip[p]);

				val_ip[p] += BBFE_elemmat_convdiff_vec_stab_source(
						fe->geo[e][p].grad_N[i], v_ip[p], tau, f_ip[p]);
			
				val_ip[p] += 1.0/DT * BBFE_elemmat_convdiff_vec_mass(
						basis->N[p][i], T_ip[p], 1.0);

				val_ip[p] += 1.0/DT * BBFE_elemmat_convdiff_vec_stab_mass(
						fe->geo[e][p].grad_N[i], 1.0, v_ip[p], T_ip[p], tau);
			}

			double integ_val = BBFE_std_integ_calc(
					np, val_ip, basis->integ_weight, Jacobian_ip);

			monolis->mat.B[ fe->conn[e][i] ] += integ_val;
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x, fe->local_num_nodes, 3);
	BB_std_free_1d_double(local_T, fe->local_num_nodes);
	
	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(T_ip, np);
	BB_std_free_1d_double(f_ip, np);
}



int main (
		int argc,
		char* argv[])
{
	printf("\n");

	FE_SYSTEM sys;

	monolis_global_initialize();
	double t1 = monolis_get_time();
	
	sys.cond.directory = BBFE_convdiff_get_directory_name(argc, argv, CODENAME);	

	BBFE_convdiff_pre(
			&(sys.fe), &(sys.basis), (&sys.bc), (&sys.monolis),
			argc, argv, sys.cond.directory,
			NUM_INTEG_POINTS_EACH_AXIS, 
			true);
	
	memory_allocation_nodal_values(
			&(sys.vals),
			sys.fe.total_num_nodes);

	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, "l2_error.txt", sys.cond.directory);
	fclose(fp);

	BBFE_elemmat_set_Jacobi_mat(
			&(sys.fe),
			&(sys.basis));
	BBFE_elemmat_set_shapefunc_derivative(
			&(sys.fe),
			&(sys.basis));

	monolis_initialize(&(sys.monolis0));
	monolis_get_nonzero_pattern(
			&(sys.monolis0),
			sys.fe.total_num_nodes,
			sys.fe.local_num_nodes,
			1,
			sys.fe.total_num_elems,
			sys.fe.conn);
	
	set_element_mat(
			&(sys.monolis0),
			&(sys.fe),
			&(sys.basis),
			&(sys.vals));

	/****************** solver ********************/
	double t = 0.0;
	int step = 0;
	int file_num = 0;
	while (t < FIN_T) {
		t += DT;
		step += 1;

		printf("%s ----------------- step %d ----------------\n", CODENAME, step);
		monolis_copy_all(&(sys.monolis0), &(sys.monolis));
		monolis_clear_rhs(&(sys.monolis));
		set_element_vec(
				&(sys.monolis),
				&(sys.fe),
				&(sys.basis),
				&(sys.vals),
				t);

		manusol_set_theo_sol(&(sys.fe), sys.vals.theo_sol, t);
		BBFE_manusol_set_bc_scalar(
				&(sys.fe),
				&(sys.bc),
				sys.vals.theo_sol,
				t);

		BBFE_sys_monowrap_set_Dirichlet_bc(
				&(sys.monolis),
				sys.fe.total_num_nodes,
				BLOCK_SIZE,
				&(sys.bc),
				sys.monolis.mat.B);

		BBFE_sys_monowrap_solve(
				&(sys.monolis),
				sys.vals.T,
				monolis_iter_BiCGSTAB,
				monolis_prec_DIAG,
				MAT_MAX_ITER,
				MAT_EPSILON);
		/**********************************************/

		if(step%OUTPUT_INTERVAL == 0) {
			output_files(&sys, file_num, t);
			file_num += 1;
		}

	}

	BBFE_convdiff_finalize(&(sys.fe), &(sys.basis), &(sys.bc));
	
	monolis_finalize(&(sys.monolis));
	monolis_finalize(&(sys.monolis0));
	
	monolis_global_finalize();

	double t2 = monolis_get_time();
	printf("** Total time: %f\n", t2 - t1);

	printf("\n");

	return 0;
}
