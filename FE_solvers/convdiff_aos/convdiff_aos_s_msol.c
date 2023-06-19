
#include "convdiff_aos_core.h"

const char* ID_NUM_IP_EACH_AXIS = "#num_ip_each_axis";
const int DVAL_NUM_IP_EACH_AXIS = 3;
const char*     ID_MAT_EPSILON   = "#mat_epsilon";
const double  DVAL_MAT_EPSILON   = 1.0e-10;
const char*    ID_MAT_MAX_ITER  = "#mat_max_iter";
const int    DVAL_MAT_MAX_ITER  = 10000;


static const double DELTA    = 1.0E-05;
static const int BUFFER_SIZE = 10000;

static const char* INPUT_FILENAME_COND          = "cond.dat";
static const char* OUTPUT_FILENAME_VTK          = "result.vtk";
static const char* OUTPUT_FILENAME_ASCII_TEMP   = "temparature.dat";
static const char* OUTPUT_FILENAME_ASCII_SOURCE = "source.dat";


typedef struct
{
	int    num_ip_each_axis;
	double mat_epsilon;
	int    mat_max_iter;

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
	MONOLIS_COM  monolis_com;

	CONDITIONS   cond;
	VALUES       vals;

} FE_SYSTEM;


double manusol_get_sol(
		double x,
		double y,
		double z,
		double t)
{
	double val = sin( 0.25*x ) * sin( 0.5*y ) * sin( 1.0*z );

	return val;
}


double manusol_get_mass_coef(
		double x[3])
{
	double val = 1.0;

	return val;
}


void manusol_get_conv_vel(
		double v[3],
		double x[3])
{
	v[0] = 0.0; //1.0 / sqrt(3.0);
	v[1] = 0.0; //1.0 / sqrt(3.0);
	v[2] = 0.0; //1.0 / sqrt(3.0);
}


double manusol_get_diff_coef(
		double x[3])
{
	double val = 1.0;
	//double val = (2.0 + sin(1.0*x[0]) * sin(0.5*x[1]) * sin(0.25*x[2]));

	return val;
}


double manusol_get_source(
		double x[3],
		double a,
		double v[3],
		double k)
{
	double invD  = 1.0/DELTA;
	double invD2 = invD*invD;

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
	dT_dx[0] =  ( manusol_get_sol(x[0]+DELTA, x[1], x[2], 0.0) -
		manusol_get_sol(x[0]-DELTA, x[1], x[2], 0.0) ) * (invD/2.0);
	dT_dx[1] =  ( manusol_get_sol(x[0], x[1]+DELTA, x[2], 0.0) -
		manusol_get_sol(x[0], x[1]-DELTA, x[2], 0.0) ) * (invD/2.0);
	dT_dx[2] =  ( manusol_get_sol(x[0], x[1], x[2]+DELTA, 0.0) -
		manusol_get_sol(x[0], x[1], x[2]-DELTA, 0.0) ) * (invD/2.0);

	double d2T_dx2[3];
	d2T_dx2[0] = (
			    manusol_get_sol(x[0]+DELTA, x[1], x[2], 0.0) +
			    manusol_get_sol(x[0]-DELTA, x[1], x[2], 0.0) -
			2.0*manusol_get_sol(x[0]      , x[1], x[2], 0.0) ) * invD2;
	d2T_dx2[1] = (
			    manusol_get_sol(x[0], x[1]+DELTA, x[2], 0.0) +
			    manusol_get_sol(x[0], x[1]-DELTA, x[2], 0.0) -
			2.0*manusol_get_sol(x[0], x[1]      , x[2], 0.0) ) * invD2;
	d2T_dx2[2] = (
			    manusol_get_sol(x[0], x[1], x[2]+DELTA, 0.0) +
			    manusol_get_sol(x[0], x[1], x[2]-DELTA, 0.0) -
			2.0*manusol_get_sol(x[0], x[1], x[2]      , 0.0) ) * invD2;

	double val =
		a * (v[0]*dT_dx[0] + v[1]*dT_dx[1] + v[2]*dT_dx[2]) -
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
		double*  source)
{
	for(int i=0; i<(fe->total_num_nodes); i++) {
		double a; double v[3];  double k;
		a = manusol_get_mass_coef(fe->x[i]);
		manusol_get_conv_vel(v, fe->x[i]);
		k = manusol_get_diff_coef(fe->x[i]);
		source[i] = manusol_get_source(fe->x[i], a, v, k);
	}
}


void assign_default_values(
		VALUES*     vals)
{
	vals->num_ip_each_axis = DVAL_NUM_IP_EACH_AXIS;
	vals->mat_epsilon      = DVAL_MAT_EPSILON;
	vals->mat_max_iter     = DVAL_MAT_MAX_ITER;
}


void print_all_values(
		VALUES*  vals)
{

	printf("\n%s ---------- Calculation condition ----------\n", CODENAME);

	printf("%s %s: %d\n", CODENAME, ID_NUM_IP_EACH_AXIS, vals->num_ip_each_axis);
	printf("%s %s: %e\n", CODENAME, ID_MAT_EPSILON,      vals->mat_epsilon);
	printf("%s %s: %d\n", CODENAME, ID_MAT_MAX_ITER,     vals->mat_max_iter);

	printf("%s -------------------------------------------\n\n", CODENAME);
}


void read_calc_conditions(
		VALUES*     vals,
		const char* directory)
{
	printf("\n");

	assign_default_values(vals);

	char filename[BUFFER_SIZE];
	snprintf(filename, BUFFER_SIZE, "%s/%s", directory, INPUT_FILENAME_COND);

	FILE* fp;
	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s Calc condition file \"%s\" is not found.\n", CODENAME, filename);
		printf("%s Default values are used in this calculation.\n", CODENAME);
	}
	else {
		printf("%s Reading conditon file \"%s\".\n", CODENAME, filename);
		int num;
		num = BB_std_read_file_get_val_int_p(
				&(vals->num_ip_each_axis), filename, ID_NUM_IP_EACH_AXIS, BUFFER_SIZE, CODENAME);

		num = BB_std_read_file_get_val_double_p(
				&(vals->mat_epsilon), filename, ID_MAT_EPSILON, BUFFER_SIZE, CODENAME);

		num = BB_std_read_file_get_val_int_p(
				&(vals->mat_max_iter), filename, ID_MAT_MAX_ITER, BUFFER_SIZE, CODENAME);

		fclose(fp);
	}

	print_all_values(vals);


	printf("\n");
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
		const char*    directory)
{
	/*
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
	*/
}


void output_files(
		FE_SYSTEM* sys)
{
	const char* filename;

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, OUTPUT_FILENAME_VTK);
	output_result_file_vtk(
			&(sys->fe),
			&(sys->vals),
			filename,
			sys->cond.directory);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, OUTPUT_FILENAME_ASCII_TEMP);
	BBFE_write_ascii_nodal_vals_scalar(
			&(sys->fe),
			sys->vals.T,
			filename,
			sys->cond.directory);

	/**** for manufactured solution ****/
	double* source;
	source = BB_std_calloc_1d_double(source, sys->fe.total_num_nodes);
	manusol_set_source(&(sys->fe), source);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, OUTPUT_FILENAME_ASCII_SOURCE);
	BBFE_write_ascii_nodal_vals_scalar(
			&(sys->fe),
			source,
			filename,
			sys->cond.directory);

	double L2_error = BBFE_convdiff_equivval_relative_L2_error_scalar(
			&(sys->fe),
			&(sys->basis),
			&(sys->monolis),
			&(sys->monolis_com),
			0.0,
			sys->vals.T,
			manusol_get_sol);

	printf("%s L2 error: %e\n", CODENAME, L2_error);
	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, "l2_error.txt", sys->cond.directory);
	fprintf(fp, "%e\n", L2_error);
	fclose(fp);

	BB_std_free_1d_double(source, sys->fe.total_num_nodes);
	/***********************************/
}


void set_element_mat_vec(
		MONOLIS*     monolis,
		BBFE_DATA*     fe,
		BBFE_BASIS*  basis)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	local_x   = BB_std_calloc_2d_double(local_x, nl, 3);

	double** x_ip;  double* a_ip;  double** v_ip;  double* k_ip;  double* f_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	a_ip = BB_std_calloc_1d_double(a_ip, np);
	k_ip = BB_std_calloc_1d_double(k_ip, np);
	f_ip = BB_std_calloc_1d_double(f_ip, np);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			manusol_get_conv_vel(v_ip[p], x_ip[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);
			f_ip[p] = manusol_get_source(x_ip[p], a_ip[p], v_ip[p], k_ip[p]);
		}

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;

					double tau = BBFE_elemmat_convdiff_stab_coef(
					 	k_ip[p], a_ip[p], v_ip[p], h_e);

					val_ip[p] += BBFE_elemmat_convdiff_mat_conv(
							basis->N[p][i], fe->geo[e][p].grad_N[j], a_ip[p], v_ip[p]);

					val_ip[p] -= BBFE_elemmat_convdiff_mat_diff(
							fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], k_ip[p]);

					val_ip[p] += BBFE_elemmat_convdiff_mat_stab_conv(
							fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], a_ip[p], v_ip[p], tau);
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				monolis_add_scalar_to_sparse_matrix_R(
						monolis,
						fe->conn[e][i], fe->conn[e][j], 0, 0,
						integ_val);
			}
		}

		for(int i=0; i<(fe->local_num_nodes); i++) {
			for(int p=0; p<(basis->num_integ_points); p++) {
				val_ip[p] = 0.0;

				double tau = BBFE_elemmat_convdiff_stab_coef(
						k_ip[p], a_ip[p], v_ip[p], h_e);

				val_ip[p] += BBFE_elemmat_convdiff_vec_source(
						basis->N[p][i], f_ip[p]);

				val_ip[p] += BBFE_elemmat_convdiff_vec_stab_source(
						fe->geo[e][p].grad_N[i], a_ip[p], v_ip[p], tau, f_ip[p]);

			}
			double integ_val = BBFE_std_integ_calc(
					np, val_ip, basis->integ_weight, Jacobian_ip);

			monolis->mat.R.B[ fe->conn[e][i] ] += integ_val;
		}
	}

	BB_std_free_1d_double(val_ip,     basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x,   fe->local_num_nodes, 3);

	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(a_ip, np);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(f_ip, np);
}


int main (
		int argc,
		char* argv[])
{
	printf("\n");

	FE_SYSTEM sys;

	monolis_global_initialize();
	double t1 = monolis_get_time_global_sync();

	sys.cond.directory = BBFE_convdiff_get_directory_name(argc, argv, CODENAME);
	read_calc_conditions(&(sys.vals), sys.cond.directory);

	BBFE_convdiff_pre(
			&(sys.fe), &(sys.basis), (&sys.bc), (&sys.monolis), (&sys.monolis_com),
			argc, argv, sys.cond.directory,
			sys.vals.num_ip_each_axis,
			true);

	memory_allocation_nodal_values(
			&(sys.vals),
			sys.fe.total_num_nodes);

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

	monolis_show_timelog (&(sys.monolis), true);
	monolis_show_iterlog (&(sys.monolis), true);
	BBFE_sys_monowrap_set_Dirichlet_bc(
			&(sys.monolis),
			sys.fe.total_num_nodes,
			BLOCK_SIZE,
			&(sys.bc),
			sys.monolis.mat.R.B);

	double t2 = monolis_get_time_global_sync();
	printf("** Pre  time: %f\n", t2 - t1);

	BBFE_sys_monowrap_solve(
			&(sys.monolis),
			&(sys.monolis_com),
			sys.vals.T,
			MONOLIS_ITER_BICGSTAB,
			MONOLIS_PREC_DIAG,
			sys.vals.mat_max_iter,
			sys.vals.mat_epsilon);
	/**********************************************/

	double t3 = monolis_get_time_global_sync();
	output_files(&sys);

	BBFE_convdiff_finalize(&(sys.fe), &(sys.basis), &(sys.bc));

	double t4 = monolis_get_time_global_sync();
	printf("** Post  time: %f\n", t4 - t3);
	printf("** Total time: %f\n", t4 - t1);
	printf("\n");

	monolis_finalize(&(sys.monolis));
	monolis_global_finalize();

	return 0;
}
