
#include "fluid_core.h"

const char* ID_NUM_IP_EACH_AXIS = "#num_ip_each_axis";
const int DVAL_NUM_IP_EACH_AXIS = 2;
const char*     ID_MAT_EPSILON  = "#mat_epsilon";
const double  DVAL_MAT_EPSILON  = 1.0e-8;
const char*    ID_MAT_MAX_ITER  = "#mat_max_iter";
const int    DVAL_MAT_MAX_ITER  = 10000;
const char*              ID_DT  = "#time_spacing";
const double           DVAL_DT  = 0.01;
const char*     ID_FINISH_TIME  = "#finish_time";
const double  DVAL_FINISH_TIME  = 10.0;
const char* ID_OUTPUT_INTERVAL  = "#output_interval";
const int DVAL_OUTPUT_INTERVAL  = 4;


const double DELTA    = 1.0E-06;
const int BUFFER_SIZE = 10000;

static const char* INPUT_FILENAME_COND    = "cond.dat";
static const char* INPUT_FILENAME_D_BC_V  = "D_bc_v.dat";
static const char* INPUT_FILENAME_D_BC_P  = "D_bc.dat";

static const char* OUTPUT_FILENAME_VTK    = "result_%06d.vtk";


double RE = 10000.0;

typedef struct
{
	int    num_ip_each_axis;
	double mat_epsilon;
	int    mat_max_iter;

	double dt;
	double finish_time;
	int    output_interval;

	double** v;
	double*  p;

} VALUES;


typedef struct
{
	const char* directory;

} CONDITIONS;


typedef struct
{
	BBFE_BASIS   basis;
	BBFE_DATA    fe;

	CONDITIONS   cond;
	VALUES       vals;

	BBFE_BC      bc_v;
	BBFE_BC      bc_p;

	MONOLIS      mono_pred;
	MONOLIS      mono_corr;
	MONOLIS      mono_corr0;
	MONOLIS      mono_ppe;
	MONOLIS      mono_ppe0;

} FE_SYSTEM;


void memory_allocation_nodal_values(
		VALUES*         vals,
		const int       total_num_nodes)
{
	vals->v = BB_std_calloc_2d_double(vals->v, total_num_nodes, 3);
	vals->p = BB_std_calloc_1d_double(vals->p, total_num_nodes);
}


void assign_default_values(
		VALUES*     vals)
{
	vals->num_ip_each_axis = DVAL_NUM_IP_EACH_AXIS;
	vals->mat_epsilon      = DVAL_MAT_EPSILON;
	vals->mat_max_iter     = DVAL_MAT_MAX_ITER;

	vals->dt               = DVAL_DT;
	vals->finish_time      = DVAL_FINISH_TIME;
	vals->output_interval  = DVAL_OUTPUT_INTERVAL;
}


void print_all_values(
		VALUES*  vals)
{
	printf("\n%s ---------- Calculation condition ----------\n", CODENAME);

	printf("%s %s: %d\n", CODENAME, ID_NUM_IP_EACH_AXIS, vals->num_ip_each_axis);
	printf("%s %s: %e\n", CODENAME, ID_MAT_EPSILON,      vals->mat_epsilon);
	printf("%s %s: %d\n", CODENAME, ID_MAT_MAX_ITER,     vals->mat_max_iter);

	printf("%s %s: %e\n", CODENAME, ID_DT,               vals->dt);
	printf("%s %s: %e\n", CODENAME, ID_FINISH_TIME,      vals->finish_time);
	printf("%s %s: %d\n", CODENAME, ID_OUTPUT_INTERVAL,  vals->output_interval);

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
		num = BB_std_read_file_get_val_double_p(
				&(vals->dt), filename, ID_DT, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->finish_time), filename, ID_FINISH_TIME, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_int_p(
				&(vals->output_interval), filename, ID_OUTPUT_INTERVAL, BUFFER_SIZE, CODENAME);
	}

	print_all_values(vals);

	fclose(fp);

	printf("\n");
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
	BB_vtk_write_point_vals_vector(fp, vals->v, fe->total_num_nodes, "Velocity");
	BB_vtk_write_point_vals_scalar(fp, vals->p, fe->total_num_nodes, "Pressure");

	fclose(fp);

}


void output_files(
		FE_SYSTEM* sys,
		int file_num,
		double t)
{
	char fname_vtk[BUFFER_SIZE];
	snprintf(fname_vtk, BUFFER_SIZE, OUTPUT_FILENAME_VTK, file_num);

	output_result_file_vtk(
			&(sys->fe),
			&(sys->vals),
			fname_vtk,
			sys->cond.directory,
			t);
}


void set_element_mat_pred(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);
	
	double** local_v;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);

	double** v_ip; 
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);
		
		BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
		}

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;
	
					double tau = BBFE_elemmat_convdiff_stab_coef_ns(
							1.0/RE, v_ip[p], 1.0, h_e, vals->dt);

					val_ip[p] += basis->N[p][i] * basis->N[p][j];
					val_ip[p] += tau * BB_calc_vec3d_dot(v_ip[p], fe->geo[e][p].grad_N[i]) * basis->N[p][j];
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				monolis_add_scalar_to_sparse_matrix(
						monolis, integ_val,
						fe->conn[e][i], fe->conn[e][j], 0, 0);
				monolis_add_scalar_to_sparse_matrix(
						monolis, integ_val,
						fe->conn[e][i], fe->conn[e][j], 1, 1);
				monolis_add_scalar_to_sparse_matrix(
						monolis, integ_val,
						fe->conn[e][i], fe->conn[e][j], 2, 2);
			}
		}
	}

	BB_std_free_1d_double(val_ip,      np);
	BB_std_free_1d_double(Jacobian_ip, np);
	
	BB_std_free_2d_double(local_v, nl, 3);

	BB_std_free_2d_double(v_ip, np, 3);
}


void set_element_vec_pred(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);

	double** v_ip; 
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);

	double*** grad_v_ip;
	grad_v_ip = BB_std_calloc_3d_double(grad_v_ip, np, 3, 3);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(
				np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
			BBFE_std_mapping_vector3d_grad(grad_v_ip[p], nl, local_v, fe->geo[e][p].grad_N);
		}

		for(int i=0; i<nl; i++) {

			for(int d=0; d<3; d++) {
				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;	
					
					double tau = BBFE_elemmat_convdiff_stab_coef_ns(
							1.0/RE, v_ip[p], 1.0, h_e, vals->dt);

					val_ip[p] += - basis->N[p][i]*(
							v_ip[p][0] * grad_v_ip[p][d][0] + 
							v_ip[p][1] * grad_v_ip[p][d][1] + 
							v_ip[p][2] * grad_v_ip[p][d][2] 
							);
					
					val_ip[p] += -1.0/RE*(
							fe->geo[e][p].grad_N[i][0] * grad_v_ip[p][d][0] +
							fe->geo[e][p].grad_N[i][1] * grad_v_ip[p][d][1] +
							fe->geo[e][p].grad_N[i][2] * grad_v_ip[p][d][2]
							);

					val_ip[p] += -tau * BB_calc_vec3d_dot(v_ip[p], fe->geo[e][p].grad_N[i]) * 
						BB_calc_vec3d_dot(v_ip[p], grad_v_ip[p][d]);
					
					val_ip[p] *= vals->dt;

					val_ip[p] += basis->N[p][i] * v_ip[p][d];
					val_ip[p] += tau * BB_calc_vec3d_dot(v_ip[p], fe->geo[e][p].grad_N[i]) * v_ip[p][d];
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				monolis->mat.B[ 3*(fe->conn[e][i]) + d ] += integ_val;
			}
		}
	}

	BB_std_free_1d_double(val_ip,      np);
	BB_std_free_1d_double(Jacobian_ip, np);

	BB_std_free_2d_double(local_v, nl, 3);

	BB_std_free_2d_double(v_ip, np, 3);

	BB_std_free_3d_double(grad_v_ip, np, 3, 3);
}


void set_element_mat_corr(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;

					val_ip[p] += basis->N[p][i] * basis->N[p][j];
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				monolis_add_scalar_to_sparse_matrix(
						monolis, integ_val,
						fe->conn[e][i], fe->conn[e][j], 0, 0);
				monolis_add_scalar_to_sparse_matrix(
						monolis, integ_val,
						fe->conn[e][i], fe->conn[e][j], 1, 1);
				monolis_add_scalar_to_sparse_matrix(
						monolis, integ_val,
						fe->conn[e][i], fe->conn[e][j], 2, 2);
			}
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);
}


void set_element_vec_corr(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v;  double* local_p;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);
	local_p = BB_std_calloc_1d_double(local_p, nl);

	double** v_ip;  double** grad_p_ip;
	v_ip      = BB_std_calloc_2d_double(v_ip     , np, 3);
	grad_p_ip = BB_std_calloc_2d_double(grad_p_ip, np, 3);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);
		BBFE_elemmat_set_local_array_scalar(local_p, fe, vals->p, e);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
			BBFE_std_mapping_scalar_grad(grad_p_ip[p], nl, local_p, fe->geo[e][p].grad_N);
		}

		for(int i=0; i<nl; i++) {

			for(int d=0; d<3; d++) {
				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;

					val_ip[p] += - basis->N[p][i]*grad_p_ip[p][d];

					val_ip[p] *= vals->dt;

					val_ip[p] += basis->N[p][i] * v_ip[p][d];
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				monolis->mat.B[ 3*(fe->conn[e][i]) + d ] += integ_val;
			}
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_v, fe->local_num_nodes, 3);
	BB_std_free_1d_double(local_p, fe->local_num_nodes);
	
	BB_std_free_2d_double(v_ip     , np, 3);
	BB_std_free_2d_double(grad_p_ip, np, 3);
}


void set_element_mat_ppe(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;

					val_ip[p] += 
						-( fe->geo[e][p].grad_N[i][0] * fe->geo[e][p].grad_N[j][0] + 
								fe->geo[e][p].grad_N[i][1] * fe->geo[e][p].grad_N[j][1] + 
								fe->geo[e][p].grad_N[i][2] * fe->geo[e][p].grad_N[j][2] );
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
}


void set_element_vec_ppe(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);
	double* div_v_ip;
	div_v_ip = BB_std_calloc_1d_double(div_v_ip, np);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);
		for(int p=0; p<np; p++) {
			div_v_ip[p] = BBFE_std_mapping_vector3d_div(nl, local_v, fe->geo[e][p].grad_N);
		}

		for(int i=0; i<nl; i++) {
			for(int p=0; p<np; p++) {
				val_ip[p] = 0.0;

				val_ip[p] += 1.0/(vals->dt) * div_v_ip[p] * basis->N[p][i];
			}

			double integ_val = BBFE_std_integ_calc(
					np, val_ip, basis->integ_weight, Jacobian_ip);

			monolis->mat.B[ fe->conn[e][i] ] += integ_val;
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_v, fe->local_num_nodes, 3);
	BB_std_free_1d_double(div_v_ip, np);
}


int main(
		int   argc,
		char* argv[])
{
	printf("\n");

	FE_SYSTEM sys;

	monolis_global_initialize();
	double t1 = monolis_get_time();

	sys.cond.directory = BBFE_fluid_get_directory_name(argc, argv, CODENAME);	
	read_calc_conditions(&(sys.vals), sys.cond.directory);

	BBFE_fluid_pre(
			&(sys.fe), &(sys.basis),
			argc, argv, sys.cond.directory,
			sys.vals.num_ip_each_axis);

	BBFE_sys_read_Dirichlet_bc(
			&(sys.bc_p),
			INPUT_FILENAME_D_BC_P,
			sys.cond.directory,
			sys.fe.total_num_nodes,
			1);

	memory_allocation_nodal_values(
			&(sys.vals),
			sys.fe.total_num_nodes);

	BBFE_elemmat_set_Jacobi_mat(&(sys.fe), &(sys.basis));
	BBFE_elemmat_set_shapefunc_derivative(&(sys.fe), &(sys.basis));

	BBFE_sys_monowrap_init_monomat(&(sys.mono_pred) , &(sys.fe), 3);
	BBFE_sys_monowrap_init_monomat(&(sys.mono_ppe)  , &(sys.fe), 1);
	BBFE_sys_monowrap_init_monomat(&(sys.mono_ppe0) , &(sys.fe), 1);
	BBFE_sys_monowrap_init_monomat(&(sys.mono_corr) , &(sys.fe), 3);
	BBFE_sys_monowrap_init_monomat(&(sys.mono_corr0), &(sys.fe), 3);


	for(int i=0; i<(sys.fe.total_num_nodes); i++) {
		double x[3];
		x[0] = sys.fe.x[i][0];
		x[1] = sys.fe.x[i][1];
		x[2] = sys.fe.x[i][2];
		if(
				4.0 < x[0] && x[0] < 6.0 &&
				2.0 < x[1] && x[1] < 4.0 &&
				4.0 < x[2] && x[2] < 6.0 
		  ) {
			sys.vals.v[i][0] = 0.0;
			sys.vals.v[i][1] = 2.0;
			sys.vals.v[i][2] = 0.0;
		}
	}

	set_element_mat_ppe(
			&(sys.mono_ppe0),
			&(sys.fe),
			&(sys.basis),
			&(sys.vals));
	set_element_mat_corr(
			&(sys.mono_corr0),
			&(sys.fe),
			&(sys.basis),
			&(sys.vals));
	/****************** solver ********************/
	double t = 0.0;
	int step = 0;
	int file_num = 0;
	while (t < sys.vals.finish_time) {
		t += sys.vals.dt;
		step += 1;

		printf("\n%s ----------------- step %d ----------------\n", CODENAME, step);
		
		BBFE_sys_monowrap_copy_mat(&(sys.mono_ppe0) , &(sys.mono_ppe));
		BBFE_sys_monowrap_copy_mat(&(sys.mono_corr0), &(sys.mono_corr));
	
		monolis_clear(&(sys.mono_pred));
		monolis_clear_rhs(&(sys.mono_ppe));
		monolis_clear_rhs(&(sys.mono_corr));

		printf("%s --- prediction step ---\n", CODENAME);
		set_element_mat_pred(
				&(sys.mono_pred),
				&(sys.fe),
				&(sys.basis),
				&(sys.vals));
		set_element_vec_pred(
				&(sys.mono_pred),
				&(sys.fe),
				&(sys.basis),
				&(sys.vals));
		BBFE_sys_monowrap_solve(
				&(sys.mono_pred),
				sys.mono_pred.mat.X,
				monolis_iter_BiCGSTAB,
				monolis_prec_SOR,
				sys.vals.mat_max_iter,
				sys.vals.mat_epsilon);

		BBFE_fluid_renew_velocity(
				sys.vals.v, 
				sys.mono_pred.mat.X,
				sys.fe.total_num_nodes);

		printf("%s --- pressure Poisson eq. ---\n", CODENAME);
		set_element_vec_ppe(
				&(sys.mono_ppe),
				&(sys.fe),
				&(sys.basis),
				&(sys.vals));
		BBFE_sys_monowrap_set_Dirichlet_bc(
				&(sys.mono_ppe),
				sys.fe.total_num_nodes,
				BLOCK_SIZE,
				&(sys.bc_p),
				sys.mono_ppe.mat.B);
		BBFE_sys_monowrap_solve(
				&(sys.mono_ppe),
				sys.vals.p,
				monolis_iter_CG,
				monolis_prec_SOR,
				sys.vals.mat_max_iter,
				sys.vals.mat_epsilon);

		printf("%s --- Correction step ---\n", CODENAME);
		set_element_vec_corr(
				&(sys.mono_corr),
				&(sys.fe),
				&(sys.basis),
				&(sys.vals));
		BBFE_sys_monowrap_solve(
				&(sys.mono_corr),
				sys.mono_corr.mat.X,
				monolis_iter_CG,
				monolis_prec_SOR,
				sys.vals.mat_max_iter,
				sys.vals.mat_epsilon);

		BBFE_fluid_renew_velocity(
				sys.vals.v, 
				sys.mono_corr.mat.X,
				sys.fe.total_num_nodes);

		/**********************************************/

		if(step%sys.vals.output_interval == 0) {
			output_files(&sys, file_num, t);
			file_num += 1;
		}

	}

	BBFE_fluid_finalize(&(sys.fe), &(sys.basis), &(sys.bc_p));

	monolis_finalize(&(sys.mono_ppe));
	monolis_global_finalize();

	double t2 = monolis_get_time();
	printf("** Total time: %f\n", t2 - t1);

	printf("\n");

}
