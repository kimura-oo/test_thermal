
#include "solid_core.h"

const char* ID_NUM_IP_EACH_AXIS  = "#num_ip_each_axis";
const int DVAL_NUM_IP_EACH_AXIS  = 2;
const char*      ID_MAT_EPSILON  = "#mat_epsilon";
const double   DVAL_MAT_EPSILON  = 1.0e-8;
const char*     ID_MAT_MAX_ITER  = "#mat_max_iter";
const int     DVAL_MAT_MAX_ITER  = 10000;
const char*          ID_DENSITY  = "#density";
const double       DVAL_DENSITY  = 1000.0;
const char*    ID_YOUNG_MODULUS  = "#young_modulus";
const double DVAL_YOUNG_MODULUS  = 1.0e09;
const char*    ID_POISSON_RATIO  = "#poisson_ratio";
const double DVAL_POISSON_RATIO  = 0.3;
const char*          ID_GRAVITY  = "#gravity";
const double       DVAL_GRAVITY  = -9.81;

const int BUFFER_SIZE = 10000;
const double DISPLACEMENT_SCALE = 1.0;

static const char* INPUT_FILENAME_COND    = "cond.dat";
static const char* INPUT_FILENAME_D_BC    = "D_bc.dat";
static const char* INPUT_FILENAME_N_BC    = "N_bc.dat";

static const char* OUTPUT_FILENAME_VTK    = "result.vtk";


typedef struct
{
	int    num_ip_each_axis;
	double mat_epsilon;
	int    mat_max_iter;

	double density;
	double young;
	double poisson;
	double g[3];

	double** u;

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
	BBFE_BC      bc;
	MONOLIS      mono;
	MONOLIS_COM  mono_com;
} FE_SYSTEM;


void memory_allocation_nodal_values(
		VALUES*         vals,
		const int       total_num_nodes)
{
	vals->u = BB_std_calloc_2d_double(vals->u, total_num_nodes, 3);
}


void assign_default_values(
		VALUES*     vals)
{
	vals->num_ip_each_axis = DVAL_NUM_IP_EACH_AXIS;
	vals->mat_epsilon      = DVAL_MAT_EPSILON;
	vals->mat_max_iter     = DVAL_MAT_MAX_ITER;

	vals->density          = DVAL_DENSITY;
	vals->young            = DVAL_YOUNG_MODULUS;
	vals->poisson          = DVAL_POISSON_RATIO;

	vals->g[0] = 0.0;
	vals->g[1] = 0.0;
	vals->g[2] = DVAL_GRAVITY;
}


void print_all_values(
		VALUES*  vals)
{
	printf("\n%s ---------- Calculation condition ----------\n", CODENAME);

	printf("%s %s: %d\n", CODENAME, ID_NUM_IP_EACH_AXIS, vals->num_ip_each_axis);
	printf("%s %s: %e\n", CODENAME, ID_MAT_EPSILON,      vals->mat_epsilon);
	printf("%s %s: %d\n", CODENAME, ID_MAT_MAX_ITER,     vals->mat_max_iter);

	printf("%s %s: %e\n", CODENAME, ID_DENSITY,          vals->density);
	printf("%s %s: %e\n", CODENAME, ID_YOUNG_MODULUS,    vals->young);
	printf("%s %s: %e\n", CODENAME, ID_POISSON_RATIO,    vals->poisson);

	printf("%s %s: %e, %e, %e\n", CODENAME, ID_GRAVITY,  vals->g[0], vals->g[1], vals->g[2]);
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
				&(vals->density), filename, ID_DENSITY, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->young), filename, ID_YOUNG_MODULUS, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->poisson), filename, ID_POISSON_RATIO, BUFFER_SIZE, CODENAME);		
		
		num = BB_std_read_file_get_val_double_p(
				vals->g, filename, ID_GRAVITY, BUFFER_SIZE, CODENAME);

		fclose(fp);
	}

	print_all_values(vals);

	printf("\n");
}


void output_result_file_vtk(
		BBFE_DATA*       fe,
		VALUES*        vals,
		const char*    filename,
		const char*    directory,
		double         t,
		double         scale)
{
	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, filename, directory);

	switch( fe->local_num_nodes ) {
		case 4:
			BBFE_sys_write_vtk_shape_with_disp(fp, fe, TYPE_VTK_TETRA, vals->u, scale);
			break;

		case 8:
			BBFE_sys_write_vtk_shape_with_disp(fp, fe, TYPE_VTK_HEXAHEDRON, vals->u, scale);
			break;
	}

	fprintf(fp, "POINT_DATA %d\n", fe->total_num_nodes);
	BB_vtk_write_point_vals_vector(fp, vals->u, fe->total_num_nodes, "Displacement");

	fclose(fp);
}


void output_files(
		FE_SYSTEM* sys,
		int file_num,
		double t,
		double scale)
{
	const char* filename;
	char fname_vtk[BUFFER_SIZE];
	snprintf(fname_vtk, BUFFER_SIZE, OUTPUT_FILENAME_VTK, file_num);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_vtk);
	output_result_file_vtk(
			&(sys->fe), &(sys->vals), filename, sys->cond.directory, t, scale);
}


void set_element_mat(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double*** val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_3d_double(val_ip, 3, 3, np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {
				for(int p=0; p<np; p++) {
					double mat[3][3];
					BBFE_elemmat_solid_mat_linear(
							mat, 
							fe->geo[e][p].grad_N[i], 
							fe->geo[e][p].grad_N[j], 
							vals->young, vals->poisson);

					for(int k=0; k<3; k++) { 
						for(int l=0; l<3; l++) {
							val_ip[k][l][p] = mat[k][l];
						}
					}
				}

				for(int k=0; k<3; k++) {
					for(int l=0; l<3; l++) {
						double integ_val = BBFE_std_integ_calc(
								np, val_ip[k][l], basis->integ_weight, Jacobian_ip);
						monolis_add_scalar_to_sparse_matrix_R(
								monolis, fe->conn[e][i], fe->conn[e][j], k, l, integ_val);
					}
				}
			}
		}
	}

	BB_std_free_3d_double(val_ip, 3, 3, np);
	BB_std_free_1d_double(Jacobian_ip, np);
}


void set_element_vec(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double** val_ip;
	double*  Jacobian_ip;
	val_ip      = BB_std_calloc_2d_double(val_ip, 3, np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		for(int i=0; i<nl; i++) {
			double integ_val[3];

			for(int p=0; p<np; p++) {
				double vec[3];
				vec[0] = vals->g[0];  vec[1] = vals->g[1];  vec[2] = vals->g[2];

				for(int d=0; d<3; d++) {
					val_ip[d][p] = vec[d] * vals->density * basis->N[p][i];
				}
			}

			for(int d=0; d<3; d++) {
				integ_val[d] = BBFE_std_integ_calc(
						np, val_ip[d], basis->integ_weight, Jacobian_ip);

				monolis->mat.R.B[ 3*fe->conn[e][i] + d ] += integ_val[d];
			}
		}
	}

	BB_std_free_2d_double(val_ip, 3, np);
	BB_std_free_1d_double(Jacobian_ip, np);
}


int main(
		int   argc,
		char* argv[])
{
	printf("\n");

	FE_SYSTEM sys;
	const char* filename;

	monolis_global_initialize();
	double t1 = monolis_get_time();

	sys.cond.directory = BBFE_solid_get_directory_name(argc, argv, CODENAME);	
	read_calc_conditions(&(sys.vals), sys.cond.directory);

	BBFE_solid_pre(
			&(sys.fe), &(sys.basis),
			argc, argv, sys.cond.directory,
			sys.vals.num_ip_each_axis);

	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_D_BC);
	BBFE_sys_read_Dirichlet_bc(
			&(sys.bc),
			filename,
			sys.cond.directory,
			sys.fe.total_num_nodes, 3);

	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_N_BC);
	BBFE_sys_read_Neumann_bc(
			&(sys.bc),
			filename,
			sys.cond.directory,
			sys.fe.total_num_nodes, 3);

	memory_allocation_nodal_values(
			&(sys.vals), sys.fe.total_num_nodes);

	BBFE_elemmat_set_Jacobi_mat(&(sys.fe), &(sys.basis));
	BBFE_elemmat_set_shapefunc_derivative(&(sys.fe), &(sys.basis));

	BBFE_sys_monowrap_init_monomat(&(sys.mono), &(sys.mono_com), &(sys.fe), 3, sys.cond.directory);

	/****************** solver ********************/
	set_element_mat(
			&(sys.mono), &(sys.fe), &(sys.basis), &(sys.vals));

	set_element_vec(
			&(sys.mono), &(sys.fe), &(sys.basis), &(sys.vals));

	BBFE_sys_monowrap_set_Dirichlet_bc(
			&(sys.mono), sys.fe.total_num_nodes, 3,
			&(sys.bc), sys.mono.mat.R.B);
	BBFE_sys_monowrap_set_Neumann_bc(
			sys.fe.total_num_nodes, 3,
			&(sys.bc), sys.mono.mat.R.B);
	BBFE_sys_monowrap_solve(
			&(sys.mono),
			&(sys.mono_com),
			sys.mono.mat.R.X,
			MONOLIS_ITER_CG,
			MONOLIS_PREC_DIAG,
			sys.vals.mat_max_iter,
			sys.vals.mat_epsilon);

	BBFE_solid_renew_vector(
			sys.vals.u, 
			sys.mono.mat.R.X,
			sys.fe.total_num_nodes);
	/**********************************************/
	
	output_files(&sys, 0, 0, DISPLACEMENT_SCALE);

	BBFE_solid_finalize(&(sys.fe), &(sys.basis), &(sys.bc));
	monolis_finalize(&(sys.mono));
	monolis_global_finalize();

	double t2 = monolis_get_time();
	printf("** Total time: %f\n", t2 - t1);

	printf("\n");

	return 0;
}
