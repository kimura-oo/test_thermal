#include "main.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>


const int DIM = 3;

const int NUM_INTEG_POINTS   = 5;
const int POL_ORDER          = 1;
const int NUM_NODES_IN_ELEM  = 4;
const int BLOCK_SIZE         = 1;

const double MAT_EPSILON  = 1.0e-10;
const double MAT_MAX_ITER = 100000;

const char* CODENAME = "test_thermal >";

const char* INPUT_FILENAME_NODE        = "node.dat";
const char* INPUT_FILENAME_ELEM        = "elem.dat";
const char* INPUT_FILENAME_D_BC        = "D_bc.dat";
const char* OUTPUT_FILENAME_VTK        = "result.vtk";
const char* OUTPUT_FILENAME_ASCII_TEMP = "temparature.dat";
const char* OUTPUT_FILENAME_ASCII_RHS  = "rhs.dat";



typedef struct
{
	double* T;
	double* error;

} NODAL_VALUES;


typedef struct
{
	char* data_directory;
	char* filename_node;
	char* filename_elem;
	char* filename_D_bc;
	char* filename_N_bc;

} CONDITIONS;

typedef struct
{
	FE_3D_BASIS  basis;
	FE_DATA      fe;
	NODAL_VALUES vals;
	BC_DATA      bc;
	MONOLIS      monolis;

} FE_SYSTEM;


void memory_allocation_nodal_values(
		NODAL_VALUES*   vals,
		const int       total_num_nodes)
{
	vals->T     = BB_std_calloc_1d_double(vals->T,     total_num_nodes);
	vals->error = BB_std_calloc_1d_double(vals->error, total_num_nodes);
}


void set_basis(
		FE_3D_BASIS*  basis)
{
	BBFE_std_integ_tet_set_5points(
			&(basis->num_integ_points),
			basis->integ_point,
			basis->integ_weight
			);

	for(int i=0; i<(basis->num_integ_points); i++) {
		BBFE_std_shapefunc_tet1st_get_val(
				basis->integ_point[i],
				basis->N[i]);

		BBFE_std_shapefunc_tet1st_get_derivative(
				basis->integ_point[i],
				basis->dN_dxi[i],
				basis->dN_det[i],
				basis->dN_dze[i]);
	}

}


void output_result_file_vtk(
		FE_DATA*       fe,
		NODAL_VALUES*  vals,
		const char*    filename)
{
	FILE* fp;
	fp = fopen(filename, "w");
	if( fp == NULL ) {
		printf("%s ERROR: File \"%s\" cannot be opened.\n",
				CODENAME, filename);
	}

	BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_TETRA);
	
	fprintf(fp, "POINT_DATA %d\n", fe->total_num_nodes);
	BB_vtk_write_point_vals_scalar(fp, vals->T, fe->total_num_nodes, "temperature");
	BB_vtk_write_point_vals_scalar(fp, vals->error, fe->total_num_nodes, "abs_error");

	fclose(fp);

}


const char* get_directory_name(
		int   argc,
		char* argv[])
{
	const char* dir_name;

	if(argc < 2) { dir_name = "."; }
	else         { dir_name = argv[1]; }

	printf("%s Main directory: %s\n", CODENAME, dir_name);

	return dir_name;
}


void set_element_mat(
		MONOLIS*     monolis,
		FE_DATA*     fe,
		FE_3D_BASIS* basis)
{
	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip      , basis->num_integ_points);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip , basis->num_integ_points);

	for(int e=0; e<(fe->total_num_elems); e++) {
		for(int i=0; i<(fe->local_num_nodes); i++) {
			for(int j=0; j<(fe->local_num_nodes); j++) {

				BBFE_elemmat_set_Jacobian_array(
						Jacobian_ip,
						basis->num_integ_points,
						e,
						fe);

				for(int p=0; p<(basis->num_integ_points); p++) {
					val_ip[p] =
						BBFE_elemmat_thermal_steady_linear(
							fe->geo[e][p].grad_N[i],
							fe->geo[e][p].grad_N[j]);
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

	BB_std_free_1d_double(val_ip,     basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);
}


int main (
		int argc,
		char* argv[])
{
	printf("\n");

	FE_SYSTEM sys;

	double t1 = monolis_get_time();
	
	monolis_global_initialize();

	const char* dir_name = get_directory_name(argc, argv);	

	BBFE_sys_memory_allocation_basis(
			&(sys.basis),
			NUM_INTEG_POINTS,
			POL_ORDER,
			NUM_NODES_IN_ELEM,
			3);

	BBFE_sys_read_node(
			&(sys.fe),
			INPUT_FILENAME_NODE,
			dir_name);
	BBFE_sys_read_elem(
			&(sys.fe),
			INPUT_FILENAME_ELEM,
			dir_name,
			sys.basis.num_integ_points);

	memory_allocation_nodal_values(
			&(sys.vals),
			sys.fe.total_num_nodes);

	// for manufactured solution
	BBFE_manusol_overwrite_bc_file(
			&(sys.fe),
			BLOCK_SIZE, 
			INPUT_FILENAME_D_BC,
			dir_name);

	BBFE_sys_read_Dirichlet_bc(
			&(sys.bc),
			INPUT_FILENAME_D_BC,
			dir_name,
			sys.fe.total_num_nodes);

	set_basis(
			&(sys.basis));

	monolis_initialize(&(sys.monolis));
	monolis_get_nonzero_pattern(
			&(sys.monolis),
			sys.fe.total_num_nodes,
			sys.fe.local_num_nodes,
			1,
			sys.fe.total_num_elems,
			sys.fe.conn);

	BBFE_elemmat_set_Jacobi_mat(
			&(sys.fe),
			&(sys.basis));
	BBFE_elemmat_set_shapefunc_derivative(
			&(sys.fe),
			&(sys.basis));
	set_element_mat(
			&(sys.monolis),
			&(sys.fe),
			&(sys.basis));

	// for manufactured solution
	BBFE_manusol_set_bc_scalar(
			&(sys.fe),
			&(sys.bc),
			0.0);
	BBFE_manusol_add_rhs_scalar(
			&(sys.fe),
			&(sys.basis),
			0.0,
			sys.monolis.mat.B);

	BBFE_sys_monowrap_set_Dirichlet_bc(
			&(sys.monolis),
			sys.fe.total_num_nodes,
			BLOCK_SIZE,
			&(sys.bc),
			sys.monolis.mat.B);

	monolis_set_method   (&(sys.monolis), monolis_iter_BiCGSTAB_noprec);
	monolis_set_precond  (&(sys.monolis), monolis_prec_DIAG);
	monolis_set_maxiter  (&(sys.monolis), MAT_MAX_ITER);
	monolis_set_tolerance(&(sys.monolis), MAT_EPSILON);

	monolis_solve(
			&(sys.monolis),
			sys.monolis.mat.B,
			sys.vals.T);

	BBFE_manusol_calc_nodal_error_scalar(
			&(sys.fe), 0.0, sys.vals.error, sys.vals.T);

	output_result_file_vtk(
			&(sys.fe),
			&(sys.vals),
			OUTPUT_FILENAME_VTK);

	BBFE_write_ascii_nodal_vals_scalar(
			&(sys.fe),
			sys.vals.T,
			OUTPUT_FILENAME_ASCII_TEMP,
			".");
	BBFE_write_ascii_nodal_vals_scalar(
			&(sys.fe),
			sys.monolis.mat.B,
			OUTPUT_FILENAME_ASCII_RHS,
			".");
	
	double L2_error = BBFE_elemmat_equivval_relative_L2_error_scalar(
			&(sys.fe),
			&(sys.basis),
			0.0,
			sys.vals.T,
			BBFE_manusol_get_sol_scalar_3d);
	printf("L2 error: %e\n", L2_error);

	monolis_finalize(&(sys.monolis));
	monolis_global_finalize();

	double t2 = monolis_get_time();
	printf("** Total time: %f\n", t2 - t1);

	printf("\n");

	return 0;
}
