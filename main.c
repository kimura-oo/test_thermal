#include "main.h"
#include "monolis.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

typedef struct
{
	FE_3D_BASIS  basis;
	FE_DATA      fe;
	NODAL_VALUES vals;
	BC_DATA      bc;

} FE_SYSTEM;

const int DIM = 3;

const int NUM_INTEG_POINTS   = 5;
const int POL_ORDER          = 1;
const int NUM_NODES_IN_ELEM  = 4;
const int BLOCK_SIZE         = 1;

const double MAT_EPSILON  = 1.0e-10;
const double MAT_MAX_ITER = 100000;

const char* CODENAME = "test_thermal >";
const int BUFFER_SIZE = 1000;

const char* INPUT_FILENAME_NODE        = "./util/meshgen/node.dat";
const char* INPUT_FILENAME_ELEM        = "./util/meshgen/elem.dat";
const char* OUTPUT_FILENAME_VTK        = "result.vtk";
const char* OUTPUT_FILENAME_ASCII_TEMP = "temparature.dat";
const char* OUTPUT_FILENAME_ASCII_RHS  = "rhs.dat";


/**********************************************************
 * equivval
 **********************************************************/

void BBFE_elemmat_equivval_volume_smooth_function(
		double* equiv_val,
		FE_DATA* fe,
		FE_3D_BASIS* basis,
		double (*func)(double, double, double)) // scalar function(x, y, z)
{

	for(int i=0; i<(fe->total_num_nodes); i++) {
		equiv_val[i] = 0.0;
	}

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip      , basis->num_integ_points);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip , basis->num_integ_points);

	double** local_x;
	local_x = BB_std_calloc_2d_double(local_x, fe->local_num_nodes, DIM);

	for(int e=0; e<(fe->total_num_elems); e++) {
		for(int i=0; i<(fe->local_num_nodes); i++) {
			local_x[i][0] = fe->x[ fe->conn[e][i] ][0];
			local_x[i][1] = fe->x[ fe->conn[e][i] ][1];
			local_x[i][2] = fe->x[ fe->conn[e][i] ][2];
		}

		for(int i=0; i<(fe->local_num_nodes); i++) {
			BBFE_elemmat_set_Jacobian_array(
					Jacobian_ip,
					basis->num_integ_points,
					e,
					fe);

			for(int p=0; p<(basis->num_integ_points); p++) {
				double x_ip[3];
				BBFE_std_mapping_integ_point(
						x_ip,
						fe->local_num_nodes,
						local_x,
						basis->N[p]);

				double rhs_ip = func(x_ip[0], x_ip[1], x_ip[2]);
				val_ip[p] = rhs_ip * basis->N[p][i];
			}

			double integ_val = BBFE_std_integ_calc(
					basis->num_integ_points,
					val_ip,
					basis->integ_weight,
					Jacobian_ip);

			equiv_val[ fe->conn[e][i] ] += integ_val;
		}
	}

	BB_std_free_1d_double(val_ip,     basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x, fe->local_num_nodes, DIM);
}

/**********************************************************
 * element matrix
 **********************************************************/



double BBFE_elemmat_thermal_steady_linear(
		double grad_N_i[3],
		double grad_N_j[3])
{
	double val = -(  
			grad_N_i[0] * grad_N_j[0] +
			grad_N_i[1] * grad_N_j[1] +
			grad_N_i[2] * grad_N_j[2]);

	return val;
}


/**********************************************************
 * main function
 **********************************************************/


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


int main (
		int argc,
		char* argv[])
{
	printf("\n");

	FE_SYSTEM sys;
	MONOLIS monolis;

	double t1 = monolis_get_time();

	monolis_global_initialize();
	monolis_initialize(&monolis);

	BBFE_sys_memory_allocation_basis(
			&(sys.basis),
			NUM_INTEG_POINTS,
			POL_ORDER,
			NUM_NODES_IN_ELEM,
			3);

	BBFE_sys_read_node(
			&(sys.fe),
			INPUT_FILENAME_NODE);
	BBFE_sys_read_elem(
			&(sys.fe),
			INPUT_FILENAME_ELEM,
			sys.basis.num_integ_points);

	memory_allocation_nodal_values(
			&(sys.vals),
			sys.fe.total_num_nodes);

	// for manufactured solution
	BBFE_manusol_overwrite_bc_file(
			&(sys.fe),
			BLOCK_SIZE);

	BBFE_sys_read_Dirichlet_bc(
			&(sys.bc),
			"bc_D.dat",
			sys.fe.total_num_nodes);

	set_basis(
			&(sys.basis));

	monolis_get_nonzero_pattern(
			&monolis,
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
	BBFE_elemmat_set_element_mat(
			&monolis,
			&(sys.fe),
			&(sys.basis));

	// for manufactured solution
	BBFE_manusol_set_bc_scalar(
			&(sys.fe),
			&(sys.bc));
	BBFE_manusol_add_rhs_scalar(
			&(sys.fe),
			&(sys.basis),
			monolis.mat.B);

	BBFE_sys_monowrap_set_Dirichlet_bc(
			&monolis,
			sys.fe.total_num_nodes,
			BLOCK_SIZE,
			&(sys.bc),
			monolis.mat.B);

	monolis_set_method   (&monolis, monolis_iter_CG);
	monolis_set_precond  (&monolis, monolis_prec_DIAG);
	monolis_set_maxiter  (&monolis, MAT_MAX_ITER);
	monolis_set_tolerance(&monolis, MAT_EPSILON);

	monolis_solve(
			&monolis,
			monolis.mat.B,
			sys.vals.T);

	BBFE_manusol_calc_nodal_error_scalar(
			&(sys.fe), sys.vals.error, sys.vals.T);

	output_result_file_vtk(
			&(sys.fe),
			&(sys.vals),
			OUTPUT_FILENAME_VTK);

	BBFE_write_ascii_nodal_vals_scalar(
			&(sys.fe),
			sys.vals.T,
			OUTPUT_FILENAME_ASCII_TEMP);

	BBFE_write_ascii_nodal_vals_scalar(
			&(sys.fe),
			monolis.mat.B,
			OUTPUT_FILENAME_ASCII_RHS);

	monolis_finalize(&monolis);
	monolis_global_finalize();

	double t2 = monolis_get_time();
	printf("** Total time: %f\n", t2 - t1);

	printf("\n");

	return 0;
}
