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
 * memory allocation
 **********************************************************/


void BBFE_std_memory_allocation_basis(
		FE_3D_BASIS*    basis,
		const int       num_integ_points,
		const int       pol_order,
		const int       num_nodes_in_elem)
{
	BBFE_std_memory_allocation_integ(
			basis,
			num_integ_points);

	BBFE_std_memory_allocation_shapefunc(
			basis,
			num_nodes_in_elem,
			pol_order,
			num_integ_points);
}


void BBFE_std_memory_allocation_integ(
		FE_3D_BASIS*    basis,
		const int       num_integ_points)
{
	int num = num_integ_points;

	basis->num_integ_points = num;

	basis->integ_point  = BB_std_calloc_2d_double(basis->integ_point , num, DIM);
	basis->integ_weight = BB_std_calloc_1d_double(basis->integ_weight, num);
}


void BBFE_std_memory_allocation_shapefunc(
		FE_3D_BASIS*    basis,
		const int       num_nodes_in_elem,
		const int       pol_order,
		const int       num_integ_points)
{
	int nn = num_nodes_in_elem;
	int ni = num_integ_points;

	basis->num_nodes = nn;
	basis->pol_order = pol_order;

	basis->N      = BB_std_calloc_2d_double(basis->N     , ni, nn);
	basis->dN_dxi = BB_std_calloc_2d_double(basis->dN_dxi, ni, nn);
	basis->dN_det = BB_std_calloc_2d_double(basis->dN_det, ni, nn);
	basis->dN_dze = BB_std_calloc_2d_double(basis->dN_dze, ni, nn);
}


void BBFE_std_memory_allocation_node(
		FE_DATA*  fe)
{
	fe->x = BB_std_calloc_2d_double(fe->x, fe->total_num_nodes, DIM);
}


void BBFE_std_memory_allocation_elem(
		FE_DATA*  fe,
		int       num_integ_points)
{
	fe->conn = BB_std_calloc_2d_int(fe->conn, fe->total_num_elems, fe->local_num_nodes);
	
	fe->geo  = (FE_3D_GEO**)calloc(fe->total_num_elems, sizeof(FE_3D_GEO*));
	for(int e=0; e<(fe->total_num_elems); e++) {
		fe->geo[e]  = (FE_3D_GEO*)calloc(num_integ_points, sizeof(FE_3D_GEO));

		for(int p=0; p<num_integ_points; p++) {
			fe->geo[e][p].grad_N = BB_std_calloc_2d_double(fe->geo[e][p].grad_N, fe->local_num_nodes, DIM);
		}
	}
}


/**********************************************************
 * read
 **********************************************************/


void BBFE_sys_read_node(
		FE_DATA*     fe,
		const char*  filename)
{
	FILE* fp;

	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s ERROR: File \"%s\" cannot be opened.\n",
				CODENAME, filename);
	}

	// read the number of nodes
	BB_std_scan_line(
			&fp, BUFFER_SIZE, "%d", &(fe->total_num_nodes));
	printf("%s Num. nodes: %d\n", CODENAME, fe->total_num_nodes);
	BBFE_std_memory_allocation_node(fe);

	// read positions of nodes
	for(int i=0; i<(fe->total_num_nodes); i++) {
		BB_std_scan_line(&fp, BUFFER_SIZE,
				"%lf %lf %lf", &(fe->x[i][0]), &(fe->x[i][1]), &(fe->x[i][2]));
	}

	fclose(fp);
}


void BBFE_sys_read_elem(
		FE_DATA*     fe,
		const char*  filename,
		int          num_integ_points)
{
	FILE* fp;

	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s ERROR: File \"%s\" cannot be opened.\n",
				CODENAME, filename);
	}

	// read the number of elements
	BB_std_scan_line(
			&fp, BUFFER_SIZE, "%d %d",&(fe->total_num_elems), &(fe->local_num_nodes));
	printf("%s Num. elements: %d\n", CODENAME, fe->total_num_elems);
	BBFE_std_memory_allocation_elem(fe, num_integ_points);

	// read the connectivities of elements
	for(int e=0; e<(fe->total_num_elems); e++) {
		for(int i=0; i<(fe->local_num_nodes); i++) {
			fscanf(fp, "%d", &(fe->conn[e][i]));
		}
	}

	fclose(fp);
}


void BBFE_sys_read_Dirichlet_bc(
		BC_DATA*     bc,
		const char*  filename,
		const int    total_num_nodes)
{
	FILE* fp;
	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s ERROR: File \"%s\" cannot be opened.\n",
				CODENAME, filename);
	}

	bc->total_num_nodes = total_num_nodes;

	BB_std_scan_line(&fp, BUFFER_SIZE,
			"%d %d", &(bc->num_D_bcs), &(bc->block_size));
	printf("%s Num. Dirichlet B.C.: %d\n", CODENAME, bc->num_D_bcs);

	int n = total_num_nodes * bc->block_size;

	bc->D_bc_exists   = BB_std_calloc_1d_bool(  bc->D_bc_exists  , n);
	bc->imposed_D_val = BB_std_calloc_1d_double(bc->imposed_D_val, n);
	for(int i=0; i<n; i++) {
		bc->D_bc_exists[i]   = false;
		bc->imposed_D_val[i] = 0.0;
	}

	for(int i=0; i<(bc->num_D_bcs); i++) {
		int node_id;  int block_id;  double val;
		BB_std_scan_line(&fp, BUFFER_SIZE,
				"%d %d %lf", &node_id, &block_id, &val);

		int index = (bc->block_size)*node_id + block_id;
		bc->D_bc_exists[ index ]   = true;
		bc->imposed_D_val[ index ] = val;
	}
}


/**********************************************************
 * write
 **********************************************************/


void BBFE_sys_write_vtk_shape(
		FILE*     fp,
		FE_DATA*  fe,
		const int cell_type)
{
	BB_vtk_write_header(fp);
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
	BB_vtk_write_points_3d(fp, fe->total_num_nodes, fe->x);
	BB_vtk_write_cells(fp, fe->total_num_elems, fe->local_num_nodes, fe->conn);
	BB_vtk_write_cell_types(fp, fe->total_num_elems, cell_type);

}




void BBFE_write_ascii_nodal_vals_scalar(
		FE_DATA*     fe,
		double*      vals,
		const char*  filename)
{
	FILE* fp;
	fp = fopen(filename, "w");
	if( fp == NULL ) {
		printf("%s ERROR: File \"%s\" cannot be opened.\n",
				CODENAME, filename);
	}

	fprintf(fp, "%d\n", fe->total_num_nodes);

	for(int i=0; i<(fe->total_num_nodes); i++) {
		fprintf(fp, "%e\n", vals[i]);
	}

	fclose(fp);

}


/**********************************************************
 * equivval
 **********************************************************/

static void get_Jacobian_array(
		double*    Jacobian_ip,
		const int  num_integ_points,
		const int  elem_num,
		FE_DATA*   fe)
{
	for(int p=0; p<num_integ_points; p++) {
		Jacobian_ip[p] = fe->geo[ elem_num ][p].Jacobian;
	}
}


void BBFE_std_equivval_volume_smooth_function(
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
			get_Jacobian_array(
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
 * monolis wrapper
 **********************************************************/


void BBFE_sys_monowrap_set_Dirichlet_bc(
		MONOLIS*      monolis,
		int           num_nodes,
		int           num_dofs_on_node,
		BC_DATA*      bc,
		double*       g_rhs)
	{
	for(int i=0; i<num_nodes; i++) {
		for(int k=0; k<num_dofs_on_node; k++) {
			if( bc->D_bc_exists[ num_dofs_on_node*i+k ] ) {
				monolis_set_Dirichlet_bc(
						monolis,
						g_rhs,
						i,
						k,
						bc->imposed_D_val[ num_dofs_on_node*i+k ]);
			}
		}
	}
}


/**********************************************************
 * element matrix
 **********************************************************/

void BBFE_elemmat_set_Jacobi_mat(
		FE_DATA*     fe,
		FE_3D_BASIS* basis)
{
	double** local_x;
	local_x = BB_std_calloc_2d_double(local_x, fe->local_num_nodes, DIM);

	for(int e=0; e<(fe->total_num_elems); e++) {

		for(int i=0; i<(fe->local_num_nodes); i++) {
			local_x[i][0] = fe->x[ fe->conn[e][i] ][0];
			local_x[i][1] = fe->x[ fe->conn[e][i] ][1];
			local_x[i][2] = fe->x[ fe->conn[e][i] ][2];
		}

		for(int p=0; p<(basis->num_integ_points); p++) {
			BBFE_std_mapping_calc_Jacobi_mat_3d(
					fe->geo[e][p].J,
					fe->local_num_nodes,
					local_x,
					basis->dN_dxi[p],
					basis->dN_det[p],
					basis->dN_dze[p]);

			fe->geo[e][p].Jacobian = BB_calc_mat3d_determinant(
					fe->geo[e][p].J);
		}
	}

	BB_std_free_2d_double(local_x, fe->local_num_nodes, DIM);
}


void BBFE_elemmat_set_shapefunc_derivative(
		FE_DATA*      fe,
		FE_3D_BASIS*  basis)
{
	double J_inv[3][3];

	for(int e=0; e<(fe->total_num_elems); e++) {
		for(int p=0; p<(basis->num_integ_points); p++) {
			BB_calc_mat3d_inverse(
					fe->geo[e][p].J, fe->geo[e][p].Jacobian, J_inv);

			for(int i=0; i<(fe->local_num_nodes); i++) {
				for(int d=0; d<3; d++) {
					fe->geo[e][p].grad_N[i][d] =
						J_inv[d][0] * basis->dN_dxi[p][i] +
						J_inv[d][1] * basis->dN_det[p][i] +
						J_inv[d][2] * basis->dN_dze[p][i];
				}
			}
		}
	}
}


void BBFE_elemmat_set_element_mat(
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

				get_Jacobian_array(
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
 * manufactured solution
 **********************************************************/

void BBFE_std_manusol_calc_nodal_error_scalar(
		FE_DATA*      fe,
		double*       error,
		const double* val)
{
	for(int i=0; i<fe->total_num_nodes; i++) {
		double x[3];
		for(int d=0; d<3; d++) {
			x[d] = fe->x[i][d];
		}
		double theo_sol = BBFE_std_manusol_get_sol_scalar_3d(x[0], x[1], x[2]);
		error[i] = val[i] - theo_sol;
	}

}


void BBFE_std_manusol_overwrite_bc_file(
		FE_DATA*   fe,
		const int  block_size)
{
	bool* node_is_on_surface;
	node_is_on_surface = BB_std_calloc_1d_bool(node_is_on_surface, fe->total_num_nodes);
	for(int i=0; i<(fe->total_num_nodes); i++) {
		node_is_on_surface[i] = false;
	}
	
	int num_bcs;
	num_bcs = BBFE_std_surface_tet1st_get_surface_node(
			node_is_on_surface, 
			fe->total_num_nodes, fe->x,
			fe->total_num_elems, fe->conn);

	const char* filename = "bc_D.dat";
	FILE* fp;
	fp = fopen(filename, "w");
	if( fp == NULL ) {
		printf("%s ERROR: File \"%s\" cannot be opened.\n",
				CODENAME, filename);
	}
	fprintf(fp, "%d %d\n", block_size*num_bcs, block_size);
	for(int i=0; i<(fe->total_num_nodes); i++) {
		if(node_is_on_surface[i]) {
			for(int b=0; b<block_size; b++) {
				fprintf(fp, "%d %d %e\n", i, b, 0.0);
			}
		}
	}

	fclose(fp);

	BB_std_free_1d_bool(node_is_on_surface, fe->total_num_nodes);
}


double BBFE_std_manusol_get_sol_scalar_3d(
		double x,
		double y,
		double z)
{
	double sol = sin( 0.5*x ) * sin( 1.0*y ) * sin( 2.0*z );
	return sol;
}


double BBFE_std_manusol_get_rhs_scalar_3d(
		double x,
		double y,
		double z)
{
	double rhs = -5.25*sin( 0.5*x ) * sin( y ) * sin( 2.0*z );
	return rhs;
}


void BBFE_std_manusol_set_bc_scalar(
		FE_DATA* fe,
		BC_DATA* bc)
{
	for(int i=0; i<(fe->total_num_nodes); i++) {
		if( bc->D_bc_exists[i] ) {
			double x[3];
			for(int d=0; d<3; d++) {
				x[d] = fe->x[i][d];
			}

			bc->imposed_D_val[i] =
				BBFE_std_manusol_get_sol_scalar_3d(x[0], x[1], x[2]);
		}
	}
}


void BBFE_std_manusol_add_rhs_scalar(
		FE_DATA* fe,
		FE_3D_BASIS* basis,
		double* rhs)
{
	double* equiv_val;
	equiv_val = BB_std_calloc_1d_double(equiv_val, fe->total_num_nodes);
	
	BBFE_std_equivval_volume_smooth_function(
			equiv_val, fe, basis, BBFE_std_manusol_get_rhs_scalar_3d);

	BB_std_free_1d_double(equiv_val, fe->total_num_nodes);
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

	BBFE_std_memory_allocation_basis(
			&(sys.basis),
			NUM_INTEG_POINTS,
			POL_ORDER,
			NUM_NODES_IN_ELEM);

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
	BBFE_std_manusol_overwrite_bc_file(
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
	BBFE_std_manusol_set_bc_scalar(
			&(sys.fe),
			&(sys.bc));
	BBFE_std_manusol_add_rhs_scalar(
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

	BBFE_std_manusol_calc_nodal_error_scalar(
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
