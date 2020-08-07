#include "main.h"
#include "define_cell_type_vtk.h"
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
 * libBB_math
 **********************************************************/
void BB_math_vec3d_cross(
		double        ans[3],
		const double  vec_1[3],
		const double  vec_2[3])
{
	ans[0] = vec_1[1]*vec_2[2] - vec_1[2]*vec_2[1];
	ans[1] = vec_1[2]*vec_2[0] - vec_1[0]*vec_2[2];
	ans[2] = vec_1[0]*vec_2[1] - vec_1[1]*vec_2[0];
}


double BB_math_vec3d_length(
		double vec[3])
{
	double len = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
	return( sqrt(len) );
}

void BB_math_vec3d_normal_vec(
		double vec[3])
{
	double len = BB_math_vec3d_length(vec);
	if(len == 0.0) {
		return;
	}
	else {
		for(int d=0; d<3; d++) {
			vec[d] = vec[d]/len;
		}
	}
}

/**********************************************************
 * libBB_vtk
 **********************************************************/
void BB_vtk_write_header(
		FILE* fp)
{
	fprintf(fp, "# vtk DataFile Version 3.0\n");
	fprintf(fp, "vtk output\n");
	fprintf(fp, "ASCII\n");
}


void BB_vtk_write_points_3d(
		FILE*    fp,
		int      num_points,
		double** x) //[num_points][3]
{
	fprintf(fp, "POINTS %d float\n", num_points);

	for(int i=0; i<num_points; i++) {
		fprintf(fp, "%e %e %e\n", x[i][0], x[i][1], x[i][2]);
	}
}


void BB_vtk_write_cells(
		FILE* fp,
		int   num_cells,
		int   num_points_in_cell,
		int** connectivity) //[num_cells][num_points_in_cell]
{
	fprintf(fp, "CELLS %d %d\n",
			num_cells, num_cells*(num_points_in_cell + 1));
	for(int e=0; e<num_cells; e++) {
		fprintf(fp, "%d ", num_points_in_cell);

		for(int i=0; i<num_points_in_cell; i++) {
			fprintf(fp, "%d ", connectivity[e][i]);
		}
		fprintf(fp, "\n");
	}
}


void BB_vtk_write_cell_types(
		FILE* fp,
		int   num_cells,
		int   elem_type)
{
	fprintf(fp, "CELL_TYPES %d\n", num_cells);

	for(int e=0; e<num_cells; e++) {
		fprintf(fp, "%d\n", elem_type);
	}
}


void BB_vtk_write_point_vals_scalar(
		FILE*        fp,
		double*      val,
		const int    num_points,
		const char*  label)
{
	fprintf(fp, "SCALARS %s float\n", label);
	fprintf(fp, "LOOKUP_TABLE default\n");

	for(int i=0; i<num_points; i++) {
		fprintf(fp, "%e\n", val[i]);
	}
}


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
 * numerical integration
 **********************************************************/
void BBFE_std_integ_set_tet_5(
		int*        num_integ_points,
		double**    integ_point,
		double*     integ_weight)
{
	(*num_integ_points) = 5;

	integ_point[0][0] = 1.0/4.0;  integ_point[0][1] = 1.0/4.0;  integ_point[0][2] = 1.0/4.0;
	integ_point[1][0] = 1.0/6.0;  integ_point[1][1] = 1.0/6.0;  integ_point[1][2] = 1.0/6.0;
	integ_point[2][0] = 1.0/6.0;  integ_point[2][1] = 1.0/6.0;  integ_point[2][2] = 1.0/2.0;
	integ_point[3][0] = 1.0/6.0;  integ_point[3][1] = 1.0/2.0;  integ_point[3][2] = 1.0/6.0;
	integ_point[4][0] = 1.0/2.0;  integ_point[4][1] = 1.0/6.0;  integ_point[4][2] = 1.0/6.0;

	integ_weight[0] = -4.0/30.0;
	integ_weight[1] =  9.0/120.0;
	integ_weight[2] =  9.0/120.0;
	integ_weight[3] =  9.0/120.0;
	integ_weight[4] =  9.0/120.0;
}


double BBFE_std_integ_calc(
		const int      num_integ_points,
		const double*  value,
		const double*  weight,
		const double*  Jacobian)
{
	double val = 0.0;

	for(int i=0; i<num_integ_points; i++) {
		val += value[i] * weight[i] * Jacobian[i];
	}

	return val;
}


/**********************************************************
 * shape function
 **********************************************************/
void BBFE_std_shapefunc_get_val_tet_1st(
		const double    xi[3],
		double*         N)
{
	N[0] = 1.0 - xi[0] - xi[1] - xi[2];
	N[1] = xi[0];
	N[2] = xi[1];
	N[3] = xi[2];
}


void BBFE_std_shapefunc_get_derivative_tet_1st(
		const double    xi[3],
		double*         dN_dxi,
		double*         dN_det,
		double*         dN_dze)
{
	dN_dxi[0] = -1.0;  dN_det[0] = -1.0;  dN_dze[0] = -1.0;
	dN_dxi[1] =  1.0;  dN_det[1] =  0.0;  dN_dze[1] =  0.0;
	dN_dxi[2] =  0.0;  dN_det[2] =  1.0;  dN_dze[2] =  0.0;
	dN_dxi[3] =  0.0;  dN_det[3] =  0.0;  dN_dze[3] =  1.0;
}


void BBFE_std_shapefunc_get_surface_tet_1st(
		int        surf_conn[3],
		const int  surf_num)
{
	switch(surf_num) {
		case 0:
			surf_conn[0] = 2;  surf_conn[1] = 1;  surf_conn[2] = 3;
			break;
		case 1:
			surf_conn[0] = 0;  surf_conn[1] = 2;  surf_conn[2] = 3;
			break;
		case 2:
			surf_conn[0] = 1;  surf_conn[1] = 0;  surf_conn[2] = 3;
			break;
		case 3:
			surf_conn[0] = 0;  surf_conn[1] = 1;  surf_conn[2] = 2;
			break;
	}
}


static double matrixDeterminant_3x3(
		double mat[3][3])
{
	return ( mat[0][0]*mat[1][1]*mat[2][2] +
			mat[0][1]*mat[1][2]*mat[2][0] +
			mat[0][2]*mat[2][1]*mat[1][0] -
			mat[0][2]*mat[1][1]*mat[2][0] -
			mat[0][1]*mat[1][0]*mat[2][2] -
			mat[0][0]*mat[2][1]*mat[1][2]   );
}


static void inverseMatrix_3x3(
		double mat[3][3], /* input matrix */
		double det_mat,   /* determinant of input matrix */
		double invMat[3][3]     /* inverse matrix */)
{
	double invDet = 1.0/det_mat;

	invMat[0][0] = invDet * (mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1]);
	invMat[0][1] = invDet * (mat[0][2]*mat[2][1] - mat[0][1]*mat[2][2]);
	invMat[0][2] = invDet * (mat[0][1]*mat[1][2] - mat[0][2]*mat[1][1]);

	invMat[1][0] = invDet * (mat[1][2]*mat[2][0] - mat[1][0]*mat[2][2]);
	invMat[1][1] = invDet * (mat[0][0]*mat[2][2] - mat[0][2]*mat[2][0]);
	invMat[1][2] = invDet * (mat[0][2]*mat[1][0] - mat[0][0]*mat[1][2]);

	invMat[2][0] = invDet * (mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0]);
	invMat[2][1] = invDet * (mat[0][1]*mat[2][0] - mat[0][0]*mat[2][1]);
	invMat[2][2] = invDet * (mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0]);
}


/**********************************************************
 * mapping
 **********************************************************/


void BBFE_std_mapping_calc_Jacobi_mat_3d(
		double     J[3][3],
		const int  local_num_nodes,
		double**   local_x,
		double*    local_dN_dxi,
		double*    local_dN_det,
		double*    local_dN_dze)
{
	for(int i=0; i<3; i++) {
		for(int j=0; j<3; j++) {
			J[i][j] = 0.0;
		}
	}

	for(int n=0; n<local_num_nodes; n++) {
		for(int i=0; i<3; i++) {
			J[0][i] += local_dN_dxi[n] * local_x[n][i];
			J[1][i] += local_dN_det[n] * local_x[n][i];
			J[2][i] += local_dN_dze[n] * local_x[n][i];
		}
	}
}


void BBFE_std_mapping_integ_point(
		double        x[3],
		const int     local_num_nodes,
		double**      local_x,
		double*       N)
{
	for(int d=0; d<3; d++) {
		x[d] = 0.0;
	}

	for(int i=0; i<local_num_nodes; i++) {
		for(int d=0; d<3; d++) {
			x[d] += N[i]*local_x[i][d];
		}
	}
}


/**********************************************************
 * surface
 **********************************************************/

int BBFE_std_surface_get_surface_node(
		bool*    node_is_on_surface,
		FE_DATA* fe)
{
	double** norm;
	norm = BB_std_calloc_2d_double(norm, fe->total_num_nodes, DIM);

	for(int e=0; e<(fe->total_num_elems); e++) {
		for(int i=0; i<4; i++) {
			int surf_conn[3];
			BBFE_std_shapefunc_get_surface_tet_1st(
					surf_conn, i);

			double ans[3];  double vec_1[3];  double vec_2[3];
			int nid_0 = fe->conn[e][ surf_conn[0] ];
			int nid_1 = fe->conn[e][ surf_conn[1] ];
			int nid_2 = fe->conn[e][ surf_conn[2] ];
			for(int d=0; d<3; d++) {
				vec_1[d] = fe->x[ nid_1 ][d] - fe->x[ nid_0 ][d];
				vec_2[d] = fe->x[ nid_2 ][d] - fe->x[ nid_0 ][d];
			}
			BB_math_vec3d_cross(ans, vec_1, vec_2);
			BB_math_vec3d_normal_vec(ans);

			for(int d=0; d<3; d++) {
				norm[ nid_0 ][d] += ans[d];
				norm[ nid_1 ][d] += ans[d];
				norm[ nid_2 ][d] += ans[d];
			}

		}
	}
	
	int num_surface_nodes = 0;
	for(int i=0; i<(fe->total_num_nodes); i++) {
		double len = BB_math_vec3d_length(norm[i]);
		if(len > 1.0) {
			num_surface_nodes++;
			node_is_on_surface[i] = true;
		}
	}
	
	BB_std_free_2d_double(norm, fe->total_num_nodes, DIM);

	return num_surface_nodes;
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

			fe->geo[e][p].Jacobian = matrixDeterminant_3x3(
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
			inverseMatrix_3x3(
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
	num_bcs = BBFE_std_surface_get_surface_node(node_is_on_surface, fe);

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
	BBFE_std_integ_set_tet_5(
			&(basis->num_integ_points),
			basis->integ_point,
			basis->integ_weight
			);

	for(int i=0; i<(basis->num_integ_points); i++) {
		BBFE_std_shapefunc_get_val_tet_1st(
				basis->integ_point[i],
				basis->N[i]);

		BBFE_std_shapefunc_get_derivative_tet_1st(
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
