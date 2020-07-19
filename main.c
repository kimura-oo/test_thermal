#include "main.h"
#include "define_cell_type_vtk.h"

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

typedef struct
{
	FE_3D_BASIS  basis;
	FE_DATA      fe;
	NODAL_VALUES vals;
	Dataset_CSR  csr;

} FE_SYSTEM;

const int DIM = 3;

const int NUM_INTEG_POINTS  = 5;
const int POL_ORDER         = 1;
const int NUM_NODES_IN_ELEM = 4;

const double MAT_EPSILON  = 1.0e-10;
const double MAT_MAX_ITER = 100000;

const char* CODENAME = "test_thermal >";
const int BUFFER_SIZE = 1000;

const char* OUTPUT_FILENAME = "result.vtk";


/**********************************************************
 * memory allocation
 **********************************************************/
void memory_allocation_basis(
		FE_3D_BASIS*    basis,
		const int       num_integ_points,
		const int       pol_order,
		const int       num_nodes_in_elem)
{
	memory_allocation_integ(
			basis, 
			num_integ_points);

	memory_allocation_shapefunc(
			basis,
			num_nodes_in_elem,
			pol_order,
			num_integ_points);
}


void memory_allocation_integ(
		FE_3D_BASIS*    basis,
		const int       num_integ_points)
{
	int num = num_integ_points;

	basis->num_integ_points = num;

	basis->integ_point  = (double**)calloc(num, sizeof(double*));
	basis->integ_weight = (double* )calloc(num, sizeof(double ));
	for(int i=0; i<num; i++) {
		basis->integ_point[i] = (double*)calloc(DIM, sizeof(double));
	}
}


void memory_allocation_shapefunc(
		FE_3D_BASIS*    basis,
		const int       num_nodes_in_elem,
		const int       pol_order,
		const int       num_integ_points)
{
	int nn = num_nodes_in_elem;
	int ni = num_integ_points;

	basis->num_nodes = nn;
	basis->pol_order = pol_order;

	basis->N      = (double**)calloc(ni, sizeof(double*));
	basis->dN_dxi = (double**)calloc(ni, sizeof(double*));
	basis->dN_det = (double**)calloc(ni, sizeof(double*));
	basis->dN_dze = (double**)calloc(ni, sizeof(double*));
	for(int i=0; i<ni; i++) {
		basis->N[i]      = (double*)calloc(nn, sizeof(double));
		basis->dN_dxi[i] = (double*)calloc(nn, sizeof(double));
		basis->dN_det[i] = (double*)calloc(nn, sizeof(double));
		basis->dN_dze[i] = (double*)calloc(nn, sizeof(double));
	}
}


void memory_allocation_nodal_values(
		NODAL_VALUES*   vals,
		const int       total_num_nodes)
{
	vals->T = (double*)calloc(total_num_nodes, sizeof(double));
}


/**********************************************************
 * initializers
 **********************************************************/
void initialize_basis(
		FE_3D_BASIS*  basis)
{
	integ_point_tet_5(
			&(basis->num_integ_points),
			basis->integ_point,
			basis->integ_weight
			);

	for(int i=0; i<(basis->num_integ_points); i++) {
		shapefunc_3d_tet_1st_value(
				basis->integ_point[i],
				basis->N[i]);

		shapefunc_3d_tet_1st_der_value(
				basis->integ_point[i],
				basis->dN_dxi[i],
				basis->dN_det[i],
				basis->dN_dze[i]);
	}

}


/**********************************************************
 * input
 **********************************************************/
static bool BEBOPS_IO_scan_line(
		FILE** fp,
		const int buffer_size,
		const char* format,
		...)
{
	char buf[buffer_size];
	if( fgets(buf, sizeof(buf), *fp) == NULL )
	{
		return false;
	}

	va_list va;
	va_start(va, format);
	vsscanf(buf, format, va);
	va_end(va);

	return true;
}


static bool BEBOPS_IO_read_file_return_char(
		char* ret_char,
		const char* filename,
		const char* identifier,
		const int buffer_size)
{
	int identical = 0;

	FILE* fp;
	fp = fopen(filename, "r");
	if(fp == NULL) {
		ret_char = NULL;
		return false;
	}

	char buf[ buffer_size];
	char buf2[buffer_size];
	char buf3[buffer_size];
	while(1) {
		if( fgets(buf, sizeof(buf), fp) == NULL) {
			break;
		}
		strcpy(buf2, buf);
		if( strstr(buf2, identifier) != NULL ) {
			sscanf(buf, "%s %s", buf3, ret_char);

			return true;
		}
	}

	return false;
}


static void memory_allocation_node(
		FE_DATA*  fe)
{
	fe->x = (double**)calloc(fe->total_num_nodes, sizeof(double*));
	
	for(int i=0; i<(fe->total_num_nodes); i++) {
		fe->x[i] = (double*)calloc(3, sizeof(double));
	}
}


static void memory_allocation_elem(
		FE_DATA*  fe,
		int       num_integ_points)
{
	fe->conn = (int**)calloc(fe->total_num_elems, sizeof(int*));
	fe->geo  = (FE_3D_GEO**)calloc(fe->total_num_elems, sizeof(FE_3D_GEO*));

	for(int e=0; e<(fe->total_num_elems); e++) {
		fe->conn[e] = (int*)calloc(fe->local_num_nodes, sizeof(int));
		fe->geo[e]  = (FE_3D_GEO*)calloc(num_integ_points, sizeof(FE_3D_GEO));
		
		for(int p=0; p<num_integ_points; p++) {
			fe->geo[e][p].grad_N = (double**)calloc(fe->local_num_nodes, sizeof(double*));

			for(int i=0; i<(fe->local_num_nodes); i++) {
				fe->geo[e][p].grad_N[i] = (double*)calloc(3, sizeof(double));
			}
		}
	}
}

void read_and_memory_allocation_FE_data(
		FE_DATA*  fe,
		char*     filename, 
		int       num_integ_points)
{
	FILE* fp;

	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s ERROR: File \"%s\" cannot be opened.\n", 
				CODENAME, filename);
	}

	// read the number of nodes
	BEBOPS_IO_scan_line(
			&fp, BUFFER_SIZE, "%d", &(fe->total_num_nodes));
	printf("%s Num. nodes: %d\n", CODENAME, fe->total_num_nodes);
	memory_allocation_node(fe);

	// read positions of nodes
	for(int i=0; i<(fe->total_num_nodes); i++) {
		BEBOPS_IO_scan_line(&fp, BUFFER_SIZE, 
				"%lf %lf %lf", &(fe->x[i][0]), &(fe->x[i][1]), &(fe->x[i][2]));
	}

	// read the number of elements
	BEBOPS_IO_scan_line(
			&fp, BUFFER_SIZE, "%d %d",&(fe->total_num_elems), &(fe->local_num_nodes));
	printf("%s Num. elements: %d\n", CODENAME, fe->total_num_elems);
	memory_allocation_elem(fe, num_integ_points);

	// read the connectivities of elements
	for(int e=0; e<(fe->total_num_elems); e++) {
		for(int i=0; i<(fe->local_num_nodes); i++) {	
			fscanf(fp, "%d", &(fe->conn[e][i]));
		}
	}

	fclose(fp);
}



/**********************************************************
 * output
 **********************************************************/
void write_vtk_shape(
		FE_DATA*  fe,
		FILE*     fp)
{
	fprintf(fp, "# vtk DataFile Version 3.0\n");
	fprintf(fp, "vtk output\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

	fprintf(fp, "POINTS %d float\n", fe->total_num_nodes);

	for(int i=0; i<(fe->total_num_nodes); i++) {
		fprintf(fp, "%e %e %e\n", fe->x[i][0], fe->x[i][1], fe->x[i][2]);
	}

	int num_nodes_in_cell = 1;
	fprintf(fp, "CELLS %d %d\n", 
			fe->total_num_elems, fe->total_num_elems*(fe->local_num_nodes + 1));
	for(int e=0; e<(fe->total_num_elems); e++) {
		fprintf(fp, "%d ", fe->local_num_nodes);

		for(int i=0; i<(fe->local_num_nodes); i++) {
			fprintf(fp, "%d ", fe->conn[e][i]);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "CELL_TYPES %d\n", fe->total_num_elems);
	int elem_type = TYPE_VTK_TETRA ;

	for(int e=0; e<(fe->total_num_elems); e++) {
		fprintf(fp, "%d\n", elem_type);
	}

}

void write_nodal_value_scalar(
		FE_DATA*     fe,
		FILE*        fp,
		double*      val,
		const char*  label)
{
	fprintf(fp, "SCALARS %s float\n", label);
	fprintf(fp, "LOOKUP_TABLE default\n");
	
	for(int i=0; i<(fe->total_num_nodes); i++) {	
		fprintf(fp, "%e\n", val[i]);
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

	write_vtk_shape(fe, fp);
	
	fprintf(fp, "POINT_DATA %d\n", fe->total_num_nodes);
	write_nodal_value_scalar(fe, fp, vals->T, "temperature");

	fclose(fp);

}


/**********************************************************
 * numerical integration
 **********************************************************/
void integ_point_tet_5(
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

	integ_weight[0] = -2.0/15.0;
	integ_weight[1] =  3.0/40.0;
	integ_weight[2] =  3.0/40.0;
	integ_weight[3] =  3.0/40.0;
	integ_weight[4] =  3.0/40.0;
}


double calc_integ(
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
void shapefunc_3d_tet_1st_value(
		const double    xi[3],
		double*         N)
{
	N[0] = 1.0 - xi[0] - xi[1] - xi[2];
	N[1] = xi[0];
	N[2] = xi[1];
	N[3] = xi[2];
}


void shapefunc_3d_tet_1st_der_value(
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


void calc_3d_Jacobi_matrix(
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


/**********************************************************
 * boundary condition
 **********************************************************/
void set_Dirichlet_bc_CSR_mat(
		Dataset_CSR*   csr,
		const bool*    node_is_Dirichlet_bc, //[num_nodes*num_dof_on_node]
		const double*  imposed_val)          //[num_nodes*num_dof_on_node]
{
	int n = csr->num_nodes;
	int s = csr->num_dofs_on_node;

	for(int i=0; i<n; i++) {
		for(int k=0; k<s; k++) {
			if( node_is_Dirichlet_bc[ s*i+k ] ) {

				for(int j=0; j<n; j++) {
					for(int l=0; l<s; l++) {
						MatrixCSR_set_nonzero_value(csr, 0.0, i, j, k, l);
						MatrixCSR_set_nonzero_value(csr, 0.0, j, i, l, k);
					}
				}

				MatrixCSR_set_nonzero_value(csr, 1.0, i, i, k, k);
			}
		}
	}
}


void set_Dirichlet_bc_CSR_vec(
		Dataset_CSR* csr,
		const bool* node_is_Dirichlet_bc,
		double* imposed_val,
		double*  g_rhs)
{
	int n = csr->num_nodes;
	int s = csr->num_dofs_on_node;

	for(int i=0; i<n; i++) {
		for(int k=0; k<s; k++) {
			if( !node_is_Dirichlet_bc[ s*i+k ] ) {
				imposed_val[ s*i+k ] = 0.0;
			}
		}
	}

	double* bc_vec;
	bc_vec = (double*)calloc(n*s, sizeof(double));

	MatrixCSR_mat_vec_multiplication(
			csr, imposed_val, bc_vec);

	for(int i=0; i<n; i++) {
		for(int k=0; k<s; k++) {
			g_rhs[ s*i+k ] -= bc_vec[ s*i+k ];

			if( node_is_Dirichlet_bc[ s*i+k ] ) {
				g_rhs[ s*i+k ] = imposed_val[ s*i+k ];
			}
		}
	}

	free(bc_vec);
}



/**********************************************************
 * element matrix
 **********************************************************/
void set_Jacobi_matrix(
		FE_DATA*     fe,
		FE_3D_BASIS* basis)
{
	double** local_x;
	local_x = (double**)calloc(fe->local_num_nodes, sizeof(double*));
	for(int i=0; i<(fe->local_num_nodes); i++) {
		local_x[i] = (double*)calloc(3, sizeof(double));
	}

	for(int e=0; e<(fe->total_num_elems); e++) {

		for(int i=0; i<(fe->local_num_nodes); i++) {
			local_x[i][0] = fe->x[ fe->conn[e][i] ][0];
			local_x[i][1] = fe->x[ fe->conn[e][i] ][1];
			local_x[i][2] = fe->x[ fe->conn[e][i] ][2];
		}

		for(int p=0; p<(basis->num_integ_points); p++) {
			calc_3d_Jacobi_matrix(
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

	for(int i=0; i<(fe->local_num_nodes); i++) {
		free(local_x[i]);
	}
	free(local_x);
}


void get_Jacobian_array(
		double*    Jacobian_ip,
		const int  num_integ_points,
		const int  elem_num,
		FE_DATA*   fe)
{
	for(int p=0; p<num_integ_points; p++) {
		Jacobian_ip[p] = fe->geo[ elem_num ][p].Jacobian;
	}
}


void set_shapefunc_derivative(
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


void set_element_matrix(
		FE_DATA*     fe,
		FE_3D_BASIS* basis,
		Dataset_CSR* csr)
{
	double* val_ip;
	double* Jacobian_ip;
	val_ip      = (double*)calloc(basis->num_integ_points, sizeof(double));
	Jacobian_ip = (double*)calloc(basis->num_integ_points, sizeof(double));

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
						fe->geo[e][p].grad_N[i][0] * fe->geo[e][p].grad_N[j][0] +
						fe->geo[e][p].grad_N[i][1] * fe->geo[e][p].grad_N[j][1] +
						fe->geo[e][p].grad_N[i][2] * fe->geo[e][p].grad_N[j][2];
				}

				double integ_val = calc_integ(
						basis->num_integ_points,
						val_ip,
						basis->integ_weight,
						Jacobian_ip);

				MatrixCSR_add_nonzero_value(
						csr,
						integ_val,
						fe->conn[e][i], fe->conn[e][j], 0, 0);

			}
		}
	}

	free(val_ip);
	free(Jacobian_ip);
}



/**********************************************************
 * main function
 **********************************************************/
int main (
		int argc, 
		char* argv[])
{

	printf("\n");

	FE_SYSTEM sys;

	memory_allocation_basis(
			&(sys.basis),
			NUM_INTEG_POINTS,
			POL_ORDER,
			NUM_NODES_IN_ELEM);

	read_and_memory_allocation_FE_data(
			&(sys.fe), 
			"./util/meshgen/mesh.txt",
			sys.basis.num_integ_points);

	memory_allocation_nodal_values(
			&(sys.vals),
			sys.fe.total_num_nodes);

	initialize_basis(&(sys.basis));

	MatrixCSR_initialize(&sys.csr, 
			sys.fe.total_num_nodes,
			sys.fe.total_num_elems,
			sys.fe.local_num_nodes,
			1,
			sys.fe.conn);

	set_Jacobi_matrix(
			&(sys.fe), 
			&(sys.basis));
	set_shapefunc_derivative(
			&(sys.fe), 
			&(sys.basis));
	set_element_matrix(
			&(sys.fe), 
			&(sys.basis), 
			&(sys.csr));

	/* test code for B.C. (tentative) */
	bool* node_is_Dirichlet_bc;  double* imposed_val;
	node_is_Dirichlet_bc = (bool*)calloc(sys.fe.total_num_nodes, sizeof(bool));
	imposed_val = (double*)calloc(sys.fe.total_num_nodes, sizeof(double));
	for(int i=0; i<sys.fe.total_num_nodes; i++) {
		node_is_Dirichlet_bc[i] = false;
		imposed_val[i] = 0.0;
	}
	int last = sys.fe.total_num_nodes - 1;
	node_is_Dirichlet_bc[0]    = true;  imposed_val[0]    = 0.0;
	node_is_Dirichlet_bc[last] = true;  imposed_val[last] = 1.0;

	set_Dirichlet_bc_CSR_vec(&(sys.csr), node_is_Dirichlet_bc, imposed_val, sys.csr.rhs);
	set_Dirichlet_bc_CSR_mat(&(sys.csr), node_is_Dirichlet_bc, imposed_val);
	free(node_is_Dirichlet_bc);  free(imposed_val);
	/************************************/

	MatrixCSR_solver_CG(
			&(sys.csr), 
			sys.csr.rhs, 
			sys.vals.T, 
			MAT_EPSILON, 
			MAT_MAX_ITER, 
			true);
	
	output_result_file_vtk(
			&(sys.fe),
			&(sys.vals),
			OUTPUT_FILENAME);

	printf("\n");

	return 0;
}
