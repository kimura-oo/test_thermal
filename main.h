#pragma once

#include <stdio.h>
#include "solve_mat.h"

// Information of polynomials and numerical integration in normalized space
// Instance of this structure is made for each types of FE 
typedef struct
{
	// information of numerical integration in normalized space
	int      num_integ_points;  // the number of integration points in an element
	double** integ_point;       // position of integration point
	                            // [num_integ_points][x, y, z]
	double*  integ_weight;      // weight at intergration points [num_integ_points]

	/* inforamtion of shape functions in normalized space */
	int      pol_order;  // polynomial order
	int      num_nodes;  // the number of nodes in an element
	double** N;          // array of shape function 
	                     // [num_integ_points][num_nodes_in_elem]
	double** dN_dxi;	 // array of derivative (xi) shape function
	                     // [num_integ_points][num_nodes_in_elem]
	double** dN_det;     // (eta)
	double** dN_dze;     // (zeta)
} FE_3D_BASIS;


// geometric 
typedef struct
{
	double** grad_N;     // [local_num_nodes][3 (dimension)]

	double J[3][3];
	double Jacobian;
} FE_3D_GEO;


typedef struct
{
	int      total_num_nodes;
	double** x;
	
	int        total_num_elems;
	int        local_num_nodes;
	int**      conn;
	FE_3D_GEO** geo;  //[total_num_elems][num_integ_points]

} FE_DATA;


typedef struct
{
	double* T;
	
} NODAL_VALUES;


/**********************************************************
 * memory allocation
 **********************************************************/
void memory_allocation_basis(
		FE_3D_BASIS*    basis,
		const int       num_integ_points,
		const int       pol_order,
		const int       num_nodes_in_elem);

void memory_allocation_integ(
		FE_3D_BASIS*    basis,
		const int       num_integ_points);

void memory_allocation_shapefunc(
		FE_3D_BASIS*    basis,
		const int       num_nodes_in_elem,
		const int       pol_order,
		const int       num_integ_points);

void memory_allocation_nodal_values(
		NODAL_VALUES*   vals,
		const int       total_num_nodes);

/**********************************************************
 * initializers
 **********************************************************/
void initialize_basis(
		FE_3D_BASIS*  basis);


/**********************************************************
 * input
 **********************************************************/
void read_and_memory_allocation_FE_node(
		FE_DATA*  fe,
		char*     filename);


void read_and_memory_allocation_FE_elem(
		FE_DATA*  fe,
		char*     filename,
		int       num_integ_points);


/**********************************************************
 * output
 **********************************************************/
void write_vtk_shape(
		FE_DATA*  fe,
		char*     filename);

void write_nodal_value_scalar(
		FE_DATA*     fe,
		FILE*        fp,
		double*      val,
		const char*  label);

void output_result_file_vtk(
		FE_DATA*       fe,
		NODAL_VALUES*  vals,
		const char*    filename);

/**********************************************************
 * numerical integration
 **********************************************************/
void integ_point_tet_5(
		int*        num_integ_points,
		double**    integ_point,
		double*     integ_weight);

double calc_integ(
		const int      num_integ_points,
		const double*  value,
		const double*  weight,
		const double*  Jacobian);

/**********************************************************
 * shape function
 **********************************************************/
void shapefunc_3d_tet_1st_value(
		const double    xi[3],
		double*         N);

void shapefunc_3d_tet_1st_der_value(
		const double    xi[3],
		double*         dN_dxi,
		double*         dN_det,
		double*         dN_dze);

void calc_3d_Jacobi_matrix(
		double     J[3][3],
		const int  local_num_nodes,
		double**   local_x,
		double*    local_dN_dxi,
		double*    local_dN_det,
		double*    local_dN_dze);

/**********************************************************
 * boundary condition
 **********************************************************/

void set_Dirichlet_bc_CSR_mat(
		Dataset_CSR*   csr,
		const bool*    node_is_Dirichlet_bc, //[num_nodes*num_dof_on_node]
		const double*  imposed_val);         //[num_nodes*num_dof_on_node]

void set_Dirichlet_bc_CSR_vec(
		Dataset_CSR* csr,
		const bool* node_is_Dirichlet_bc,
		double* imposed_val,
		double*  g_rhs);

/**********************************************************
 * element matrix
 **********************************************************/
void set_Jacobi_matrix(
		FE_DATA*     fe,
		FE_3D_BASIS* basis);

void get_Jacobian_array(
		double*    Jacobian_ip,
		const int  num_integ_points,
		const int  elem_num,
		FE_DATA*   fe);

void set_shapefunc_derivative(
		FE_DATA*      fe,
		FE_3D_BASIS*  basis);

void set_element_matrix(
		FE_DATA*     fe,
		FE_3D_BASIS* basis,
		Dataset_CSR* csr);

