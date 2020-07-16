#pragma once

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


typedef struct
{
	int      total_num_nodes;
	double** x;
	
	int   total_num_elems;
	int   local_num_nodes;
	int** conn;
} FE_3D_DATA;


typedef struct
{
	double* T;
	
} NODAL_VALUES;


typedef struct
{
	FE_3D_BASIS  basis;
	FE_3D_DATA   fe;
	NODAL_VALUES vals;

} FE_SYSTEM;


/**********************************************************
 * memory allocation
 **********************************************************/
void memory_allocation_basis(
		FE_3D_BASIS* 	basis,
		const int 		num_integ_points,
		const int 		pol_order,
		const int 		num_nodes_in_elem);

void memory_allocation_integ(
		FE_3D_BASIS* 	basis,
		const int 		num_integ_points);

void memory_allocation_shapefunc(
		FE_3D_BASIS* 	basis,
		const int 		num_nodes_in_elem,
		const int 		pol_order,
		const int 		num_integ_points);


/**********************************************************
 * initializers
 **********************************************************/
void initialize_basis(
		FE_3D_BASIS* basis);


/**********************************************************
 * numerical integration
 **********************************************************/
void integ_point_tet_5(
		int* 		num_integ_points,
		double**	integ_point,
		double*		integ_weight);


/**********************************************************
 * shape function
 **********************************************************/
void shapefunc_3d_tet_1st_value(
		const double 	xi[3],
		double* 		N);

void shapefunc_3d_tet_1st_der_value(
		const double 	xi[3],
		double* 		dN_dxi,
		double* 		dN_det,
		double* 		dN_dze);
