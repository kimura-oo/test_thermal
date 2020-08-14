#pragma once

#include <stdbool.h>

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
	                     // [num_integ_points][local_num_nodes]
	double** dN_dxi;	 // array of derivative (xi) shape function
	                     // [num_integ_points][local_num_nodes]
	double** dN_det;     // (eta)
	double** dN_dze;     // (zeta)
} BBFE_BASIS;


// geometric
typedef struct
{
	double** grad_N;     // [local_num_nodes][3 (dimension)]

	double J[3][3];
	double Jacobian;
} BBFE_GEO;


typedef struct
{
	int      total_num_nodes;
	double** x;

	int         total_num_elems;
	int         local_num_nodes;
	int**       conn;
	BBFE_GEO**  geo;  //[total_num_elems][num_integ_points]

} BBFE_DATA;



typedef struct
{
	int total_num_nodes;
	int block_size;

	int num_D_bcs;
	bool*   D_bc_exists;
	double* imposed_D_val;

	int num_N_bcs;
	bool*   N_bc_exists;
	double* imposed_N_val;
} BBFE_BC;
