#pragma once

#include <stdio.h>
#include <stdbool.h>

#include "monolis.h"

#include "libBB/std.h"
#include "libBB/calc.h"
#include "libBB/vtk.h"

#include "FE_std/integ.h"
#include "FE_std/shapefunc.h"

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

	int         total_num_elems;
	int         local_num_nodes;
	int**       conn;
	FE_3D_GEO** geo;  //[total_num_elems][num_integ_points]

} FE_DATA;


typedef struct
{
	double* T;
	double* error;

} NODAL_VALUES;


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
} BC_DATA;


/**********************************************************
 * memory allocation
 **********************************************************/

void BBFE_std_memory_allocation_basis(
		FE_3D_BASIS*    basis,
		const int       num_integ_points,
		const int       pol_order,
		const int       num_nodes_in_elem);

void BBFE_std_memory_allocation_integ(
		FE_3D_BASIS*    basis,
		const int       num_integ_points);

void BBFE_std_memory_allocation_shapefunc(
		FE_3D_BASIS*    basis,
		const int       num_nodes_in_elem,
		const int       pol_order,
		const int       num_integ_points);

void BBFE_std_memory_allocation_node(
		FE_DATA*  fe);

void BBFE_std_memory_allocation_elem(
		FE_DATA*  fe,
		int       num_integ_points);

/**********************************************************
 * read
 **********************************************************/
void BBFE_sys_read_node(
		FE_DATA*     fe,
		const char*  filename);

void BBFE_sys_read_elem(
		FE_DATA*     fe,
		const char*  filename,
		int          num_integ_points);

void BBFE_sys_read_Dirichlet_bc(
		BC_DATA*     bc,
		const char*  filename,
		const int    total_num_nodes);

/**********************************************************
 * write
 **********************************************************/

void BBFE_sys_write_vtk_shape(
		FILE*     fp,
		FE_DATA*  fe,
		const int cell_type);

void BBFE_write_ascii_nodal_vals_scalar(
		FE_DATA*     fe,
		double*      vals,
		const char*  filename);


/**********************************************************
 * shape function
 **********************************************************/
void BBFE_std_shapefunc_get_val_tet_1st(
		const double    xi[3],
		double*         N);

void BBFE_std_shapefunc_get_derivative_tet_1st(
		const double    xi[3],
		double*         dN_dxi,
		double*         dN_det,
		double*         dN_dze);

void BBFE_std_shapefunc_get_surface_tet_1st(
		int        surf_conn[3],
		const int  surf_num);

/**********************************************************
 * mapping
 **********************************************************/

void BBFE_std_mapping_calc_Jacobi_mat_3d(
		double     J[3][3],
		const int  local_num_nodes,
		double**   local_x,
		double*    local_dN_dxi,
		double*    local_dN_det,
		double*    local_dN_dze);


void BBFE_std_mapping_integ_point(
		double        x[3],
		const int     local_num_nodes,
		double**      local_x,
		double*       N);


/**********************************************************
 * surface
 **********************************************************/

int BBFE_std_surface_get_surface_node(
		bool*    node_is_on_surface,
		FE_DATA* fe);


/**********************************************************
 * equivval
 **********************************************************/

void BBFE_std_equivval_volume_smooth_function(
		double* equiv_val,
		FE_DATA* fe,
		FE_3D_BASIS* basis,
		double (*func)(double, double, double)); // scalar function(x, y, z)


/**********************************************************
 * monolis wrapper
 **********************************************************/

void BBFE_sys_monowrap_set_Dirichlet_bc(
		MONOLIS*      monolis,
		int           num_nodes,
		int           num_dofs_on_node,
		BC_DATA*      bc,
		double*       g_rhs);

/**********************************************************
 * element matrix
 **********************************************************/
void BBFE_elemmat_set_Jacobi_mat(
		FE_DATA*     fe,
		FE_3D_BASIS* basis);

void BBFE_elemmat_set_shapefunc_derivative(
		FE_DATA*      fe,
		FE_3D_BASIS*  basis);

void BBFE_elemmat_set_element_mat(
		MONOLIS*     monolis,
		FE_DATA*     fe,
		FE_3D_BASIS* basis);

double BBFE_elemmat_thermal_steady_linear(
		double grad_N_i[3],
		double grad_N_j[3]);

/**********************************************************
 * manufactured solution
 **********************************************************/
void BBFE_std_manusol_calc_nodal_error_scalar(
		FE_DATA*      fe,
		double*       error,
		const double* val);

void BBFE_std_manusol_overwrite_bc_file(
		FE_DATA*   fe,
		const int  block_size);

double BBFE_std_manusol_get_sol_scalar_3d(
		double x,
		double y,
		double z);

double BBFE_std_manusol_get_rhs_scalar_3d(
		double x,
		double y,
		double z);

void BBFE_std_manusol_set_bc_scalar(
		FE_DATA* fe,
		BC_DATA* bc);


void BBFE_std_manusol_add_rhs_scalar(
		FE_DATA* fe,
		FE_3D_BASIS* basis,
		double* rhs);
