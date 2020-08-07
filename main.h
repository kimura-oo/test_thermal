#pragma once

#include <stdio.h>
#include <stdbool.h>

#include "monolis.h"

#include "libBB/std.h"
#include "libBB/calc.h"
#include "libBB/vtk.h"

#include "FE_std/integ.h"
#include "FE_std/shapefunc.h"
#include "FE_std/mapping.h"
#include "FE_std/surface.h"

#include "FE_sys/FE_dataset.h"
#include "FE_sys/memory.h"
#include "FE_sys/read.h"



typedef struct
{
	double* T;
	double* error;

} NODAL_VALUES;

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
 * equivval
 **********************************************************/

void BBFE_elemmat_equivval_volume_smooth_function(
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
void BBFE_sys_manusol_calc_nodal_error_scalar(
		FE_DATA*      fe,
		double*       error,
		const double* val);

void BBFE_sys_manusol_overwrite_bc_file(
		FE_DATA*   fe,
		const int  block_size);

double BBFE_sys_manusol_get_sol_scalar_3d(
		double x,
		double y,
		double z);

double BBFE_sys_manusol_get_rhs_scalar_3d(
		double x,
		double y,
		double z);

void BBFE_sys_manusol_set_bc_scalar(
		FE_DATA* fe,
		BC_DATA* bc);


void BBFE_sys_manusol_add_rhs_scalar(
		FE_DATA* fe,
		FE_3D_BASIS* basis,
		double* rhs);
