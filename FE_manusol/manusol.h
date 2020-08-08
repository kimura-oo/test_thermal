#pragma once

#include "../FE_sys/FE_dataset.h"

void BBFE_manusol_calc_nodal_error_scalar(
		FE_DATA*      fe,
		double        t,
		double*       error,
		const double* val);

void BBFE_manusol_overwrite_bc_file(
		FE_DATA*   fe,
		const int  block_size);

double BBFE_manusol_get_sol_scalar_3d(
		double x,
		double y,
		double z, 
		double t);

double BBFE_manusol_get_rhs_scalar_3d(
		double x,
		double y,
		double z,
		double t);

void BBFE_manusol_set_bc_scalar(
		FE_DATA* fe,
		BC_DATA* bc,
		double   t);


void BBFE_manusol_add_rhs_scalar(
		FE_DATA*     fe,
		FE_3D_BASIS* basis,
		double       t,
		double*      rhs);
