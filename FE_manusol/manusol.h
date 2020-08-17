#pragma once

#include "BBFE/sys/FE_dataset.h"

void BBFE_manusol_calc_nodal_error_scalar(
		BBFE_DATA*    fe,
		double*       error,
		double*       theo_sol,
		const double* val);

void BBFE_manusol_set_bc_scalar(
		BBFE_DATA* fe,
		BBFE_BC*   bc,
		double*    theo_sol,
		double     t);

