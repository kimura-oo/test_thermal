#pragma once

#include "../FE_sys/FE_dataset.h"

void BBFE_manusol_calc_nodal_error_scalar(
		BBFE_DATA*    fe,
		double*       error,
		double*       theo_sol,
		const double* val);

void BBFE_manusol_overwrite_bc_file_hex(
		BBFE_DATA*  fe,
		const int   block_size,
		const char* filename,
		const char* directory);

void BBFE_manusol_overwrite_bc_file_tet(
		BBFE_DATA*  fe,
		const int   block_size,
		const char* filename,
		const char* directory);

void BBFE_manusol_set_bc_scalar(
		BBFE_DATA* fe,
		BBFE_BC*   bc,
		double*    theo_sol,
		double     t);

