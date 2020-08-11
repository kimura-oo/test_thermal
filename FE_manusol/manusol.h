#pragma once

#include "../FE_sys/FE_dataset.h"

void BBFE_manusol_calc_nodal_error_scalar(
		FE_DATA*      fe,
		double*       error,
		double*       theo_sol,
		const double* val);

void BBFE_manusol_overwrite_bc_file_hex(
		FE_DATA*    fe,
		const int   block_size,
		const char* filename,
		const char* directory);

void BBFE_manusol_overwrite_bc_file_tet(
		FE_DATA*    fe,
		const int   block_size,
		const char* filename,
		const char* directory);

void BBFE_manusol_set_bc_scalar(
		FE_DATA* fe,
		BC_DATA* bc,
		double*  theo_sol,
		double   t);

