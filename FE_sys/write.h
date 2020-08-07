#pragma once

#include <stdio.h>
#include "FE_dataset.h"


void BBFE_sys_write_vtk_shape(
		FILE*     fp,
		FE_DATA*  fe,
		const int cell_type);

void BBFE_write_ascii_nodal_vals_scalar(
		FE_DATA*     fe,
		double*      vals,
		const char*  filename);
