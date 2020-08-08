#pragma once

#include "../FE_sys/FE_dataset.h"


void BBFE_elemmat_equivval_volume_smooth_function(
		double*      equiv_val,
		FE_DATA*     fe,
		FE_3D_BASIS* basis,
		double       t,
		double       (*func)(double, double, double, double)); // scalar function(x, y, z, t)

