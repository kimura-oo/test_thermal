#pragma once

#include "BBFE/sys/FE_dataset.h"
#include "monolis.h"

void BBFE_elemmat_equivval_volume_smooth_function(
		double*      equiv_val,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		double       t,
		double       (*func)(double, double, double, double)); // scalar function(x, y, z, t)

double BBFE_elemmat_equivval_relative_L2_error_scalar(
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
		MONOLIS_COM*  monolis_com,
		double        t,
		const double* comp_vec, // [total_num_nodes]
		double        (*func)(double, double, double, double)); // scalar function(x, y, z, t)
