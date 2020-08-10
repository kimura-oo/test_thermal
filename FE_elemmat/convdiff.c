
#include "convdiff.h"



double BBFE_elemmat_thermal_steady_linear(
		double grad_N_i[3],
		double grad_N_j[3])
{
	double val = -(  
			grad_N_i[0] * grad_N_j[0] +
			grad_N_i[1] * grad_N_j[1] +
			grad_N_i[2] * grad_N_j[2]);

	return val;
}
