
#include "convdiff_core.h"


void BBFE_convdiff_pre_surface(
		BBFE_DATA*    surf,
		BBFE_BASIS*   basis,
		const char*   directory,
		int           num_integ_points_each_axis);

void BBFE_convdiff_set_basis_surface(
		BBFE_BASIS*   basis,
		int           local_num_nodes,
		int           num_integ_points_each_axis);

void BBFE_convdiff_set_equiv_val_Neumann(
		double*     equiv_val,
		BBFE_DATA*  surf,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		double      t,
		double      delta,
		double      (*func)(double, double, double, double)); // scalar function(x, y, z, t);
