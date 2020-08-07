#pragma once

/**********************************************************
 * 3D tetrahedron
 **********************************************************/
// rule: BBFE_std_integ_tet_
void BBFE_std_integ_tet_set_5points(
		int*        num_integ_points,
		double**    integ_point,
		double*     integ_weight);

double BBFE_std_integ_calc(
		const int      num_integ_points,
		const double*  value,
		const double*  weight,
		const double*  Jacobian);

