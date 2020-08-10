#pragma once


double BBFE_std_integ_calc(
		const int      num_integ_points,
		const double*  value,
		const double*  weight,
		const double*  Jacobian);

/**********************************************************
 * 1D line (generalized)
 **********************************************************/
// set the integ points and weights for Gauss-Legendre quadrature
int BBFE_std_integ_line_set_arbitrary_points(
		int     num_points,  // the number of integration points
		double* integ_point,       // array of integration points
		double* integ_weight );    // array of integration weight

/**********************************************************
 * 3D hexahedron
 **********************************************************/
int BBFE_std_integ_hex_set_arbitrary_points(
		int      num_points_in_each_axis,  // the number of integration points
		double** integ_point,       
		double*  integ_weight );    

/**********************************************************
 * 3D tetrahedron
 **********************************************************/
// rule: BBFE_std_integ_tet_
int BBFE_std_integ_tet_set_arbitrary_points(
		int      num_points_in_each_axis,  // the number of integration points
		double** integ_point,       
		double*  integ_weight );    

int BBFE_std_integ_tet_set_1points(
		double**    integ_point,
		double*     integ_weight);

int BBFE_std_integ_tet_set_4points(
		double**    integ_point,
		double*     integ_weight);

int BBFE_std_integ_tet_set_5points(
		double**    integ_point,
		double*     integ_weight);

