
/**********************************************************
 * 3D tetrahedron
 **********************************************************/

void BBFE_std_integ_tet_set_5points(
		int*        num_integ_points,
		double**    integ_point,
		double*     integ_weight)
{
	(*num_integ_points) = 5;

	integ_point[0][0] = 1.0/4.0;  integ_point[0][1] = 1.0/4.0;  integ_point[0][2] = 1.0/4.0;
	integ_point[1][0] = 1.0/6.0;  integ_point[1][1] = 1.0/6.0;  integ_point[1][2] = 1.0/6.0;
	integ_point[2][0] = 1.0/6.0;  integ_point[2][1] = 1.0/6.0;  integ_point[2][2] = 1.0/2.0;
	integ_point[3][0] = 1.0/6.0;  integ_point[3][1] = 1.0/2.0;  integ_point[3][2] = 1.0/6.0;
	integ_point[4][0] = 1.0/2.0;  integ_point[4][1] = 1.0/6.0;  integ_point[4][2] = 1.0/6.0;

	integ_weight[0] = -4.0/30.0;
	integ_weight[1] =  9.0/120.0;
	integ_weight[2] =  9.0/120.0;
	integ_weight[3] =  9.0/120.0;
	integ_weight[4] =  9.0/120.0;
}


double BBFE_std_integ_calc(
		const int      num_integ_points,
		const double*  value,
		const double*  weight,
		const double*  Jacobian)
{
	double val = 0.0;

	for(int i=0; i<num_integ_points; i++) {
		val += value[i] * weight[i] * Jacobian[i];
	}

	return val;
}

