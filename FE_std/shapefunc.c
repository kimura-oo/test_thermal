
void BBFE_std_shapefunc_tet1st_get_val(
		const double    xi[3],
		double*         N)
{
	N[0] = 1.0 - xi[0] - xi[1] - xi[2];
	N[1] = xi[0];
	N[2] = xi[1];
	N[3] = xi[2];
}


void BBFE_std_shapefunc_tet1st_get_derivative(
		const double    xi[3],
		double*         dN_dxi,
		double*         dN_det,
		double*         dN_dze)
{
	dN_dxi[0] = -1.0;  dN_det[0] = -1.0;  dN_dze[0] = -1.0;
	dN_dxi[1] =  1.0;  dN_det[1] =  0.0;  dN_dze[1] =  0.0;
	dN_dxi[2] =  0.0;  dN_det[2] =  1.0;  dN_dze[2] =  0.0;
	dN_dxi[3] =  0.0;  dN_det[3] =  0.0;  dN_dze[3] =  1.0;
}


void BBFE_std_shapefunc_tet1st_get_surface(
		int        surf_conn[3],
		const int  surf_num)
{
	switch(surf_num) {
		case 0:
			surf_conn[0] = 2;  surf_conn[1] = 1;  surf_conn[2] = 3;
			break;
		case 1:
			surf_conn[0] = 0;  surf_conn[1] = 2;  surf_conn[2] = 3;
			break;
		case 2:
			surf_conn[0] = 1;  surf_conn[1] = 0;  surf_conn[2] = 3;
			break;
		case 3:
			surf_conn[0] = 0;  surf_conn[1] = 1;  surf_conn[2] = 2;
			break;
	}
}

