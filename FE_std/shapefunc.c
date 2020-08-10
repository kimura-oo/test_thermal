

/**********************************************************
 * 3D hexahedron, 1st order
 **********************************************************/
void BBFE_std_shapefunc_hex1st_get_val(
		const double    xi[3],
		double*         N)
{
	double coef = 1.0/8.0;
	N[0] = coef * (1.0-xi[0]) * (1.0-xi[1]) * (1.0-xi[2]);
	N[1] = coef * (1.0+xi[0]) * (1.0-xi[1]) * (1.0-xi[2]);
	N[2] = coef * (1.0+xi[0]) * (1.0+xi[1]) * (1.0-xi[2]);
	N[3] = coef * (1.0-xi[0]) * (1.0+xi[1]) * (1.0-xi[2]);
	N[4] = coef * (1.0-xi[0]) * (1.0-xi[1]) * (1.0+xi[2]);
	N[5] = coef * (1.0+xi[0]) * (1.0-xi[1]) * (1.0+xi[2]);
	N[6] = coef * (1.0+xi[0]) * (1.0+xi[1]) * (1.0+xi[2]);
	N[7] = coef * (1.0-xi[0]) * (1.0+xi[1]) * (1.0+xi[2]);
}


void BBFE_std_shapefunc_hex1st_get_derivative(
		const double    xi[3],
		double*         dN_dxi,
		double*         dN_det,
		double*         dN_dze)
{	
	double coef = 1.0/8.0;
	dN_dxi[0] = -coef *               (1.0-xi[1]) * (1.0-xi[2]);
	dN_dxi[1] =  coef *               (1.0-xi[1]) * (1.0-xi[2]);
	dN_dxi[2] =  coef *               (1.0+xi[1]) * (1.0-xi[2]);
	dN_dxi[3] = -coef *               (1.0+xi[1]) * (1.0-xi[2]);
	dN_dxi[4] = -coef *               (1.0-xi[1]) * (1.0+xi[2]);
	dN_dxi[5] =  coef *               (1.0-xi[1]) * (1.0+xi[2]);
	dN_dxi[6] =  coef *               (1.0+xi[1]) * (1.0+xi[2]);
	dN_dxi[7] = -coef *               (1.0+xi[1]) * (1.0+xi[2]);

	dN_det[0] = -coef * (1.0-xi[0])               * (1.0-xi[2]);
	dN_det[1] = -coef * (1.0+xi[0])               * (1.0-xi[2]);
	dN_det[2] =  coef * (1.0+xi[0])               * (1.0-xi[2]);
	dN_det[3] =  coef * (1.0-xi[0])               * (1.0-xi[2]);
	dN_det[4] = -coef * (1.0-xi[0])               * (1.0+xi[2]);
	dN_det[5] = -coef * (1.0+xi[0])               * (1.0+xi[2]);
	dN_det[6] =  coef * (1.0+xi[0])               * (1.0+xi[2]);
	dN_det[7] =  coef * (1.0-xi[0])               * (1.0+xi[2]);

	dN_dze[0] = -coef * (1.0-xi[0]) * (1.0-xi[1])              ;
	dN_dze[1] = -coef * (1.0+xi[0]) * (1.0-xi[1])              ;
	dN_dze[2] = -coef * (1.0+xi[0]) * (1.0+xi[1])              ;
	dN_dze[3] = -coef * (1.0-xi[0]) * (1.0+xi[1])              ;
	dN_dze[4] =  coef * (1.0-xi[0]) * (1.0-xi[1])              ;
	dN_dze[5] =  coef * (1.0+xi[0]) * (1.0-xi[1])              ;
	dN_dze[6] =  coef * (1.0+xi[0]) * (1.0+xi[1])              ;
	dN_dze[7] =  coef * (1.0-xi[0]) * (1.0+xi[1])              ;
	
}


void BBFE_std_shapefunc_hex1st_get_surface(
		int        surf_conn[4],
		const int  surf_num)
{
	switch(surf_num) {
		case 0:
			surf_conn[0] = 0;  surf_conn[1] = 1;  surf_conn[2] = 2;  surf_conn[3] = 3;
			break;
		case 1:
			surf_conn[0] = 4;  surf_conn[1] = 7;  surf_conn[2] = 6;  surf_conn[3] = 5;
			break;
		case 2:
			surf_conn[0] = 0;  surf_conn[1] = 4;  surf_conn[2] = 5;  surf_conn[3] = 1;
			break;
		case 3:
			surf_conn[0] = 3;  surf_conn[1] = 2;  surf_conn[2] = 6;  surf_conn[3] = 7;
			break;
		case 4:
			surf_conn[0] = 0;  surf_conn[1] = 3;  surf_conn[2] = 7;  surf_conn[3] = 4;
			break;
		case 5:
			surf_conn[0] = 1;  surf_conn[1] = 5;  surf_conn[2] = 6;  surf_conn[3] = 2;
			break;
	}
}


/**********************************************************
 * 3D tetrahedron, 1st order
 **********************************************************/
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

