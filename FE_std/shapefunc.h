#pragma once


/**********************************************************
 * 2D rectangle, 1st order
 **********************************************************/
void BBFE_std_shapefunc_rec1st_get_val(
		const double    xi[3],
		double*         N);

void BBFE_std_shapefunc_rec1st_get_derivative(
		const double    xi[3],
		double*         dN_dxi,
		double*         dN_det);

/**********************************************************
 * 2D triangle, 1st order
 **********************************************************/
void BBFE_std_shapefunc_tri1st_get_val(
		const double    xi[2],
		double*         N);

void BBFE_std_shapefunc_tri1st_get_derivative(
		const double    xi[2],
		double*         dN_dxi,
		double*         dN_det);

/**********************************************************
 * 2D triangle, 2nd order
 **********************************************************/
void BBFE_std_shapefunc_tri2nd_get_val(
		const double    xi[2],
		double*         N);

void BBFE_std_shapefunc_tri2nd_get_derivative(
		const double    xi[2],
		double*         dN_dxi,
		double*         dN_det);
/**********************************************************
 * 3D hexahedron, 1st order
 **********************************************************/
void BBFE_std_shapefunc_hex1st_get_val(
		const double    xi[3],
		double*         N);

void BBFE_std_shapefunc_hex1st_get_derivative(
		const double    xi[3],
		double*         dN_dxi,
		double*         dN_det,
		double*         dN_dze);

void BBFE_std_shapefunc_hex1st_get_surface(
		int        surf_conn[4],
		const int  surf_num);

/**********************************************************
 * 3D tetrahedron, 1st order
 **********************************************************/
void BBFE_std_shapefunc_tet1st_get_val(
		const double    xi[3],
		double*         N);

void BBFE_std_shapefunc_tet1st_get_derivative(
		const double    xi[3],
		double*         dN_dxi,
		double*         dN_det,
		double*         dN_dze);

void BBFE_std_shapefunc_tet1st_get_surface(
		int        surf_conn[3],
		const int  surf_num);

/**********************************************************
 * 3D tetrahedron, 2nd order
 **********************************************************/
void BBFE_std_shapefunc_tet2nd_get_val(
		const double xi[3], 	
		double       N[10]);

void BBFE_std_shapefunc_tet2nd_get_derivative(
		const double xi[3], 
		double*      dN_dxi,	
		double*      dN_det,
		double*      dN_dze);

void BBFE_std_shapefunc_tet2nd_get_surface(
		int        surf_conn[6],
		const int  surf_num);
