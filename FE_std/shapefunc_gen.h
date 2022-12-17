#pragma once


/**********************************************************
 * 1D, arbitrary order
 **********************************************************/
void BBFE_std_shapefunc_gen_1d_get_val(
		const int 		p_order, 
		const double 	xi,
		double* 	 	N);

void BBFE_std_shapefunc_gen_1d_get_derivative(
		const int 		p_order, 
		const double 	xi,
		double* 	 	dN_dxi);


/**********************************************************
 * 2D, arbitrary order
 **********************************************************/
void BBFE_std_shapefunc_gen_2d_get_val(
		const int 		p_order, 
		const double 	xi[2],
		double* 	 	N);


void BBFE_std_shapefunc_gen_2d_get_derivative(
		const int 		p_order, 
		const double 	xi[2],
		double* 	 	dN_dxi,
		double* 	 	dN_det);


/**********************************************************
 * 3D, arbitrary order
 **********************************************************/
void BBFE_std_shapefunc_gen_3d_get_val(
		const int 		p_order, 
		const double 	xi[3],
		double* 	 	N);


void BBFE_std_shapefunc_gen_3d_get_derivative(
		const int 		p_order, 
		const double 	xi[2],
		double* 	 	dN_dxi,
		double* 	 	dN_det,
		double* 	 	dN_dze);
