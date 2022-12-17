
#include "BB/std.h"

/**********************************************************
 * 1D, arbitrary order
 **********************************************************/
static double calc_ith_Lagrange_polynomial(
		const int 		p_order, // polynomial order
		const int 		i,       // local node number in element
		const double 	xi)      // position in natural coordinates
{
	double N_i = 1.0;
	double denom = 1.0;

	int num_nodes = p_order+1;
	double dxi = 2.0/(double)(p_order);
	double xi_i = dxi * i - 1.0;

	for(int j=0; j<num_nodes; j++) {
		double xi_j = dxi * j - 1.0;

		if(i != j) {
			N_i *= (xi - xi_j);
			denom *= (xi_i - xi_j); 
		}
	}

	N_i /= denom;

	return N_i;
}


static double calc_ith_Lagrange_polynomial_derivative(
		const int 		p_order, // polynomial order
		const int 		i,       // local node number in element
		const double 	xi)      // position in natural coordinates
{
	double dNi_dxi = 0.0;
	double denom = 1.0;

	int num_nodes = p_order+1;
	double dxi = 2.0/(double)(p_order);
	double xi_i = dxi * i - 1.0;

	for(int k=0; k<num_nodes; k++) {

		double term = 1.0;

		if(i != k) {
			for(int j=0; j<num_nodes; j++) {
				double xi_j = dxi * j - 1.0;

				if(i != j && j != k) {
					term *= (xi - xi_j);
				}
			}

			dNi_dxi += term;
		}
	}

	for(int j=0; j<num_nodes; j++) {
		double xi_j = dxi * j - 1.0;
		if(i != j) {
			denom *= (xi_i - xi_j); 
		}
	}

	dNi_dxi /= denom;

	return dNi_dxi;
}


void BBFE_std_shapefunc_gen_1d_get_val(
		const int 		p_order, 
		const double 	xi,
		double* 	 	N)
{
	int n = p_order + 1;

	for(int i=0; i<n; i++) {
		N[i] = calc_ith_Lagrange_polynomial(p_order, i, xi);
	}
}


void BBFE_std_shapefunc_gen_1d_get_derivative(
		const int 		p_order, 
		const double 	xi,
		double* 	 	dN_dxi)
{
		int n = p_order + 1;

	for(int i=0; i<n; i++) {
		dN_dxi[i] = calc_ith_Lagrange_polynomial_derivative(p_order, i, xi);
	}
}


/**********************************************************
 * 2D, arbitrary order
 **********************************************************/
static int set_index_2d(
		const int 	ix,
		const int 	iy,
		const int 	nx,
		const int 	ny)
{
	return( nx*iy + ix );
}


void BBFE_std_shapefunc_gen_2d_get_val(
		const int 		p_order, 
		const double 	xi[2],
		double* 	 	N)
{
	int n = p_order + 1;

	double* L_xi;
	L_xi = BB_std_calloc_1d_double(L_xi, n);
	double* L_et;
	L_et = BB_std_calloc_1d_double(L_et, n);

	BBFE_std_shapefunc_gen_1d_get_val(p_order, xi[0], L_xi);
	BBFE_std_shapefunc_gen_1d_get_val(p_order, xi[1], L_et);

	for(int i=0; i<n; i++) {
		for(int j=0; j<n; j++) {
			N[ set_index_2d(i, j, n, n) ] = L_xi[i]*L_et[j];
		}
	}

	BB_std_free_1d_double(L_xi, n);
	BB_std_free_1d_double(L_et, n);
}


void BBFE_std_shapefunc_gen_2d_get_derivative(
		const int 		p_order, 
		const double 	xi[2],
		double* 	 	dN_dxi,
		double* 	 	dN_det)
{
	int n = p_order + 1;

	double* L_xi;  double* dL_xi;
	 L_xi = BB_std_calloc_1d_double( L_xi, n);
	dL_xi = BB_std_calloc_1d_double(dL_xi, n);
	double* L_et;  double* dL_et;
	 L_et = BB_std_calloc_1d_double( L_et, n);
	dL_et = BB_std_calloc_1d_double(dL_et, n);

	BBFE_std_shapefunc_gen_1d_get_val(       p_order, xi[0],  L_xi);
	BBFE_std_shapefunc_gen_1d_get_derivative(p_order, xi[0], dL_xi);
	BBFE_std_shapefunc_gen_1d_get_val(       p_order, xi[1],  L_et);
	BBFE_std_shapefunc_gen_1d_get_derivative(p_order, xi[1], dL_et);

	for(int i=0; i<n; i++) {
		for(int j=0; j<n; j++) {
			dN_dxi[ set_index_2d(i, j, n, n) ] = dL_xi[i]* L_et[j];
			dN_det[ set_index_2d(i, j, n, n) ] =  L_xi[i]*dL_et[j];
		}
	}

	BB_std_free_1d_double( L_xi, n);
	BB_std_free_1d_double(dL_xi, n);
	BB_std_free_1d_double( L_et, n);
	BB_std_free_1d_double(dL_et, n);

}


/**********************************************************
 * 3D, arbitrary order
 **********************************************************/
static int set_index_3d(
		const int 	ix,
		const int 	iy,
		const int 	iz,
		const int 	nx,
		const int 	ny,
		const int 	nz)
{
	return( nx*ny*iz + nx*iy + ix );
}


void BBFE_std_shapefunc_gen_3d_get_val(
		const int 		p_order, 
		const double 	xi[3],
		double* 	 	N)
{
	int n = p_order + 1;

	double* L_xi;
	L_xi = BB_std_calloc_1d_double(L_xi, n);
	double* L_et;
	L_et = BB_std_calloc_1d_double(L_et, n);
	double* L_ze;
	L_ze = BB_std_calloc_1d_double(L_ze, n);

	BBFE_std_shapefunc_gen_1d_get_val(p_order, xi[0], L_xi);
	BBFE_std_shapefunc_gen_1d_get_val(p_order, xi[1], L_et);
	BBFE_std_shapefunc_gen_1d_get_val(p_order, xi[2], L_ze);

	for(int i=0; i<n; i++) {
		for(int j=0; j<n; j++) {
			for(int k=0; k<n; k++) {
				N[ set_index_3d(i, j, k, n, n, n) ] = L_xi[i]*L_et[j]*L_ze[k];
			}
		}
	}

	BB_std_free_1d_double(L_xi, n);
	BB_std_free_1d_double(L_et, n);
	BB_std_free_1d_double(L_ze, n);
}


void BBFE_std_shapefunc_gen_3d_get_derivative(
		const int 		p_order, 
		const double 	xi[3],
		double* 	 	dN_dxi,
		double* 	 	dN_det,
		double* 	 	dN_dze)
{
	int n = p_order + 1;

	double* L_xi;  double* dL_xi;
	 L_xi = BB_std_calloc_1d_double( L_xi, n);
	dL_xi = BB_std_calloc_1d_double(dL_xi, n);
	double* L_et;  double* dL_et;
	 L_et = BB_std_calloc_1d_double( L_et, n);
	dL_et = BB_std_calloc_1d_double(dL_et, n);
	double* L_ze;  double* dL_ze;
	 L_ze = BB_std_calloc_1d_double( L_ze, n);
	dL_ze = BB_std_calloc_1d_double(dL_ze, n);

	BBFE_std_shapefunc_gen_1d_get_val(       p_order, xi[0],  L_xi);
	BBFE_std_shapefunc_gen_1d_get_derivative(p_order, xi[0], dL_xi);
	BBFE_std_shapefunc_gen_1d_get_val(       p_order, xi[1],  L_et);
	BBFE_std_shapefunc_gen_1d_get_derivative(p_order, xi[1], dL_et);
	BBFE_std_shapefunc_gen_1d_get_val(       p_order, xi[2],  L_ze);
	BBFE_std_shapefunc_gen_1d_get_derivative(p_order, xi[2], dL_ze);

	for(int i=0; i<n; i++) {
		for(int j=0; j<n; j++) {
			for(int k=0; k<n; k++) {
				dN_dxi[ set_index_3d(i, j, k, n, n, n) ] = dL_xi[i]* L_et[j]* L_ze[k];
				dN_det[ set_index_3d(i, j, k, n, n, n) ] =  L_xi[i]*dL_et[j]* L_ze[k];
				dN_dze[ set_index_3d(i, j, k, n, n, n) ] =  L_xi[i]* L_et[j]*dL_ze[k];
			}
		}
	}

	BB_std_free_1d_double( L_xi, n);
	BB_std_free_1d_double(dL_xi, n);
	BB_std_free_1d_double( L_et, n);
	BB_std_free_1d_double(dL_et, n);
	BB_std_free_1d_double( L_ze, n);
	BB_std_free_1d_double(dL_ze, n);
}
