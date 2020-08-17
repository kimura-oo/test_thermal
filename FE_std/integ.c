
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "BB/std.h"

static const char* CODENAME = "FE_std/integ >";
static const double EPSILON_INTEG = 1.0e-12;


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


double BBFE_std_integ_calc_volume(
		const int      num_integ_points,
		const double*  weight,
		const double*  Jacobian)
{
	double val = 0.0;

	for(int i=0; i<num_integ_points; i++) {
		val += weight[i] * Jacobian[i];
	}

	return val;
}


/**********************************************************
 * 1D line (generalized)
 **********************************************************/
static double calc_Legendre_polynominal(
		const int n,
		const double x)
{
	if(n < 0) {
		printf("%s ERROR: %d-th order Legendre polynominal cannot be calculated\n", CODENAME, n);
		exit(EXIT_FAILURE);
	}
	else if(n == 0) {
		return 1.0;
	}
	else if(n == 1) {
		return x;
	}
	else {
		int i;
		double* p;
		p = (double*)calloc(n+1, sizeof(double) );
		p[0] = 1.0;
		p[1] = x;

		for(i=2; i<=n; i++) {
			p[i] = (2.0*i - 1.0) * x * p[i-1] - (i-1.0) * p[i-2];
			p[i] /= (double)(i);
		}

		free(p);

		return p[n];
	}
}


static double calc_Legendre_polynominal_derivative(
		const int n,
		const double x)
{
	double P_der = 0.0;
	double P_n = calc_Legendre_polynominal(n, x);
	double P_n_1 = calc_Legendre_polynominal(n-1, x);

	P_der = (double)(n)/(x*x - 1.0) * ( x*P_n - P_n_1 );

	return P_der;

}


static double calc_kth_integration_point(
		const int n,
		const int k) 
{
	double x;
	double x_k = cos( M_PI * ( (double)(k) + 1.0 - 0.25 )/( (double)(n) + 0.5 ) );

	int iter = 0;
	while(1) {
		x = x_k;
		x_k -= calc_Legendre_polynominal(n, x) / calc_Legendre_polynominal_derivative(n, x);
		iter++;

		if(fabs(x_k - x) < EPSILON_INTEG) break;
	}

	return x;
}


static double calc_kth_weight(
		const int n,
		const double x_k)
{
	double P_dash_n = calc_Legendre_polynominal_derivative(n, x_k);
	double weight = 2.0/( (1.0 - x_k*x_k) * (P_dash_n*P_dash_n) );

	return weight;
}


// set the integ points and weights for Gauss-Legendre quadrature
int BBFE_std_integ_line_set_arbitrary_points(
		int     num_points,  // the number of integration points
		double* integ_point,   // array of integration points
		double* integ_weight ) // array of integration weight
{
	int return_num_points = num_points;

	for(int i=0; i<num_points; i++) {
		integ_point[i]  = calc_kth_integration_point(num_points, i);
		integ_weight[i] = calc_kth_weight(num_points, integ_point[i]);
	}	

	return return_num_points;
}


/**********************************************************
 * 2D rectangle
 **********************************************************/
static int index_rec2d(
		const int i_x,
		const int i_y,
		const int num_points)
{
	return (num_points*i_y + i_x);
}


int BBFE_std_integ_rec_set_arbitrary_points(
		int     num_points_in_each_axis,  // the number of integration points
		double** integ_point,       
		double*  integ_weight)
{
	int n = num_points_in_each_axis;
	int num_integ_points = n*n;

	double* point_1d;  double* weight_1d;
	point_1d  = BB_std_calloc_1d_double(point_1d , n);
	weight_1d = BB_std_calloc_1d_double(weight_1d, n);

	BBFE_std_integ_line_set_arbitrary_points(n, point_1d, weight_1d);

	for(int i_y=0; i_y<n; i_y++) {
		for(int i_x=0; i_x<n; i_x++) {
			int g = index_rec2d(i_x, i_y, n);
			integ_point[g][0] = point_1d[i_x];
			integ_point[g][1] = point_1d[i_y];

			integ_weight[g] = weight_1d[i_x] * weight_1d[i_y];
		}
	}

	BB_std_free_1d_double(point_1d , n);
	BB_std_free_1d_double(weight_1d, n);

	return num_integ_points;
}


/**********************************************************
 * 2D triangle
 **********************************************************/
static void set_shapefunc_rec1st(
		double       N[4],
		const double xi[2])
{
	double coef = 1.0/4.0;
	N[0] = coef * (1.0-xi[0]) * (1.0-xi[1]);
	N[1] = coef * (1.0+xi[0]) * (1.0-xi[1]);
	N[2] = coef * (1.0+xi[0]) * (1.0+xi[1]);
	N[3] = coef * (1.0-xi[0]) * (1.0+xi[1]);
}


static void set_degeneration_nodes_2d(
		double     x[4][2])
{
	// Degeneration rules:
	//   tri3 = rec3 = rec4
	x[0][0] = 0.0;  x[0][1] = 0.0; 
	x[1][0] = 1.0;  x[1][1] = 0.0; 
	x[2][0] = 0.0;  x[2][1] = 1.0; 
	x[3][0] = 0.0;  x[3][1] = 1.0; 
}


static void mapping_rec_to_tri(
		double       xi_tri[2],
		const double xi_rec[2])
{
	double x[4][2];
	set_degeneration_nodes_2d(x);

	double N[4];
	set_shapefunc_rec1st(N, xi_rec);

	for(int d=0; d<2; d++) {
		xi_tri[d] = 0.0;
	}

	for(int i=0; i<4; i++) {
		for(int d=0; d<2; d++) {
			xi_tri[d] += N[i] * x[i][d];
		}
	}	
}


int BBFE_std_integ_tri_set_arbitrary_points(
		int      num_points_in_each_axis,  // the number of integration points
		double** integ_point,       
		double*  integ_weight )
{
	int n = num_points_in_each_axis;

	int num_integ_points = BBFE_std_integ_rec_set_arbitrary_points(
			n, integ_point, integ_weight);

	double w=0.0;
	for(int g=0; g<num_integ_points; g++) {
		double tri_ip[2];
		mapping_rec_to_tri(tri_ip, integ_point[g]);

		for(int d=0; d<2; d++) {
			integ_point[g][d] = tri_ip[d];
		}

		integ_weight[g] = integ_weight[g]/(4.0*2.0);
	}

	return num_integ_points;

}


/**********************************************************
 * 3D hexahedron
 **********************************************************/
static int index_hex3d(
		const int i_x,
		const int i_y,
		const int i_z,
		const int num_points)
{
	return (num_points*num_points*i_z + num_points*i_y + i_x);
}


int BBFE_std_integ_hex_set_arbitrary_points(
		int     num_points_in_each_axis,  // the number of integration points
		double** integ_point,       
		double*  integ_weight)
{
	int n = num_points_in_each_axis;
	int num_integ_points = n*n*n;

	double* point_1d;  double* weight_1d;
	point_1d  = BB_std_calloc_1d_double(point_1d , n);
	weight_1d = BB_std_calloc_1d_double(weight_1d, n);

	BBFE_std_integ_line_set_arbitrary_points(n, point_1d, weight_1d);

	for(int i_z=0; i_z<n; i_z++) {
		for(int i_y=0; i_y<n; i_y++) {
			for(int i_x=0; i_x<n; i_x++) {
				int g = index_hex3d(i_x, i_y, i_z, n);
				integ_point[g][0] = point_1d[i_x];
				integ_point[g][1] = point_1d[i_y];
				integ_point[g][2] = point_1d[i_z];

				integ_weight[g] = weight_1d[i_x] * weight_1d[i_y] * weight_1d[i_z];
			}
		}
	}

	BB_std_free_1d_double(point_1d , n);
	BB_std_free_1d_double(weight_1d, n);

	return num_integ_points;
}


/**********************************************************
 * 3D tetrahedron
 **********************************************************/
static void set_shapefunc_hex1st(
		double       N[8],
		const double xi[3])
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


static void set_degeneration_nodes_3d(
		double     x[8][3])
{
	// Degeneration rules:
	//   tet3 = hex3 = hex4
	//   tet5 = hex5 = hex6 = hex7 = hex8
	x[0][0] = 0.0;  x[0][1] = 0.0;  x[0][2] = 0.0; 
	x[1][0] = 1.0;  x[1][1] = 0.0;  x[1][2] = 0.0; 
	x[2][0] = 0.0;  x[2][1] = 1.0;  x[2][2] = 0.0; 
	x[3][0] = 0.0;  x[3][1] = 1.0;  x[3][2] = 0.0; 
	x[4][0] = 0.0;  x[4][1] = 0.0;  x[4][2] = 1.0; 
	x[5][0] = 0.0;  x[5][1] = 0.0;  x[5][2] = 1.0; 
	x[6][0] = 0.0;  x[6][1] = 0.0;  x[6][2] = 1.0; 
	x[7][0] = 0.0;  x[7][1] = 0.0;  x[7][2] = 1.0; 
}


static void mapping_hex_to_tet(
		double       xi_tet[3],
		const double xi_hex[3])
{
	double x[8][3];
	set_degeneration_nodes_3d(x);

	double N[8];
	set_shapefunc_hex1st(N, xi_hex);

	for(int d=0; d<3; d++) {
		xi_tet[d] = 0.0;
	}

	for(int i=0; i<8; i++) {
		for(int d=0; d<3; d++) {
			xi_tet[d] += N[i] * x[i][d];
		}
	}	
}


int BBFE_std_integ_tet_set_arbitrary_points(
		int      num_points_in_each_axis,  // the number of integration points
		double** integ_point,       
		double*  integ_weight )
{
	int n = num_points_in_each_axis;

	int num_integ_points = BBFE_std_integ_hex_set_arbitrary_points(
			n, integ_point, integ_weight);

	double w=0.0;
	for(int g=0; g<num_integ_points; g++) {
		double tet_ip[3];
		mapping_hex_to_tet(tet_ip, integ_point[g]);

		for(int d=0; d<3; d++) {
			integ_point[g][d] = tet_ip[d];
		}

		integ_weight[g] = integ_weight[g]/(8.0*6.0);
	}

	return num_integ_points;

}


int BBFE_std_integ_tet_set_1points(
		double**    integ_point,
		double*     integ_weight)
{
	integ_point[0][0] = 1.0/4.0;  integ_point[0][1] = 1.0/4.0;  integ_point[0][2] = 1.0/4.0;

	integ_weight[0] = 1.0 * 1.0/6.0;

	return 1;
}


int BBFE_std_integ_tet_set_4points(
		double**    integ_point,
		double*     integ_weight)
{
	double a = (5.0 + 3.0*sqrt(5.0))/20.0;
	double b = (5.0 - sqrt(5.0))/20.0;

	integ_point[0][0] = a;  integ_point[0][1] = b;  integ_point[0][2] = b;
	integ_point[1][0] = b;  integ_point[1][1] = a;  integ_point[1][2] = b;
	integ_point[2][0] = a;  integ_point[2][1] = b;  integ_point[2][2] = a;
	integ_point[3][0] = a;  integ_point[3][1] = b;  integ_point[3][2] = b;

	integ_weight[0] = 1.0/4.0 * 1.0/6.0;
	integ_weight[1] = 1.0/4.0 * 1.0/6.0;
	integ_weight[2] = 1.0/4.0 * 1.0/6.0;
	integ_weight[3] = 1.0/4.0 * 1.0/6.0;

	return 4;
}


int BBFE_std_integ_tet_set_5points(
		double**    integ_point,
		double*     integ_weight)
{
	integ_point[0][0] = 1.0/4.0;  integ_point[0][1] = 1.0/4.0;  integ_point[0][2] = 1.0/4.0;
	integ_point[1][0] = 1.0/6.0;  integ_point[1][1] = 1.0/6.0;  integ_point[1][2] = 1.0/6.0;
	integ_point[2][0] = 1.0/6.0;  integ_point[2][1] = 1.0/6.0;  integ_point[2][2] = 1.0/2.0;
	integ_point[3][0] = 1.0/6.0;  integ_point[3][1] = 1.0/2.0;  integ_point[3][2] = 1.0/6.0;
	integ_point[4][0] = 1.0/2.0;  integ_point[4][1] = 1.0/6.0;  integ_point[4][2] = 1.0/6.0;

	integ_weight[0] = -4.0/5.0  * 1.0/6.0;
	integ_weight[1] =  9.0/20.0 * 1.0/6.0;
	integ_weight[2] =  9.0/20.0 * 1.0/6.0;
	integ_weight[3] =  9.0/20.0 * 1.0/6.0;
	integ_weight[4] =  9.0/20.0 * 1.0/6.0;

	return 5;
}
