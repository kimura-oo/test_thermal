
#include "solid.h"
#include "BB/calc.h"

#include <math.h>

/**********************************************************
 * solid linear 
 **********************************************************/
void BBFE_elemmat_solid_mat_dispstr_linear(
		double       mat[6][3],
		const double grad_N[3])
{
	mat[0][0] = grad_N[0];  mat[0][1] = 0.0;        mat[0][2] = 0.0;
	mat[1][0] = 0.0;        mat[1][1] = grad_N[1];  mat[1][2] = 0.0;
	mat[2][0] = 0.0;        mat[2][1] = 0.0;        mat[2][2] = grad_N[2];
	mat[3][0] = grad_N[1];  mat[3][1] = grad_N[0];  mat[3][2] = 0.0;
	mat[4][0] = 0.0;        mat[4][1] = grad_N[2];  mat[4][2] = grad_N[1];
	mat[5][0] = grad_N[2];  mat[5][1] = 0.0;        mat[5][2] = grad_N[0];
}


void BBFE_elemmat_solid_mat_Hooke(
		double       mat[6][6],
		const double e, /* Young's mudulus */
		const double v) /* Poisson's ratio */
{
	double coef  = e/( (1.0 + v)*(1.0 - 2.0*v) );
	double val1  = coef * (1.0 - v);
	double val2  = coef * v;
	double val3  = coef * (1.0 - 2.0*v)/2.0;

	for(int i=0; i<6; i++) {
		for(int j=0; j<6; j++) {
			mat[i][j] = 0.0;
		}
	}

	mat[0][0] = val1;  mat[0][1] = val2;  mat[0][2] = val2;
	mat[1][0] = val2;  mat[1][1] = val1;  mat[1][2] = val2;
	mat[2][0] = val2;  mat[2][1] = val2;  mat[2][2] = val1;

	mat[3][3] = val3;                                 
	mat[4][4] = val3;                
	mat[5][5] = val3;

}


void BBFE_elemmat_solid_mat_linear(
		double       mat[3][3],
		const double grad_N_i[3],
		const double grad_N_j[3],
		const double e,
		const double v)
{
	double B_i[6][3];
	double B_j[6][3];
	double D[6][6];

	BBFE_elemmat_solid_mat_dispstr_linear(B_i, grad_N_i);
	BBFE_elemmat_solid_mat_dispstr_linear(B_j, grad_N_j);
	BBFE_elemmat_solid_mat_Hooke(D, e, v);

	for(int i=0; i<3; i++) {
		for(int j=0; j<3; j++) {
			mat[i][j] = 0.0;

			for(int k=0; k<6; k++) {
				double c = 0.0;
				for(int l=0; l<6; l++) {
					c += B_i[l][i] * D[l][k];
				}
				mat[i][j] += c * B_j[k][j];
			}
		}
	}

}


