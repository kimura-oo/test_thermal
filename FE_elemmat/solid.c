
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


/**********************************************************
 * solid nonlinear 
 **********************************************************/
void BBFE_elemmat_solid_tensor_defgrad(
		double       mat[3][3],
		double**     grad_u)
{
	for(int i=0; i<3; i++) {
		for(int j=0; j<3; j++) {
			mat[i][j] = grad_u[i][j];
		}
		mat[i][i] += 1.0;
	}
}


void BBFE_elemmat_solid_mat_dispstr_tl(
		double       mat[6][3],
		double**     grad_u,
		const double grad_N[3])
{
	double F[3][3];
	BBFE_elemmat_solid_tensor_defgrad(F, grad_u);

	mat[0][0] = F[0][0]*grad_N[0];  mat[0][1] = F[1][0]*grad_N[0];  mat[0][2] = F[2][0]*grad_N[0];
	mat[1][0] = F[0][1]*grad_N[1];  mat[1][1] = F[1][1]*grad_N[1];  mat[1][2] = F[2][1]*grad_N[1];
	mat[2][0] = F[0][2]*grad_N[2];  mat[2][1] = F[1][2]*grad_N[2];  mat[2][2] = F[2][2]*grad_N[2];

	mat[3][0] = F[0][0]*grad_N[1] + F[0][1]*grad_N[0];
	mat[3][1] = F[1][0]*grad_N[1] + F[1][1]*grad_N[0];
	mat[3][2] = F[2][0]*grad_N[1] + F[2][1]*grad_N[0];

	mat[4][0] = F[0][1]*grad_N[2] + F[0][2]*grad_N[1];
	mat[4][1] = F[1][1]*grad_N[2] + F[1][2]*grad_N[1];
	mat[4][2] = F[2][1]*grad_N[2] + F[2][2]*grad_N[1];

	mat[5][0] = F[0][2]*grad_N[0] + F[0][0]*grad_N[2];
	mat[5][1] = F[1][2]*grad_N[0] + F[1][0]*grad_N[2];
	mat[5][2] = F[2][2]*grad_N[0] + F[2][0]*grad_N[2];
}


void BBFE_elemmat_solid_tensor_Green_Lagrange_mat_notation(
		double    mat[6],
		double    F[3][3])
{
	double a[3][3];

	for(int i=0; i<3; i++) {
		for(int j=0; j<3; j++) {
			a[i][j] = 0.0;
			for(int k=0; k<3; k++) {
				a[i][j] += F[k][i]*F[k][j];
			}
		}
	}

	mat[0] = 0.5*a[0][0] - 0.5;  mat[1] = 0.5*a[1][1] - 0.5;  mat[2] = 0.5*a[2][2] - 0.5;
	mat[3] = a[0][1];            mat[4] = a[1][2];            mat[5] = a[2][0];
}


void BBFE_elemmat_solid_tensor_second_Piora_Kirchhoff_mat_notation(
		double       mat[6],
		double       D[6][6],
		const double E[6])
{
	for(int i=0; i<6; i++) {
		mat[i] = 0.0;
		for(int j=0; j<6; j++) {
			mat[i] += D[i][j]*E[j];
		}
	}
}


void BBFE_elemmat_solid_mat_tl(
		double       mat[3][3],
		double       D[6][6],
		const double grad_N_i[3],
		const double grad_N_j[3],
		double**     grad_u )
{
	double B_i[6][3];  double B_j[6][3];
	BBFE_elemmat_solid_mat_dispstr_tl(B_i, grad_u, grad_N_i);
	BBFE_elemmat_solid_mat_dispstr_tl(B_j, grad_u, grad_N_j);
	
	double F[3][3];  double E[6];  double S[6];
	BBFE_elemmat_solid_tensor_defgrad(F, grad_u);
	BBFE_elemmat_solid_tensor_Green_Lagrange_mat_notation(E, F);
	BBFE_elemmat_solid_tensor_second_Piora_Kirchhoff_mat_notation(S, D, E);

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

	double val = 
		S[0] *   grad_N_i[0]*grad_N_j[0]                              +
		S[1] *   grad_N_i[1]*grad_N_j[1]                              +
		S[2] *   grad_N_i[2]*grad_N_j[2]                              +
		S[3] * ( grad_N_i[0]*grad_N_j[1] +  grad_N_i[1]*grad_N_j[0] ) +
		S[4] * ( grad_N_i[1]*grad_N_j[2] +  grad_N_i[2]*grad_N_j[1] ) +
		S[5] * ( grad_N_i[2]*grad_N_j[0] +  grad_N_i[0]*grad_N_j[2] );

	for(int i=0; i<3; i++) {
		mat[i][i] += val;
	}
}


void BBFE_elemmat_solid_vec_inner_force_tl(
		double       vec[3], 
		double       D[6][6],
		const double grad_N_i[3],
		double**     grad_u)
{
	double B_i[6][3];
	BBFE_elemmat_solid_mat_dispstr_tl(B_i, grad_u, grad_N_i);
	
	double F[3][3];  double E[6];  double S[6];
	BBFE_elemmat_solid_tensor_defgrad(F, grad_u);
	BBFE_elemmat_solid_tensor_Green_Lagrange_mat_notation(E, F);
	BBFE_elemmat_solid_tensor_second_Piora_Kirchhoff_mat_notation(S, D, E);

	for(int i=0; i<3; i++) {
		vec[i] = 0.0;
		for(int j=0; j<6; j++) {
			vec[i] += B_i[j][i]*S[j];
		}
	}
}
