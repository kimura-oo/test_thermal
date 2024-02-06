
#include "fluid.h"
#include "BB/calc.h"

#include <math.h>

const double ZERO_CRITERION  = 1.0e-10;

/**********************************************************
 * Stabilized FEM related
 **********************************************************/
double BBFE_elemmat_fluid_supg_coef(
		const double density,
		const double viscosity,
		const double v[3],
		const double h_e,
		const double dt)
{

	double l_v = BB_calc_vec3d_length(v);
	if( l_v < ZERO_CRITERION) { return 0.0; }

	double nu = viscosity/density;
	double denom = (2.0/dt)*(2.0/dt) + (2.0*l_v/h_e)*(2.0*l_v/h_e) 
		+ (4.0*nu/(h_e*h_e))*(4.0*nu/(h_e*h_e));

	double val = sqrt(1.0/denom);

	return (val);
}


//sups = supg + pspg
double BBFE_elemmat_fluid_sups_coef(
		const double density,
		const double viscosity,
		const double v[3],
		const double h_e,
		const double dt)
{
	double l_v = BB_calc_vec3d_length(v);

	double nu = viscosity/density;
	double denom = (2.0/dt)*(2.0/dt) + (2.0*l_v/h_e)*(2.0*l_v/h_e) 
		+ (4.0*nu/(h_e*h_e))*(4.0*nu/(h_e*h_e));

	double val = sqrt(1.0/denom);

	return (val);
}


/**********************************************************
 * fractional step method
 **********************************************************/

double BBFE_elemmat_fluid_fs_mat_pred_expl(
		const double N_i,
		const double N_j,
		const double grad_N_i[3],
		const double v[3],
		const double tau)
{
	double val = 
		N_i * N_j + 
		tau * BB_calc_vec3d_dot(v, grad_N_i) * N_j;

	return val;
}


void BBFE_elemmat_fluid_fs_vec_pred_expl(
		double         vec[3],
		const double   N_i,
		const double   grad_N_i[3],
		const double   v[3],
		double**       grad_v,
		const double   density,
		const double   viscosity,
		const double   tau,
		const double   dt)
{
	double dyn_vis = viscosity/density;

	for(int d=0; d<3; d++) {
		double val = 0.0;
		
		val += -N_i * ( 
				 v[0]*grad_v[d][0] + 
				 v[1]*grad_v[d][1] + 
				 v[2]*grad_v[d][2] 
				 );
	   
		val += -dyn_vis * ( 
				 grad_N_i[0]*grad_v[d][0] +
				 grad_N_i[1]*grad_v[d][1] +
				 grad_N_i[2]*grad_v[d][2] 
				 );

	   	val += -tau * BB_calc_vec3d_dot(v, grad_N_i) * BB_calc_vec3d_dot(v, grad_v[d]);

		val *= dt;

		val += N_i * v[d];

		val += tau * BB_calc_vec3d_dot(v, grad_N_i) * v[d];

		vec[d] = val;
	}
}


double BBFE_elemmat_fluid_fs_vec_ppe(
		const double N_i,
		const double div_v,
		const double density,
		const double dt)
{
	double val = density/dt * div_v * N_i;
	
	return val;
}


void BBFE_elemmat_fluid_fs_vec_corr(
		double       vec[3],
		const double N_i,
		const double grad_p[3],
		const double v[3],
		const double density,
		const double dt)
{
	for(int d=0; d<3; d++) {
		vec[d] = -1.0/density * N_i * grad_p[d] * dt + v[d] * N_i;
	}
}


/**********************************************************
 * supg + pspg method
 **********************************************************/
void BBFE_elemmat_fluid_sups_mat(
		double         mat[4][4],
		const double   N_i,
		const double   N_j,
		const double   grad_N_i[3],
		const double   grad_N_j[3],
		const double   v[3],
		const double   density,
		const double   viscosity,
		const double   tau,
		const double   dt)
{
	double M = density * N_i * N_j;
	double A = dt * density * N_i * BB_calc_vec3d_dot(v, grad_N_j);

	double G1 = - dt * grad_N_i[0] * N_j;
	double G2 = - dt * grad_N_i[1] * N_j;
	double G3 = - dt * grad_N_i[2] * N_j;

	double D_11 = dt * viscosity * ( BB_calc_vec3d_dot(grad_N_i, grad_N_j) + grad_N_i[0]*grad_N_j[0] );
	double D_12 = dt * viscosity * grad_N_i[1] * grad_N_j[0];
	double D_13 = dt * viscosity * grad_N_i[2] * grad_N_j[0];
	double D_21 = dt * viscosity * grad_N_i[0] * grad_N_j[1];
	double D_22 = dt * viscosity * ( BB_calc_vec3d_dot(grad_N_i, grad_N_j) + grad_N_i[1]*grad_N_j[1] );
	double D_23 = dt * viscosity * grad_N_i[2] * grad_N_j[1];
	double D_31 = dt * viscosity * grad_N_i[0] * grad_N_j[2];
	double D_32 = dt * viscosity * grad_N_i[1] * grad_N_j[2];
	double D_33 = dt * viscosity * ( BB_calc_vec3d_dot(grad_N_i, grad_N_j) + grad_N_i[2]*grad_N_j[2] );

	double C1 = dt * N_i * grad_N_j[0];
	double C2 = dt * N_i * grad_N_j[1];
	double C3 = dt * N_i * grad_N_j[2];

	//SUPG 項
	double M_s = density * tau * BB_calc_vec3d_dot(v, grad_N_i) * N_j;
	double A_s = dt * density * tau * BB_calc_vec3d_dot(v, grad_N_i) * BB_calc_vec3d_dot(v, grad_N_j);

	double G_s1 = dt * tau  * BB_calc_vec3d_dot(v, grad_N_i) * grad_N_j[0];
	double G_s2 = dt * tau  * BB_calc_vec3d_dot(v, grad_N_i) * grad_N_j[1];
	double G_s3 = dt * tau  * BB_calc_vec3d_dot(v, grad_N_i) * grad_N_j[2];

	//PSPG 項
	double M_p1 = tau * grad_N_i[0] * N_j;
	double M_p2 = tau * grad_N_i[1] * N_j;
	double M_p3 = tau * grad_N_i[2] * N_j;

	double A_p1 = dt * tau * grad_N_i[0] * BB_calc_vec3d_dot(v, grad_N_j);
	double A_p2 = dt * tau * grad_N_i[1] * BB_calc_vec3d_dot(v, grad_N_j);
	double A_p3 = dt * tau * grad_N_i[2] * BB_calc_vec3d_dot(v, grad_N_j);

	double G_p = dt * tau * (
			grad_N_i[0] * grad_N_j[0] +
			grad_N_i[1] * grad_N_j[1] +
			grad_N_i[2] * grad_N_j[2]) / density;

	mat[0][0] = M + M_s + A + A_s + D_11;
	mat[0][1] = D_12;
	mat[0][2] = D_13;
	mat[0][3] = G1 + G_s1;
	mat[1][0] = D_21;
	mat[1][1] = M + M_s + A + A_s + D_22;
	mat[1][2] = D_23;
	mat[1][3] = G2 + G_s2;
	mat[2][0] = D_31;
	mat[2][1] = D_32;
	mat[2][2] = M + M_s + A + A_s + D_33;
	mat[2][3] = G3 + G_s3;
	mat[3][0] = C1 + M_p1 + A_p1;
	mat[3][1] = C2 + M_p2 + A_p2;
	mat[3][2] = C3 + M_p3 + A_p3;
	mat[3][3] = G_p;
}


void BBFE_elemmat_fluid_sups_vec(
		double         vec[4],
		const double   N_i,
		const double   grad_N_i[3],
		const double   v[3],
		const double   density,
		const double   tau,
		const double   dt)
{
	for(int d=0; d<3; d++) {
		double val = 0.0;

		val += density * N_i * v[d];

		val += density * tau * BB_calc_vec3d_dot(v, grad_N_i) * v[d];

		vec[d] = val;
	}

	vec[3] = tau * BB_calc_vec3d_dot(grad_N_i, v);
}
