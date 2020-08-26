
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

