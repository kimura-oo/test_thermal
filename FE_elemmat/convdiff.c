
#include "convdiff.h"
#include "BB/calc.h"

#include <math.h>

const double ZERO_CRITERION  = 1.0e-10;


double BBFE_elemmat_convdiff_mat_conv(
		const double  N_i,
		const double  grad_N_j[3],
		const double  a,
		const double  v[3])
{
	double val = a*(  
			v[0] * N_i * grad_N_j[0] +
			v[1] * N_i * grad_N_j[1] +
			v[2] * N_i * grad_N_j[2]);
	
	return val;
}


double BBFE_elemmat_convdiff_mat_diff(
		const double  grad_N_i[3],
		const double  grad_N_j[3],
		const double  k)
{
	double val = -k*(  
			grad_N_i[0] * grad_N_j[0] +
			grad_N_i[1] * grad_N_j[1] +
			grad_N_i[2] * grad_N_j[2]);

	return val;
}


double BBFE_elemmat_convdiff_vec_source(
		const double  N_i,
		const double  f)
{
	double val = N_i * f;

	return val;
}


double BBFE_elemmat_convdiff_stab_coef(
		const double k,
		const double a,
		const double v[3],
		const double h_e)
{
	
	double l_v = a * BB_calc_vec3d_length(v);
	
	if( l_v < ZERO_CRITERION) { return 0.0; }

	double alpha = l_v * h_e/(2.0*k);
	double tilxi = 1.0/tanh(alpha) - 1.0/alpha;

	return (h_e/(2.0*l_v) * tilxi);
}


double BBFE_elemmat_convdiff_mat_stab_conv(
		const double  grad_N_i[3],
		const double  grad_N_j[3],
		const double  a,
		const double  v[3],
		const double  tau)
{
	double val = tau*a*a*(
			v[0]*grad_N_i[0] * (v[0]*grad_N_j[0] + v[1]*grad_N_j[1] + v[2]*grad_N_j[2]) +
			v[1]*grad_N_i[1] * (v[0]*grad_N_j[0] + v[1]*grad_N_j[1] + v[2]*grad_N_j[2]) +
			v[2]*grad_N_i[2] * (v[0]*grad_N_j[0] + v[1]*grad_N_j[1] + v[2]*grad_N_j[2]) 
			);

	return val;
}


double BBFE_elemmat_convdiff_vec_stab_source(
		const double  grad_N_i[3],
		const double  a,
		const double  v[3],
		const double  tau,
		const double  f)
{
	double val = a*tau*(v[0]*grad_N_i[0] + v[1]*grad_N_i[1] + v[2]*grad_N_i[2]) * f;

	return val;
}


/**********************************************************
 * non-steady
 **********************************************************/
double BBFE_elemmat_convdiff_mat_mass(
		const double  N_i,
		const double  N_j,
		const double  a)
{
	double val = a * N_i * N_j;
	
	return val;
}


double BBFE_elemmat_convdiff_vec_mass(
		const double  N_i,
		const double  T,
		const double  a)
{
	double val = a * N_i * T;
	
	return val;
}


double BBFE_elemmat_convdiff_stab_coef_ns(
		const double k,
		const double v[3],
		const double a, 
		const double h_e,
		const double dt)
{
	
	double l_v = BB_calc_vec3d_length(v);
	
	if( l_v < ZERO_CRITERION) { return 0.0; }

	double sub_k = k/a;
		

	double denom = (2.0/dt)*(2.0/dt) + (2.0*l_v/h_e)*(2.0*l_v/h_e) + (4.0*sub_k/(h_e*h_e))*(4.0*sub_k/(h_e*h_e));

	double val = sqrt(1.0/denom);

	return (val);
}


double BBFE_elemmat_convdiff_mat_stab_mass(
		const double  grad_N_i[3],
		const double  N_j,
		const double  a,
		const double  v[3],
		const double  tau)
{
	double val = tau * a * a * (v[0]*grad_N_i[0] + v[1]*grad_N_i[1] + v[2]*grad_N_i[2]) *  N_j;
	
	return val;
}


double BBFE_elemmat_convdiff_vec_stab_mass(
		const double  grad_N_i[3],
		const double  a,
		const double  v[3],
		const double  T,
		const double  tau)
{
	double val = tau * a * a * (v[0]*grad_N_i[0] + v[1]*grad_N_i[1] + v[2]*grad_N_i[2]) *  T;
	
	return val;
}
