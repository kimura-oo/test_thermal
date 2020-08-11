
#include "convdiff.h"
#include "../libBB/calc.h"

#include <math.h>

const double ZERO_CRITERION  = 1.0e-8;


double BBFE_elemmat_convdiff_mat_conv(
		const double  N_i,
		const double  grad_N_j[3],
		const double  a[3])
{
	double val = (  
			a[0] * N_i * grad_N_j[0] +
			a[1] * N_i * grad_N_j[1] +
			a[2] * N_i * grad_N_j[2]);
	
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
		const double a[3],
		const double h_e)
{
	
	double l_a = BB_calc_vec3d_length(a);
	
	if( l_a < ZERO_CRITERION) { return 0.0; }

	double alpha = l_a * h_e/(2.0*k);
	double tilxi = 1.0/tanh(alpha) - 1.0/alpha;

	return (h_e/(2.0*l_a) * tilxi);
}


double BBFE_elemmat_convdiff_mat_stab_conv(
		const double  grad_N_i[3],
		const double  grad_N_j[3],
		const double  a[3],
		const double  tau)
{
	double val = tau*(
			a[0]*grad_N_i[0] * (a[0]*grad_N_j[0] + a[1]*grad_N_j[1] + a[2]*grad_N_j[2]) +
			a[1]*grad_N_i[1] * (a[0]*grad_N_j[0] + a[1]*grad_N_j[1] + a[2]*grad_N_j[2]) +
			a[2]*grad_N_i[2] * (a[0]*grad_N_j[0] + a[1]*grad_N_j[1] + a[2]*grad_N_j[2]) 
			);

	return val;
}


double BBFE_elemmat_convdiff_vec_stab_source(
		const double  grad_N_i[3],
		const double  a[3],
		const double  tau,
		const double  f)
{
	double val = tau * (a[0]*grad_N_i[0] + a[1]*grad_N_i[1] + a[2]*grad_N_i[2]) * f;

	return val;
}
