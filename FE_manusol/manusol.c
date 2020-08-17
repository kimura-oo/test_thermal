
#include "manusol.h"
#include "BB/std.h"
#include "BBFE/std/surface.h"
#include "BBFE/sys/write.h"
#include "BBFE/elemmat/equivval.h"

#include <stdio.h>
#include <math.h>

static const char* CODENAME = "FE_manusol/manusol >";


void BBFE_manusol_calc_nodal_error_scalar(
		BBFE_DATA*    fe,
		double*       error,
		double*       theo_sol,
		const double* val)
{
	for(int i=0; i<fe->total_num_nodes; i++) {
		error[i] = val[i] - theo_sol[i];
	}

}


void BBFE_manusol_set_bc_scalar(
		BBFE_DATA* fe,
		BBFE_BC*   bc,
		double*    theo_sol,
		double     t)
{
	for(int i=0; i<(fe->total_num_nodes); i++) {
		if( bc->D_bc_exists[i] ) {
			double x[3];
			for(int d=0; d<3; d++) {
				x[d] = fe->x[i][d];
			}

			bc->imposed_D_val[i] = theo_sol[i];
		}
	}
}

