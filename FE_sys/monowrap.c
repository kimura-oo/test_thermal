
#include "monowrap.h"


void BBFE_sys_monowrap_solve(
		MONOLIS*      monolis,
		double*       ans_vec,
		const int     solver_type,
		const int     precond_type,
		const int     num_max_iters,
		const double  epsilon)
{
	monolis_set_method   (monolis, solver_type);
	monolis_set_precond  (monolis, precond_type);
	monolis_set_maxiter  (monolis, num_max_iters);
	monolis_set_tolerance(monolis, epsilon);

	monolis_solve(
			monolis,
			monolis->mat.B,
			ans_vec);
}


void BBFE_sys_monowrap_set_Dirichlet_bc(
		MONOLIS*      monolis,
		int           num_nodes,
		int           num_dofs_on_node,
		BBFE_BC*      bc,
		double*       g_rhs)
	{
	for(int i=0; i<num_nodes; i++) {
		for(int k=0; k<num_dofs_on_node; k++) {
			if( bc->D_bc_exists[ num_dofs_on_node*i+k ] ) {
				monolis_set_Dirichlet_bc(
						monolis,
						g_rhs,
						i,
						k,
						bc->imposed_D_val[ num_dofs_on_node*i+k ]);
			}
		}
	}
}
