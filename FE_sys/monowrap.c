
#include "monowrap.h"
#include <math.h>


void BBFE_sys_monowrap_init_monomat(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		const int   block_size,
		const char* dirname)
{
	monolis_initialize(monolis, dirname);
	
	monolis_get_nonzero_pattern(
			monolis,
			fe->total_num_nodes,
			fe->local_num_nodes,
			block_size,
			fe->total_num_elems,
			fe->conn);
}


void BBFE_sys_monowrap_copy_mat(
		MONOLIS* in,
		MONOLIS* out)
{
	int num_nonzeros = (in->mat.NDOF) * (in->mat.NDOF) * (in->mat.NZ);

	for(int i=0; i<num_nonzeros; i++) {
		out->mat.A[i] = in->mat.A[i];
	}
}


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


void BBFE_sys_monowrap_set_Neumann_bc(
		int           num_nodes,
		int           num_dofs_on_node,
		BBFE_BC*      bc,
		double*       g_rhs)
	{
	for(int i=0; i<num_nodes; i++) {
		for(int k=0; k<num_dofs_on_node; k++) {
			if( bc->N_bc_exists[ num_dofs_on_node*i+k ] ) {
				g_rhs[ num_dofs_on_node*i+k ] += 
					bc->imposed_N_val[ num_dofs_on_node*i+k ];
			}
		}
	}
}


double BBFE_sys_monowrap_calc_error_norm(
		int        num_nodes,
		int        num_dofs_on_node,
		double*    vec)
{
	double norm = 0.0;
	
	for(int i=0; i<num_nodes; i++) {
		for(int k=0; k<num_dofs_on_node; k++) {
			norm += vec[ num_dofs_on_node*i+k ]*vec[ num_dofs_on_node*i+k ];
		}
	}

	norm = sqrt(norm);

	return norm;
}
