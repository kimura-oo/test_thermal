
#include "monowrap.h"
#include <math.h>


void BBFE_sys_monowrap_init_monomat(
		MONOLIS*    monolis,
		MONOLIS_COM* monolis_com,
		BBFE_DATA*  fe,
		const int   block_size,
		const char* dirname)
{
	monolis_initialize(monolis);

	monolis_com_initialize_by_parted_files(
		monolis_com,
		monolis_mpi_get_global_comm(),
		MONOLIS_DEFAULT_TOP_DIR,
		MONOLIS_DEFAULT_PART_DIR,
		"node.dat");

	monolis_get_nonzero_pattern_by_simple_mesh_R(
			monolis,
			fe->total_num_nodes,
			fe->local_num_nodes,
			block_size,
			fe->total_num_elems,
			fe->conn);
}


void BBFE_sys_monowrap_init_monomat_C(
		MONOLIS*    monolis,
		MONOLIS_COM* monolis_com,
		BBFE_DATA*  fe,
		const int   block_size,
		const char* dirname)
{
	monolis_initialize(monolis);

	monolis_com_initialize_by_parted_files(
		monolis_com,
		monolis_mpi_get_global_comm(),
		MONOLIS_DEFAULT_TOP_DIR,
		MONOLIS_DEFAULT_PART_DIR,
		"node.dat");

	monolis_get_nonzero_pattern_by_simple_mesh_C(
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
	monolis_copy_mat_value_R(in, out);
}


void BBFE_sys_monowrap_copy_mat_C(
		MONOLIS* in,
		MONOLIS* out)
{
	monolis_copy_mat_value_C(in, out);
}


void BBFE_sys_monowrap_solve(
		MONOLIS*      monolis,
		MONOLIS_COM*  monolis_com,
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
	monolis_show_iterlog (monolis, false);

	monolis_solve_R(
			monolis,
			monolis_com,
			monolis->mat.R.B,
			ans_vec);
}


void BBFE_sys_monowrap_solve_C(
		MONOLIS*      monolis,
		MONOLIS_COM* monolis_com,
		double _Complex *       ans_vec,
		const int     solver_type,
		const int     precond_type,
		const int     num_max_iters,
		const double  epsilon)
{
	monolis_set_method   (monolis, solver_type);
	monolis_set_precond  (monolis, precond_type);
	monolis_set_maxiter  (monolis, num_max_iters);
	monolis_set_tolerance(monolis, epsilon);

	monolis_solve_C(
			monolis,
			monolis_com,
			monolis->mat.C.B,
			ans_vec);
}


void BBFE_sys_monowrap_set_Dirichlet_bc(
		MONOLIS*      monolis,
		int           num_nodes,
		int           num_dofs_on_node,
		BBFE_BC*      bc,
		double*        g_rhs)
	{
	for(int i=0; i<num_nodes; i++) {
		for(int k=0; k<num_dofs_on_node; k++) {
			if( bc->D_bc_exists[ num_dofs_on_node*i+k ] ) {
				monolis_set_Dirichlet_bc_R(
						monolis,
						g_rhs,
						i,
						k,
						bc->imposed_D_val[ num_dofs_on_node*i+k ]);
			}
		}
	}
}

void BBFE_sys_monowrap_set_Dirichlet_bc_C(
		MONOLIS*      monolis,
		int           num_nodes,
		int           num_dofs_on_node,
		BBFE_BC*      bc,
		double _Complex*       g_rhs)
	{
	for(int i=0; i<num_nodes; i++) {
		for(int k=0; k<num_dofs_on_node; k++) {
			if( bc->D_bc_exists[ num_dofs_on_node*i+k ] ) {
				monolis_set_Dirichlet_bc_C(
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
		MONOLIS*     monolis,
		MONOLIS_COM* monolis_com,
		int        num_nodes,
		int        num_dofs_on_node,
		double*    vec)
{
	double norm = 0.0;

	monolis_inner_product_R(
		monolis,
		monolis_com,
		3,
		vec,
		vec,
		&norm);

	norm = sqrt(norm);

	return norm;
}
