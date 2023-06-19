#pragma once

#include "FE_dataset.h"
#include "monolis.h"

void BBFE_sys_monowrap_init_monomat(
		MONOLIS*     monolis,
		MONOLIS_COM* monolis_com,
		BBFE_DATA*   fe,
		const int    block_size,
		const char*  dirname);

void BBFE_sys_monowrap_init_monomat_C(
		MONOLIS*     monolis,
		MONOLIS_COM* monolis_com,
		BBFE_DATA*   fe,
		const int    block_size,
		const char*  dirname);

void BBFE_sys_monowrap_copy_mat(
		MONOLIS* in,
		MONOLIS* out);

void BBFE_sys_monowrap_copy_mat_C(
		MONOLIS* in,
		MONOLIS* out);

void BBFE_sys_monowrap_solve(
		MONOLIS*      monolis,
		MONOLIS_COM*  monolis_com,
		double*       ans_vec,
		const int     solver_type,
		const int     precond_type,
		const int     num_max_iters,
		const double  epsilon);

void BBFE_sys_monowrap_solve_C(
		MONOLIS*      monolis,
		MONOLIS_COM*  monolis_com,
		double _Complex *       ans_vec,
		const int     solver_type,
		const int     precond_type,
		const int     num_max_iters,
		const double  epsilon);

void BBFE_sys_monowrap_set_Dirichlet_bc(
		MONOLIS*      monolis,
		int           num_nodes,
		int           num_dofs_on_node,
		BBFE_BC*      bc,
		double*       g_rhs);

void BBFE_sys_monowrap_set_Dirichlet_bc_C(
		MONOLIS*      monolis,
		int           num_nodes,
		int           num_dofs_on_node,
		BBFE_BC*      bc,
		double _Complex*       g_rhs);

void BBFE_sys_monowrap_set_Neumann_bc(
		int           num_nodes,
		int           num_dofs_on_node,
		BBFE_BC*      bc,
		double*       g_rhs);

double BBFE_sys_monowrap_calc_error_norm(
		MONOLIS*     monolis,
		MONOLIS_COM* monolis_com,
		int        num_nodes,
		int        num_dofs_on_node,
		double*    vec);
