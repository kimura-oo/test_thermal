#pragma once

#include "FE_dataset.h"
#include "monolis.h"


void BBFE_sys_monowrap_solve(
		MONOLIS*      monolis,
		double*       ans_vec,
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
