#pragma once

#include "FE_dataset.h"
#include "monolis.h"


void BBFE_sys_monowrap_set_Dirichlet_bc(
		MONOLIS*      monolis,
		int           num_nodes,
		int           num_dofs_on_node,
		BC_DATA*      bc,
		double*       g_rhs);
