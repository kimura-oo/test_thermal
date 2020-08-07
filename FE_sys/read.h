#pragma once

#include "FE_dataset.h"

void BBFE_sys_read_node(
		FE_DATA*     fe,
		const char*  filename);

void BBFE_sys_read_elem(
		FE_DATA*     fe,
		const char*  filename,
		int          num_integ_points);

void BBFE_sys_read_Dirichlet_bc(
		BC_DATA*     bc,
		const char*  filename,
		const int    total_num_nodes);
