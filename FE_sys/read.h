#pragma once

#include <stdio.h>
#include "FE_dataset.h"


FILE* BBFE_sys_read_fopen(
		FILE*        fp,
		const char*  filename,
		const char*  directory);

FILE* BBFE_sys_read_fopen_without_error(
		FILE*        fp,
		const char*  filename,
		const char*  directory);

void BBFE_sys_read_node(
		BBFE_DATA*   fe,
		const char*  filename,
		const char*  directory);

void BBFE_sys_read_elem(
		BBFE_DATA*   fe,
		const char*  filename,
		const char*  directory,
		int          num_integ_points);

void BBFE_sys_read_Dirichlet_bc(
		BBFE_BC*     bc,
		const char*  filename,
		const char*  directory,
		const int    total_num_nodes,
		const int    block_size);

void BBFE_sys_read_Neumann_bc(
		BBFE_BC*     bc,
		const char*  filename,
		const char*  directory,
		const int    total_num_nodes,
		const int    block_size);
