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
		FE_DATA*     fe,
		const char*  filename,
		const char*  directory);

void BBFE_sys_read_elem(
		FE_DATA*     fe,
		const char*  filename,
		const char*  directory,
		int          num_integ_points);

void BBFE_sys_read_Dirichlet_bc(
		BC_DATA*     bc,
		const char*  filename,
		const char*  directory,
		const int    total_num_nodes);
