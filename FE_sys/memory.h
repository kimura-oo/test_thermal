#pragma once


#include "FE_dataset.h"

void BBFE_sys_memory_allocation_basis(
		FE_3D_BASIS*    basis,
		const int       num_integ_points,
		const int       pol_order,
		const int       num_nodes_in_elem,
		const int       dimension);

void BBFE_sys_memory_allocation_integ(
		FE_3D_BASIS*    basis,
		const int       num_integ_points,
		const int       dimension);

void BBFE_sys_memory_allocation_shapefunc(
		FE_3D_BASIS*    basis,
		const int       num_nodes_in_elem,
		const int       pol_order,
		const int       num_integ_points);

void BBFE_sys_memory_allocation_node(
		FE_DATA*  fe,
		const int dimension);

void BBFE_sys_memory_allocation_elem(
		FE_DATA*  fe,
		const int num_integ_points,
		const int dimension);
