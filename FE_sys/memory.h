#pragma once


#include "FE_dataset.h"

void BBFE_sys_memory_allocation_integ(
		FE_3D_BASIS*    basis,
		const int       num_integ_points,
		const int       dimension);

void BBFE_sys_memory_free_integ(
		FE_3D_BASIS*    basis,
		const int       dimension);

void BBFE_sys_memory_allocation_shapefunc(
		FE_3D_BASIS*    basis,
		const int       num_nodes_in_elem,
		const int       pol_order,
		const int       num_integ_points);

void BBFE_sys_memory_free_shapefunc(
		FE_3D_BASIS*    basis);

void BBFE_sys_memory_allocation_node(
		FE_DATA*  fe,
		const int dimension);

void BBFE_sys_memory_free_node(
		FE_DATA*  fe,
		const int dimension);

void BBFE_sys_memory_allocation_elem(
		FE_DATA*  fe,
		const int num_integ_points,
		const int dimension);

void BBFE_sys_memory_free_elem(
		FE_DATA*  fe,
		const int num_integ_points,
		const int dimension);

void BBFE_sys_memory_allocation_Dirichlet_bc(
		BC_DATA*   bc,
		const int  total_num_nodes,
		const int  block_size);

void BBFE_sys_memory_free_Dirichlet_bc(
		BC_DATA*   bc,
		const int  total_num_nodes,
		const int  block_size);
