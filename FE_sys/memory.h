#pragma once


#include "FE_dataset.h"

void BBFE_sys_memory_allocation_integ(
		BBFE_BASIS*     basis,
		const int       num_integ_points,
		const int       dimension);

void BBFE_sys_memory_free_integ(
		BBFE_BASIS*     basis,
		const int       dimension);

void BBFE_sys_memory_allocation_shapefunc(
		BBFE_BASIS*     basis,
		const int       num_nodes_in_elem,
		const int       pol_order,
		const int       num_integ_points);

void BBFE_sys_memory_free_shapefunc(
		BBFE_BASIS*    basis);

void BBFE_sys_memory_allocation_node(
		BBFE_DATA*  fe,
		const int   dimension);

void BBFE_sys_memory_free_node(
		BBFE_DATA*  fe,
		const int   dimension);

void BBFE_sys_memory_allocation_elem(
		BBFE_DATA*  fe,
		const int   num_integ_points,
		const int   dimension);

void BBFE_sys_memory_free_elem(
		BBFE_DATA*  fe,
		const int   num_integ_points,
		const int   dimension);

void BBFE_sys_memory_allocation_Dirichlet_bc(
		BBFE_BC*   bc,
		const int  total_num_nodes,
		const int  block_size);

void BBFE_sys_memory_free_Dirichlet_bc(
		BBFE_BC*   bc,
		const int  total_num_nodes,
		const int  block_size);

void BBFE_sys_memory_allocation_Neumann_bc(
		BBFE_BC*   bc,
		const int  total_num_nodes,
		const int  block_size);

void BBFE_sys_memory_free_Neumann_bc(
		BBFE_BC*   bc,
		const int  total_num_nodes,
		const int  block_size);
