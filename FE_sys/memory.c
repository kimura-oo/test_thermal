
#include <stdio.h>
#include <stdlib.h>

#include "memory.h"
#include "BB/std.h"


void BBFE_sys_memory_allocation_integ(
		BBFE_BASIS*     basis,
		const int       num_integ_points,
		const int       dimension)
{
	int num = num_integ_points;

	basis->num_integ_points = num;

	basis->integ_point  = BB_std_calloc_2d_double(basis->integ_point , num, dimension);
	basis->integ_weight = BB_std_calloc_1d_double(basis->integ_weight, num);
}


void BBFE_sys_memory_free_integ(
		BBFE_BASIS*     basis,
		const int       dimension)
{
	int num = basis->num_integ_points;

	BB_std_free_2d_double(basis->integ_point , num, dimension);
	BB_std_free_1d_double(basis->integ_weight, num);
}


void BBFE_sys_memory_allocation_shapefunc(
		BBFE_BASIS*     basis,
		const int       num_nodes_in_elem,
		const int       pol_order,
		const int       num_integ_points)
{
	int nn = num_nodes_in_elem;
	int ni = num_integ_points;

	basis->num_nodes = nn;
	basis->pol_order = pol_order;

	basis->N      = BB_std_calloc_2d_double(basis->N     , ni, nn);
	basis->dN_dxi = BB_std_calloc_2d_double(basis->dN_dxi, ni, nn);
	basis->dN_det = BB_std_calloc_2d_double(basis->dN_det, ni, nn);
	basis->dN_dze = BB_std_calloc_2d_double(basis->dN_dze, ni, nn);
}


void BBFE_sys_memory_free_shapefunc(
		BBFE_BASIS*    basis)
{
	int nn = basis->num_nodes;
	int ni = basis->num_integ_points;

	BB_std_free_2d_double(basis->N     , ni, nn);
	BB_std_free_2d_double(basis->dN_dxi, ni, nn);
	BB_std_free_2d_double(basis->dN_det, ni, nn);
	BB_std_free_2d_double(basis->dN_dze, ni, nn);
}


void BBFE_sys_memory_allocation_node(
		BBFE_DATA*  fe,
		const int   dimension)
{
	fe->x = BB_std_calloc_2d_double(fe->x, fe->total_num_nodes, dimension);
}


void BBFE_sys_memory_free_node(
		BBFE_DATA*  fe,
		const int   dimension)
{
	BB_std_free_2d_double(fe->x, fe->total_num_nodes, dimension);
}


void BBFE_sys_memory_allocation_elem(
		BBFE_DATA*  fe,
		const int   num_integ_points,
		const int   dimension)
{
	fe->conn = BB_std_calloc_2d_int(fe->conn, fe->total_num_elems, fe->local_num_nodes);
	
	fe->geo  = (BBFE_GEO**)calloc(fe->total_num_elems, sizeof(BBFE_GEO*));
	for(int e=0; e<(fe->total_num_elems); e++) {
		fe->geo[e]  = (BBFE_GEO*)calloc(num_integ_points, sizeof(BBFE_GEO));

		for(int p=0; p<num_integ_points; p++) {
			fe->geo[e][p].grad_N = BB_std_calloc_2d_double(fe->geo[e][p].grad_N, fe->local_num_nodes, dimension);
		}
	}
}


void BBFE_sys_memory_free_elem(
		BBFE_DATA*  fe,
		const int   num_integ_points,
		const int   dimension)
{
	BB_std_free_2d_int(fe->conn, fe->total_num_elems, fe->local_num_nodes);
	
	for(int e=0; e<(fe->total_num_elems); e++) {
		for(int p=0; p<num_integ_points; p++) {
			BB_std_calloc_2d_double(fe->geo[e][p].grad_N, fe->local_num_nodes, dimension);
		}
		free(fe->geo[e]);
	}
	free(fe->geo);
}


void BBFE_sys_memory_allocation_Dirichlet_bc(
		BBFE_BC*   bc,
		const int  total_num_nodes,
		const int  block_size)
{
	int n = total_num_nodes * block_size;

	bc->D_bc_exists   = BB_std_calloc_1d_bool(  bc->D_bc_exists  , n);
	bc->imposed_D_val = BB_std_calloc_1d_double(bc->imposed_D_val, n);
}


void BBFE_sys_memory_free_Dirichlet_bc(
		BBFE_BC*   bc,
		const int  total_num_nodes,
		const int  block_size)
{
	int n = total_num_nodes * block_size;

	BB_std_free_1d_bool(  bc->D_bc_exists  , n);
	BB_std_free_1d_double(bc->imposed_D_val, n);
}


void BBFE_sys_memory_allocation_Neumann_bc(
		BBFE_BC*   bc,
		const int  total_num_nodes,
		const int  block_size)
{
	int n = total_num_nodes * block_size;

	bc->N_bc_exists   = BB_std_calloc_1d_bool(  bc->N_bc_exists  , n);
	bc->imposed_N_val = BB_std_calloc_1d_double(bc->imposed_N_val, n);
}


void BBFE_sys_memory_free_Neumann_bc(
		BBFE_BC*   bc,
		const int  total_num_nodes,
		const int  block_size)
{
	int n = total_num_nodes * block_size;

	BB_std_free_1d_bool(  bc->N_bc_exists  , n);
	BB_std_free_1d_double(bc->imposed_N_val, n);
}
