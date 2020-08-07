
#include <stdio.h>
#include <stdlib.h>

#include "memory.h"
#include "../libBB/std.h"


void BBFE_sys_memory_allocation_basis(
		FE_3D_BASIS*    basis,
		const int       num_integ_points,
		const int       pol_order,
		const int       num_nodes_in_elem,
		const int       dimension)
{
	BBFE_sys_memory_allocation_integ(
			basis,
			num_integ_points,
			dimension);

	BBFE_sys_memory_allocation_shapefunc(
			basis,
			num_nodes_in_elem,
			pol_order,
			num_integ_points);
}


void BBFE_sys_memory_allocation_integ(
		FE_3D_BASIS*    basis,
		const int       num_integ_points,
		const int       dimension)
{
	int num = num_integ_points;

	basis->num_integ_points = num;

	basis->integ_point  = BB_std_calloc_2d_double(basis->integ_point , num, dimension);
	basis->integ_weight = BB_std_calloc_1d_double(basis->integ_weight, num);
}


void BBFE_sys_memory_allocation_shapefunc(
		FE_3D_BASIS*    basis,
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


void BBFE_sys_memory_allocation_node(
		FE_DATA*  fe,
		const int dimension)
{
	fe->x = BB_std_calloc_2d_double(fe->x, fe->total_num_nodes, dimension);
}


void BBFE_sys_memory_allocation_elem(
		FE_DATA*  fe,
		const int num_integ_points,
		const int dimension)
{
	fe->conn = BB_std_calloc_2d_int(fe->conn, fe->total_num_elems, fe->local_num_nodes);
	
	fe->geo  = (FE_3D_GEO**)calloc(fe->total_num_elems, sizeof(FE_3D_GEO*));
	for(int e=0; e<(fe->total_num_elems); e++) {
		fe->geo[e]  = (FE_3D_GEO*)calloc(num_integ_points, sizeof(FE_3D_GEO));

		for(int p=0; p<num_integ_points; p++) {
			fe->geo[e][p].grad_N = BB_std_calloc_2d_double(fe->geo[e][p].grad_N, fe->local_num_nodes, dimension);
		}
	}
}

