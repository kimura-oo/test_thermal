#include "main.h"

#include <stdio.h>
#include <stdlib.h>


const int DIM = 3;

const int NUM_INTEG_POINTS  = 5;
const int POL_ORDER         = 1;
const int NUM_NODES_IN_ELEM = 4;


/**********************************************************
 * memory allocation
 **********************************************************/
void memory_allocation_basis(
		FE_3D_BASIS* basis,
		const int 		num_integ_points,
		const int 		pol_order,
		const int 		num_nodes_in_elem)
{
	memory_allocation_integ(
			basis, 
			num_integ_points);

	memory_allocation_shapefunc(
			basis,
			num_nodes_in_elem,
			pol_order,
			num_integ_points);
}


void memory_allocation_integ(
		FE_3D_BASIS* 	basis,
		const int 		num_integ_points)
{
	int num = num_integ_points;

	basis->num_integ_points = num;
	
	basis->integ_point  = (double**)calloc(num, sizeof(double*));
	basis->integ_weight = (double* )calloc(num, sizeof(double ));
	for(int i=0; i<num; i++) {
		basis->integ_point[i] = (double*)calloc(DIM, sizeof(double));
	}
}


void memory_allocation_shapefunc(
		FE_3D_BASIS*	basis,
		const int 		num_nodes_in_elem,
		const int 		pol_order,
		const int 		num_integ_points)
{
	int nn = num_nodes_in_elem;
	int ni = num_integ_points;

	basis->num_nodes = nn;
	basis->pol_order = pol_order;

	basis->N      = (double**)calloc(ni, sizeof(double*));
	basis->dN_dxi = (double**)calloc(ni, sizeof(double*));
	basis->dN_det = (double**)calloc(ni, sizeof(double*));
	basis->dN_dze = (double**)calloc(ni, sizeof(double*));
	for(int i=0; i<ni; i++) {
		basis->N[i]      = (double*)calloc(nn, sizeof(double));
		basis->dN_dxi[i] = (double*)calloc(nn, sizeof(double));
		basis->dN_det[i] = (double*)calloc(nn, sizeof(double));
		basis->dN_dze[i] = (double*)calloc(nn, sizeof(double));
	}
}

/**********************************************************
 * initializers
 **********************************************************/
void initialize_basis(
		FE_3D_BASIS* basis)
{
	integ_point_tet_5(
			&(basis->num_integ_points),
			basis->integ_point,
			basis->integ_weight
			);

	for(int i=0; i<(basis->num_integ_points); i++) {
		shapefunc_3d_tet_1st_value(
				basis->integ_point[i],
				basis->N[i]);
	
		for(int j=0; j<(basis->num_nodes); j++) {
			printf("%e, ", basis->N[i][j]);
		}
		printf("\n");
	}

}


/**********************************************************
 * numerical integration
 **********************************************************/
void integ_point_tet_5(
		int* 		num_integ_points,
		double**	integ_point,
		double*		integ_weight)
{
	(*num_integ_points) = 5;
	
	integ_point[0][0] = 1.0/4.0;  integ_point[0][1] = 1.0/4.0;  integ_point[0][2] = 1.0/4.0;
	integ_point[1][0] = 1.0/6.0;  integ_point[1][1] = 1.0/6.0;  integ_point[1][2] = 1.0/6.0;
	integ_point[2][0] = 1.0/6.0;  integ_point[2][1] = 1.0/6.0;  integ_point[2][2] = 1.0/2.0;
	integ_point[3][0] = 1.0/6.0;  integ_point[3][1] = 1.0/2.0;  integ_point[3][2] = 1.0/6.0;
	integ_point[4][0] = 1.0/2.0;  integ_point[4][1] = 1.0/6.0;  integ_point[4][2] = 1.0/6.0;

	integ_weight[0] = -2.0/15.0;
	integ_weight[1] =  3.0/40.0;
	integ_weight[2] =  3.0/40.0;
	integ_weight[3] =  3.0/40.0;
	integ_weight[4] =  3.0/40.0;
}


/**********************************************************
 * shape function
 **********************************************************/
void shapefunc_3d_tet_1st_value(
		const double 	xi[3],
		double* 		N)
{
	N[0] = 1.0 - xi[0] - xi[1] - xi[2];
	N[1] = xi[0];
	N[2] = xi[1];
	N[3] = xi[2];
}


void shapefunc_3d_tet_1st_der_value(
		const double 	xi[3],
		double* 		dN_dxi,
		double* 		dN_det,
		double* 		dN_dze)
{
	dN_dxi[0] = -1.0;  dN_det[0] = -1.0;  dN_dze[0] = -1.0; 
	dN_dxi[1] =  1.0;  dN_det[1] =  0.0;  dN_dze[1] =  0.0; 
	dN_dxi[2] =  0.0;  dN_det[2] =  1.0;  dN_dze[2] =  0.0; 
	dN_dxi[3] =  0.0;  dN_det[3] =  0.0;  dN_dze[3] =  1.0; 
}


int main (
		int argc, 
		char* argv[])
{
 
	printf("\n Test Thermal \n");
	
	FE_SYSTEM sys;

	memory_allocation_basis(
			&(sys.basis),
			NUM_INTEG_POINTS,
			POL_ORDER,
			NUM_NODES_IN_ELEM);

	initialize_basis(
			&(sys.basis));
	
	printf("\n");

	return 0;
}
