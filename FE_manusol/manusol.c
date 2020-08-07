
#include "manusol.h"
#include "../libBB/std.h"
#include "../FE_std/surface.h"
#include "../main.h"

#include <stdio.h>
#include <math.h>

static const char* CODENAME = "FE_sys/manusol >";


void BBFE_manusol_calc_nodal_error_scalar(
		FE_DATA*      fe,
		double*       error,
		const double* val)
{
	for(int i=0; i<fe->total_num_nodes; i++) {
		double x[3];
		for(int d=0; d<3; d++) {
			x[d] = fe->x[i][d];
		}
		double theo_sol = BBFE_manusol_get_sol_scalar_3d(x[0], x[1], x[2]);
		error[i] = val[i] - theo_sol;
	}

}


void BBFE_manusol_overwrite_bc_file(
		FE_DATA*   fe,
		const int  block_size)
{
	bool* node_is_on_surface;
	node_is_on_surface = BB_std_calloc_1d_bool(node_is_on_surface, fe->total_num_nodes);
	for(int i=0; i<(fe->total_num_nodes); i++) {
		node_is_on_surface[i] = false;
	}
	
	int num_bcs;
	num_bcs = BBFE_std_surface_tet1st_get_surface_node(
			node_is_on_surface, 
			fe->total_num_nodes, fe->x,
			fe->total_num_elems, fe->conn);

	const char* filename = "bc_D.dat";
	FILE* fp;
	fp = fopen(filename, "w");
	if( fp == NULL ) {
		printf("%s ERROR: File \"%s\" cannot be opened.\n",
				CODENAME, filename);
	}
	fprintf(fp, "%d %d\n", block_size*num_bcs, block_size);
	for(int i=0; i<(fe->total_num_nodes); i++) {
		if(node_is_on_surface[i]) {
			for(int b=0; b<block_size; b++) {
				fprintf(fp, "%d %d %e\n", i, b, 0.0);
			}
		}
	}

	fclose(fp);

	BB_std_free_1d_bool(node_is_on_surface, fe->total_num_nodes);
}


double BBFE_manusol_get_sol_scalar_3d(
		double x,
		double y,
		double z)
{
	double sol = sin( 0.5*x ) * sin( 1.0*y ) * sin( 2.0*z );
	return sol;
}


double BBFE_manusol_get_rhs_scalar_3d(
		double x,
		double y,
		double z)
{
	double rhs = -5.25*sin( 0.5*x ) * sin( y ) * sin( 2.0*z );
	return rhs;
}


void BBFE_manusol_set_bc_scalar(
		FE_DATA* fe,
		BC_DATA* bc)
{
	for(int i=0; i<(fe->total_num_nodes); i++) {
		if( bc->D_bc_exists[i] ) {
			double x[3];
			for(int d=0; d<3; d++) {
				x[d] = fe->x[i][d];
			}

			bc->imposed_D_val[i] =
				BBFE_manusol_get_sol_scalar_3d(x[0], x[1], x[2]);
		}
	}
}


void BBFE_manusol_add_rhs_scalar(
		FE_DATA* fe,
		FE_3D_BASIS* basis,
		double* rhs)
{
	double* equiv_val;
	equiv_val = BB_std_calloc_1d_double(equiv_val, fe->total_num_nodes);
	
	BBFE_elemmat_equivval_volume_smooth_function(
			equiv_val, fe, basis, BBFE_manusol_get_rhs_scalar_3d);

	BB_std_free_1d_double(equiv_val, fe->total_num_nodes);
}
