
#include "manusol.h"
#include "../libBB/std.h"
#include "../FE_std/surface.h"
#include "../FE_sys/write.h"
#include "../FE_elemmat/equivval.h"

#include <stdio.h>
#include <math.h>

static const char* CODENAME = "FE_manusol/manusol >";


void BBFE_manusol_calc_nodal_error_scalar(
		BBFE_DATA*    fe,
		double*       error,
		double*       theo_sol,
		const double* val)
{
	for(int i=0; i<fe->total_num_nodes; i++) {
		error[i] = val[i] - theo_sol[i];
	}

}


void BBFE_manusol_overwrite_bc_file_hex(
		BBFE_DATA*  fe,
		const int   block_size,
		const char* filename,
		const char* directory)
{
	bool* node_is_on_surface;
	node_is_on_surface = BB_std_calloc_1d_bool(node_is_on_surface, fe->total_num_nodes);
	for(int i=0; i<(fe->total_num_nodes); i++) {
		node_is_on_surface[i] = false;
	}
	
	int num_bcs;
	num_bcs = BBFE_std_surface_hex1st_get_surface_node(
			node_is_on_surface, 
			fe->total_num_nodes, fe->x,
			fe->total_num_elems, fe->conn);

	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, filename, directory);

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


void BBFE_manusol_overwrite_bc_file_tet(
		BBFE_DATA*  fe,
		const int   block_size,
		const char* filename,
		const char* directory)
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

	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, filename, directory);

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


void BBFE_manusol_set_bc_scalar(
		BBFE_DATA* fe,
		BBFE_BC*   bc,
		double*    theo_sol,
		double     t)
{
	for(int i=0; i<(fe->total_num_nodes); i++) {
		if( bc->D_bc_exists[i] ) {
			double x[3];
			for(int d=0; d<3; d++) {
				x[d] = fe->x[i][d];
			}

			bc->imposed_D_val[i] = theo_sol[i];
		}
	}
}

