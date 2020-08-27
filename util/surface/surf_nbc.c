#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <BB/std.h>
#include <BB/calc.h>
#include <BB/vtk.h>
#include <BBFE/std/integ.h>
#include <BBFE/std/shapefunc.h>
#include <BBFE/std/mapping.h>
#include <BBFE/std/surface.h>
#include <BBFE/sys/FE_dataset.h>
#include <BBFE/sys/read.h>
#include <BBFE/sys/write.h>
#include <BBFE/sys/memory.h>

#include "surf_core.h"

static const char* CODENAME            = "surf_nbc >";
static const char* VOIDNAME            = "          ";

static const char* FILENAME_N_BC = "N_bc.dat";

void set_basis(
		BBFE_BASIS*   basis,
		int           local_num_nodes,
		int           num_integ_points_each_axis)
{
	switch( local_num_nodes ) {
		case 3:
			basis->num_integ_points = 
				BBFE_std_integ_tri_set_arbitrary_points(
						num_integ_points_each_axis,
						basis->integ_point,
						basis->integ_weight);

			for(int i=0; i<(basis->num_integ_points); i++) {
				BBFE_std_shapefunc_tri1st_get_val(
						basis->integ_point[i],
						basis->N[i]);

				BBFE_std_shapefunc_tri1st_get_derivative(
						basis->integ_point[i],
						basis->dN_dxi[i],
						basis->dN_det[i]);
			}
			printf("%s Element type: 1st-order triangle.\n", CODENAME);
			break;

		case 4:
			basis->num_integ_points = 
				BBFE_std_integ_rec_set_arbitrary_points(
						num_integ_points_each_axis,
						basis->integ_point,
						basis->integ_weight);

			for(int i=0; i<(basis->num_integ_points); i++) {
				BBFE_std_shapefunc_rec1st_get_val(
						basis->integ_point[i],
						basis->N[i]);

				BBFE_std_shapefunc_rec1st_get_derivative(
						basis->integ_point[i],
						basis->dN_dxi[i],
						basis->dN_det[i]);
			}
			printf("%s Element type: 1st-order rectangle.\n", CODENAME);
			break;

		case 6:
		case 9:
			// should be implemented for higher order elements
			break;
	}
	printf("%s The number of integration points: %d\n", CODENAME, basis->num_integ_points);
}


void set_local_array_vector(
		double**       local_val,
		BBFE_DATA*     fe,
		double**       val,
		const int      elem_num,
		const int      dimension)
{
	for(int i=0; i<(fe->local_num_nodes); i++) {
		for(int d=0; d<dimension; d++) {
			local_val[i][d] = val[ fe->conn[elem_num][i] ][d];
		}
	}
}


void equivval_surface_ss(
		double*     equiv_val,
		BBFE_DATA*  surf,
		BBFE_BASIS* basis,
		double      n_val)
{
	double* val_ip;
	double* cross_ip;
	val_ip   = BB_std_calloc_1d_double(val_ip  , basis->num_integ_points);
	cross_ip = BB_std_calloc_1d_double(cross_ip, basis->num_integ_points);

	double** local_x;
	local_x = BB_std_calloc_2d_double(local_x, surf->local_num_nodes, 3);

	for(int e=0; e<(surf->total_num_elems); e++) {
		set_local_array_vector(local_x, surf, surf->x, e, 3);

		for(int p=0; p<(basis->num_integ_points); p++) {
			double dx_dxi[3];  double dx_det[3];
			BBFE_std_mapping_vector3d(
					dx_dxi, surf->local_num_nodes, local_x, basis->dN_dxi[p]);
			BBFE_std_mapping_vector3d(
					dx_det, surf->local_num_nodes, local_x, basis->dN_det[p]);

			double cross[3];
			BB_calc_vec3d_cross(cross, dx_dxi, dx_det);
			cross_ip[p] = BB_calc_vec3d_length(cross);
		}

		for(int i=0; i<(surf->local_num_nodes); i++) {
			for(int p=0; p<(basis->num_integ_points); p++) {
				val_ip[p] = basis->N[p][i] * n_val;
			}

			double integ_val = BBFE_std_integ_calc(
					basis->num_integ_points,
					val_ip,
					basis->integ_weight,
					cross_ip);

			equiv_val[ surf->conn[e][i] ] += integ_val;
		}

	}
	BB_std_free_1d_double(val_ip,   basis->num_integ_points);
	BB_std_free_1d_double(cross_ip, basis->num_integ_points);
	BB_std_free_2d_double(local_x,  surf->local_num_nodes, 3);
}


void equivval_surface_sv(
		double**    equiv_val,
		BBFE_DATA*  surf,
		BBFE_BASIS* basis,
		double      n_val)
{
	double** val_ip;
	double**  cross_ip;
	val_ip   = BB_std_calloc_2d_double(val_ip  , 3, basis->num_integ_points);
	cross_ip = BB_std_calloc_2d_double(cross_ip, 3, basis->num_integ_points);

	double** local_x;
	local_x = BB_std_calloc_2d_double(local_x, surf->local_num_nodes, 3);

	for(int e=0; e<(surf->total_num_elems); e++) {
		set_local_array_vector(local_x, surf, surf->x, e, 3);

		for(int p=0; p<(basis->num_integ_points); p++) {
			double dx_dxi[3];  double dx_det[3];
			BBFE_std_mapping_vector3d(
					dx_dxi, surf->local_num_nodes, local_x, basis->dN_dxi[p]);
			BBFE_std_mapping_vector3d(
					dx_det, surf->local_num_nodes, local_x, basis->dN_det[p]);

			double cross[3];
			BB_calc_vec3d_cross(cross, dx_dxi, dx_det);
			
			BB_calc_vec3d_normal_vec(cross);
			for(int d=0; d<3; d++) {
				cross_ip[d][p] = cross[d];
			}
		}

		for(int i=0; i<(surf->local_num_nodes); i++) {
			for(int p=0; p<(basis->num_integ_points); p++) {
				for(int d=0; d<3; d++) {
					val_ip[d][p] = basis->N[p][i] * n_val;
				}
			}

			double integ_val[3];
			for(int d=0; d<3; d++) {
				integ_val[d] = BBFE_std_integ_calc(
					basis->num_integ_points,
					val_ip[d],
					basis->integ_weight,
					cross_ip[d]);
				equiv_val[ surf->conn[e][i] ][d] += integ_val[d];
			}

		}

	}
	BB_std_free_2d_double(val_ip,   3, basis->num_integ_points);
	BB_std_free_2d_double(cross_ip, 3, basis->num_integ_points);

	BB_std_free_2d_double(local_x,  surf->local_num_nodes, 3);
}


void write_surface_data_vtk(
		BBFE_DATA*  surf,
		double*     equiv_val,
		double**    norm_vec,
		const char* filename,
		const char* directory)
{
	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, filename, directory);

	BB_vtk_write_header(fp);
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

	BB_vtk_write_points_3d(fp, surf->total_num_nodes, surf->x);
	BB_vtk_write_cells(fp, surf->total_num_elems, 
			surf->local_num_nodes, surf->conn);

	int cell_type;
	switch(surf->local_num_nodes) {
		case 4:
		case 10:
			cell_type = TYPE_VTK_TRIANGLE;
			break;
		case 8:
		case 27:
			cell_type = TYPE_VTK_QUAD;
			break;
	}

	BB_vtk_write_cell_types(fp, surf->total_num_elems, cell_type);

	fprintf(fp, "POINT_DATA %d\n", surf->total_num_nodes);	

	BB_vtk_write_point_vals_vector(fp, norm_vec, surf->total_num_nodes, "Normal_vector");

	fclose(fp);
}


void write_dbc_data(
		double**    equivval,
		const char* filename,
		const char* directory)
{
	FILE* fp;

	fp = BBFE_sys_write_fopen(fp, filename, directory);

}


int main(
		int   argc,
		char* argv[])
{
	printf("\n");

	int n_axis   = 2;

	BBFE_BASIS basis;
	BBFE_DATA  surf;
	SETTINGS   sets;

	cmd_args_reader_bc(&sets, argc, argv, CODENAME, VOIDNAME, FILENAME_N_BC);

	BBFE_sys_read_node(
			&surf,
			sets.infile_node,
			sets.directory);
	BBFE_sys_read_elem(
			&surf,
			sets.infile_elem,
			sets.directory,
			1);

	BBFE_sys_memory_allocation_integ(
			&basis,
			n_axis*n_axis,
			2);
	BBFE_sys_memory_allocation_shapefunc(
			&basis,
			surf.local_num_nodes,
			1,
			n_axis*n_axis);

	set_basis(&basis, surf.local_num_nodes, n_axis);

	double** equiv_val;
	equiv_val = BB_std_calloc_2d_double(equiv_val, sets.block_size, surf.total_num_nodes);
	double** norm_vec;
	norm_vec = BB_std_calloc_2d_double(norm_vec, surf.total_num_nodes, 3);
	
	//equivval_surface_sv(norm_vec, &surf, &basis, 1.0);
	
	for(int b=0; b<sets.block_size; b++) {
		equivval_surface_ss(equiv_val[b], &surf, &basis, sets.bc_value[b]);
	}
	
	bool* node_has_bc;
	node_has_bc = BB_std_calloc_1d_bool(node_has_bc, surf.total_num_nodes);
	
	int num_bc_nodes;
	num_bc_nodes = get_bc_node_list(node_has_bc, &surf);

	write_bc_file_nonconst(&surf, node_has_bc, num_bc_nodes, sets.block_size,
			equiv_val, sets.outfile_bc, sets.directory);

	BB_std_free_2d_double(equiv_val, sets.block_size, surf.total_num_nodes);
	BB_std_free_2d_double(norm_vec,  surf.total_num_nodes, 3);
	BB_std_free_1d_bool(node_has_bc, surf.total_num_nodes);

	printf("\n");

	return 0;

}
