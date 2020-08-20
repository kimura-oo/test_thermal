#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <BB/std.h>
#include <BB/calc.h>
#include <BB/vtk.h>
#include <BBFE/std/shapefunc.h>
#include <BBFE/std/surface.h>
#include <BBFE/sys/FE_dataset.h>
#include <BBFE/sys/read.h>
#include <BBFE/sys/write.h>

#include "surf_core.h"


void read_fe_data(
		BBFE_DATA*  fe,
		const char* directory)
{
	BBFE_sys_read_node(
			fe,
			FILENAME_NODE,
			directory);
	BBFE_sys_read_elem(
			fe,
			FILENAME_ELEM,
			directory,
			1);
}


void memory_allocation_surface(
		SURFACE*    surf,
		const int   total_num_nodes,
		const int   total_num_elems,
		const int   local_num_nodes,
		const char* codename)
{
	surf->node_is_on_surface = 
		BB_std_calloc_1d_bool(surf->node_is_on_surface, total_num_nodes);

	for(int i=0; i<total_num_nodes; i++) {
		surf->node_is_on_surface[i] = false;
	}

	surf->num_nodes_on_surf = 
		BBFE_std_surface_get_num_nodes_on_surf(local_num_nodes);
	surf->num_surfs_in_elem = 
		BBFE_std_surface_get_num_surfs_in_elem(local_num_nodes);

	if( surf->num_nodes_on_surf == 0 ) {
		printf("%s ERROR: unknown element type (num. nodes in element: %d\n)", 
				codename, local_num_nodes);
		exit(EXIT_FAILURE);
	}

	surf->surf_is_on_surface = BB_std_calloc_2d_bool(
			surf->surf_is_on_surface, total_num_elems, surf->num_surfs_in_elem);
}


void get_surface_nodes(
		BBFE_DATA*  fe,
		SURFACE*    surf,
		const char* codename)
{
	int num_bcs = 0;

	num_bcs = BBFE_std_surface_get_surface_node_3d(
			surf->node_is_on_surface,
			fe->total_num_nodes,
			fe->x,
			fe->total_num_elems,
			fe->local_num_nodes,
			fe->conn);

	switch(fe->local_num_nodes) {
		case 4:
		case 10:
			printf("%s Element type: tetrahedron\n", codename);
			break;

		case 8:
		case 27:
			printf("%s Element type: hexahedron\n", codename);
			break;

		default:
			printf("%s ERROR: unknown element type (num. nodes in element: %d\n)", 
					codename, fe->local_num_nodes);
			exit(EXIT_FAILURE);
	}

	surf->num_bc_nodes = num_bcs;
}


void get_surface_info(
		BBFE_DATA*  fe,
		SURFACE*    surf,
		const char* codename)
{
	int num_surfs = 0;

	num_surfs = BBFE_std_surface_get_surface(
			surf->surf_is_on_surface,
			surf->node_is_on_surface,
			fe->total_num_elems,
			fe->local_num_nodes,
			fe->conn);

	surf->num_bc_surfs = num_surfs;
}


void memory_allocation_surface_conn(
		SURFACE* surf)
{
	surf->conn_surf = BB_std_calloc_2d_int(
			surf->conn_surf, surf->num_bc_surfs, surf->num_nodes_on_surf);
	surf->orig_elem_num = BB_std_calloc_1d_int(
			surf->orig_elem_num, surf->num_bc_surfs);
	surf->orig_surf_num = BB_std_calloc_1d_int(
			surf->orig_surf_num, surf->num_bc_surfs);
}


void set_surface_conn(
		BBFE_DATA* fe,
		SURFACE*   surf)
{
	int c = 0;

	for(int e=0; e<(fe->total_num_elems); e++) {
		for(int s=0; s<(surf->num_surfs_in_elem); s++) {
			if( surf->surf_is_on_surface[e][s] ) {
				surf->orig_elem_num[c] = e;
				surf->orig_surf_num[c] = s;

				int loc_conn[ surf->num_nodes_on_surf ];

				switch(fe->local_num_nodes) {
					case 4:
						BBFE_std_shapefunc_tet1st_get_surface(loc_conn, s);
						break;
					case 10:
						// not implemented...
						//BBFE_std_shapefunc_tet2nd_get_surface(loc_conn, s);
						break;
					case 8:
						BBFE_std_shapefunc_hex1st_get_surface(loc_conn, s);
						break;
					case 27:
						// not implemented...
						//BBFE_std_shapefunc_tet2nd_get_surface(loc_conn, s);
						break;
				}

				for(int i=0; i<(surf->num_nodes_on_surf); i++) {
					surf->conn_surf[c][i] = fe->conn[e][ loc_conn[i] ];
				}
				c++;
			}
		}
	}
}


void write_surface_vtk(
		BBFE_DATA*  fe,
		SURFACE*    surf,
		const char* directory)
{
	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, OUTPUT_FILENAME_SURF_VTK, directory);

	BB_vtk_write_header(fp);
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

	BB_vtk_write_points_3d(fp, fe->total_num_nodes, fe->x);
	BB_vtk_write_cells(fp, surf->num_bc_surfs, 
			surf->num_nodes_on_surf, surf->conn_surf);

	int cell_type;
	switch(fe->local_num_nodes) {
		case 4:
		case 10:
			cell_type = TYPE_VTK_TRIANGLE;
			break;
		case 8:
		case 27:
			cell_type = TYPE_VTK_QUAD;
			break;
	}

	BB_vtk_write_cell_types(fp, surf->num_bc_surfs, cell_type);

	fclose(fp);
}
