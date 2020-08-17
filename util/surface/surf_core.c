#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <BB/std.h>
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
			INPUT_FILENAME_NODE,
			directory);
	BBFE_sys_read_elem(
			fe,
			INPUT_FILENAME_ELEM,
			directory,
			1);
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
