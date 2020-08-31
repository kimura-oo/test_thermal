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


void cmd_args_reader_bc(
		SETTINGS* set,
		int       argc,
		char*     argv[],
		const char* codename,
		const char* voidname,
		const char* def_filename_bc)
{
	if(argc < 2) {
		printf("%s Please specify parameters.\n", codename);
		printf("%s Format: \n", voidname);
		printf("%s     %s [block size: n] [bc value 1] [bc value 2] ... [bc value n]\n\n", voidname, argv[0]);
		printf("%s Options: \n", voidname);
		printf("%s     %s [input & output directory]\n", voidname, OPTION_DIRECTORY);
		printf("%s     %s [input filename (nodes)]\n",   voidname, OPTION_INFILE_NODE);
		printf("%s     %s [input filename (surface elements)]\n", voidname, OPTION_INFILE_SURF);
		printf("%s     %s [output filename for B.C.]\n", voidname, OPTION_OUTFILE);
		printf("\n");

		exit(0);
	}

	set->block_size = atoi(argv[1]);
	printf("%s Block size: %d\n", codename, set->block_size);

	if(argc < 2 + set->block_size) {
		printf("%s Please specify parameters.\n", codename);
		printf("%s Format: \n", voidname);
		printf("%s     %s [block size: n] [bc value 1] [bc value 2] ... [bc value n]\n\n", voidname, argv[0]);

		exit(0);
	}

	set->bc_value = BB_std_calloc_1d_double(set->bc_value, set->block_size);
	printf("%s B.C. values: ", codename);
	for(int i=0; i<(set->block_size); i++) {
		set->bc_value[i] = atof( argv[i+2] );
		printf("%e, ", set->bc_value[i]);
	}
	printf("\n");


	int num; 
	num = BB_std_read_args_return_char_num(
			argc, argv, OPTION_DIRECTORY);
	if(num == -1) {
		set->directory = DEF_DIRECTORY;
		printf("%s Input & output directory: %s (default)\n", codename, set->directory);
	}
	else {
		set->directory = argv[num+1];
		printf("%s Input & output directory: %s\n", codename, set->directory);
	}

	num = BB_std_read_args_return_char_num(
			argc, argv, OPTION_INFILE_NODE);
	if(num == -1) {
		set->infile_node = FILENAME_NODE;
		printf("%s Input filename (nodes): %s (default)\n", codename, set->infile_node);
	}
	else {
		set->infile_node = argv[num+1];
		printf("%s Input filename (nodes): %s\n", codename, set->infile_node);
	}

	num = BB_std_read_args_return_char_num(
			argc, argv, OPTION_INFILE_SURF);
	if(num == -1) {
		set->infile_elem = FILENAME_SURF;
		printf("%s Input filename (surface elements): %s (default)\n", codename, set->infile_elem);
	}
	else {
		set->infile_elem = argv[num+1];
		printf("%s Input filename (surface elements): %s\n", codename, set->infile_elem);
	}

	num = BB_std_read_args_return_char_num(
			argc, argv, OPTION_OUTFILE);
	if(num == -1) {
		set->outfile_bc = def_filename_bc;
		printf("%s Output filename for B.C.: %s (default)\n", codename, set->outfile_bc);
	}
	else {
		set->outfile_bc = argv[num+1];
		printf("%s Output filename for B.C.: %s\n", codename, set->outfile_bc);
	}

	printf("\n");
}


void read_fe_data(
		BBFE_DATA*  fe,
		const char* directory,
		const char* infile_node,
		const char* infile_elem)
{
	BB_calc_void();

	BBFE_sys_read_node(
			fe,
			infile_node,
			directory);
	BBFE_sys_read_elem(
			fe,
			infile_elem,
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
			printf("%s Element type: 1st-order tetrahedron\n", codename);
			break;
		case 10:
			printf("%s Element type: 2nd-order tetrahedron\n", codename);
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
						BBFE_std_shapefunc_tet2nd_get_surface(loc_conn, s);
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


int get_bc_node_list(
		bool* node_has_bc,
		BBFE_DATA* fe)
{
	for(int i=0; i<(fe->total_num_nodes); i++) {
		node_has_bc[i] = false;
	}

	for(int e=0; e<(fe->total_num_elems); e++) {
		for(int i=0; i<(fe->local_num_nodes); i++) {
			node_has_bc[ fe->conn[e][i] ] = true;
		}
	}

	int count = 0;
	for(int i=0; i<(fe->total_num_nodes); i++) {
		if( node_has_bc[i] ){ count++; }
	}

	return count;
}


void write_bc_file_const(
		BBFE_DATA*  fe,
		bool*       node_has_bc,
		int         num_bc_nodes,
		int         block_size,
		double*     val,
		const char* filename,
		const char* directory)
{
	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, filename, directory);

	int num_bcs = num_bc_nodes * block_size;

	fprintf(fp, "%d %d\n", num_bcs, block_size);

	for(int i=0; i<(fe->total_num_nodes); i++) {
		if(node_has_bc[i]) {
			for(int b=0; b<block_size; b++) {
				fprintf(fp, "%d %d %.15e\n", i, b, val[b]);
			}
		}
	}

	fclose(fp);
}


void write_bc_file_nonconst(
		BBFE_DATA*  fe,
		bool*       node_has_bc,
		int         num_bc_nodes,
		int         block_size,
		double**    val,
		const char* filename,
		const char* directory)
{
	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, filename, directory);

	int num_bcs = num_bc_nodes * block_size;

	fprintf(fp, "%d %d\n", num_bcs, block_size);

	for(int i=0; i<(fe->total_num_nodes); i++) {
		if(node_has_bc[i]) {
			for(int b=0; b<block_size; b++) {
				fprintf(fp, "%d %d %.15e\n", i, b, val[b][i]);
			}
		}
	}

	fclose(fp);
}
