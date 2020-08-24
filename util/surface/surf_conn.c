#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <BB/std.h>
#include <BBFE/std/shapefunc.h>
#include <BBFE/std/surface.h>
#include <BBFE/sys/FE_dataset.h>
#include <BBFE/sys/read.h>
#include <BBFE/sys/write.h>

#include "surf_core.h"

static const char* CODENAME            = "surf_conn >";
static const char* VOIDNAME            = "           ";


void cmd_args_reader(
		SETTINGS* set,
		int       argc,
		char*     argv[])
{
	if(argc < 2) {
		printf("%s Please specify parameters.\n", CODENAME);
		printf("%s Format: \n", VOIDNAME);
		printf("%s     %s [block size: n]\n\n", VOIDNAME, argv[0]);
		printf("%s Options: \n", VOIDNAME);
		printf("%s     %s [input & output directory]\n", VOIDNAME, OPTION_DIRECTORY);
		printf("%s     %s [input filename (nodes)]\n",   VOIDNAME, OPTION_INFILE_NODE);
		printf("%s     %s [input filename (elements)]\n", VOIDNAME, OPTION_INFILE_SURF);
		printf("%s     %s [output filename (surface elements)]\n", VOIDNAME, OPTION_OUTFILE);
		printf("\n");

		exit(0);
	}

	set->block_size = atoi(argv[1]);
	printf("%s Block size: %d\n", CODENAME, set->block_size);

	int num; 
	num = BB_std_read_args_return_char_num(
			argc, argv, OPTION_DIRECTORY);
	if(num == -1) {
		set->directory = DEF_DIRECTORY;
		printf("%s Input & output directory: %s (default)\n", CODENAME, set->directory);
	}
	else {
		set->directory = argv[num+1];
		printf("%s Input & output directory: %s\n", CODENAME, set->directory);
	}

	num = BB_std_read_args_return_char_num(
			argc, argv, OPTION_INFILE_NODE);
	if(num == -1) {
		set->infile_node = FILENAME_NODE;
		printf("%s Input filename (nodes): %s (default)\n", CODENAME, set->infile_node);
	}
	else {
		set->infile_node = argv[num+1];
		printf("%s Input filename (nodes): %s\n", CODENAME, set->infile_node);
	}

	num = BB_std_read_args_return_char_num(
			argc, argv, OPTION_INFILE_SURF);
	if(num == -1) {
		set->infile_elem = FILENAME_ELEM;
		printf("%s Input filename (elements): %s (default)\n", CODENAME, set->infile_elem);
	}
	else {
		set->infile_elem = argv[num+1];
		printf("%s Input filename (elements): %s\n", CODENAME, set->infile_elem);
	}

	num = BB_std_read_args_return_char_num(
			argc, argv, OPTION_OUTFILE);
	if(num == -1) {
		set->outfile_bc = FILENAME_SURF;
		printf("%s Output filename (surface elements): %s (default)\n", CODENAME, set->outfile_bc);
	}
	else {
		set->outfile_bc = argv[num+1];
		printf("%s Output filename (surface elements): %s\n", CODENAME, set->outfile_bc);
	}

	printf("\n");
}




void write_surface_conn(
		BBFE_DATA*  fe,
		SURFACE*    surf,
		const char* filename,
		const char* directory)
{
	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, filename, directory);

	fprintf(fp, "%d %d\n", surf->num_bc_surfs, surf->num_nodes_on_surf);
	for(int s=0; s<(surf->num_bc_surfs); s++) {
		for(int i=0; i<(surf->num_nodes_on_surf); i++) {
			fprintf(fp, "%d ", surf->conn_surf[s][i]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}


int main(
		int   argc,
		char* argv[])
{
	printf("\n");

	BBFE_DATA fe;
	SETTINGS  sets;
	SURFACE   surf;

	cmd_args_reader(&sets, argc, argv);

	read_fe_data(&fe, sets.directory, sets.infile_node, sets.infile_elem);
	memory_allocation_surface(&surf, 
			fe.total_num_nodes, fe.total_num_elems, fe.local_num_nodes, CODENAME);

	get_surface_nodes(&fe, &surf, CODENAME);

	get_surface_info( &fe, &surf, CODENAME);

	memory_allocation_surface_conn(&surf);

	set_surface_conn(&fe, &surf);

	write_surface_vtk(&fe, &surf, sets.directory);

	write_surface_conn(&fe, &surf, sets.outfile_bc, sets.directory);
	
	printf("\n");

	return 0;

}
