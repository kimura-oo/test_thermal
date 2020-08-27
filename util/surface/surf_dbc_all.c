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

static const char* CODENAME            = "surf_dbc_all >";
static const char* VOIDNAME            = "              ";

static const char* FILENAME_D_BC = "D_bc.dat";

void cmd_args_reader(
		SETTINGS* set,
		int       argc,
		char*     argv[],
		const char* codename,
		const char* voidname)
{
	if(argc < 2) {
		printf("%s Please specify parameters.\n", codename);
		printf("%s Format: \n", voidname);
		printf("%s     %s [block size: n] [bc value 1] [bc value 2] ... [bc value n]\n\n", voidname, argv[0]);
		printf("%s Options: \n", voidname);
		printf("%s     %s [input & output directory]\n", voidname, OPTION_DIRECTORY);
		printf("%s     %s [input filename (nodes)]\n",   voidname, OPTION_INFILE_NODE);
		printf("%s     %s [input filename (surface elements)]\n", voidname, OPTION_INFILE_SURF);
		printf("%s     %s [output filename for Dirichlet B.C.]\n", voidname, OPTION_OUTFILE);
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
		set->infile_elem = FILENAME_ELEM;
		printf("%s Input filename (elements): %s (default)\n", codename, set->infile_elem);
	}
	else {
		set->infile_elem = argv[num+1];
		printf("%s Input filename (elements): %s\n", codename, set->infile_elem);
	}

	num = BB_std_read_args_return_char_num(
			argc, argv, OPTION_OUTFILE);
	if(num == -1) {
		set->outfile_bc = FILENAME_D_BC;
		printf("%s Output filename for Dirichlet B.C.: %s (default)\n", codename, set->outfile_bc);
	}
	else {
		set->outfile_bc = argv[num+1];
		printf("%s Output filename for Dirichlet B.C.: %s\n", codename, set->outfile_bc);
	}

	printf("\n");
}


int main(
		int   argc,
		char* argv[])
{
	printf("\n");

	BBFE_DATA fe;
	SETTINGS  sets;
	SURFACE   surf;

	cmd_args_reader(&sets, argc, argv, CODENAME, VOIDNAME);

	read_fe_data(&fe, sets.directory, sets.infile_node, sets.infile_elem);
	memory_allocation_surface(&surf, 
			fe.total_num_nodes, fe.total_num_elems, fe.local_num_nodes, CODENAME);

	get_surface_nodes(&fe, &surf, CODENAME);

	// function for higher order elements should be implemented!!!
	
	write_bc_file_const(&fe, surf.node_is_on_surface, 
			surf.num_bc_nodes, sets.block_size, sets.bc_value, sets.outfile_bc, sets.directory);

	printf("\n");

	return 0;
}
