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

static const char* CODENAME            = "surf_dbc >";
static const char* VOIDNAME            = "          ";

static const char* OPTION_DIRECTORY    = "-o";
static const char* DEFAULT_DIRECTORY   = ".";


void write_Dirichlet_bc_data(
		BBFE_DATA*  fe,
		SURFACE*    surf,
		const int   block_size,
		const char* directory)
{
	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, FILENAME_D_BC, directory);

	printf("\n%s Writing Dirichlet B.C. data\n", CODENAME);
	printf("%s     Num. B.C. nodes: %d\n", CODENAME, surf->num_bc_nodes);
	printf("%s     Block size     : %d\n", CODENAME, block_size);

	fprintf(fp, "%d %d\n", block_size*(surf->num_bc_nodes), block_size);
	for(int i=0; i<(fe->total_num_nodes); i++) {
		if( surf->node_is_on_surface[i] ) {
			for(int b=0; b<block_size; b++) {
				fprintf(fp, "%d %d %e\n", i, b, 0.0);
			}
		}
	}

	fclose(fp);
}


void cmd_args_reader(
		SETTINGS* sets,
		int       argc,
		char*     argv[])
{
	if(argc < 2) {
		printf("%s Please specify parameters.\n", CODENAME);
		printf("%s Format: \n", VOIDNAME);
		printf("%s     ./surf_dbc [block size]\n\n", VOIDNAME);
		printf("%s Options: \n", VOIDNAME);
		printf("%s     %s [input & output directory]\n", VOIDNAME, OPTION_DIRECTORY);
		printf("\n");

		exit(0);
	}

	sets->block_size = atoi(argv[1]);
	printf("%s Block size: %d\n", CODENAME, sets->block_size);

	int num = BB_std_read_args_return_char_num(
			argc, argv, OPTION_DIRECTORY);
	if(num == -1) {
		printf("%s Input & output directory is not specified.\n", CODENAME);
		sets->directory = DEFAULT_DIRECTORY;
		printf("%s Input & output directory: %s (default)\n", CODENAME, sets->directory);
	}
	else {
		sets->directory = argv[num+1];
		printf("%s Input & output directory: %s\n", CODENAME, sets->directory);
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

	cmd_args_reader(&sets, argc, argv);

	read_fe_data(&fe, sets.directory);
	memory_allocation_surface(&surf, 
			fe.total_num_nodes, fe.total_num_elems, fe.local_num_nodes, CODENAME);

	get_surface_nodes(&fe, &surf, CODENAME);

	// function for higher order elements should be implemented!!!
	
	write_Dirichlet_bc_data(&fe, &surf, sets.block_size, sets.directory);

	printf("\n");

	return 0;
}
