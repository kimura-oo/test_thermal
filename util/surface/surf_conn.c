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

static const char* CODENAME            = "surf_nbc_eq >";
static const char* VOIDNAME            = "             ";

static const char* OPTION_DIRECTORY    = "-o";
static const char* DEFAULT_DIRECTORY   = ".";


void memory_allocation_surface(
		SURFACE*  surf,
		const int total_num_nodes,
		const int total_num_elems,
		const int local_num_nodes)
{
	surf->node_is_on_surface = 
		BB_std_calloc_1d_bool(surf->node_is_on_surface, total_num_nodes);

	for(int i=0; i<total_num_nodes; i++) {
		surf->node_is_on_surface[i] = false;
	}
	
	int num_surfs;
	switch(local_num_nodes) {
		case 4:
		case 10:
			num_surfs = 4;
			break;

		case 8:
		case 27:
			num_surfs = 6;
			break;

		default:
			printf("%s ERROR: unknown element type (num. nodes in element: %d\n)", 
					CODENAME, local_num_nodes);
			exit(EXIT_FAILURE);
	}
	surf->surf_is_on_surface = BB_std_calloc_2d_bool(
			surf->surf_is_on_surface, total_num_elems, num_surfs);
}


void cmd_args_reader(
		SETTINGS* sets,
		int       argc,
		char*     argv[])
{
	if(argc < 2) {
		printf("%s Please specify parameters.\n", CODENAME);
		printf("%s Format: \n", VOIDNAME);
		printf("%s     ./surf_nbc_eq [block size]\n\n", VOIDNAME);
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
			fe.total_num_nodes, fe.total_num_elems, fe.local_num_nodes);

	get_surface_nodes(&fe, &surf, CODENAME);

	get_surface_info( &fe, &surf, CODENAME);

	printf("%d\n", surf.num_bc_surfs);
	
	printf("\n");

	return 0;

}
