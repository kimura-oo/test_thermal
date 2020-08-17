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


void write_surface_conn(
		BBFE_DATA*  fe,
		SURFACE*    surf,
		const char* directory)
{
	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, FILENAME_SURF, directory);

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

	read_fe_data(&fe, sets.directory);
	memory_allocation_surface(&surf, 
			fe.total_num_nodes, fe.total_num_elems, fe.local_num_nodes, CODENAME);

	get_surface_nodes(&fe, &surf, CODENAME);

	get_surface_info( &fe, &surf, CODENAME);

	memory_allocation_surface_conn(&surf);

	set_surface_conn(&fe, &surf);

	write_surface_vtk(&fe, &surf, sets.directory);

	write_surface_conn(&fe, &surf, sets.directory);
	
	printf("\n");

	return 0;

}
