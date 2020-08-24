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

static const char* CODENAME            = "surf_dbc >";
static const char* VOIDNAME            = "          ";

static const char* OPTION_DIRECTORY    = "-d";

static const char* OPTION_INFILE_NODE = "-in";
static const char* OPTION_INFILE_SURF = "-ie";
static const char* OPTION_OUTFILE_DBC = "-o";

static const char* DEF_DIRECTORY   = ".";

void cmd_args_reader(
		SETTINGS* set,
		int       argc,
		char*     argv[])
{
	if(argc < 2) {
		printf("%s Please specify parameters.\n", CODENAME);
		printf("%s Format: \n", VOIDNAME);
		printf("%s     %s [block size]\n\n", VOIDNAME, argv[0]);
		printf("%s Options: \n", VOIDNAME);
		printf("%s     %s [input & output directory]\n", VOIDNAME, OPTION_DIRECTORY);
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
		set->directory = argv[num+1];
		printf("%s Input filename (nodes): %s\n", CODENAME, set->infile_node);
	}

	printf("\n");
}


int main(
		int   argc,
		char* argv[])
{
	SETTINGS set;
	BBFE_DATA fe;
	
	cmd_args_reader(&set, argc, argv);
	
}
