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

static const char* FILENAME_D_BC = "D_bc.dat";


int main(
		int   argc,
		char* argv[])
{
	SETTINGS set;
	BBFE_DATA fe;
	
	cmd_args_reader_bc(&set, argc, argv, CODENAME, VOIDNAME, FILENAME_D_BC);

	BBFE_sys_read_node(&fe, set.infile_node, set.directory);
	BBFE_sys_read_elem(&fe, set.infile_elem, set.directory, 1);

	bool* node_has_bc;
	node_has_bc = BB_std_calloc_1d_bool(node_has_bc, fe.total_num_nodes);
	
	int num_bc_nodes;
	num_bc_nodes = get_bc_node_list(node_has_bc, &fe);
	
	write_bc_file_const(&fe, node_has_bc, num_bc_nodes, 
			set.block_size, set.bc_value, set.outfile_bc, set.directory);

	BB_std_free_1d_bool(node_has_bc, fe.total_num_nodes);
	
	return 0;
}
