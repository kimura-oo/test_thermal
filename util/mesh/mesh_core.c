#include <stdio.h>
#include <stdbool.h>

#include "BB/std.h"
#include "BB/vtk.h"
#include "BBFE/sys/memory.h"
#include "BBFE/sys/read.h"
#include "BBFE/sys/write.h"
#include "mesh_core.h"


static const char* CODENAME = "mesh_core >";
static const char* VOIDNAME = "           ";


void args_manager(
		SETTINGS*   set,
		int         argc,
		char*       argv[],
		bool        surface_mesh)
{
	if(argc < 7) {
		printf("%s Please input parameters.\n", CODENAME);
		printf("%s Format: \n", CODENAME);
		printf("%s    %s [x_min] [y_min] [z_min] [x_max] [y_max] [z_max]\n", 
				VOIDNAME, argv[0]);

		printf("%s Options: \n", CODENAME);
		printf("%s     %s [input/output directory]\n", 
				VOIDNAME, OPTION_DIRECTORY);
		printf("%s     %s [input filename (nodes)]\n", 
				VOIDNAME, OPTION_INFILE_NODE);
		printf("%s     %s [input filename (elements)]\n", 
				VOIDNAME, OPTION_INFILE_ELEM);
		printf("%s     %s [output filename (nodes)]\n", 
				VOIDNAME, OPTION_OUTFILE_NODE);
		printf("%s     %s [output filename (elements)]\n", 
				VOIDNAME, OPTION_OUTFILE_ELEM);
		printf("%s     %s [output filename (vtk)]\n", 
				VOIDNAME, OPTION_OUTFILE_VTK);

		printf("\n");
		exit(0);
	}

	set->x_min[0] = atof(argv[1]);
	set->x_min[1] = atof(argv[2]);
	set->x_min[2] = atof(argv[3]);

	set->x_max[0] = atof(argv[4]);
	set->x_max[1] = atof(argv[5]);
	set->x_max[2] = atof(argv[6]);

	int num;
	num = BB_std_read_args_return_char_num(
				argc, argv, OPTION_DIRECTORY);
	if(num == -1) {
		set->directory = DEF_DIRECTORY;
	}
	else {
		set->directory = argv[num+1];
	}

	num = BB_std_read_args_return_char_num(
				argc, argv, OPTION_INFILE_NODE);
	if(num == -1) {
		set->infile_node = DEF_FILENAME_NODE;
	}
	else {
		set->infile_node = argv[num+1];
	}

	num = BB_std_read_args_return_char_num(
				argc, argv, OPTION_INFILE_ELEM);
	if(num == -1) {
		if(surface_mesh) {
			set->infile_elem = DEF_FILENAME_SURF;
		}
		else{
			set->infile_elem = DEF_FILENAME_ELEM;
		}
	}
	else {
		set->infile_elem = argv[num+1];
	}

	num = BB_std_read_args_return_char_num(
				argc, argv, OPTION_OUTFILE_NODE);
	if(num == -1) {
		set->outfile_node = DEF_FILENAME_NODE;
	}
	else {
		set->outfile_node = argv[num+1];
	}

	num = BB_std_read_args_return_char_num(
				argc, argv, OPTION_OUTFILE_ELEM);
	if(num == -1) {
		if(surface_mesh) {
			set->outfile_elem = DEF_FILENAME_SURF;
		}
		else{
			set->outfile_elem = DEF_FILENAME_ELEM;
		}
	}
	else {
		set->outfile_elem = argv[num+1];
	}

	num = BB_std_read_args_return_char_num(
				argc, argv, OPTION_OUTFILE_VTK);
	if(num == -1) {
		if(surface_mesh) {
			set->outfile_vtk = DEF_FILENAME_SVTK;
		}
		else{
			set->outfile_vtk = DEF_FILENAME_VTK;
		}
	}
	else {
		set->outfile_vtk = argv[num+1];
	}
	
	printf("%s x_min: (%e, %e, %e)\n", CODENAME, set->x_min[0], set->x_min[1], set->x_min[2]);
	printf("%s x_max: (%e, %e, %e)\n", CODENAME, set->x_max[0], set->x_max[1], set->x_max[2]);
	printf("%s Input/output directory    : %s\n", CODENAME, set->directory);
	printf("%s Input filename (nodes)    : %s\n", CODENAME, set->infile_node);
	printf("%s Input filename (elements) : %s\n", CODENAME, set->infile_elem);
	printf("%s Output filename (nodes)   : %s\n", CODENAME, set->outfile_node);
	printf("%s Output filename (elements): %s\n", CODENAME, set->outfile_elem);
	printf("%s Output filename (vtk)     : %s\n", CODENAME, set->outfile_vtk);
	
	printf("\n");
}


void open_fe_files(
		BBFE_DATA*  fe,
		const char* filename_node,
		const char* filename_elem,
		const char* directory)
{
	BBFE_sys_read_node(
			fe,
			filename_node,
			directory);
	BBFE_sys_read_elem(
			fe,
			filename_elem,
			directory,
			1);
}


bool point_is_inside_3d(
		double       x[3],
		const double x_min[3],
		const double x_max[3])
{
	if( 
			x_min[0] <= x[0] && x[0] <= x_max[0] && 
			x_min[1] <= x[1] && x[1] <= x_max[1] && 
			x_min[2] <= x[2] && x[2] <= x_max[2] 
			)
	{
		return true;
	}
	else {
		return false;
	}
}


bool point_cloud_is_inside_3d(
		double**     x,
		const int    num_points,
		const double x_min[3],
		const double x_max[3])
{
	for(int i=0; i<num_points; i++) {
		if( point_is_inside_3d(x[i], x_min, x_max) ) {
			// do nothing
		}
		else {
			return false;
		}
	}

	return true;
}


bool element_is_inside_3d(
		BBFE_DATA*   fe,
		const int    elem_num,
		const double x_min[3],
		const double x_max[3])
{
	bool is_inside = false;

	double** local_x;
	local_x = BB_std_calloc_2d_double(local_x, fe->local_num_nodes, 3);
	
	for(int i=0; i<(fe->local_num_nodes); i++) {
		for(int d=0; d<3; d++) {	
			local_x[i][d] = fe->x[ fe->conn[elem_num][i] ][d];
		}
	}

	is_inside = point_cloud_is_inside_3d(
				local_x, fe->local_num_nodes, 
				x_min, x_max);

	BB_std_free_2d_double(local_x, fe->local_num_nodes, 3);

	return is_inside;
}


void extract_elements(
		BBFE_DATA* fe_ext,
		BBFE_DATA* fe_orig,
		const double x_min[3],
		const double x_max[3])
{
	int count = 0;
	for(int e=0; e<(fe_orig->total_num_elems); e++) {
		if( !element_is_inside_3d(
				fe_orig, e, x_min, x_max) )
		{
			count++;
		}
	}
	
	fe_ext->total_num_elems = fe_orig->total_num_elems - count;
	fe_ext->local_num_nodes = fe_orig->local_num_nodes;
	BBFE_sys_memory_allocation_elem(fe_ext, 1, 3);
	printf("%s The number of new elements: %d\n", 
			CODENAME, fe_ext->total_num_elems);

	int elem_num = 0;
	for(int e=0; e<(fe_orig->total_num_elems); e++) {
		if( element_is_inside_3d(
				fe_orig, e, x_min, x_max) )
		{
			for(int i=0; i<(fe_ext->local_num_nodes); i++) {
				fe_ext->conn[elem_num][i] = fe_orig->conn[e][i];
			}
			
			elem_num++;
		}
	}
}


void remove_elements(
		BBFE_DATA* fe_remove,
		BBFE_DATA* fe_orig,
		const double x_min[3],
		const double x_max[3])
{
	int count = 0;
	for(int e=0; e<(fe_orig->total_num_elems); e++) {
		if( element_is_inside_3d(
				fe_orig, e, x_min, x_max) )
		{
			count++;
		}
	}
	
	fe_remove->total_num_elems = fe_orig->total_num_elems - count;
	fe_remove->local_num_nodes = fe_orig->local_num_nodes;
	BBFE_sys_memory_allocation_elem(fe_remove, 1, 3);
	printf("%s The number of new elements: %d\n", 
			CODENAME, fe_remove->total_num_elems);

	int elem_num = 0;
	for(int e=0; e<(fe_orig->total_num_elems); e++) {
		if( !element_is_inside_3d(
				fe_orig, e, x_min, x_max) )
		{
			for(int i=0; i<(fe_remove->local_num_nodes); i++) {
				fe_remove->conn[elem_num][i] = fe_orig->conn[e][i];
			}
			
			elem_num++;
		}
	}
}


void remove_floating_nodes(
		BBFE_DATA* fe_remove,
		BBFE_DATA* fe_orig)
{
	bool* floating_node;
	floating_node = BB_std_calloc_1d_bool(
			floating_node, fe_orig->total_num_nodes);
	int* new_node_num;
	new_node_num  = BB_std_calloc_1d_int(
			new_node_num, fe_orig->total_num_nodes);

	for(int i=0; i<(fe_orig->total_num_nodes); i++) {
		floating_node[i] = true;
	}
	
	for(int e=0; e<(fe_remove->total_num_elems); e++) {
		for(int i=0; i<(fe_remove->local_num_nodes); i++) {
			floating_node[ fe_remove->conn[e][i] ] = false;
		}
	}

	int count = 0;
	for(int i=0; i<(fe_orig->total_num_nodes); i++) {
		if( floating_node[i] ) { count++; }
	}

	fe_remove->total_num_nodes = fe_orig->total_num_nodes - count;
	BBFE_sys_memory_allocation_node(fe_remove, 3);
	printf("%s The number of new nodes: %d\n", 
			CODENAME, fe_remove->total_num_nodes);
	
	int node_num = 0;
	for(int i=0; i<(fe_orig->total_num_nodes); i++) {
		new_node_num[i] = node_num;

		if( !floating_node[i] ) { 
			for(int d=0; d<3; d++) {
				fe_remove->x[node_num][d] = fe_orig->x[i][d];
			}

			node_num++; 
		}
	}
	
	for(int e=0; e<(fe_remove->total_num_elems); e++) {
		for(int i=0; i<(fe_remove->local_num_nodes); i++) {
			fe_remove->conn[e][i] = new_node_num[ fe_remove->conn[e][i] ];
		}
	}

 	BB_std_free_1d_bool(
			floating_node, fe_orig->total_num_nodes);
	 BB_std_calloc_1d_int(
			new_node_num, fe_orig->total_num_nodes);
}


void allocate_and_copy_node_data(
		BBFE_DATA* fe_copied,
		BBFE_DATA* fe_orig)
{
	fe_copied->total_num_nodes = fe_orig->total_num_nodes;
	BBFE_sys_memory_allocation_node(fe_copied, 3);

	for(int i=0; i<(fe_copied->total_num_nodes); i++) {
		for(int d=0; d<3; d++) {
			fe_copied->x[i][d] = fe_orig->x[i][d];
		}
	}
}


void write_node_file(
		BBFE_DATA*  fe,
		const char* filename,
		const char* directory)
{
	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, filename, directory);
	
	fprintf(fp, "%d\n", fe->total_num_nodes);
	
	for(int i=0; i<(fe->total_num_nodes); i++) {
		fprintf(fp, "%.15e %.15e %.15e\n", 
				fe->x[i][0], fe->x[i][1], fe->x[i][2]);
	}

	fclose(fp);
}


void write_element_file(
		BBFE_DATA*  fe,
		const char* filename,
		const char* directory)
{
	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, filename, directory);
	
	fprintf(fp, "%d %d\n", 
			fe->total_num_elems, fe->local_num_nodes);
	
	for(int e=0; e<(fe->total_num_elems); e++) {
		for(int i=0; i<(fe->local_num_nodes); i++) {
			fprintf(fp, "%d ", fe->conn[e][i]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}


int get_cell_type_vtk_3d(
		int num_points_in_cell)
{
	int cell_type = 0;

	switch(num_points_in_cell) {
		case 4: 
		case 10: 
			cell_type = TYPE_VTK_TETRA;      
			break;
		case 8:
		case 27:
			cell_type = TYPE_VTK_HEXAHEDRON; 
			break;
	}

	return cell_type;
}


int get_cell_type_vtk_2d(
		int num_points_in_cell)
{
	int cell_type = 0;

	switch(num_points_in_cell) {
		case 3: 
		case 10: 
			cell_type = TYPE_VTK_TRIANGLE;      
			break;
		case 4: 
		case 27:
			cell_type = TYPE_VTK_QUAD; 
			break;
	}

	return cell_type;
}


void write_vtk_shape(
		BBFE_DATA*  fe,
		const char* filename,
		const char* directory,
		const int   cell_type)
{
	BB_vtk_void();

	FILE* fp;
	
	fp = BBFE_sys_write_fopen(fp, filename, directory);
	BBFE_sys_write_vtk_shape(fp, fe, cell_type);

	fclose(fp);
}


