
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "BB/std.h"

static const int BUFFER_SIZE = 10000;
static const char* CODENAME     = "vtk2dat >";
static const char* ID_NODE      = "POINTS";
static const char* ID_ELEM      = "CELLS";
static const char* OP_DIRECTORY = "-d";
static const char* DEFAULT_DIR  = ".";

static const int L_NUM_NODES = 4;

int main(
		int   argc,
		char* argv[])
{
	printf("\n");

	char buf[BUFFER_SIZE], buf2[BUFFER_SIZE];

	const char* filename = argv[1];

	if( argc < 2 ) {
		printf("%s ERROR: Input file should be specified.\n", CODENAME);
		printf("%s Format:  ./vtk2dat [input filename]\n", CODENAME);
		printf("%s Options: -d [output directory] (default: .)\n\n", CODENAME);
		
		exit(EXIT_FAILURE);
	}

	const char* directory = 
		BB_std_read_args_return_next_arg(argc, argv, OP_DIRECTORY);
	if(directory == NULL) {
		directory = DEFAULT_DIR;
		printf("%s Output directory: %s (default)\n", CODENAME, directory);
	}
	else {
		printf("%s Output directory: %s\n", CODENAME, directory);
	}

	FILE* fp;
	fp = fopen(filename, "r");
	if( fp == NULL ) {		
		printf("%s ERROR: Input file \"%s\" cannot be opened.\n\n", CODENAME, filename);
		exit(EXIT_FAILURE);
	}

	/********************** nodal data ***********************/
	char fname_node[BUFFER_SIZE];
	snprintf(fname_node, BUFFER_SIZE, "%s/node.dat", directory);
	printf("%s Output file (node): %s\n", CODENAME, fname_node);
	
	FILE* fp_node;
	fp_node = fopen(fname_node, "w");
	if( fp_node == NULL ) {		
		printf("%s ERROR: Output file \"%s\" cannot be opened.\n\n", CODENAME, fname_node);
		exit(EXIT_FAILURE);
	}

	fp = BB_std_read_file_search_line(fp, buf, ID_NODE, BUFFER_SIZE);
	if( fp == NULL ) {
		printf("%s ERROR: Identifier \"%s\" cannot be found.\n\n", CODENAME, ID_NODE);
		exit(EXIT_FAILURE);
	}

	int total_num_nodes;
	sscanf(buf, "%s %d %s", buf2, &total_num_nodes, buf2);
	printf("%s The total number of nodes: %d\n", CODENAME, total_num_nodes);
	fprintf(fp_node, "%d\n", total_num_nodes);
	
	for(int i=0; i<total_num_nodes; i++) {
		double x[3];
		BB_std_scan_line(&fp, BUFFER_SIZE, "%lf %lf %lf", &x[0], &x[1], &x[2]);
		fprintf(fp_node, "%.15e %.15e %.15e\n", x[0], x[1], x[2]);
	}
	fclose(fp_node);
	/*********************************************************/
	
	/********************* element data **********************/
	fp = BB_std_read_file_search_line(fp, buf, ID_ELEM, BUFFER_SIZE);
	if( fp == NULL ) {
		printf("%s ERROR: Identifier \"%s\" cannot be found.\n\n", CODENAME, ID_ELEM);
		exit(EXIT_FAILURE);
	}

	int total_num_cells;
	int total_num_elems = 0;
	sscanf(buf, "%s %d %s", buf2, &total_num_cells, buf2);
	for(int c=0; c<total_num_cells; c++) {
		int nl;
		fscanf(fp, "%d", &nl);
		if(nl == L_NUM_NODES) {  total_num_elems++;  }

		for(int l=0; l<nl; l++) {
			int nn;
			fscanf(fp, "%d", &nn);
		}
	}

	fclose(fp);
	fp = fopen(filename, "r");

	char fname_elem[BUFFER_SIZE];
	snprintf(fname_elem, BUFFER_SIZE, "%s/elem.dat", directory);
	printf("%s Output file (elem): %s\n", CODENAME, fname_elem);

	FILE* fp_elem;
	fp_elem = fopen(fname_elem, "w");
	if( fp_elem == NULL ) {		
		printf("%s ERROR: Output file \"%s\" cannot be opened.\n\n", CODENAME, fname_elem);
		exit(EXIT_FAILURE);
	}
	printf("%s The total number of elements: %d\n", CODENAME, total_num_elems);
	fprintf(fp_elem, "%d %d\n", total_num_elems, L_NUM_NODES);

	fp = BB_std_read_file_search_line(fp, buf, ID_ELEM, BUFFER_SIZE);
	sscanf(buf, "%s %d %s", buf2, &total_num_cells, buf2);
	for(int c=0; c<total_num_cells; c++) {
		int nl;
		fscanf(fp, "%d", &nl);

		int conn[nl];
		for(int l=0; l<nl; l++) {
			int nn;
			fscanf(fp, "%d", &nn);
			if(nl == L_NUM_NODES) {	fprintf(fp_elem, "%d ", nn); }
		}
		if(nl == L_NUM_NODES) {	 fprintf(fp_elem, "\n"); }
	}
	
	fclose(fp_elem);

	/*********************************************************/
	
	fclose(fp);

	printf("\n");

	return 0;
}
