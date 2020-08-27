
#include "read.h"
#include "memory.h"
#include "BB/std.h"

#include <stdlib.h>

static const char* CODENAME = "FE_sys/read >";
static const int BUFFER_SIZE = 10000;


FILE* BBFE_sys_read_fopen(
		FILE*        fp,
		const char*  filename,
		const char*  directory)
{
	char fname[BUFFER_SIZE];
	snprintf(fname, BUFFER_SIZE, "%s/%s", directory, filename);

	fp = fopen(fname, "r");
	if( fp == NULL ) {
		printf("%s ERROR: File \"%s\" cannot be opened.\n",
				CODENAME, fname);
		exit(EXIT_FAILURE);
	}
	else {
		printf("%s Reading file \"%s\".\n", CODENAME, fname);
	}

	return fp;
}


FILE* BBFE_sys_read_fopen_without_error(
		FILE*        fp,
		const char*  filename,
		const char*  directory)
{
	char fname[BUFFER_SIZE];
	snprintf(fname, BUFFER_SIZE, "%s/%s", directory, filename);
	fp = fopen(fname, "r");
	printf("%s Reading file \"%s\".\n", CODENAME, fname);

	return fp;
}


void BBFE_sys_read_node(
		BBFE_DATA*   fe,
		const char*  filename,
		const char*  directory)
{
	FILE* fp;
	fp = BBFE_sys_read_fopen(fp, filename, directory);

	// read the number of nodes
	BB_std_scan_line(
			&fp, BUFFER_SIZE, "%d", &(fe->total_num_nodes));
	printf("%s Num. nodes: %d\n", CODENAME, fe->total_num_nodes);
	BBFE_sys_memory_allocation_node(fe, 3);

	// read positions of nodes
	for(int i=0; i<(fe->total_num_nodes); i++) {
		BB_std_scan_line(&fp, BUFFER_SIZE,
				"%lf %lf %lf", &(fe->x[i][0]), &(fe->x[i][1]), &(fe->x[i][2]));
	}

	fclose(fp);
}


void BBFE_sys_read_elem(
		BBFE_DATA*   fe,
		const char*  filename,
		const char*  directory,
		int          num_integ_points)
{
	FILE* fp;
	fp = BBFE_sys_read_fopen(fp, filename, directory);

	// read the number of elements
	BB_std_scan_line(
			&fp, BUFFER_SIZE, "%d %d",&(fe->total_num_elems), &(fe->local_num_nodes));
	printf("%s Num. elements: %d\n", CODENAME, fe->total_num_elems);
	BBFE_sys_memory_allocation_elem(fe, num_integ_points, 3);

	// read the connectivities of elements
	for(int e=0; e<(fe->total_num_elems); e++) {
		for(int i=0; i<(fe->local_num_nodes); i++) {
			fscanf(fp, "%d", &(fe->conn[e][i]));
		}
	}

	fclose(fp);
}


void BBFE_sys_read_Dirichlet_bc(
		BBFE_BC*     bc,
		const char*  filename,
		const char*  directory,
		const int    total_num_nodes,
		const int    block_size)
{

	bc->total_num_nodes = total_num_nodes;
	bc->block_size      = block_size;

	BBFE_sys_memory_allocation_Dirichlet_bc(bc, total_num_nodes, bc->block_size);
	int n = total_num_nodes * bc->block_size;

	for(int i=0; i<n; i++) {
		bc->D_bc_exists[i]   = false;
		bc->imposed_D_val[i] = 0.0;
	}

	FILE* fp;
	fp = BBFE_sys_read_fopen_without_error(fp, filename, directory);
	if( fp == NULL ) {
		printf("%s WARNING: Dirichlet B.C. file, \"%s\", is not found.\n",
				CODENAME, filename);
		return;
	}


	BB_std_scan_line(&fp, BUFFER_SIZE,
			"%d %d", &(bc->num_D_bcs), &(bc->block_size));
	printf("%s Num. Dirichlet B.C.: %d\n", CODENAME, bc->num_D_bcs);

	for(int i=0; i<(bc->num_D_bcs); i++) {
		int node_id;  int block_id;  double val;
		BB_std_scan_line(&fp, BUFFER_SIZE,
				"%d %d %lf", &node_id, &block_id, &val);

		int index = (bc->block_size)*node_id + block_id;
		bc->D_bc_exists[ index ]   = true;
		bc->imposed_D_val[ index ] = val;
	}
	
	fclose(fp);
}


void BBFE_sys_read_Neumann_bc(
		BBFE_BC*     bc,
		const char*  filename,
		const char*  directory,
		const int    total_num_nodes,
		const int    block_size)
{

	bc->total_num_nodes = total_num_nodes;
	bc->block_size      = block_size;

	BBFE_sys_memory_allocation_Neumann_bc(bc, total_num_nodes, bc->block_size);
	int n = total_num_nodes * bc->block_size;

	for(int i=0; i<n; i++) {
		bc->N_bc_exists[i]   = false;
		bc->imposed_N_val[i] = 0.0;
	}

	FILE* fp;
	fp = BBFE_sys_read_fopen_without_error(fp, filename, directory);
	if( fp == NULL ) {
		return;
	}


	BB_std_scan_line(&fp, BUFFER_SIZE,
			"%d %d", &(bc->num_N_bcs), &(bc->block_size));
	printf("%s Num. Neumann B.C.: %d\n", CODENAME, bc->num_N_bcs);

	for(int i=0; i<(bc->num_N_bcs); i++) {
		int node_id;  int block_id;  double val;
		BB_std_scan_line(&fp, BUFFER_SIZE,
				"%d %d %lf", &node_id, &block_id, &val);

		int index = (bc->block_size)*node_id + block_id;
		bc->N_bc_exists[ index ]    = true;
		bc->imposed_N_val[ index ] += val;
	}

	fclose(fp);
}

