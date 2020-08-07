
#include "read.h"
#include "memory.h"
#include "../libBB/std.h"

static const char* CODENAME = "FE_sys/read >";
static const int BUFFER_SIZE = 10000;


void BBFE_sys_read_node(
		FE_DATA*     fe,
		const char*  filename)
{
	FILE* fp;

	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s ERROR: File \"%s\" cannot be opened.\n",
				CODENAME, filename);
	}

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
		FE_DATA*     fe,
		const char*  filename,
		int          num_integ_points)
{
	FILE* fp;

	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s ERROR: File \"%s\" cannot be opened.\n",
				CODENAME, filename);
	}

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
		BC_DATA*     bc,
		const char*  filename,
		const int    total_num_nodes)
{
	FILE* fp;
	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s ERROR: File \"%s\" cannot be opened.\n",
				CODENAME, filename);
	}

	bc->total_num_nodes = total_num_nodes;

	BB_std_scan_line(&fp, BUFFER_SIZE,
			"%d %d", &(bc->num_D_bcs), &(bc->block_size));
	printf("%s Num. Dirichlet B.C.: %d\n", CODENAME, bc->num_D_bcs);

	int n = total_num_nodes * bc->block_size;

	bc->D_bc_exists   = BB_std_calloc_1d_bool(  bc->D_bc_exists  , n);
	bc->imposed_D_val = BB_std_calloc_1d_double(bc->imposed_D_val, n);
	for(int i=0; i<n; i++) {
		bc->D_bc_exists[i]   = false;
		bc->imposed_D_val[i] = 0.0;
	}

	for(int i=0; i<(bc->num_D_bcs); i++) {
		int node_id;  int block_id;  double val;
		BB_std_scan_line(&fp, BUFFER_SIZE,
				"%d %d %lf", &node_id, &block_id, &val);

		int index = (bc->block_size)*node_id + block_id;
		bc->D_bc_exists[ index ]   = true;
		bc->imposed_D_val[ index ] = val;
	}
}
