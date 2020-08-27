
#include "write.h"
#include "BB/vtk.h"

#include <stdlib.h>

static const char* CODENAME = "FE_sys/write >";
static const int BUFFER_SIZE = 10000;


FILE* BBFE_sys_write_fopen(
		FILE*        fp,
		const char*  filename,
		const char*  directory)
{
	char fname[BUFFER_SIZE];
	snprintf(fname, BUFFER_SIZE, "%s/%s", directory, filename);

	fp = fopen(fname, "w");
	if( fp == NULL ) {
		printf("%s ERROR: File \"%s\" cannot be opened.\n",
				CODENAME, fname);
		exit(EXIT_FAILURE);
	}
	else {
		printf("%s Writing file \"%s\".\n", CODENAME, fname);
	}

	return fp;
}


FILE* BBFE_sys_write_fopen_without_error(
		FILE*        fp,
		const char*  filename,
		const char*  directory)
{
	char fname[BUFFER_SIZE];
	snprintf(fname, BUFFER_SIZE, "%s/%s", directory, filename);

	fp = fopen(fname, "w");
	printf("%s Writing file \"%s\".\n", CODENAME, fname);

	return fp;
}


FILE* BBFE_sys_write_add_fopen(
		FILE*        fp,
		const char*  filename,
		const char*  directory)
{
	char fname[BUFFER_SIZE];
	snprintf(fname, BUFFER_SIZE, "%s/%s", directory, filename);

	fp = fopen(fname, "a");
	if( fp == NULL ) {
		printf("%s ERROR: File \"%s\" cannot be opened.\n",
				CODENAME, fname);
		exit(EXIT_FAILURE);
	}
	else {
		printf("%s Writing (adding) file \"%s\".\n", CODENAME, fname);
	}

	return fp;
}


void BBFE_sys_write_vtk_shape(
		FILE*       fp,
		BBFE_DATA*  fe,
		const int cell_type)
{
	BB_vtk_write_header(fp);
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
	BB_vtk_write_points_3d(fp, fe->total_num_nodes, fe->x);
	BB_vtk_write_cells(fp, fe->total_num_elems, fe->local_num_nodes, fe->conn);
	BB_vtk_write_cell_types(fp, fe->total_num_elems, cell_type);

}


void BBFE_sys_write_vtk_shape_with_disp(
		FILE*       fp,
		BBFE_DATA*  fe,
		const int cell_type,
		double**    u,  // displacement
		double      scale)
{
	BB_vtk_write_header(fp);
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
	BB_vtk_write_points_3d_with_disp(fp, fe->total_num_nodes, fe->x, u, scale);
	BB_vtk_write_cells(fp, fe->total_num_elems, fe->local_num_nodes, fe->conn);
	BB_vtk_write_cell_types(fp, fe->total_num_elems, cell_type);

}


void BBFE_write_ascii_nodal_vals_scalar(
		BBFE_DATA*   fe,
		double*      vals,
		const char*  filename,
		const char*  directory)
{
	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, filename, directory);

	fprintf(fp, "%d\n", fe->total_num_nodes);

	for(int i=0; i<(fe->total_num_nodes); i++) {
		fprintf(fp, "%e\n", vals[i]);
	}

	fclose(fp);

}
