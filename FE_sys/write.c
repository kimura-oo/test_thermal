
#include "write.h"
#include "../libBB/vtk.h"

static const char* CODENAME = "FE_sys/write >";


void BBFE_sys_write_vtk_shape(
		FILE*     fp,
		FE_DATA*  fe,
		const int cell_type)
{
	BB_vtk_write_header(fp);
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
	BB_vtk_write_points_3d(fp, fe->total_num_nodes, fe->x);
	BB_vtk_write_cells(fp, fe->total_num_elems, fe->local_num_nodes, fe->conn);
	BB_vtk_write_cell_types(fp, fe->total_num_elems, cell_type);

}


void BBFE_write_ascii_nodal_vals_scalar(
		FE_DATA*     fe,
		double*      vals,
		const char*  filename)
{
	FILE* fp;
	fp = fopen(filename, "w");
	if( fp == NULL ) {
		printf("%s ERROR: File \"%s\" cannot be opened.\n",
				CODENAME, filename);
	}

	fprintf(fp, "%d\n", fe->total_num_nodes);

	for(int i=0; i<(fe->total_num_nodes); i++) {
		fprintf(fp, "%e\n", vals[i]);
	}

	fclose(fp);

}
