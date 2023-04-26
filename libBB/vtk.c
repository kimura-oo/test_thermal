
#include <stdio.h>

void BB_vtk_void(){}


void BB_vtk_write_header(
		FILE* fp)
{
	fprintf(fp, "# vtk DataFile Version 3.0\n");
	fprintf(fp, "vtk output\n");
	fprintf(fp, "ASCII\n");
}


void BB_vtk_write_points_3d(
		FILE*    fp,
		int      num_points,
		double** x) //[num_points][3]
{
	fprintf(fp, "POINTS %d float\n", num_points);

	for(int i=0; i<num_points; i++) {
		fprintf(fp, "%e %e %e\n", x[i][0], x[i][1], x[i][2]);
	}
}


void BB_vtk_write_points_3d_with_disp(
		FILE*    fp,
		int      num_points,
		double** x,
		double** u,     //[num_points][3] displacement
		double   scale) 

{
	fprintf(fp, "POINTS %d float\n", num_points);

	for(int i=0; i<num_points; i++) {
		fprintf(fp, "%e %e %e\n", 
				x[i][0] + scale*u[i][0], 
				x[i][1] + scale*u[i][1], 
				x[i][2] + scale*u[i][2]);
	}
}


void BB_vtk_write_cells(
		FILE* fp,
		int   num_cells,
		int   num_points_in_cell,
		int** connectivity) //[num_cells][num_points_in_cell]
{
	fprintf(fp, "CELLS %d %d\n",
			num_cells, num_cells*(num_points_in_cell + 1));
	for(int e=0; e<num_cells; e++) {
		fprintf(fp, "%d ", num_points_in_cell);

		for(int i=0; i<num_points_in_cell; i++) {
			fprintf(fp, "%d ", connectivity[e][i]);
		}
		fprintf(fp, "\n");
	}
}


void BB_vtk_write_cell_types(
		FILE* fp,
		int   num_cells,
		int   elem_type)
{
	fprintf(fp, "CELL_TYPES %d\n", num_cells);

	for(int e=0; e<num_cells; e++) {
		fprintf(fp, "%d\n", elem_type);
	}
}


void BB_vtk_write_point_vals_scalar(
		FILE*        fp,
		double*      val,
		const int    num_points,
		const char*  label)
{
	fprintf(fp, "SCALARS %s float\n", label);
	fprintf(fp, "LOOKUP_TABLE default\n");

	for(int i=0; i<num_points; i++) {
		fprintf(fp, "%e\n", val[i]);
	}
}


void BB_vtk_write_point_vals_vector(
		FILE*        fp,
		double**     val,
		const int    num_points,
		const char*  label)
{
	fprintf(fp, "VECTORS %s float\n", label);

	for(int i=0; i<num_points; i++) {
		fprintf(fp, "%e %e %e\n", val[i][0], val[i][1], val[i][2]);
	}
}

