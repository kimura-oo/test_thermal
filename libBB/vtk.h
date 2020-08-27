#pragma once


#define TYPE_VTK_VERTEX               1
#define TYPE_VTK_POLY_VERTEX          2
#define TYPE_VTK_LINE                 3
#define TYPE_VTK_POLY_LINE            4

#define TYPE_VTK_TRIANGLE             5
#define TYPE_VTK_TRIANGLE_STRIP       6
#define TYPE_VTK_POLYGON              7
#define TYPE_VTK_PIXEL                8
#define TYPE_VTK_QUAD                 9
                                      
#define TYPE_VTK_TETRA                10
#define TYPE_VTK_VOXEL                11
#define TYPE_VTK_HEXAHEDRON           12
#define TYPE_VTK_WEDGE                13
#define TYPE_VTK_PYRAMID              14
                                      
#define TYPE_VTK_QUADRATIC_EDGE       21
#define TYPE_VTK_QUADRATIC_TRIANGLE   22
#define TYPE_VTK_QUADRATIC_QUAD       23
#define TYPE_VTK_QUADRATIC_TETRA      24
#define TYPE_VTK_QUADRATIC_HEXAHEDRON 25



void BB_vtk_void();

void BB_vtk_write_header(
		FILE* fp);

void BB_vtk_write_points_3d(
		FILE*    fp,
		int      num_points,
		double** x); //[num_points][3]

void BB_vtk_write_points_3d_with_disp(
		FILE*    fp,
		int      num_points,
		double** x,
		double** u,     //[num_points][3] displacement
		double   scale);

void BB_vtk_write_cells(
		FILE* fp,
		int   num_cells,
		int   num_points_in_cell,
		int** connectivity); //[num_cells][num_points_in_cell]

void BB_vtk_write_cell_types(
		FILE* fp,
		int   num_cells,
		int   elem_type);

void BB_vtk_write_point_vals_scalar(
		FILE*        fp,
		double*      val,
		const int    num_points,
		const char*  label);

void BB_vtk_write_point_vals_vector(
		FILE*        fp,
		double**     val,
		const int    num_points,
		const char*  label);
