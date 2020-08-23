#include "BBFE/sys/FE_dataset.h"

static const char* DEF_DIRECTORY     = ".";
static const char* DEF_FILENAME_NODE = "node.dat";
static const char* DEF_FILENAME_ELEM = "elem.dat";
static const char* DEF_FILENAME_SURF = "surf.dat";
static const char* DEF_FILENAME_VTK  = "elem.vtk";
static const char* DEF_FILENAME_SVTK = "surf.vtk";

static const char* OPTION_INFILE_NODE = "-in";
static const char* OPTION_INFILE_ELEM = "-ie";
static const char* OPTION_OUTFILE_NODE = "-on";
static const char* OPTION_OUTFILE_ELEM = "-oe";
static const char* OPTION_OUTFILE_VTK  = "-ov";

static const char* OPTION_DIRECTORY = "-d";

typedef struct {
	double x_min[3];
	double x_max[3];

	const char* infile_node;
	const char* infile_elem;
	const char* outfile_node;
	const char* outfile_elem;
	const char* outfile_vtk;
	const char* directory;
} SETTINGS;

void args_manager(
		SETTINGS*   set,
		int         argc,
		char*       argv[],
		bool        surface_mesh);

void open_fe_files(
		BBFE_DATA*  fe,
		const char* filename_node,
		const char* filename_elem,
		const char* directory);

bool point_is_inside_3d(
		double       x[3],
		const double x_min[3],
		const double x_max[3]);

bool point_cloud_is_inside_3d(
		double**     x,
		const int    num_points,
		const double x_min[3],
		const double x_max[3]);

bool element_is_inside_3d(
		BBFE_DATA*   fe,
		const int    elem_num,
		const double x_min[3],
		const double x_max[3]);

void extract_elements(
		BBFE_DATA* fe_ext,
		BBFE_DATA* fe_orig,
		const double x_min[3],
		const double x_max[3]);

void remove_elements(
		BBFE_DATA* fe_remove,
		BBFE_DATA* fe_orig,
		const double x_min[3],
		const double x_max[3]);

void remove_floating_nodes(
		BBFE_DATA* fe_remove,
		BBFE_DATA* fe_orig);

void allocate_and_copy_node_data(
		BBFE_DATA* fe_copied,
		BBFE_DATA* fe_orig);

void write_node_file(
		BBFE_DATA*  fe,
		const char* filename,
		const char* directory);

void write_element_file(
		BBFE_DATA*  fe,
		const char* filename,
		const char* directory);

int get_cell_type_vtk_3d(
		int num_points_in_cell);

int get_cell_type_vtk_2d(
		int num_points_in_cell);

void write_vtk_shape(
		BBFE_DATA*  fe,
		const char* filename,
		const char* directory,
		const int   cell_type);
