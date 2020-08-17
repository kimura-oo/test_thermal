
#include <BBFE/sys/FE_dataset.h>

static const char* INPUT_FILENAME_NODE = "node.dat";
static const char* INPUT_FILENAME_ELEM = "elem.dat";
static const char* INPUT_FILENAME_D_BC = "D_bc.dat";
static const char* INPUT_FILENAME_N_BC = "N_bc.dat";

static const char* OUTPUT_FILENAME_SURF_VTK = "surf.vtk";

typedef struct {
	const char* directory;
	int         block_size;

} SETTINGS;


typedef struct {
	int    num_bc_nodes;
	bool*  node_is_on_surface;

	int    num_bc_surfs;
	int    num_nodes_on_surf;
	int*   orig_elem_num;
	int*   orig_surf_num;
	int**  conn_surf;
	bool** surf_is_on_surface;
	int    num_surfs_in_elem;

} SURFACE;


void read_fe_data(
		BBFE_DATA*  fe,
		const char* directory);

void memory_allocation_surface(
		SURFACE*    surf,
		const int   total_num_nodes,
		const int   total_num_elems,
		const int   local_num_nodes,
		const char* codename);

void get_surface_nodes(
		BBFE_DATA*  fe,
		SURFACE*    surf,
		const char* codename);

void get_surface_info(
		BBFE_DATA*  fe,
		SURFACE*    surf,
		const char* codename);

void memory_allocation_surface_conn(
		SURFACE* surf);

void set_surface_conn(
		BBFE_DATA* fe,
		SURFACE*   surf);

void write_surface_vtk(
		BBFE_DATA*  fe,
		SURFACE*    surf,
		const char* directory);
