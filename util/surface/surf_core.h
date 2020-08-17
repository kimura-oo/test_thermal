
#include <BBFE/sys/FE_dataset.h>

static const char* INPUT_FILENAME_NODE = "node.dat";
static const char* INPUT_FILENAME_ELEM = "elem.dat";
static const char* INPUT_FILENAME_D_BC = "D_bc.dat";
static const char* INPUT_FILENAME_N_BC = "N_bc.dat";


typedef struct {
	const char* directory;
	int         block_size;

} SETTINGS;


typedef struct {
	int    num_bc_nodes;
	bool*  node_is_on_surface;

	int    num_bc_surfs;
	int**  conn_surf;
	bool** surf_is_on_surface;

} SURFACE;


void read_fe_data(
		BBFE_DATA*  fe,
		const char* directory);

void get_surface_nodes(
		BBFE_DATA*  fe,
		SURFACE*    surf,
		const char* codename);

void get_surface_info(
		BBFE_DATA*  fe,
		SURFACE*    surf,
		const char* codename);
