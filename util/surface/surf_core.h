
#include <BBFE/sys/FE_dataset.h>

static const char* FILENAME_NODE = "node.dat";
static const char* FILENAME_ELEM = "elem.dat";
static const char* FILENAME_SURF = "surf.dat";

static const char* OPTION_DIRECTORY    = "-d";
static const char* OPTION_INFILE_NODE = "-in";
static const char* OPTION_INFILE_SURF = "-ie";
static const char* OPTION_OUTFILE     = "-o";

static const char* DEF_DIRECTORY   = ".";

static const char* OUTPUT_FILENAME_SURF_VTK = "surf.vtk";

static int BUFFER_SIZE = 10000;


typedef struct {
	const char* directory;

	const char* infile_node;
	const char* infile_elem;
	const char* outfile_bc;

	int         block_size;

	double*     bc_value;

} SETTINGS;


typedef struct {
	// node related data
	int    num_bc_nodes;
	bool*  node_is_on_surface;

	// surface related data (basic)
	int    num_bc_surfs;
	int    num_nodes_on_surf;
	int**  conn_surf;
	// surface related data (additional)
	int*   orig_elem_num;
	int*   orig_surf_num;
	bool** surf_is_on_surface;
	int    num_surfs_in_elem;

} SURFACE;

void cmd_args_reader_bc(
		SETTINGS* set,
		int       argc,
		char*     argv[],
		const char* codename,
		const char* voidname,
		const char* def_filename_bc);

void read_fe_data(
		BBFE_DATA*  fe,
		const char* directory,
		const char* infile_node,
		const char* infile_elem);


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

int get_bc_node_list(
		bool* node_has_bc,
		BBFE_DATA* fe);

void write_bc_file_const(
		BBFE_DATA*  fe,
		bool*       node_has_bc,
		int         num_bc_nodes,
		int         block_size,
		double*     val,
		const char* filename,
		const char* directory);

void write_bc_file_nonconst(
		BBFE_DATA*  fe,
		bool*       node_has_bc,
		int         num_bc_nodes,
		int         block_size,
		double**    val,
		const char* filename,
		const char* directory);
