#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>


typedef struct
{
	int      total_num_nodes;
	double** x;
	
	int   total_num_elems;
	int   local_num_nodes;
	int** conn;
} FE_DATA;


typedef struct
{
	int div_x;   int div_y;   int div_z;
	double l_x;  double l_y;  double l_z;
	double x0;   double y0;   double z0;

	int  order;

	bool dbc_is_imposed;
} CONDITION;


static const char* CODENAME = "meshgen_aos >";
static const char* VOIDNAME = "            >";
static const int BUFFER_SIZE = 10000;


static const char* read_args_return_next_arg(
		int argc,
		char* argv[],
		const char* c_option)
{
	int num = 0;

	for(int i=1; i<argc-1; i++) {
		if(strcmp(argv[i], c_option) == 0 ) {
			num = i;
		}
	}

	if(num == 0) {
		return NULL;
	}
	else {
		return argv[num+1];
	}

}


static bool read_args_find_option(
		int argc,
		char* argv[],
		const char* c_option)
{

	for(int i=1; i<argc; i++) {
		if(strcmp(argv[i], c_option) == 0 ) {
			return true;
		}
	}

	return false;

}


const char* arg_manager(
		int argc,
		char* argv[])
{
	if(argc < 8) {
		printf("%s Please specify parameters for meshing\n", CODENAME);
		printf("%s Format: \n", VOIDNAME);
		printf("%s ./meshgen_aos [1: num. elements (x)] [2: num. elements (y)] [3: num. elements (z)] [4: total length (x)] [5: total length (y)] [6: total length (z)] [7: polynomial order]\n\n", VOIDNAME);
		printf("%s Options: \n", VOIDNAME);
		printf("%s -d [output directory] (Default: ./) \n", VOIDNAME);
		printf("%s --dbc (Output 'D_bc.dat' with zero Dirichlet B.C. for all surfaces) \n\n", VOIDNAME);

		exit(0);
	}

	const char* dir_name;
	dir_name = read_args_return_next_arg(argc, argv, "-d");
	if(dir_name == NULL) {
		return "./";
	}
	else {
		return dir_name;
	}
}


void set_condition(
		int argc,
		char* argv[],
		CONDITION* cond)
{
	cond->div_x = atoi( argv[1] );
	cond->div_y = atoi( argv[2] );
	cond->div_z = atoi( argv[3] );
	cond->l_x   = atof( argv[4] );
	cond->l_y   = atof( argv[5] );
	cond->l_z   = atof( argv[6] );
	cond->order = atoi( argv[7] );

	cond->x0   = 0.0;
	cond->y0   = 0.0;
	cond->z0   = 0.0;
	if(argc >= 10) {
		cond->x0   = atof( argv[8] );
		cond->y0   = atof( argv[9] );
		cond->z0   = atof( argv[10] );
	}

	printf("%s Num. elements (%d %d %d)\n", CODENAME, 
			cond->div_x, cond->div_y, cond->div_z);
	printf("%s Total length  (%e %e %e)\n", CODENAME, 
			cond->l_x, cond->l_y, cond->l_z);
	printf("%s Reference point  (%e, %e, %e)\n", CODENAME, 
			cond->x0, cond->y0, cond->z0);
	printf("%s Polynomial order  %d\n", CODENAME, 
			cond->order);
}


void calc_num_nodes_and_elems(
		FE_DATA* fe,
		const int num_elems_x,
		const int num_elems_y,
		const int num_elems_z,
		const int pol_order)
{
	fe->total_num_nodes = 
		(num_elems_x*pol_order + 1) * 
		(num_elems_y*pol_order + 1) * 
		(num_elems_z*pol_order + 1);
	printf("%s Total num. nodes %d\n", CODENAME, fe->total_num_nodes);


	fe->total_num_elems = num_elems_x*num_elems_y*num_elems_z;
	printf("%s Total num. elems %d\n", CODENAME, fe->total_num_elems);

	fe->local_num_nodes = (pol_order+1)*(pol_order+1)*(pol_order+1);
	printf("%s Local num. nodes %d\n", CODENAME, fe->local_num_nodes);
}


void init_fe_data(
		FE_DATA* fe)
{
	fe->x = (double**)calloc(fe->total_num_nodes, sizeof(double*));
	for(int i=0; i<(fe->total_num_nodes); i++) {
		fe->x[i] = (double*)calloc(3, sizeof(double));
	}

	fe->conn = (int**)calloc(fe->total_num_elems, sizeof(int*));
	for(int e=0; e<(fe->total_num_elems); e++) {
		fe->conn[e] = (int*)calloc(fe->local_num_nodes, sizeof(int));
	}
}


static int set_index_elem(
		CONDITION* cond,
		const int ex,
		const int ey,
		const int ez)
{
	int num_elems_x = cond->div_x;
	int num_elems_y = cond->div_y;
	int num_elems_z = cond->div_z;

	int index = 
		num_elems_x*num_elems_y*ez + 
		num_elems_x*            ey + 
		                        ex;
	
	return index;
}


static int set_index_global_node(
		CONDITION* cond,
		const int ix,
		const int iy,
		const int iz,
		const int order)
{
	int num_nodes_x = order*cond->div_x + 1;
	int num_nodes_y = order*cond->div_y + 1;
	int num_nodes_z = order*cond->div_z + 1;
	
	int index = 
		num_nodes_x*num_nodes_y*iz + 
		num_nodes_x            *iy + 
		                        ix;

	return index;
}


static int set_index_local_node(
		const int ix,
		const int iy,
		const int iz,
		const int order)
{
	int n = order+1;
	int index = n*n*iz + n*iy + ix;

	return index;
}


static int set_index_start_node_num_in_elem(
		CONDITION* cond,
		const int ex,
		const int ey,
		const int ez,
		const int order)
{
	int n = order+1;
	
	int num_nodes_x = order*cond->div_x + 1;
	int num_nodes_y = order*cond->div_y + 1;
	int num_nodes_z = order*cond->div_z + 1;

	int index = 
		num_nodes_x*num_nodes_y*order*ez + 
		num_nodes_x            *order*ey + 
		                        order*ex;

	return index;
}


static int set_index_global_node_corr_local_index(
		CONDITION* cond,
		const int ix,
		const int iy,
		const int iz,
		const int order,
		const int start_num)
{
	int num_nodes_x = order*cond->div_x + 1;
	int num_nodes_y = order*cond->div_y + 1;
	int num_nodes_z = order*cond->div_z + 1;

	int index = start_num + 
		num_nodes_x*num_nodes_y*iz + 
		num_nodes_x            *iy + 
		                        ix;

	return index;
}


void set_nodes(
		FE_DATA* fe,
		CONDITION* cond)
{
	int order = cond->order;
	int num_nodes_x = order*cond->div_x + 1;
	int num_nodes_y = order*cond->div_y + 1;
	int num_nodes_z = order*cond->div_z + 1;

	double dx = cond->l_x/(order*cond->div_x);
	double dy = cond->l_y/(order*cond->div_y);
	double dz = cond->l_z/(order*cond->div_z);

	for(int k=0; k<num_nodes_z; k++) {
		for(int j=0; j<num_nodes_y; j++) {
			for(int i=0; i<num_nodes_x; i++) {
				int n = set_index_global_node(cond, i, j, k, order);

				fe->x[n][0] = cond->x0 + dx*i;
				fe->x[n][1] = cond->y0 + dy*j;
				fe->x[n][2] = cond->z0 + dz*k;
					
			}
		}
	}
}


void set_elems(
		FE_DATA* fe,
		CONDITION* cond)
{
	int order = cond->order;
	int num_elems_x = cond->div_x;
	int num_elems_y = cond->div_y;
	int num_elems_z = cond->div_z;

	int nl = order+1;

	for(int ek=0; ek<num_elems_z; ek++) {
		for(int ej=0; ej<num_elems_y; ej++) {
			for(int ei=0; ei<num_elems_x; ei++) {
				int elem_num  = set_index_elem(cond, ei, ej, ek);
				int start_num = set_index_start_node_num_in_elem(
						cond, ei, ej, ek, order);
				
				for(int k=0; k<nl; k++) {
					for(int j=0; j<nl; j++) {
						for(int i=0; i<nl; i++) {
							int n  = set_index_local_node(i, j, k, order);
							int gn = set_index_global_node_corr_local_index(
									cond, i, j, k, order, start_num);
							
							fe->conn[elem_num][n] = gn;
						}
					}
				}

			}
		}
	}
}


void output_data(
		FE_DATA* fe,
		const char* filename_node,
		const char* filename_elem,
		const char* directory)
{
	char fname_node[BUFFER_SIZE];
	snprintf(fname_node, BUFFER_SIZE, "%s/%s", directory, filename_node);
	FILE* fp_node;
	fp_node = fopen(fname_node, "w");
	if( fp_node == NULL ) {
		printf("%s ERROR: File \"%s\" cannot be opened.\n", 
				CODENAME, fname_node);
	}

	fprintf(fp_node, "%d\n", fe->total_num_nodes);
	for(int i=0; i<fe->total_num_nodes; i++) {
		fprintf(fp_node, "%e %e %e\n", fe->x[i][0], fe->x[i][1], fe->x[i][2]);
	}

	char fname_elem[BUFFER_SIZE];
	snprintf(fname_elem, BUFFER_SIZE, "%s/%s", directory, filename_elem);
	FILE* fp_elem;
	fp_elem = fopen(fname_elem, "w");
	if( fp_elem == NULL ) {
		printf("%s ERROR: File \"%s\" cannot be opened.\n", 
				CODENAME, fname_elem);
	}

	fprintf(fp_elem, "%d %d\n", fe->total_num_elems, fe->local_num_nodes);
	for(int i=0; i<fe->total_num_elems; i++) {
		for(int j=0; j<fe->local_num_nodes; j++) {
			fprintf(fp_elem, "%d ", fe->conn[i][j]);
		}
		fprintf(fp_elem, "\n");
	}

	fclose(fp_node);
	fclose(fp_elem);
}


static bool node_is_inside_plane(
		const double 	x,
		const double 	val,
		const double 	epsilon) 
{
	double min_val = val - epsilon;
	double max_val = val + epsilon;

	if ( min_val < x && x < max_val ) {
		return true;
	}
	else {
		return false;
	}
}


void detect_and_output_Dirichlet_bc(
		FE_DATA* fe,
		CONDITION* cond,
		const char* filename_dbc,
		const char* directory)
{
	double x_min = cond->x0;
	double y_min = cond->y0;
	double z_min = cond->z0;
	double x_max = cond->x0 + cond->l_x;
	double y_max = cond->y0 + cond->l_y;
	double z_max = cond->z0 + cond->l_z;

	double ratio = 100.0;
	double epsilon_x =  cond->l_x/(ratio*cond->order*cond->div_x);
	double epsilon_y =  cond->l_y/(ratio*cond->order*cond->div_y);
	double epsilon_z =  cond->l_z/(ratio*cond->order*cond->div_z);
	
	// Count the num. of Dirichlet B.C.
	int num_dbc = 0;
	for(int i=0; i<(fe->total_num_nodes); i++) {
		bool dbc_flag = false;

		if( node_is_inside_plane(fe->x[i][0], x_min, epsilon_x) ) { dbc_flag = true; }
		if( node_is_inside_plane(fe->x[i][0], x_max, epsilon_x) ) { dbc_flag = true; }
		if( node_is_inside_plane(fe->x[i][1], y_min, epsilon_y) ) { dbc_flag = true; }
		if( node_is_inside_plane(fe->x[i][1], y_max, epsilon_y) ) { dbc_flag = true; }
		if( node_is_inside_plane(fe->x[i][2], z_min, epsilon_z) ) { dbc_flag = true; }
		if( node_is_inside_plane(fe->x[i][2], z_max, epsilon_z) ) { dbc_flag = true; }
		
		if( dbc_flag ) { 
			num_dbc++;
		}
	}

	// Outout Dirichlet B.C. file
	char fname_dbc[BUFFER_SIZE];
	snprintf(fname_dbc, BUFFER_SIZE, "%s/%s", directory, filename_dbc);
	FILE* fp_dbc;
	fp_dbc = fopen(fname_dbc, "w");
	if( fp_dbc == NULL ) {
		printf("%s ERROR: File \"%s\" cannot be opened.\n", 
				CODENAME, fname_dbc);
	}

	int block_length = 1;
	fprintf(fp_dbc, "%d %d\n", num_dbc, block_length);
	for(int i=0; i<(fe->total_num_nodes); i++) {
		bool dbc_flag = false;

		if( node_is_inside_plane(fe->x[i][0], x_min, epsilon_x) ) { dbc_flag = true; }
		if( node_is_inside_plane(fe->x[i][0], x_max, epsilon_x) ) { dbc_flag = true; }
		if( node_is_inside_plane(fe->x[i][1], y_min, epsilon_y) ) { dbc_flag = true; }
		if( node_is_inside_plane(fe->x[i][1], y_max, epsilon_y) ) { dbc_flag = true; }
		if( node_is_inside_plane(fe->x[i][2], z_min, epsilon_z) ) { dbc_flag = true; }
		if( node_is_inside_plane(fe->x[i][2], z_max, epsilon_z) ) { dbc_flag = true; }
		
		if( dbc_flag ) { 
			fprintf(fp_dbc, "%d 0 0.0\n", i);
		}
	}

	fclose(fp_dbc);
}


int main(
		int 	argc,
		char* 	argv[])
{
	printf("\n");

	const char* dir_name = arg_manager(argc, argv);
	printf("%s Output directory: %s\n", CODENAME, dir_name);

	FE_DATA   fe;
	CONDITION cond;

	set_condition(argc, argv, &cond);

	calc_num_nodes_and_elems(
			&fe, cond.div_x, cond.div_y, cond.div_z, cond.order);
	init_fe_data(&fe);

	set_nodes(&fe, &cond);
	set_elems(&fe, &cond);

	output_data(&fe, "node.dat", "elem.dat", dir_name);

	if(read_args_find_option(argc, argv, "--dbc")) {
		printf("%s '--dbc' option is selected.\n", CODENAME);
		detect_and_output_Dirichlet_bc(&fe, &cond, "D_bc.dat", dir_name);
	}	

	printf("\n");
	return 0;
}
