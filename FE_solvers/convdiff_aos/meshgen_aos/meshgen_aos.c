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


	if(read_args_find_option(argc, argv, "--dbc")) {
		printf("%s '--dbc' option is selected.\n", CODENAME);
		// implement Dirichlet BC routine here
		
	}	

	printf("\n");
	return 0;
}
