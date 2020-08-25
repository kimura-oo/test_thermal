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

} CONDITION;


static const char* CODENAME = "meshgen_tet >";
static const char* VOIDNAME = "            >";
static const int BUFFER_SIZE = 10000;


void init_fe_data(
		FE_DATA* fe,
		const int total_num_nodes,
		const int total_num_elems,
		const int local_num_nodes)
{
	fe->total_num_nodes = total_num_nodes;
	fe->total_num_elems = total_num_elems;
	fe->local_num_nodes = local_num_nodes;

	fe->x = (double**)calloc(fe->total_num_nodes, sizeof(double*));
	for(int i=0; i<(fe->total_num_nodes); i++) {
		fe->x[i] = (double*)calloc(3, sizeof(double));
	}

	fe->conn = (int**)calloc(fe->total_num_elems, sizeof(int*));
	for(int e=0; e<(fe->total_num_elems); e++) {
		fe->conn[e] = (int*)calloc(fe->local_num_nodes, sizeof(int));
	}
}


void set_nodes(
		FE_DATA* fe,
		CONDITION* cond)
{
	int n_xy = (cond->div_x+1)*(cond->div_y+1);
	int n_x  = cond->div_x+1;

	double dx = cond->l_x/cond->div_x;
	double dy = cond->l_y/cond->div_y;
	double dz = cond->l_z/cond->div_z;


	for(int k=0; k<(cond->div_z)+1; k++) {
		for(int j=0; j<(cond->div_y)+1; j++) {
			for(int i=0; i<(cond->div_x)+1; i++) {
				int n = n_xy*k + n_x*j + i;
				fe->x[n][0] = cond->x0 + dx*i;
				fe->x[n][1] = cond->y0 + dy*j;
				fe->x[n][2] = cond->z0 + dz*k;
					
			}
		}
	}
}


void set_hex_elems(
		FE_DATA* fe,
		CONDITION* cond)
{
	int e_xy = (cond->div_x)*(cond->div_y);
	int e_x  = cond->div_x;
	int n_xy = (cond->div_x+1)*(cond->div_y+1);
	int n_x  = cond->div_x+1;

	int counter = 0;
	for(int k=0; k<(cond->div_z); k++) {
		for(int j=0; j<(cond->div_y); j++) {
			for(int i=0; i<(cond->div_x); i++) {
				int e = e_xy*k + e_x*j + i;
				int n = e_xy*k + e_x*j + i + counter;
				fe->conn[e][0] = n;
				fe->conn[e][1] = n + 1;
				fe->conn[e][2] = n + n_x + 1;
				fe->conn[e][3] = n + n_x;
				fe->conn[e][4] = n_xy + n;
				fe->conn[e][5] = n_xy + n + 1;
				fe->conn[e][6] = n_xy + n + n_x + 1;
				fe->conn[e][7] = n_xy + n + n_x;

			}
			counter++;
		}
		counter += n_x;
	}
}


void set_tet(
		FE_DATA* fe_tet,
		FE_DATA* fe_hex)
{
	for(int i=0; i<fe_hex->total_num_nodes; i++) {
		fe_tet->x[i][0] = fe_hex->x[i][0];
		fe_tet->x[i][1] = fe_hex->x[i][1];
		fe_tet->x[i][2] = fe_hex->x[i][2];
	}
	
	for(int e=0; e<fe_hex->total_num_elems; e++) {
		int n[8];
		for(int i=0; i<8; i++) {
			n[i] = fe_hex->conn[e][i];
		}

		fe_tet->conn[6*e  ][0] = n[0];  fe_tet->conn[6*e  ][1] = n[4];  fe_tet->conn[6*e  ][2] = n[1];  fe_tet->conn[6*e  ][3] = n[7];
		fe_tet->conn[6*e+1][0] = n[0];  fe_tet->conn[6*e+1][1] = n[1];  fe_tet->conn[6*e+1][2] = n[3];  fe_tet->conn[6*e+1][3] = n[7];
		fe_tet->conn[6*e+2][0] = n[4];  fe_tet->conn[6*e+2][1] = n[7];  fe_tet->conn[6*e+2][2] = n[5];  fe_tet->conn[6*e+2][3] = n[1];
		fe_tet->conn[6*e+3][0] = n[5];  fe_tet->conn[6*e+3][1] = n[7];  fe_tet->conn[6*e+3][2] = n[6];  fe_tet->conn[6*e+3][3] = n[1];
		fe_tet->conn[6*e+4][0] = n[2];  fe_tet->conn[6*e+4][1] = n[3];  fe_tet->conn[6*e+4][2] = n[1];  fe_tet->conn[6*e+4][3] = n[7];
		fe_tet->conn[6*e+5][0] = n[2];  fe_tet->conn[6*e+5][1] = n[1];  fe_tet->conn[6*e+5][2] = n[6];  fe_tet->conn[6*e+5][3] = n[7];
	}
}


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



const char* arg_manager(
		int argc,
		char* argv[])
{
	if(argc < 7) {
		printf("%s Please specify parameters for meshing\n", CODENAME);
		printf("%s Format: \n", VOIDNAME);
		printf("%s ./meshgen_tet [1: num. elements (x)] [2: num. elements (y)] [3: num. elements (z)] [4: total length (x)] [5: total length (y)] [6: total length (z)] \n\n", VOIDNAME);
		printf("%s Options: \n", VOIDNAME);
		printf("%s -d [output directory] (Default: ./) \n\n", VOIDNAME);

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

	cond->x0   = 0.0;
	cond->y0   = 0.0;
	cond->z0   = 0.0;
	if(argc >= 10) {
		cond->x0   = atof( argv[7] );
		cond->y0   = atof( argv[8] );
		cond->z0   = atof( argv[9] );
	}

	printf("%s Num. elements (%d %d %d)\n", CODENAME, 
			cond->div_x, cond->div_y, cond->div_z);
	printf("%s Total length  (%e %e %e)\n", CODENAME, 
			cond->l_x, cond->l_y, cond->l_z);
	printf("%s Reference point  (%e, %e, %e)\n", CODENAME, 
			cond->x0, cond->y0, cond->z0);
}


void print_data(
		FE_DATA* fe)
{
	for(int i=0; i<fe->total_num_nodes; i++) {
		printf("%e %e %e\n", fe->x[i][0], fe->x[i][1], fe->x[i][2]);
	}

	for(int i=0; i<fe->total_num_elems; i++) {
		printf("%d: ", i);
		for(int j=0; j<fe->local_num_nodes; j++) {
			printf("%d ", fe->conn[i][j]);
		}
		printf("\n");
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



int main(
		int   argc,
		char* argv[])
{
	printf("\n");

	const char* dir_name = arg_manager(argc, argv);
	printf("%s Output directory: %s\n", CODENAME, dir_name);

	FE_DATA fe_hex;
	FE_DATA fe_tet;

	CONDITION cond;

	set_condition(argc, argv, &cond);

	init_fe_data(
			&fe_hex,
			(cond.div_x+1)*(cond.div_y+1)*(cond.div_z+1),
			cond.div_x*cond.div_y*cond.div_z,
			8);


	set_nodes(&fe_hex, &cond);
	set_hex_elems(&fe_hex, &cond);

	//print_data(&fe_tet);
	output_data(&fe_hex, "node.dat", "elem.dat", dir_name);

	printf("\n");
	return 0;
	
}
