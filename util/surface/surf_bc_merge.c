#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <BB/std.h>
#include <BB/calc.h>
#include <BB/vtk.h>
#include <BBFE/std/shapefunc.h>
#include <BBFE/std/surface.h>
#include <BBFE/sys/FE_dataset.h>
#include <BBFE/sys/read.h>
#include <BBFE/sys/write.h>

#include "surf_core.h"

static const char* CODENAME            = "surf_bc_merge >";
static const char* VOIDNAME            = "               ";

static const char* FILENAME_D_BC = "D_bc.dat";


void cmd_args_reader(
		SETTINGS* set,
		int       argc,
		char*     argv[],
		const char* codename,
		const char* voidname)
{
	if(argc < 3) {
		printf("%s Please specify parameters.\n", codename);
		printf("%s Format: \n", voidname);
		printf("%s     %s [B.C. filename 1] [B.C. filename 2]\n\n", voidname, argv[0]);
		printf("%s Options: \n", voidname);
		printf("%s     %s [input & output directory]\n", voidname, OPTION_DIRECTORY);
		printf("%s     %s [output filename for Dirichlet B.C.]\n", voidname, OPTION_OUTFILE);
		printf("\n");

		exit(0);
	}

	int num; 
	num = BB_std_read_args_return_char_num(
			argc, argv, OPTION_DIRECTORY);
	if(num == -1) {
		set->directory = DEF_DIRECTORY;
		printf("%s Input & output directory: %s (default)\n", codename, set->directory);
	}
	else {
		set->directory = argv[num+1];
		printf("%s Input & output directory: %s\n", codename, set->directory);
	}

	num = BB_std_read_args_return_char_num(
			argc, argv, OPTION_OUTFILE);
	if(num == -1) {
		set->outfile_bc = FILENAME_D_BC;
		printf("%s Output filename for Dirichlet B.C.: %s (default)\n", codename, set->outfile_bc);
	}
	else {
		set->outfile_bc = argv[num+1];
		printf("%s Output filename for Dirichlet B.C.: %s\n", codename, set->outfile_bc);
	}

	printf("\n");
}


void write_merged_bc_file(
		const char* infile_1,
		const char* infile_2,
		const char* outfile,
		const char* directory)
{
	FILE* fp1;
	fp1 = BBFE_sys_read_fopen(fp1, infile_1, directory);
	int nbc1, b1;
	BB_std_scan_line(&fp1, BUFFER_SIZE, "%d %d", &nbc1, &b1);
	printf("%s The number of BCs in \"%s\": %d\n", CODENAME, infile_1, nbc1);

	FILE* fp2;
	fp2 = BBFE_sys_read_fopen(fp2, infile_2, directory);
	int nbc2, b2;
	BB_std_scan_line(&fp2, BUFFER_SIZE, "%d %d", &nbc2, &b2);
	printf("%s The number of BCs in \"%s\": %d\n", CODENAME, infile_2, nbc2);

	FILE* fpo;
	fpo = BBFE_sys_write_fopen(fpo, outfile, directory);
	int nbc = nbc1 + nbc2;
	fprintf(fpo, "%d %d\n", nbc, b1);

	for(int i=0; i<nbc1; i++) {
		int nnum, bnum; double val;
		BB_std_scan_line(&fp1, BUFFER_SIZE, "%d %d %lf", &nnum, &bnum, &val);
		fprintf(fpo, "%d %d %.15e\n", nnum, bnum, val);
	}

	for(int i=0; i<nbc2; i++) {
		int nnum, bnum; double val;
		BB_std_scan_line(&fp2, BUFFER_SIZE, "%d %d %lf", &nnum, &bnum, &val);
		fprintf(fpo, "%d %d %.15e\n", nnum, bnum, val);
	}

	fclose(fp1);
	fclose(fp2);
	fclose(fpo);
}


int main(
		int   argc,
		char* argv[])
{
	printf("\n");

	SETTINGS set;

	cmd_args_reader(&set, argc, argv, CODENAME, VOIDNAME);
	
	const char* infile_1;
	const char* infile_2;
	infile_1 = argv[1];
	infile_2 = argv[2];

	write_merged_bc_file(infile_1, infile_2, set.outfile_bc, set.directory);

	printf("\n");
	
	return 0;
}
