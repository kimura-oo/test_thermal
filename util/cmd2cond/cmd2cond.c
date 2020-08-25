#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "BB/std.h"

const char* CODENAME           = "cmd2cond >";
const char* VOIDNAME           = "          ";
const char* FILENAME           = "cond.dat";
const char* OPTION_OUTPUT_DIR  = "-d";
const char* DEFAULT_OUTPUT_DIR = ".";
const char* OPTION_ADD_MODE    = "--add";

const int BUFFER_SIZE = 10000;


FILE* output_int_data(
		FILE* fp,
		char* id,
		int   val_num,
		int*  val)
{
	fprintf(fp, "%s %d\n", id, val_num);
	
	for(int i=0; i<val_num; i++) {
		fprintf(fp, "%d\n", val[i]);
	}

	fprintf(fp, "\n");

	return fp;
}


FILE* output_double_data(
		FILE*   fp,
		char*   id,
		int     val_num,
		double* val)
{
	fprintf(fp, "%s %d\n", id, val_num);
	
	for(int i=0; i<val_num; i++) {
		fprintf(fp, "%.15e\n", val[i]);
	}

	fprintf(fp, "\n");

	return fp;
}


int main(
		int   argc,
		char* argv[])
{
	printf("\n");
	
	if(argc < 2) {
		printf("%s Please imput parameters.\n", CODENAME);
		printf("%s Format: \n", CODENAME);
		printf("%s     ./cmd2cond \"[identifier like (#id)]\" [data type (int or double)] [num. data: n] [data 1] ... [data n] ...\n", VOIDNAME);

		printf("\n%s Options: \n", CODENAME);
		printf("%s     %s [output directory] (default: .)\n", VOIDNAME, OPTION_OUTPUT_DIR);
		printf("%s     %s (add contents to existing file)\n", VOIDNAME, OPTION_ADD_MODE);
		printf("\n");
		exit(0);
	}

	const char* directory = BB_std_read_args_return_next_arg(
			argc, argv, OPTION_OUTPUT_DIR);
	if( directory == NULL ) {
		directory = DEFAULT_OUTPUT_DIR;
	}
	printf("%s Output directory: %s\n", CODENAME, directory);

	char filename[BUFFER_SIZE];
	snprintf(filename, BUFFER_SIZE, "%s/%s", directory, FILENAME);

	FILE* fp;
	if( BB_std_read_args_return_boolean(argc, argv, OPTION_ADD_MODE) ) {
		printf("%s Output file \"%s\" (add mode)\n\n", CODENAME, filename);
		fp = fopen(filename, "a");
	}
	else {
		printf("%s Output file \"%s\" (write mode)\n\n", CODENAME, filename);
		fp = fopen(filename, "w");
	}
	if(fp == NULL) {
		printf("%s ERROR: file \"%s\" cannot be opened.\n\n", CODENAME, filename);
		exit(EXIT_FAILURE);
	}

	int start_num = 1;
	while(1) {
		int num = BB_std_read_args_search_num(argc, argv, start_num, "#");
		if(num == -1 || num >= argc-2) {break;}
		
		printf("%s ID: %s", CODENAME, argv[num]);

		if( strcmp(argv[num+1], "int") == 0 ) {
			int val_num = atoi( argv[num+2] );
			
			if(val_num < 0) {
				printf("\n%s ERROR: the number of variable %d should be larger than 0.\n\n", CODENAME, val_num);
				fclose(fp);
				exit(EXIT_FAILURE);
			}
			printf(" (%d int data)\n", val_num);

			int val[val_num];
			for(int i=0; i<val_num; i++) {
				val[i] = atoi( argv[num+3+i] );
				printf("%s %d: %d\n", CODENAME, i+1, val[i]);				
			}

			fp = output_int_data(fp, argv[num], val_num, val);
		}
		else if ( strcmp(argv[num+1], "double") == 0 ) {
			int val_num = atoi( argv[num+2] );

			if(val_num < 0) {
				printf("%s ERROR: the number of variable %d should be larger than 0.\n\n", CODENAME, val_num);
				fclose(fp);
				exit(EXIT_FAILURE);
			}
			printf(" (%d double data)\n", val_num);

			double val[val_num];
			for(int i=0; i<val_num; i++) {
				val[i] = atof( argv[num+3+i] );
				printf("%s %d: %.15e\n", CODENAME, i+1, val[i]);				
			}

			fp = output_double_data(fp, argv[num], val_num, val);
		}
		else {
			printf("%s ERROR: type of variable should be specified. (int or double)\n\n", CODENAME);
			fclose(fp);
			exit(EXIT_FAILURE);
		}

		start_num = num+1;
		printf("\n");
	}

	fclose(fp);

	return 0;
}
