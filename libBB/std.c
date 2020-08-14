
#include "std.h"

/**********************************************************
 * functions for memory allocation
 **********************************************************/
double* BB_std_calloc_1d_double(
		double*   array,
		const int size)
{
	array = (double*)calloc(size, sizeof(double));

	return array;
}


double** BB_std_calloc_2d_double(
		double**   array,
		const int  size1,
		const int  size2)
{
	array = (double**)calloc(size1, sizeof(double*));
	for(int i=0; i<size1; i++) {
		array[i] = (double*)calloc(size2, sizeof(double));
	}

	return array;
}


int* BB_std_calloc_1d_int(
		int*      array,
		const int size)
{
	array = (int*)calloc(size, sizeof(int));

	return array;
}


int** BB_std_calloc_2d_int(
		int**      array,
		const int  size1,
		const int  size2)
{
	array = (int**)calloc(size1, sizeof(int*));
	for(int i=0; i<size1; i++) {
		array[i] = (int*)calloc(size2, sizeof(int));
	}

	return array;
}


bool* BB_std_calloc_1d_bool(
		bool*     array,
		const int size)
{
	array = (bool*)calloc(size, sizeof(bool));

	return array;
}


void BB_std_free_1d_double(
		double*   array,
		const int size)
{
	free(array);
	array = NULL;
}


void BB_std_free_2d_double(
		double**   array,
		const int  size1,
		const int  size2)
{
	for(int i=0; i<size1; i++) {
		free(array[i]);
		array[i] = NULL;
	}

	free(array);
	array = NULL;
}


void BB_std_free_1d_int(
		int*      array,
		const int size)
{
	free(array);
	array = NULL;
}


void BB_std_free_2d_int(
		int**   array,
		const int  size1,
		const int  size2)
{
	for(int i=0; i<size1; i++) {
		free(array[i]);
		array[i] = NULL;
	}

	free(array);
	array = NULL;
}


void BB_std_free_1d_bool(
		bool*     array,
		const int size)
{
	free(array);
	array = NULL;
}


/**********************************************************
 * functions for file IO
 **********************************************************/
bool BB_std_scan_line(
		FILE** fp,
		const int buffer_size,
		const char* format,
		...)
{
	char buf[buffer_size];
	if( fgets(buf, sizeof(buf), *fp) == NULL )
	{
		return false;
	}

	va_list va;
	va_start(va, format);
	vsscanf(buf, format, va);
	va_end(va);

	return true;
}


bool BB_std_read_file_return_char(
		char* ret_char,
		const char* filename,
		const char* identifier,
		const int buffer_size)
{
	int identical = 0;

	FILE* fp;
	fp = fopen(filename, "r");
	if(fp == NULL) {
		ret_char = NULL;
		return false;
	}

	char buf[ buffer_size];
	char buf2[buffer_size];
	char buf3[buffer_size];
	while(1) {
		if( fgets(buf, sizeof(buf), fp) == NULL) {
			break;
		}
		strcpy(buf2, buf);
		if( strstr(buf2, identifier) != NULL ) {
			sscanf(buf, "%s %s", buf3, ret_char);

			return true;
		}
	}

	return false;
}


FILE* BB_std_read_file_search_line(
		FILE*       fp,
		char*       char_line,
		const char* filename, 
		const char* identifier,
		const int   buffer_size)
{
	char buf[ buffer_size];
	char buf2[buffer_size];
	
	while(1) {
		if( fgets(buf, sizeof(buf), fp) == NULL) {
			return NULL;
		}
		
		if( strstr(buf, identifier) != NULL ) {
			strcpy(char_line, buf);
			return fp;
		}
	}
}


/**********************************************************
 * functions for treating commandline arguments
 **********************************************************/
bool BB_std_read_args_return_boolean(
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


char* BB_std_read_args_return_next_arg(
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


int BB_std_read_args_return_char_num(
		int argc,
		char* argv[],
		const char* c_option)
{
	int num = 0;

	for(int i=1; i<argc; i++) {
		if(strcmp(argv[i], c_option) == 0 ) {
			num = i;
		}
	}

	if(num == 0) {
		return -1;
	}
	else {
		return num;
	}
}
