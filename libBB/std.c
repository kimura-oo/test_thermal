
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


double*** BB_std_calloc_3d_double(
		double***  array,
		const int  size1,
		const int  size2,
		const int  size3)
{
	array = (double***)calloc(size1, sizeof(double**));
	for(int i=0; i<size1; i++) {
		array[i] = (double**)calloc(size2, sizeof(double*));
		
		for(int j=0; j<size2; j++) {
			array[i][j] = (double*)calloc(size3, sizeof(double));
		}
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


bool** BB_std_calloc_2d_bool(
		bool**     array,
		const int  size1,
		const int  size2)
{
	array = (bool**)calloc(size1, sizeof(bool*));
	for(int i=0; i<size1; i++) {
		array[i] = (bool*)calloc(size2, sizeof(bool));
	}

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


void BB_std_free_3d_double(
		double***  array,
		const int  size1,
		const int  size2,
		const int  size3)
{
	for(int i=0; i<size1; i++) {
		for(int j=0; j<size2; j++) {
			free(array[i][j]);
			array[i][j] = NULL;
		}
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


void BB_std_free_2d_bool(
		bool**     array,
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


int BB_std_read_file_get_val_double(
		double*     val,
		const char* filename,
		const char* identifier,
		const int   buffer_size)
{
	FILE* fp = fopen(filename, "r");

	char buf[ buffer_size];
	char buf2[buffer_size];
	fp = BB_std_read_file_search_line(
			fp, 
			buf,
			identifier,
			buffer_size);
	if(fp == NULL) { return 0; }

	int val_num = 0;
	sscanf(buf, "%s %d", buf2, &val_num);
	if( strcmp(buf2, identifier) != 0 ) { return 0;}
	
	if(val_num <= 0) { return val_num; }

	for(int i=0; i<val_num; i++) {
		BB_std_scan_line(&fp, buffer_size, "%lf", &val[i]);
	}
	
	return val_num;
}


int BB_std_read_file_get_val_int(
		int*        val,
		const char* filename,
		const char* identifier,
		const int   buffer_size)
{
	FILE* fp = fopen(filename, "r");

	char buf[ buffer_size];
	char buf2[buffer_size];
	fp = BB_std_read_file_search_line(
			fp, 
			buf,
			identifier,
			buffer_size);
	if(fp == NULL) { return 0; }

	int val_num = 0;
	sscanf(buf, "%s %d", buf2, &val_num);
	if( strcmp(buf2, identifier) != 0 ) { return 0;}
	
	if(val_num <= 0) { return val_num; }

	for(int i=0; i<val_num; i++) {
		BB_std_scan_line(&fp, buffer_size, "%d", &val[i]);
	}
	
	return val_num;
}


int BB_std_read_file_get_val_double_p(
		double*     val,
		const char* filename,
		const char* identifier,
		const int   buffer_size,
		const char* codename)
{
	int val_num;
	val_num = BB_std_read_file_get_val_double(
			val, filename, identifier, buffer_size);

	printf("%s %s", codename, identifier);
	if(val_num <= 0 ) {
		printf(" is not specified.\n");
	}
	else {
		printf(": ");
		for(int i=0; i<val_num; i++) {
			printf("%e ", val[i]);
		}
		printf("\n");
	}

	return val_num;
}


int BB_std_read_file_get_val_int_p(
		int*        val,
		const char* filename,
		const char* identifier,
		const int   buffer_size,
		const char* codename)
{
	int val_num;
	val_num = BB_std_read_file_get_val_int(
			val, filename, identifier, buffer_size);

	printf("%s %s", codename, identifier);
	if(val_num <= 0 ) {
		printf(" is not specified.\n");
	}
	else {
		printf(": ");
		for(int i=0; i<val_num; i++) {
			printf("%d ", val[i]);
		}
		printf("\n");
	}

	return val_num;
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


int BB_std_read_args_search_num(
		int         argc,
		char*       argv[],
		int         start_num,
		const char* identifier)
{	
	for(int i=start_num; i<argc; i++) {
		if( strstr(argv[i], identifier) != NULL ) {
			return i;
		}
	}

	return -1;
}

/**********************************************************
 complex memory allocation
 **********************************************************/

double _Complex* BB_std_calloc_1d_double_C(
	double _Complex*   array,
	const int 			size)
{
	array = (double _Complex*)calloc(size, sizeof(double _Complex));

	return array;
}


void BB_std_free_1d_double_C(
		double _Complex*   array,
		const int           size)
{
	free(array);
	array = NULL;
}


double _Complex ** BB_std_calloc_2d_double_C(
		double _Complex**   array,
		const int           size1,
		const int           size2)
{
	array = (double _Complex**)calloc(size1, sizeof(double _Complex*));
	for(int i=0; i<size1; i++) {
		array[i] = (double _Complex*)calloc(size2, sizeof(double _Complex));
	}

	return array;
}

void BB_std_free_2d_double_C(
		double _Complex **   array,
		const int            size1,
		const int            size2)
{
	for(int i=0; i<size1; i++) {
		free(array[i]);
		array[i] = NULL;
	}

	free(array);
	array = NULL;
}