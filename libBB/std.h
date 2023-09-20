#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdbool.h>


/**********************************************************
 * functions for memory allocation
 **********************************************************/
// rule: BB_std_calloc_[dimension of array]_[type of array]
double* BB_std_calloc_1d_double(
		double*   array,
		const int size);

double** BB_std_calloc_2d_double(
		double**   array,
		const int  size1,
		const int  size2);

double*** BB_std_calloc_3d_double(
		double***  array,
		const int  size1,
		const int  size2,
		const int  size3);

int* BB_std_calloc_1d_int(
		int*      array,
		const int size);

int** BB_std_calloc_2d_int(
		int**      array,
		const int  size1,
		const int  size2);

bool* BB_std_calloc_1d_bool(
		bool*     array,
		const int size);

bool** BB_std_calloc_2d_bool(
		bool**     array,
		const int  size1,
		const int  size2);

void BB_std_free_1d_double(
		double*   array,
		const int size);

void BB_std_free_2d_double(
		double**   array,
		const int  size1,
		const int  size2);

void BB_std_free_3d_double(
		double***  array,
		const int  size1,
		const int  size2,
		const int  size3);

void BB_std_free_1d_int(
		int*      array,
		const int size);

void BB_std_free_2d_int(
		int**   array,
		const int  size1,
		const int  size2);

void BB_std_free_1d_bool(
		bool*     array,
		const int size);

void BB_std_free_2d_bool(
		bool**     array,
		const int  size1,
		const int  size2);

/**********************************************************
 * functions for file IO
 **********************************************************/
bool BB_std_scan_line(
		FILE** fp,
		const int buffer_size,
		const char* format,
		...);

bool BB_std_read_file_return_char(
		char* ret_char,
		const char* filename,
		const char* identifier,
		const int buffer_size);

FILE* BB_std_read_file_search_line(
		FILE*       fp,
		char*       char_line,
		const char* identifier,
		const int   buffer_size);

int BB_std_read_file_get_val_double(
		double*     val,
		const char* filename,
		const char* identifier,
		const int   buffer_size);

int BB_std_read_file_get_val_int(
		int*        val,
		const char* filename,
		const char* identifier,
		const int   buffer_size);

// p means print
int BB_std_read_file_get_val_double_p(
		double*     val,
		const char* filename,
		const char* identifier,
		const int   buffer_size,
		const char* codename);

int BB_std_read_file_get_val_int_p(
		int*        val,
		const char* filename,
		const char* identifier,
		const int   buffer_size,
		const char* codename);

/**********************************************************
 * functions for treating commandline arguments
 **********************************************************/
bool BB_std_read_args_return_boolean(
		int argc,
		char* argv[],
		const char* c_option);

char* BB_std_read_args_return_next_arg(
		int argc,
		char* argv[],
		const char* c_option);

int BB_std_read_args_return_char_num(
		int argc,
		char* argv[],
		const char* c_option);

int BB_std_read_args_search_num(
		int         argc,
		char*       argv[],
		int         start_num,
		const char* identifier);

/**********************************************************
 complex memory allocation
 **********************************************************/

double _Complex* BB_std_calloc_1d_double_C(
	double _Complex*   array,
	const int 			size);

void BB_std_free_1d_double_C(
		double _Complex*   array,
		const int           size);

double _Complex ** BB_std_calloc_2d_double_C(
		double _Complex**   array,
		const int           size1,
		const int           size2);

void BB_std_free_2d_double_C(
		double _Complex **   array,
		const int            size1,
		const int            size2);
