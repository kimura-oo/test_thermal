
#pragma once


#ifdef WITH_MONOLIS
#include <monolis.h>
#endif

typedef struct
{
	int      total_num_nodes;
	double** x;
	
	int   total_num_elems;
	int   local_num_nodes;
	int** conn;
} FE_DATA;


typedef struct {
	int num_nodes;
	int num_dofs_on_node;
	int num_nonzero_vals;
	int vec_length;       // = num_nodes * num_dofs_on_node

	int* index;           // csr index [num_nodes+1]
	int* item;            // csr item  [num_nonzero_vals]
	double* mat;          // csr matrix [num_nonzero_vals]
	double* rhs;          // right-hand-side vector [vec_length]
	double* sol;          // solution vector [vec_length]

#ifdef WITH_MONOLIS
	monolis_mat monomat;
#endif

} Dataset_CSR;


int index_CSR_mat(
		int index_num,
		int num_dofs_on_node,
		int submat_i,
		int submat_j);


void init_dataset_CSR(
		Dataset_CSR* csr,
		FE_DATA* fe,
		int num_dofs_on_node);


void free_dataset_CSR(
		Dataset_CSR* csr);


// copy csr2 to csr1
void copy_dataset_CSR(
		Dataset_CSR* csr1,
		Dataset_CSR* csr2);


// copy csr2 to csr1 (matrix data only)
void copy_CSR_matrix(
		Dataset_CSR* csr1,
		Dataset_CSR* csr2);


void set_nonzero_value_to_CSR_matrix(
		Dataset_CSR* csr,
		double val,      
		int g_node_num_i,
		int g_node_num_j,
		int submat_i,
		int submat_j); 



void add_nonzero_value_to_CSR_matrix(
		Dataset_CSR* csr,
		double val,      
		int g_node_num_i,
		int g_node_num_j,
		int submat_i,
		int submat_j); 


double get_nonzero_value_from_CSR_matrix(
		Dataset_CSR* csr,
		int g_node_num_i,
		int g_node_num_j,
		int submat_i,
		int submat_j); 


void mat_vec_multiplication_CSR(
		const Dataset_CSR* csr,
		const double* vec,
		double* ans);


int solve_mat_CG_CSR(
		Dataset_CSR* csr,
		double* b,
		double* x,
		double epsilon,
		int max_iters,
		int show_iterations);


#ifdef WITH_MONOLIS
void set_CSR_matrix_value_to_monolis(
		Dataset_CSR* csr);


void solve_CSR_matrix_by_monolis(
		Dataset_CSR* csr,
		double* rhs,
		int method_num,
		int precond_num,
		int max_num_iters,
		double tol,
		int is_scaling,
		int is_reordering,
		int is_init_x,
		int show_iterations);
#endif
