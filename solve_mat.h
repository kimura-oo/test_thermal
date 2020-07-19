#pragma once


#ifdef WITH_MONOLIS
#include <monolis.h>
#endif


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


int MatrixCSR_get_index(
		int index_num,
		int num_dofs_on_node,
		int submat_i,
		int submat_j);


void MatrixCSR_initialize(
		Dataset_CSR* csr,
		const int total_num_nodes,
		const int total_num_elems,
		const int local_num_nodes,
		const int num_dofs_on_node,
		int** conn);


void MatrixCSR_free(
		Dataset_CSR* csr);


// copy csr2 to csr1
void MatrixCSR_copy_dataset(
		Dataset_CSR* csr1,
		Dataset_CSR* csr2);


// copy csr2 to csr1 (matrix data only)
void MatrixCSR_copy_matrix(
		Dataset_CSR* csr1,
		Dataset_CSR* csr2);


void MatrixCSR_set_nonzero_value(
		Dataset_CSR* csr,
		double val,      
		int g_node_num_i,
		int g_node_num_j,
		int submat_i,
		int submat_j); 



void MatrixCSR_add_nonzero_value(
		Dataset_CSR* csr,
		double val,      
		int g_node_num_i,
		int g_node_num_j,
		int submat_i,
		int submat_j); 


double MatrixCSR_get_nonzero_value(
		Dataset_CSR* csr,
		int g_node_num_i,
		int g_node_num_j,
		int submat_i,
		int submat_j); 


void MatrixCSR_mat_vec_multiplication(
		const Dataset_CSR* csr,
		const double* vec,
		double* ans);


int MatrixCSR_solver_CG(
		Dataset_CSR* csr,
		double* b,
		double* x,
		double epsilon,
		int max_iters,
		int show_iterations);


#ifdef WITH_MONOLIS
void MatrixCSR_set_matrix_to_monolis(
		Dataset_CSR* csr);


void MatrixCSR_solver_monolis(
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
