
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>


#include "solve_mat.h"


static const int MAX_NUM_CONNECRED_ELEMS = 15;


int index_CSR_mat(
		int index_num,
		int num_dofs_on_node,
		int submat_i,
		int submat_j)
{
	return (index_num*num_dofs_on_node*num_dofs_on_node
			+ submat_i*num_dofs_on_node
			+ submat_j);
}


static bool same_index_set_exists(
		int* item_in_index,
		int  item_val,
		int num_sets)
{
	for(int s=0; s<num_sets; s++) {
		if( item_in_index[s] == item_val) {
			return true;
		}
	}

	return false;
}


static bool sort_in_ascending_order(
		int* item_in_index,  
		int num_sets)
{
	for (int i=0; i<num_sets; ++i) {
		for (int j=i+1; j<num_sets; ++j) {
			if ( item_in_index[i] > item_in_index[j] ) {
				int tmp =  item_in_index[i];
				item_in_index[i] = item_in_index[j];
				item_in_index[j] = tmp;
			}
		}
	}
}


static void set_csr_index_and_item(
		Dataset_CSR* csr,
		const int total_num_nodes,
		const int total_num_elems,
		const int local_num_nodes,
		/*const*/ int** conn)       // connectivity info of finite elements
{
	int counter = 0;

	int n = local_num_nodes;

	int num_max_buf =  MAX_NUM_CONNECRED_ELEMS*n*n;
	if(num_max_buf > total_num_nodes) {
		num_max_buf = total_num_nodes;
	}

	int** buf;
	int*  num_buf;
	int*  num_buf2;
	buf      = (int**)calloc(total_num_nodes, sizeof(int*));
	num_buf  = (int* )calloc(total_num_nodes, sizeof(int ));
	num_buf2 = (int* )calloc(total_num_nodes, sizeof(int ));
	for(int i=0; i<total_num_nodes; i++) {
		buf[i] = (int*)calloc(num_max_buf, sizeof(int));

		for(int j=0; j<num_max_buf; j++) {
			buf[i][j] = -1;
		}
	}

	// counting the number of nonzero values in upper triangle mat
	int g_num_i, g_num_j;
	for(int e=0; e<total_num_elems; e++) {

		for(int i=0; i<local_num_nodes; i++) {
			for(int j=i+1; j<local_num_nodes; j++) {
				g_num_i = conn[e][i];
				g_num_j = conn[e][j];

				int index_set[2];
				if(g_num_i < g_num_j) {
					index_set[0] = g_num_i;
					index_set[1] = g_num_j;
				}
				else {
					index_set[1] = g_num_i;
					index_set[0] = g_num_j;
				}

				if( 
						!same_index_set_exists
						(buf[ index_set[0] ], index_set[1], num_buf[ index_set[0] ] ) 
				  ) {
					buf[ index_set[0] ][ num_buf[ index_set[0] ] ] = index_set[1];

					num_buf[  index_set[0] ]++;
					num_buf2[ index_set[0] ]++;
					counter++;
				}
			}
		}
	}

	// substituting value into lower triangle mat
	for(int i=0; i<total_num_nodes; i++) {

		for(int s=0; s<num_buf[i]; s++) {
			int index = buf[i][s];
			buf[ index ][ num_buf2[index] ] = i;

			num_buf2[ index ]++;
			counter++;
		}
	}

	// adding diagonal elems
	for(int i=0; i<total_num_nodes; i++) {
		num_buf[i] = num_buf2[i];

		buf[i][ num_buf[i] ] = i;

		num_buf[i]++;
		counter++;
	}

	csr->num_nonzero_vals = counter;
	csr->index = (int*)calloc(total_num_nodes + 1,   sizeof(int));
	csr->item  = (int*)calloc(csr->num_nonzero_vals, sizeof(int));

	// set CSR index and item
	csr->index[0] = 0;
	int sum = 0;
	for(int i=0; i<total_num_nodes; i++) {
		sort_in_ascending_order(buf[i], num_buf[i]);

		sum += num_buf[i];
		csr->index[i+1] = sum;

		int start = csr->index[i];
		for(int j=0; j<num_buf[i]; j++) {
			csr->item[start + j] = buf[i][j];
		}
	}

	for(int i=0; i<total_num_nodes; i++) {
		free(buf[i]);
	}
	free(buf);
	free(num_buf);
	free(num_buf2);
}


#ifdef WITH_MONOLIS
static void init_monolis_matrix(
		monolis_mat* monomat,
		int num_nodes,
		int num_nonzero_vals,
		int num_dofs_on_node,
		int* csr_index,
		int* csr_item)
{
	monolis_convert_csr_get_size(
			num_nodes,
			num_nonzero_vals,
			csr_index,
			csr_item,
			&(monomat->NPU),
			&(monomat->NPL));

	int N    = num_nodes;
	int NPU  = monomat->NPU;
	int NPL  = monomat->NPL;
	int NDOF = num_dofs_on_node;
	monomat->N    = N;
	monomat->NPU  = NPU;
	monomat->NPL  = NPL;
	monomat->NDOF = NDOF;

	monomat->indexU = (int*)calloc(N+1, sizeof(int));
	monomat->indexL = (int*)calloc(N+1, sizeof(int));	
	monomat->itemU  = (int*)calloc(NPU, sizeof(int));
	monomat->itemL  = (int*)calloc(NPL, sizeof(int));

	monomat->D  = (double*)calloc(N  *NDOF*NDOF, sizeof(double));
	monomat->AU = (double*)calloc(NPU*NDOF*NDOF, sizeof(double));
	monomat->AL = (double*)calloc(NPL*NDOF*NDOF, sizeof(double));

	monomat->X = (double*)calloc(N*NDOF, sizeof(double));
	monomat->B = (double*)calloc(N*NDOF, sizeof(double));

	monolis_convert_csr_get_index(
			num_nodes, 
			num_nonzero_vals, 
			csr_index, 
			csr_item,
			monomat->NPU,
			monomat->NPL,
			monomat->indexU,
			monomat->itemU,
			monomat->indexL,
			monomat->itemL);
}


static void free_monolis_matrix(
		monolis_mat* monomat)
{

	free( monomat->indexU );
	free( monomat->indexL );
	free( monomat->itemU  );
	free( monomat->itemL  );

	free( monomat->D  );
	free( monomat->AU );
	free( monomat->AL );

	free( monomat->X );
	free( monomat->B );
}
#endif


void init_dataset_CSR(
		Dataset_CSR* csr,
		FE_DATA* fe,
		int num_dofs_on_node)
{

	csr->num_nodes        = fe->total_num_nodes;
	csr->num_dofs_on_node = num_dofs_on_node;
	csr->vec_length       = csr->num_nodes * csr->num_dofs_on_node;

	set_csr_index_and_item(
			csr,
			fe->total_num_nodes,
			fe->total_num_elems,
			fe->local_num_nodes,
			fe->conn);

	int nn   = csr->num_nodes;
	int ndof = csr->num_dofs_on_node;
	int nnz  = csr->num_nonzero_vals;
	csr->mat = (double*)calloc(nnz*ndof*ndof, sizeof(double));
	csr->rhs = (double*)calloc(nn*ndof      , sizeof(double));
	csr->sol = (double*)calloc(nn*ndof      , sizeof(double));

	set_csr_index_and_item(
			csr, 
			fe->total_num_nodes,
			fe->total_num_elems,
			fe->local_num_nodes,
			fe->conn);

#ifdef WITH_MONOLIS
	init_monolis_matrix(
			&(csr->monomat),
			csr->num_nodes,
			csr->num_nonzero_vals,
			csr->num_dofs_on_node,
			csr->index,
			csr->item);
#endif
}


void free_dataset_CSR(
		Dataset_CSR* csr)
{
	free(csr->index);
	free(csr->item);

	free(csr->mat);
	free(csr->rhs);
	free(csr->sol);

#ifdef WITH_MONOLIS
	free_monolis_matrix(&(csr->monomat));
#endif
}


void copy_dataset_CSR(
		Dataset_CSR* csr1,
		Dataset_CSR* csr2)
{
	csr1->num_nodes        = csr2->num_nodes;
	csr1->num_dofs_on_node = csr2->num_dofs_on_node;
	csr1->num_nonzero_vals = csr2->num_nonzero_vals;
	csr1->vec_length       = csr2->vec_length;

	int ndof = csr1->num_dofs_on_node;
	int nn   = csr1->num_nodes;
	int nnz  = csr1->num_nonzero_vals;

	csr1->index  = (int*)calloc(nn+1,          sizeof(int));
	csr1->item   = (int*)calloc(nnz,           sizeof(int));
	csr1->mat = (double*)calloc(nnz*ndof*ndof, sizeof(double));
	csr1->rhs = (double*)calloc(nn*ndof,       sizeof(double));
	csr1->sol = (double*)calloc(nn*ndof,       sizeof(double));

	for(int i=0; i<nn+1; i++) {
		csr1->index[i] = csr2->index[i];
	}
	for(int i=0; i<nnz; i++) {
		csr1->item[i] = csr2->item[i];
	}
	for(int i=0; i<nnz*ndof*ndof; i++) {
		csr1->mat[i] = csr2->mat[i];
	}
	for(int i=0; i<nn*ndof; i++) {
		csr1->rhs[i] = csr2->rhs[i];
	}
	for(int i=0; i<nn*ndof; i++) {
		csr1->sol[i] = csr2->sol[i];
	}

#ifdef WITH_MONOLIS
	init_monolis_matrix(
			&(csr1->monomat),
			csr1->num_nodes,
			csr1->num_nonzero_vals,
			csr1->num_dofs_on_node,
			csr1->index,
			csr1->item);

	set_CSR_matrix_value_to_monolis(csr1);
#endif
}


void copy_CSR_matrix(
		Dataset_CSR* csr1,
		Dataset_CSR* csr2)
{

	int ndof = csr1->num_dofs_on_node;
	int nnz  = csr1->num_nonzero_vals;

	for(int i=0; i<nnz*ndof*ndof; i++) {
		csr1->mat[i] = csr2->mat[i];
	}

#ifdef WITH_MONOLIS
	set_CSR_matrix_value_to_monolis(csr1);
#endif
}


void set_nonzero_value_to_CSR_matrix(
		Dataset_CSR* csr,
		double val,      // block_length x block_length submatrix
		int g_node_num_i,
		int g_node_num_j,
		int submat_i,
		int submat_j)  
{
	int n_start = csr->index[ g_node_num_i  ];
	int n_end   = csr->index[ g_node_num_i+1];

	int csr_num = -1;
	for(int n=n_start; n<n_end; n++) {
		if( csr->item[n] == g_node_num_j ) {
			csr_num = n;
			break;
		}
	}

	if(csr_num == -1) {
		return;	
	}

	int n = index_CSR_mat(
			csr_num,
			csr->num_dofs_on_node,
			submat_i,
			submat_j);

	csr->mat[n] = val;

}


void add_nonzero_value_to_CSR_matrix(
		Dataset_CSR* csr,
		double val,      // block_length x block_length submatrix
		int g_node_num_i,
		int g_node_num_j,
		int submat_i,
		int submat_j)  
{
	int n_start = csr->index[ g_node_num_i   ];
	int n_end   = csr->index[ g_node_num_i+1 ];

	int csr_num = -1;
	for(int n=n_start; n<n_end; n++) {
		if( csr->item[n] == g_node_num_j ) {
			csr_num = n;
			break;
		}
	}

	if(csr_num == -1) {
		return;	
	}

	int n = index_CSR_mat(
			csr_num,
			csr->num_dofs_on_node,
			submat_i,
			submat_j);

	csr->mat[n] += val;

}


double get_nonzero_value_from_CSR_matrix(
		Dataset_CSR* csr,
		int g_node_num_i,
		int g_node_num_j,
		int submat_i,
		int submat_j)  
{
	int n_start = csr->index[ g_node_num_i   ];
	int n_end   = csr->index[ g_node_num_i+1 ];

	int csr_num = -1;
	for(int n=n_start; n<n_end; n++) {
		if( csr->item[n] == g_node_num_j ) {
			csr_num = n;
			break;
		}
	}

	if(csr_num == -1) {
		return 0.0;	
	}

	int n = index_CSR_mat(
			csr_num,
			csr->num_dofs_on_node,
			submat_i,
			submat_j);

	return( csr->mat[n] );

}


void mat_vec_multiplication_CSR(
		const Dataset_CSR* csr,
		const double* vec,
		double* ans)
{
	int nn   = csr->num_nodes; 
	int ndof = csr->num_dofs_on_node;

	for(int i=0; i<nn; i++) {
		int start = csr->index[i  ];
		int end   = csr->index[i+1];

		for(int k=0; k<ndof; k++) {
			ans[ndof*i + k] = 0.0;
			
			for(int j=start; j<end; j++) {

				for(int l=0; l<ndof; l++) {
					ans[ndof*i + k] += 
						csr->mat[ index_CSR_mat(j, ndof, k, l) ] * vec[ ndof*(csr->item[j]) + l ];
				}
			}
		}
	}
}


int solve_mat_CG_CSR(
		Dataset_CSR* csr,
		double* b,
		double* x,
		double epsilon,
		int max_iters,
		int show_iterations )
{
	int n = csr->vec_length;

	int k = 0;
	double* r;  r =  (double*)calloc(n, sizeof(double));
	double* p;  p =  (double*)calloc(n, sizeof(double));
	double* ap; ap = (double*)calloc(n, sizeof(double));
	double alpha, beta;

	mat_vec_multiplication_CSR(csr, x, ap);

	double norm0 = 0.0;
	for(int i=0; i<n; i++) {
		r[i] = b[i] - ap[i];
		p[i] = r[i];
		norm0 += r[i]*r[i];
	}
	norm0 = sqrt(norm0);

	//iteration
	for(k=0; k<max_iters; k++) {

		//calc alpha
		mat_vec_multiplication_CSR(csr, p, ap);

		double pap = 0.0;  double rr_p = 0.0;
		for(int i=0; i<n; i++) {
			pap += p[i]*ap[i];
			rr_p += r[i]*r[i];
		}
		alpha = rr_p/pap;

		//update r and x
		double rr = 0.0;  double norm =	0.0;
		for(int i=0; i<n; i++) {
			r[i] = r[i] - alpha*ap[i];
			x[i] = x[i] + alpha*p[i];
			rr  += r[i]*r[i];
		}
		norm = sqrt(rr);

		if(show_iterations) {
			printf("CGloop: %d: %e\n", k, norm/norm0);
		}

		//convergence critria
		if(norm/norm0 < epsilon) break;

		//calc beta
		beta = rr/rr_p;

		//update p
		for(int i=0; i<n; i++) {
			p[i] = r[i] + beta*p[i];
		}	
	}

	free(r);
	free(p);
	free(ap);

	return k;
}


#ifdef WITH_MONOLIS
void set_CSR_matrix_value_to_monolis(
		Dataset_CSR* csr)
{

	monolis_convert_csr_update_matrix_entry(
			csr->num_nodes, 
			csr->num_nonzero_vals, 
			csr->num_dofs_on_node, 
			csr->mat, 
			csr->index,
			csr->item, 
			csr->monomat.NPU, 
			csr->monomat.NPL, 
			csr->monomat.D, 
			csr->monomat.AU, 
			csr->monomat.AL, 
			csr->monomat.indexU, 
			csr->monomat.itemU, 
			csr->monomat.indexL, 
			csr->monomat.itemL);

	int num_mat_elems = csr->num_nodes * csr->num_dofs_on_node;
	for(int i=0; i<num_mat_elems; i++) {
		csr->monomat.B[i] = csr->rhs[i];
		csr->monomat.X[i] = csr->sol[i];
	}

}


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
		int show_iterations)
{

	monolis_serial(
			csr->monomat.N, 
			csr->monomat.NDOF, 
			csr->monomat.NPU, 
			csr->monomat.NPL, 
			csr->monomat.D, 
			csr->monomat.AU, 
			csr->monomat.AL, 
			csr->monomat.X,
			rhs,
			csr->monomat.indexU, 
			csr->monomat.itemU, 
			csr->monomat.indexL,
			csr->monomat.itemL, 
			method_num,
			precond_num, 
			max_num_iters, 
			tol, 
			is_scaling,
			is_reordering,
			is_init_x,
			show_iterations);

	int num_mat_elems = csr->num_nodes * csr->num_dofs_on_node;
	for(int i=0; i<num_mat_elems; i++) {
		csr->sol[i] = csr->monomat.X[i]; 
	}
}
#endif


int main()
{

	FE_DATA fe;
	fe.total_num_nodes = 3;
	fe.local_num_nodes = 2;
	fe.total_num_elems = 2;

	fe.conn = (int**)calloc(fe.total_num_elems, sizeof(int*));
	for(int i=0; i<fe.total_num_elems; i++) {
		fe.conn[i] = (int*)calloc(fe.local_num_nodes, sizeof(int));
	}
	fe.conn[0][0] = 0; fe.conn[0][1] = 1;
	fe.conn[1][0] = 1; fe.conn[1][1] = 2;

	Dataset_CSR csr;
	init_dataset_CSR(
			&csr, &fe, 4);

	int ndof = csr.num_dofs_on_node;

	int count = 1;
	double cmat[csr.vec_length][csr.vec_length];
	double vec[csr.vec_length];
	double cans[csr.vec_length];
	double ans[csr.vec_length];
	printf("%d\n", csr.vec_length);
	
	for(int i=0; i<csr.vec_length; i++) {
		for(int j=0; j<csr.vec_length; j++) {
			cmat[i][j] = 0.0;
		}
		vec[i] = 1.0;
		cans[i] = 0.0;
		ans[i] = 0.0;
	}

	for(int i=0; i<csr.num_nodes; i++) {
		for(int j=csr.index[i]; j<csr.index[i+1]; j++) {
			for(int k=0; k<ndof; k++) {
				for(int l=0; l<ndof; l++) {
					add_nonzero_value_to_CSR_matrix(
							&csr, count, i, csr.item[j], k, l);
					cmat[ ndof*i + k ][ ndof*csr.item[j] + l ] = count;
					count++;
				}
			}
		}
	}

	for(int i=0; i<csr.vec_length; i++) {
		for(int j=0; j<csr.vec_length; j++) {
			printf("%03.1e ", cmat[i][j]);
		}
		printf("\n");
	}
	
	for(int i=0; i<csr.vec_length; i++) {
		for(int j=0; j<csr.vec_length; j++) {
			cans[i] += cmat[i][j]*vec[j];
		}
	}
	
	mat_vec_multiplication_CSR(
			&csr, vec, ans);

	for(int i=0; i<csr.vec_length; i++) {
		printf("%e %e\n", cans[i], ans[i]);
	}

	return 0;
}
