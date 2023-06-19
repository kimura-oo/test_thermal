
#include "set.h"
#include "BB/std.h"
#include "BB/calc.h"
#include "BBFE/std/integ.h"
#include "BBFE/std/mapping.h"

/**********************************************************
 * basic utils
 **********************************************************/

void BBFE_elemmat_set_Jacobian_array(
		double*    Jacobian_ip,
		const int  num_integ_points,
		const int  elem_num,
		BBFE_DATA* fe)
{
	for(int p=0; p<num_integ_points; p++) {
		Jacobian_ip[p] = fe->geo[ elem_num ][p].Jacobian;
	}
}


void BBFE_elemmat_set_local_array_scalar(
		double*       local_val,
		BBFE_DATA*    fe,
		const double* val,
		const int     elem_num)
{
	for(int i=0; i<(fe->local_num_nodes); i++) {
		local_val[i] = val[ fe->conn[elem_num][i] ];
	}
}


void BBFE_elemmat_set_local_array_vector(
		double**       local_val,
		BBFE_DATA*     fe,
		double**       val,
		const int      elem_num,
		const int      dimension)
{
	for(int i=0; i<(fe->local_num_nodes); i++) {
		for(int d=0; d<dimension; d++) {
			local_val[i][d] = val[ fe->conn[elem_num][i] ][d];
		}
	}
}


void BBFE_elemmat_set_Jacobi_mat(
		BBFE_DATA*  fe,
		BBFE_BASIS* basis)
{
	double** local_x;
	local_x = BB_std_calloc_2d_double(local_x, fe->local_num_nodes, 3);

	for(int e=0; e<(fe->total_num_elems); e++) {

		for(int i=0; i<(fe->local_num_nodes); i++) {
			local_x[i][0] = fe->x[ fe->conn[e][i] ][0];
			local_x[i][1] = fe->x[ fe->conn[e][i] ][1];
			local_x[i][2] = fe->x[ fe->conn[e][i] ][2];
		}

		for(int p=0; p<(basis->num_integ_points); p++) {
			BBFE_std_mapping_calc_Jacobi_mat_3d(
					fe->geo[e][p].J,
					fe->local_num_nodes,
					local_x,
					basis->dN_dxi[p],
					basis->dN_det[p],
					basis->dN_dze[p]);

			fe->geo[e][p].Jacobian = BB_calc_mat3d_determinant(
					fe->geo[e][p].J);
		}
	}

	BB_std_free_2d_double(local_x, fe->local_num_nodes, 3);
}


void BBFE_elemmat_set_shapefunc_derivative(
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis)
{
	double J_inv[3][3];

	for(int e=0; e<(fe->total_num_elems); e++) {
		for(int p=0; p<(basis->num_integ_points); p++) {
			BB_calc_mat3d_inverse(
					fe->geo[e][p].J, fe->geo[e][p].Jacobian, J_inv);

			for(int i=0; i<(fe->local_num_nodes); i++) {
				for(int d=0; d<3; d++) {
					fe->geo[e][p].grad_N[i][d] =
						J_inv[d][0] * basis->dN_dxi[p][i] +
						J_inv[d][1] * basis->dN_det[p][i] +
						J_inv[d][2] * basis->dN_dze[p][i];
				}
			}
		}
	}
}

/**********************************************************
 * setter of simple global matrix
 **********************************************************/
void BBFE_elemmat_set_global_mat_cmass_const(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		double      coef,
		int         block_size)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;

					val_ip[p] += coef * basis->N[p][i] * basis->N[p][j];
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				for(int b=0; b<block_size; b++) {
					monolis_add_scalar_to_sparse_matrix_R(
							monolis,
							fe->conn[e][i], fe->conn[e][j], b, b,
							integ_val);
				}
			}
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);
}


void BBFE_elemmat_set_global_mat_cmass_const_C(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		double      coef,
		int         block_size)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;

					val_ip[p] += coef * basis->N[p][i] * basis->N[p][j];
				}

				double _Complex integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				for(int b=0; b<block_size; b++) {
					monolis_add_scalar_to_sparse_matrix_C(
							monolis,
							fe->conn[e][i], fe->conn[e][j], b, b,
							integ_val);
				}
			}
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);
}


void BBFE_elemmat_set_global_mat_Laplacian_const(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		double      coef)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;

					val_ip[p] += 
						- coef * ( fe->geo[e][p].grad_N[i][0] * fe->geo[e][p].grad_N[j][0] + 
								fe->geo[e][p].grad_N[i][1] * fe->geo[e][p].grad_N[j][1] + 
								fe->geo[e][p].grad_N[i][2] * fe->geo[e][p].grad_N[j][2] );
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				monolis_add_scalar_to_sparse_matrix_R(
						monolis,
						fe->conn[e][i], fe->conn[e][j], 0, 0,
						integ_val);
			}
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);
}

void BBFE_elemmat_set_global_mat_Laplacian_const_C(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		double      coef)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;

					val_ip[p] += 
						- coef * ( fe->geo[e][p].grad_N[i][0] * fe->geo[e][p].grad_N[j][0] + 
								fe->geo[e][p].grad_N[i][1] * fe->geo[e][p].grad_N[j][1] + 
								fe->geo[e][p].grad_N[i][2] * fe->geo[e][p].grad_N[j][2] );
				}

				double _Complex integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				monolis_add_scalar_to_sparse_matrix_C(
						monolis,
						fe->conn[e][i], fe->conn[e][j], 0, 0,
						integ_val);
			}
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);
}