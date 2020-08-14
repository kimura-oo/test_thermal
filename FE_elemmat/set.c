
#include "set.h"
#include "../libBB/std.h"
#include "../libBB/calc.h"
#include "../FE_std/integ.h"
#include "../FE_std/mapping.h"

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
