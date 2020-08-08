
#include "equivval.h"
#include "../libBB/std.h"
#include "set.h"
#include "../FE_std/mapping.h"
#include "../FE_std/integ.h"

void BBFE_elemmat_equivval_volume_smooth_function(
		double*      equiv_val,
		FE_DATA*     fe,
		FE_3D_BASIS* basis,
		double       t,
		double       (*func)(double, double, double, double)) // scalar function(x, y, z, t)
{

	for(int i=0; i<(fe->total_num_nodes); i++) {
		equiv_val[i] = 0.0;
	}

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip      , basis->num_integ_points);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip , basis->num_integ_points);

	double** local_x;
	local_x = BB_std_calloc_2d_double(local_x, fe->local_num_nodes, 3);

	for(int e=0; e<(fe->total_num_elems); e++) {
		for(int i=0; i<(fe->local_num_nodes); i++) {
			local_x[i][0] = fe->x[ fe->conn[e][i] ][0];
			local_x[i][1] = fe->x[ fe->conn[e][i] ][1];
			local_x[i][2] = fe->x[ fe->conn[e][i] ][2];
		}

		for(int i=0; i<(fe->local_num_nodes); i++) {
			BBFE_elemmat_set_Jacobian_array(
					Jacobian_ip,
					basis->num_integ_points,
					e,
					fe);

			for(int p=0; p<(basis->num_integ_points); p++) {
				double x_ip[3];
				BBFE_std_mapping_integ_point(
						x_ip,
						fe->local_num_nodes,
						local_x,
						basis->N[p]);

				double rhs_ip = func(x_ip[0], x_ip[1], x_ip[2], t);
				val_ip[p] = rhs_ip * basis->N[p][i];
			}

			double integ_val = BBFE_std_integ_calc(
					basis->num_integ_points,
					val_ip,
					basis->integ_weight,
					Jacobian_ip);

			equiv_val[ fe->conn[e][i] ] += integ_val;
		}
	}

	BB_std_free_1d_double(val_ip,     basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x, fe->local_num_nodes, 3);
}


