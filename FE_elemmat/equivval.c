
#include "equivval.h"
#include "BB/std.h"
#include "set.h"
#include "BBFE/std/mapping.h"
#include "BBFE/std/integ.h"

#include <math.h>


void BBFE_elemmat_equivval_volume_smooth_function(
		double*      equiv_val,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
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
				BBFE_std_mapping_vector3d(
						x_ip,
						fe->local_num_nodes,
						local_x,
						basis->N[p]);

				double func_ip = func(x_ip[0], x_ip[1], x_ip[2], t);
				val_ip[p] = func_ip * basis->N[p][i];
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


double BBFE_elemmat_equivval_relative_L2_error_scalar(
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
		MONOLIS_COM*  monolis_com,
		double        t,
		const double* comp_vec, // [total_num_nodes]
		double        (*func)(double, double, double, double)) // scalar function(x, y, z, t)
{
	double L2_abs_error = 0.0;
	double L2_abs_theo  = 0.0;

	double* val_ip_error;
	double* val_ip_theo;
	double* Jacobian_ip;
	bool* is_internal_elem;
	val_ip_error = BB_std_calloc_1d_double(val_ip_error, basis->num_integ_points);
	val_ip_theo  = BB_std_calloc_1d_double(val_ip_theo , basis->num_integ_points);
	Jacobian_ip  = BB_std_calloc_1d_double(Jacobian_ip , basis->num_integ_points);
	is_internal_elem = BB_std_calloc_1d_bool(is_internal_elem , fe->total_num_elems);

	double** local_x;
	local_x   = BB_std_calloc_2d_double(local_x  , fe->local_num_nodes, 3);

	double* local_val;
	local_val = BB_std_calloc_1d_double(local_val, fe->local_num_nodes);

	monolis_get_bool_list_of_internal_simple_mesh(monolis_com, fe->total_num_nodes, fe->total_num_elems,
		fe->local_num_nodes, fe->conn, is_internal_elem);

	for(int e=0; e<(fe->total_num_elems); e++) {
		if (!is_internal_elem[e]) continue;

		for(int i=0; i<(fe->local_num_nodes); i++) {

			local_val[i] = comp_vec[ fe->conn[e][i] ];

			local_x[i][0] = fe->x[ fe->conn[e][i] ][0];
			local_x[i][1] = fe->x[ fe->conn[e][i] ][1];
			local_x[i][2] = fe->x[ fe->conn[e][i] ][2];
		}

		BBFE_elemmat_set_Jacobian_array(
				Jacobian_ip,
				basis->num_integ_points,
				e,
				fe);

		for(int p=0; p<(basis->num_integ_points); p++) {
			double val_ip;
			val_ip = BBFE_std_mapping_scalar(
					fe->local_num_nodes,
					local_val,
					basis->N[p]);

			double x_ip[3];
			BBFE_std_mapping_vector3d(
					x_ip,
					fe->local_num_nodes,
					local_x,
					basis->N[p]);

			val_ip_error[p] = 
				pow( fabs( val_ip - func(x_ip[0], x_ip[1], x_ip[2], t) ), 2 );
			val_ip_theo[p] =
				pow( fabs( func(x_ip[0], x_ip[1], x_ip[2], t) ), 2 );
		}

		double integ_val_error = BBFE_std_integ_calc(
				basis->num_integ_points,
				val_ip_error,
				basis->integ_weight,
				Jacobian_ip);

		double integ_val_theo  = BBFE_std_integ_calc(
				basis->num_integ_points,
				val_ip_theo,
				basis->integ_weight,
				Jacobian_ip);

		L2_abs_error += integ_val_error;
		L2_abs_theo  += integ_val_theo;
	}

	BB_std_free_1d_double(val_ip_error, basis->num_integ_points);
	BB_std_free_1d_double(val_ip_theo,  basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip,  basis->num_integ_points);

	BB_std_free_2d_double(local_x,   fe->local_num_nodes, 3);
	BB_std_free_1d_double(local_val, fe->local_num_nodes);

	monolis_allreduce_R(
			1,
			&L2_abs_error,
			MONOLIS_MPI_SUM,
			monolis_com->comm);

	monolis_allreduce_R(
			1,
			&L2_abs_theo,
			MONOLIS_MPI_SUM,
			monolis_com->comm);

	return ( sqrt(L2_abs_error)/sqrt(L2_abs_theo) );
}

