
#include "convdiff_surface.h"


static const char* INPUT_FILENAME_SURF = "surf.dat";


void BBFE_convdiff_pre_surface(
		BBFE_DATA*    surf,
		BBFE_BASIS*   basis,
		const char*   directory,
		int           num_integ_points_each_axis)
{
	int n_axis = num_integ_points_each_axis;
	const char* filename;

	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_SURF);
	BBFE_sys_read_elem(
			surf,
			filename,
			directory,
			n_axis*n_axis);

	BBFE_sys_memory_allocation_integ(
			basis,
			n_axis*n_axis*n_axis,
			2);
	BBFE_sys_memory_allocation_shapefunc(
			basis,
			surf->local_num_nodes,
			1,
			n_axis*n_axis*n_axis);

	BBFE_convdiff_set_basis_surface(
			basis,
			surf->local_num_nodes,
			n_axis);
}


void BBFE_convdiff_set_basis_surface(
		BBFE_BASIS*   basis,
		int           local_num_nodes,
		int           num_integ_points_each_axis)
{
	switch( local_num_nodes ) {
		case 3:
			basis->num_integ_points = 
				BBFE_std_integ_tri_set_arbitrary_points(
						num_integ_points_each_axis,
						basis->integ_point,
						basis->integ_weight);

			for(int i=0; i<(basis->num_integ_points); i++) {
				BBFE_std_shapefunc_tri1st_get_val(
						basis->integ_point[i],
						basis->N[i]);

				BBFE_std_shapefunc_tri1st_get_derivative(
						basis->integ_point[i],
						basis->dN_dxi[i],
						basis->dN_det[i]);
			}
			printf("%s Surface element type: 1st-order triangle.\n", CODENAME);
			break;

		case 4:
			basis->num_integ_points = 
				BBFE_std_integ_rec_set_arbitrary_points(
						num_integ_points_each_axis,
						basis->integ_point,
						basis->integ_weight);

			for(int i=0; i<(basis->num_integ_points); i++) {
				BBFE_std_shapefunc_rec1st_get_val(
						basis->integ_point[i],
						basis->N[i]);

				BBFE_std_shapefunc_rec1st_get_derivative(
						basis->integ_point[i],
						basis->dN_dxi[i],
						basis->dN_det[i]);
			}
			printf("%s Surface element type: 1st-order rectangle.\n", CODENAME);
			break;

		case 6:
		case 9:
			// should be implemented for higher order elements
			break;
	}
}


static void calc_solution_gradient(
		const double x[3],
		const double t,
		const double delta,
		double       (*func)(double, double, double, double),
		double       grad[3]) 
{
	double invD  = 1.0/delta;

	grad[0] =  ( func(x[0]+delta, x[1], x[2], t) -
		func(x[0]-delta, x[1], x[2], t) ) * (invD/2.0);
	grad[1] =  ( func(x[0], x[1]+delta, x[2], t) -
		func(x[0], x[1]-delta, x[2], t) ) * (invD/2.0);
	grad[2] =  ( func(x[0], x[1], x[2]+delta, t) -
		func(x[0], x[1], x[2]-delta, t) ) * (invD/2.0);
}


void BBFE_convdiff_set_equiv_val_Neumann(
		double*     equiv_val,
		BBFE_DATA*  surf,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		double      t,
		double      delta,
		double      (*func)(double, double, double, double)) // scalar function(x, y, z, t)
{
	double* val_ip;
	double* grad_normal_ip;
	double* Jacobian_ip;
	val_ip         = BB_std_calloc_1d_double(val_ip,         basis->num_integ_points);
	grad_normal_ip = BB_std_calloc_1d_double(grad_normal_ip, basis->num_integ_points);
	Jacobian_ip    = BB_std_calloc_1d_double(Jacobian_ip,    basis->num_integ_points);
	
	// Jacobian values here are dummy values (all values are 1.0) 
	// because they are not needed for surface integration
	for(int p=0; p<(basis->num_integ_points); p++) {
		Jacobian_ip[p] = 1.0;
	}
	
	double** local_x;
	local_x = BB_std_calloc_2d_double(local_x, surf->local_num_nodes, 3);

	for(int e=0; e<(surf->total_num_elems); e++) {
		BBFE_elemmat_set_local_array_vector(local_x, surf, fe->x, e, 3);

		for(int p=0; p<(basis->num_integ_points); p++) {
			double normal[3];
			BBFE_set_surface_get_outward_normal_vector(
					surf->local_num_nodes,
					local_x,
					basis->dN_dxi[p],
					basis->dN_det[p],
					normal);

  			double grad[3];
			double x_ip[3];
			BBFE_std_mapping_vector3d(
					x_ip,
					surf->local_num_nodes,
					local_x,
					basis->N[p]);
			calc_solution_gradient(x_ip, t, delta, func, grad);
			
			grad_normal_ip[p] = BB_calc_vec3d_dot(grad, normal);
		}

		for(int i=0; i<(surf->local_num_nodes); i++) {
			for(int p=0; p<(basis->num_integ_points); p++) {
				val_ip[p] = basis->N[p][i] * grad_normal_ip[p];
			}

			double integ_val;
			integ_val = BBFE_std_integ_calc(
				basis->num_integ_points,
				val_ip,
				basis->integ_weight,
				Jacobian_ip);
			
			equiv_val[ surf->conn[e][i] ] -= integ_val;
		}
	}

	BB_std_free_1d_double(val_ip,         basis->num_integ_points);
	BB_std_free_1d_double(grad_normal_ip, basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip,    basis->num_integ_points);
	
	BB_std_free_2d_double(local_x,        surf->local_num_nodes, 3);
}
