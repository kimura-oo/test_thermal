#pragma once

// rule: BBFE_std_mapping_

void BBFE_std_mapping_void();

void BBFE_std_mapping_calc_Jacobi_mat_3d(
		double     J[3][3],
		const int  local_num_nodes,
		double**   local_x,
		double*    local_dN_dxi,
		double*    local_dN_det,
		double*    local_dN_dze);

double BBFE_std_mapping_scalar(
		const int     local_num_nodes,
		double*       local_val,
		double*       N);

void BBFE_std_mapping_vector3d(
		double        val_ip[3],
		const int     local_num_nodes,
		double**      local_val,
		double*       N);

void BBFE_std_mapping_scalar_grad(
		double        val_ip[3],
		const int     local_num_nodes,
		double*       local_val,
		double**      grad_N);

void BBFE_std_mapping_vector3d_grad(
		double**      val_ip,   //[3][3]
		const int     local_num_nodes,
		double**      local_val,
		double**      grad_N);

double BBFE_std_mapping_vector3d_div(
		const int     local_num_nodes,
		double**      local_val,
		double**      grad_N);
