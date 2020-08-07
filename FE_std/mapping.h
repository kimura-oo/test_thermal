#pragma once

// rule: BBFE_std_mapping_
void BBFE_std_mapping_calc_Jacobi_mat_3d(
		double     J[3][3],
		const int  local_num_nodes,
		double**   local_x,
		double*    local_dN_dxi,
		double*    local_dN_det,
		double*    local_dN_dze);


void BBFE_std_mapping_integ_point(
		double        x[3],
		const int     local_num_nodes,
		double**      local_x,
		double*       N);
