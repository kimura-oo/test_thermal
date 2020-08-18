
void BBFE_std_mapping_calc_Jacobi_mat_3d(
		double     J[3][3],
		const int  local_num_nodes,
		double**   local_x,
		double*    local_dN_dxi,
		double*    local_dN_det,
		double*    local_dN_dze)
{
	for(int i=0; i<3; i++) {
		for(int j=0; j<3; j++) {
			J[i][j] = 0.0;
		}
	}

	for(int n=0; n<local_num_nodes; n++) {
		for(int i=0; i<3; i++) {
			J[0][i] += local_dN_dxi[n] * local_x[n][i];
			J[1][i] += local_dN_det[n] * local_x[n][i];
			J[2][i] += local_dN_dze[n] * local_x[n][i];
		}
	}
}


double BBFE_std_mapping_scalar(
		const int     local_num_nodes,
		double*       local_val,
		double*       N)
{
	double val_ip = 0.0;

	for(int i=0; i<local_num_nodes; i++) {
		val_ip += N[i]*local_val[i];
	}

	return val_ip;
}


void BBFE_std_mapping_vector3d(
		double        val_ip[3],
		const int     local_num_nodes,
		double**      local_val,
		double*       N)
{
	for(int d=0; d<3; d++) {
		val_ip[d] = 0.0;
	}

	for(int i=0; i<local_num_nodes; i++) {
		for(int d=0; d<3; d++) {
			val_ip[d] += N[i]*local_val[i][d];
		}
	}
}
