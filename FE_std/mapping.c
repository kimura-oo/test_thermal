
void BBFE_std_mapping_void(){}

		
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


void BBFE_std_mapping_scalar_grad(
		double        val_ip[3],
		const int     local_num_nodes,
		double*       local_val,
		double**      grad_N)
{

	for(int i=0; i<3; i++) {
		val_ip[i] = 0.0;

		for(int k=0; k<local_num_nodes; k++) {
			val_ip[i] += local_val[k] * grad_N[k][i];
		}
	}
}


void BBFE_std_mapping_vector3d_grad(
		double**      val_ip,
		const int     local_num_nodes,
		double**      local_val,
		double**      grad_N)
{

	for(int i=0; i<3; i++) {
		for(int j=0; j<3; j++) {
			val_ip[i][j] = 0.0;

			for(int k=0; k<local_num_nodes; k++) {
				val_ip[i][j] += local_val[k][i] * grad_N[k][j];
			}
		}
	}
}


double BBFE_std_mapping_vector3d_div(
		const int     local_num_nodes,
		double**      local_val,
		double**      grad_N)
{
	double val_ip;

	double dval_dx[3];

	for(int d=0; d<3; d++) {
		dval_dx[d] = 0.0;
		for(int i=0; i<local_num_nodes; i++) {
			dval_dx[d] += local_val[i][d] * grad_N[i][d];
		}
	}

	val_ip = dval_dx[0] + dval_dx[1] + dval_dx[2];

	return val_ip;
}


