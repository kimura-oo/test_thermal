#pragma once

/**********************************************************
 * solid linear 
 **********************************************************/
void BBFE_elemmat_solid_mat_dispstr_linear(
		double       mat[6][3],
		const double grad_N[3]);

void BBFE_elemmat_solid_mat_Hooke(
		double       mat[6][6],
		const double e,  /* Young's mudulus */
		const double v); /* Poisson's ratio */

void BBFE_elemmat_solid_mat_linear(
		double       mat[3][3],
		const double grad_N_i[3],
		const double grad_N_j[3],
		const double e,
		const double v);

/**********************************************************
 * solid nonlinear 
 **********************************************************/
void BBFE_elemmat_solid_tensor_defgrad(
		double       mat[3][3],
		double**     grad_u);

void BBFE_elemmat_solid_mat_dispstr_tl(
		double       mat[6][3],
		double**     grad_u,
		const double grad_N[3]);

void BBFE_elemmat_solid_tensor_Green_Lagrange_mat_notation(
		double    mat[6],
		double    F[3][3]);

void BBFE_elemmat_solid_tensor_second_Piora_Kirchhoff_mat_notation(
		double       mat[6],
		double       D[6][6],
		const double E[6]);

void BBFE_elemmat_solid_mat_tl(
		double       mat[3][3],
		double       D[6][6],
		const double grad_N_i[3],
		const double grad_N_j[3],
		double**     grad_u );

void BBFE_elemmat_solid_vec_inner_force_tl(
		double       vec[3], 
		double       D[6][6],
		const double grad_N_i[3],
		double**     grad_u);

