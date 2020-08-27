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
