#pragma once

#include "../FE_sys/FE_dataset.h"

/**********************************************************
 * steady
 **********************************************************/
double BBFE_elemmat_convdiff_mat_conv(
		const double  N_i,
		const double  grad_N_j[3],
		const double  a[3]);

double BBFE_elemmat_convdiff_mat_diff(
		const double  grad_N_i[3],
		const double  grad_N_j[3],
		const double  k);

double BBFE_elemmat_convdiff_vec_source(
		const double  N_i,
		const double  f);

double BBFE_elemmat_convdiff_stab_coef(
		const double k,
		const double a[3],
		const double h_e);

double BBFE_elemmat_convdiff_mat_stab_conv(
		const double  grad_N_i[3],
		const double  grad_N_j[3],
		const double  a[3],
		const double  tau);

double BBFE_elemmat_convdiff_vec_stab_source(
		const double  grad_N_i[3],
		const double  a[3],
		const double  tau,
		const double  f);

/**********************************************************
 * non-steady
 **********************************************************/
double BBFE_elemmat_convdiff_mat_mass(
		const double  N_i,
		const double  N_j,
		const double  a);

double BBFE_elemmat_convdiff_vec_mass(
		const double  N_i,
		const double  T,
		const double  a);
