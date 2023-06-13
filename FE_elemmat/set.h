#pragma once

#include "BBFE/sys/FE_dataset.h"
#include "monolis.h"

/**********************************************************
 * basic utils
 **********************************************************/
void BBFE_elemmat_set_Jacobian_array(
		double*    Jacobian_ip,
		const int  num_integ_points,
		const int  elem_num,
		BBFE_DATA* fe);

void BBFE_elemmat_set_local_array_scalar(
		double*       local_val,
		BBFE_DATA*    fe,
		const double* val,
		const int     elem_num);

void BBFE_elemmat_set_local_array_vector(
		double**       local_val,
		BBFE_DATA*     fe,
		double**       val,
		const int      elem_num,
		const int      dimension);

void BBFE_elemmat_set_Jacobi_mat(
		BBFE_DATA*  fe,
		BBFE_BASIS* basis);

void BBFE_elemmat_set_shapefunc_derivative(
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis);

/**********************************************************
 * setter of simple global matrix
 **********************************************************/
void BBFE_elemmat_set_global_mat_cmass_const(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		double      coef,
		int         block_size);

void BBFE_elemmat_set_global_mat_cmass_const_C(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		double      coef,
		int         block_size);

void BBFE_elemmat_set_global_mat_Laplacian_const(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		double      coef);

void BBFE_elemmat_set_global_mat_Laplacian_const_C(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		double      coef);
