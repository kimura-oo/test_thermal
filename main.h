#pragma once

#include <stdio.h>
#include <stdbool.h>

#include "monolis.h"

#include "libBB/std.h"
#include "libBB/calc.h"
#include "libBB/vtk.h"

#include "FE_std/integ.h"
#include "FE_std/shapefunc.h"
#include "FE_std/mapping.h"
#include "FE_std/surface.h"

#include "FE_sys/FE_dataset.h"
#include "FE_sys/memory.h"
#include "FE_sys/read.h"
#include "FE_sys/write.h"
#include "FE_sys/monowrap.h"

#include "FE_manusol/manusol.h"


typedef struct
{
	double* T;
	double* error;

} NODAL_VALUES;


/**********************************************************
 * equivval
 **********************************************************/

void BBFE_elemmat_equivval_volume_smooth_function(
		double* equiv_val,
		FE_DATA* fe,
		FE_3D_BASIS* basis,
		double (*func)(double, double, double)); // scalar function(x, y, z)


/**********************************************************
 * element matrix
 **********************************************************/
void BBFE_elemmat_set_Jacobi_mat(
		FE_DATA*     fe,
		FE_3D_BASIS* basis);

void BBFE_elemmat_set_shapefunc_derivative(
		FE_DATA*      fe,
		FE_3D_BASIS*  basis);

void BBFE_elemmat_set_element_mat(
		MONOLIS*     monolis,
		FE_DATA*     fe,
		FE_3D_BASIS* basis);

double BBFE_elemmat_thermal_steady_linear(
		double grad_N_i[3],
		double grad_N_j[3]);

