#pragma once

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

#include "FE_elemmat/set.h"
#include "FE_elemmat/equivval.h"
#include "FE_elemmat/convdiff.h"

#include "FE_manusol/manusol.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

static const int BLOCK_SIZE                   = 1;
static const char* CODENAME                   = "convdiff >";

static const char* INPUT_FILENAME_NODE        = "node.dat";
static const char* INPUT_FILENAME_ELEM        = "elem.dat";
static const char* INPUT_FILENAME_D_BC        = "D_bc.dat";
static const char* OUTPUT_FILENAME_VTK        = "result.vtk";
static const char* OUTPUT_FILENAME_ASCII_TEMP = "temparature.dat";
static const char* OUTPUT_FILENAME_ASCII_RHS  = "rhs.dat";


typedef struct
{
	double* T;
	double* error;
	double* theo_sol;

} NODAL_VALUES;


typedef struct
{
	const char* directory;

} CONDITIONS;


typedef struct
{
	FE_3D_BASIS  basis;
	FE_DATA      fe;
	NODAL_VALUES vals;
	BC_DATA      bc;
	MONOLIS      monolis;
	CONDITIONS   cond;

} FE_SYSTEM;


const char* BBFE_convdiff_get_directory_name(
		int         argc,
		char*       argv[],
		const char* codename);

void BBFE_convdiff_pre(
		FE_SYSTEM* sys,
		int        argc,
		char*      argv[],
		int        num_integ_points_each_axis,
		bool       manufactured_solution);

void BBFE_convdiff_memory_allocation_nodal_values(
		NODAL_VALUES*   vals,
		const int       total_num_nodes);

void BBFE_convdiff_set_basis(
		FE_3D_BASIS*  basis,
		int           local_num_nodes,
		int           num_integ_points_each_axis);
