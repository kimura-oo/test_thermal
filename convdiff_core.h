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


const char* BBFE_convdiff_get_directory_name(
		int         argc,
		char*       argv[],
		const char* codename);

void BBFE_convdiff_pre(
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
		BBFE_BC*      bc,
		MONOLIS*      monolis,
		int           argc,
		char*         argv[],
		const char*   directory,
		int           num_integ_points_each_axis,
		bool          manufactured_solution);

void BBFE_convdiff_set_basis(
		BBFE_BASIS*   basis,
		int           local_num_nodes,
		int           num_integ_points_each_axis);

void BBFE_convdiff_finalize(
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		BBFE_BC*     bc);
