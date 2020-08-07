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

#include "FE_elemmat/set.h"
#include "FE_elemmat/equivval.h"
#include "FE_elemmat/thermal.h"

#include "FE_manusol/manusol.h"
