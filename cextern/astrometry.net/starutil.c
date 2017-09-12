/*
# This file is part of the Astrometry.net suite.
# Licensed under a 3-clause BSD style license - see LICENSE
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>

#include "os-features.h"
#include "keywords.h"
#include "mathutil.h"
#include "starutil.h"

#define POGSON 2.51188643150958
#define LOGP   0.92103403719762

#define InlineDefine InlineDefineC
#include "starutil.inc"
#undef InlineDefine
