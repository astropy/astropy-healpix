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
#include "errors.h"

#define POGSON 2.51188643150958
#define LOGP   0.92103403719762

#define InlineDefine InlineDefineC
#include "starutil.inc"
#undef InlineDefine

void radec_derivatives(double ra, double dec, double* dra, double* ddec) {
	double cosd = cos(deg2rad(dec));
	double cosra = cos(deg2rad(ra));
	double sinra = sin(deg2rad(ra));
	if (dra) {
		dra[0] = cosd * -sinra;
		dra[1] = cosd *  cosra;
		dra[2] = 0.0;
		normalize_3(dra);
	}
	if (ddec) {
		double sind = sin(deg2rad(dec));
		ddec[0] = -sind * cosra;
		ddec[1] = -sind * sinra;
		ddec[2] =  cosd;
		normalize_3(ddec);
	}
}

void radecrange2xyzrange(double ralo, double declo, double rahi, double dechi,
						 double* minxyz, double* maxxyz) {
	double minmult, maxmult;
	double uxlo, uxhi, uylo, uyhi;
	// Dec only affects z, and is monotonic (z = sin(dec))
	minxyz[2] = radec2z(0, declo);
	maxxyz[2] = radec2z(0, dechi);

	// min,max of cos(dec).  cos(dec) is concave down.
	minmult = MIN(cos(deg2rad(declo)), cos(deg2rad(dechi)));
	maxmult = MAX(cos(deg2rad(declo)), cos(deg2rad(dechi)));
	if (declo < 0 && dechi > 0)
		maxmult = 1.0;
	// unscaled x (ie, cos(ra))
	uxlo = MIN(cos(deg2rad(ralo)), cos(deg2rad(rahi)));
	if (ralo < 180 && rahi > 180)
		uxlo = -1.0;
	uxhi = MAX(cos(deg2rad(ralo)), cos(deg2rad(rahi)));
	minxyz[0] = MIN(uxlo * minmult, uxlo * maxmult);
	maxxyz[0] = MAX(uxhi * minmult, uxhi * maxmult);
	// unscaled y (ie, sin(ra))
	uylo = MIN(sin(deg2rad(ralo)), sin(deg2rad(rahi)));
	if (ralo < 270 && rahi > 270)
		uylo = -1.0;
	uyhi = MAX(sin(deg2rad(ralo)), sin(deg2rad(rahi)));
	if (ralo < 90 && rahi > 90)
		uyhi = -1.0;
	minxyz[1] = MIN(uylo * minmult, uylo * maxmult);
	maxxyz[1] = MAX(uyhi * minmult, uyhi * maxmult);
}


