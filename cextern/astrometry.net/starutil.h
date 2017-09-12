/*
# This file is part of the Astrometry.net suite.
# Licensed under a 3-clause BSD style license - see LICENSE
*/

#ifndef STARUTIL_H
#define STARUTIL_H

#include <math.h>
#include "an-bool.h"
#include "keywords.h"

#define DIM_STARS 3
#define DIM_XY 2

// upper bound of dimquads value
#define DQMAX 5
// upper bound of dimcodes value
#define DCMAX 6

InlineDeclare int dimquad2dimcode(int dimquad);

typedef unsigned char uchar;

#define ONE_OVER_SIXTY 0.016666666666666666

// pi / 180.
#define RAD_PER_DEG 0.017453292519943295
// pi / (180. * 60.)
#define RAD_PER_ARCMIN 0.00029088820866572158
// pi / (180. * 60. * 60.)
#define RAD_PER_ARCSEC 4.8481368110953598e-06

// 180. / pi
#define DEG_PER_RAD 57.295779513082323
#define DEG_PER_ARCMIN ONE_OVER_SIXTY
// 1./3600.
#define DEG_PER_ARCSEC 0.00027777777777777778

// 60. * 180. / pi
#define ARCMIN_PER_RAD 3437.7467707849396
#define ARCMIN_PER_DEG 60.0
#define ARCMIN_PER_ARCSEC ONE_OVER_SIXTY

// 60. * 60. * 180. / pi
#define ARCSEC_PER_RAD 206264.80624709636
#define ARCSEC_PER_DEG 3600.0
#define ARCSEC_PER_ARCMIN 60.0

InlineDeclare Const double rad2deg(double x);
InlineDeclare Const double rad2arcmin(double x);
InlineDeclare Const double rad2arcsec(double x);

InlineDeclare Const double deg2rad(double x);
InlineDeclare Const double deg2arcmin(double x);
InlineDeclare Const double deg2arcsec(double x);

InlineDeclare Const double arcmin2rad(double x);
InlineDeclare Const double arcmin2deg(double x);
InlineDeclare Const double arcmin2arcsec(double x);

InlineDeclare Const double arcsec2rad(double x);
InlineDeclare Const double arcsec2deg(double x);
InlineDeclare Const double arcsec2arcmin(double x);

#define MJD_JD_OFFSET 2400000.5

InlineDeclare Const double mjdtojd(double mjd);
InlineDeclare Const double jdtomjd(double jd);

// RA,Dec in radians:
#define radec2x(r,d) (cos(d)*cos(r))
#define radec2y(r,d) (cos(d)*sin(r))
#define radec2z(r,d) (sin(d))
InlineDeclare Const double xy2ra(double x, double y);
InlineDeclare Const double z2dec(double z);

double atora(const char* str);
double atodec(const char* str);

double mag2flux(double mag);

/*
 RA,Dec in degrees.  RAs in range [0, 360], Decs in range [-90, 90].
 */
void radecrange2xyzrange(double ralow, double declow, double rahigh, double dechigh,
						 double* xyzlow, double* xyzhigh);

// RA,Dec in radians:
InlineDeclare void radec2xyz(double ra, double dec, double* x, double* y, double* z);
InlineDeclare Flatten void xyz2radec(double x, double y, double z, double *ra, double *dec);
InlineDeclare Flatten void xyzarr2radec(const double* xyz, double *ra, double *dec);
InlineDeclare void xyzarr2radecarr(const double* xyz, double *radec);
InlineDeclare void radec2xyzarr(double ra, double dec, double* p_xyz);
InlineDeclare void radec2xyzarrmany(double *ra, double *dec, double* xyz, int n);

// RA,Dec in degrees:
InlineDeclare void radecdeg2xyz(double ra, double dec, double* x, double* y, double* z);
InlineDeclare Flatten void xyzarr2radecdeg(const double* xyz, double *ra, double *dec);
InlineDeclare Flatten void xyzarr2radecdegarr(double* xyz, double *radec);
InlineDeclare void radecdeg2xyzarr(double ra, double dec, double* p_xyz);
InlineDeclare void radecdegarr2xyzarr(double* radec, double* xyz);
InlineDeclare void radecdeg2xyzarrmany(double *ra, double *dec, double* xyz, int n);

#ifdef INCLUDE_INLINE_SOURCE
#define InlineDefine InlineDefineH
#include "starutil.inc"
#undef InlineDefine
#endif

#endif
