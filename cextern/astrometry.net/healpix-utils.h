/*
# This file is part of the Astrometry.net suite.
# Licensed under a 3-clause BSD style license - see LICENSE
 */

#ifndef HEALPIX_UTILS_H
#define HEALPIX_UTILS_H

#include "bl.h"

/**
 Returns healpixes that are / may be within range of the given point, resp.
 */
il* healpix_rangesearch_xyz(const double* xyz, double radius, int64_t Nside, il* hps);
il* healpix_rangesearch_xyz_approx(const double* xyz, double radius, int64_t Nside, il* hps);
il* healpix_rangesearch_radec_approx(double ra, double dec, double radius, int64_t Nside, il* hps);
il* healpix_rangesearch_radec(double ra, double dec, double radius, int64_t Nside, il* hps);
int64_t healpix_rangesearch_radec_simple(double ra, double dec, double radius, int64_t Nside, int64_t **indices);

/**
 Starting from a "seed" or list of "seeds" healpixes, grows a region
 by looking at healpix neighbours.  Accepts healpixes for which the
 "accept" function returns 1.  Returns the healpixes that are
 accepted.  The accepted results are placed in "accepted", if
 non-NULL, or in a newly-allocated list.

 If "rejected" is non-NULL, the healpixes that are rejected will be
 put there.

 If "depth" is non-zero, that number of neighbour steps will be taken.
 Zero means no limit.

 NOTE that any existing entries in the "accepted" list will be treated
 as having already been accepted: when the search reaches them, their
 neighbours will not be added to the frontier to explore.
 */
il* healpix_region_search(int64_t seed, il* seeds, int64_t Nside,
						  il* accepted, il* rejected,
						  int64_t (*accept)(int64_t hp, void* token),
						  void* token,
						  int64_t depth);


#endif
