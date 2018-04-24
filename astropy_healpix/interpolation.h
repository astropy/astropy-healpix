#ifndef ASTROPY_HEALPIX_INTERPOLATION_INCL
#define ASTROPY_HEALPIX_INTERPOLATION_INCL

#ifdef _MSC_VER
#if _MSC_VER >= 1600
#include <stdint.h>
#else
#include <stdint_msc.h>
#endif
#else
#include <stdint.h>
#endif

void interpolate_weights(double lon, double lat, int64_t *ring_indices,
                         double *weights, int Nside);

#endif
