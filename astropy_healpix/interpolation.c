#include <math.h>

#include "interpolation.h"
#include "healpix.h"

// Old versions of MSVC do not support C99 and therefore
// do not define NAN in math.h.
#ifndef NAN
static const union {
    unsigned long integer;
    float value;
} type_punned_nan = {0xFFFFFFFFFFFFFFFFul};
#define NAN (type_punned_nan.value)
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void interpolate_weights(double lon, double lat, int64_t *ring_indices,
                         double *weights, int Nside) {

  // Given a longitude and a latitude, Nside, and pre-allocated arrays of 4
  // elements ring_indices and weights, find the ring index of the four nearest
  // neighbours and the weights to use for each neighbour to interpolate.

  int64_t xy_index, npix;
  int64_t ring1, ring2, ring3, ring4;
  double lon1, lat1, lon2, lat2;
  double lon3, lat3, lon4, lat4;
  double xfrac1, xfrac2, yfrac, lon_frac;
  int ring_number, longitude_index, n_in_ring;

  // Find the xy index of the pixel in which the coordinates fall
  xy_index = radec_to_healpixl(lon, lat, Nside);

  // Find the lon/lat of the center of that pixel
  healpixl_to_radec(xy_index, Nside, 0.5, 0.5, &lon1, &lat1);

  // Take into account possible wrapping so that the pixel longitude/latitude
  // are close to the requested longitude/latitude
  if (lon - lon1 > M_PI)
    lon1 += 2 * M_PI;
  if (lon1 - lon > M_PI)
    lon1 -= 2 * M_PI;

  // Convert to a ring index and decompose into ring number and longitude index
  ring1 = healpixl_xy_to_ring(xy_index, Nside);
  if (ring1 < 0)
  {
      int i;
      for (i = 0; i < 4; i ++)
      {
          ring_indices[i] = -1;
          weights[i] = NAN;
      }
      return;
  }
  healpixl_decompose_ring(ring1, Nside, &ring_number, &longitude_index);

  // Figure out how many pixels are in the ring
  if (ring_number < Nside) {
    n_in_ring = 4 * ring_number;
  } else if (ring_number < 3 * Nside) {
    n_in_ring = 4 * Nside;
  } else {
    n_in_ring = (int)(4 * (4 * (int64_t)Nside - (int64_t)ring_number));
  }

  // We now want to find the next index in the ring so that the point to
  // interpolate is between the two. First we check what direction to look in by
  // finding the longitude/latitude of the center of the HEALPix pixel.

  if (lon < lon1) { // Go to the left
    if (longitude_index == 0) {
      ring2 = ring1 + n_in_ring - 1;
    } else {
      ring2 = ring1 - 1;
    }
  } else { // Go to the right
    if (longitude_index == n_in_ring - 1) {
      ring2 = ring1 - n_in_ring + 1;
    } else {
      ring2 = ring1 + 1;
    }
  }

  // Find the lon/lat of the new pixel
  xy_index = healpixl_ring_to_xy(ring2, Nside);
  healpixl_to_radec(xy_index, Nside, 0.5, 0.5, &lon2, &lat2);

  // Take into account possible wrapping so that the pixel longitude/latitude
  // are close to the requested longitude/latitude
  if (lon - lon2 > M_PI)
    lon2 += 2 * M_PI;
  if (lon2 - lon > M_PI)
    lon2 -= 2 * M_PI;

  // Now check whether we are moving up or down in terms of ring index

  if (lat > lat1) { // Move up (0 index is at the top)
    ring_number -= 1;
  } else { // Move down
    ring_number += 1;
  }

  if (ring_number > 0 && ring_number < 4 * Nside) {

    // Now figure out again how many pixels are in the ring

    if (ring_number < Nside) {
      n_in_ring = 4 * ring_number;
    } else if (ring_number < 3 * Nside) {
      n_in_ring = 4 * Nside;
    } else {
      n_in_ring = (int)(4 * (4 * (int64_t)Nside - (int64_t)ring_number));
    }

    // Now determine the longitude index in which the requested longitude falls.

    // In all regions, the longitude elements are spaced by 360 / n_in_ring. For
    // convenience we convert the longitude index so that the spacing is 1.
    lon_frac = lon * n_in_ring / (2 * M_PI);

    // In the equatorial region, the first ring starts at 0.5 and the second at
    // 0 (in lon_frac space). The ring number is 1-based and the first ring in
    // the equatorial region is even. In this ring we can simply take
    // int(lon_frac) to get the longitude index but in the odd rings we need to
    // adjust lon_frac
    if (n_in_ring == 4 * Nside && ring_number % 2 == 1) { // Equatorial region
      lon_frac += 0.5;
    }

    // Find the longitude index of the closest pixel
    longitude_index = (int)lon_frac;
    if (longitude_index == n_in_ring) {
      longitude_index -= 1;
    }

    // Find the longitude/latitude and ring index of this pixel
    ring3 = healpixl_compose_ring(ring_number, longitude_index, Nside);
    xy_index = healpixl_ring_to_xy(ring3, Nside);
    healpixl_to_radec(xy_index, Nside, 0.5, 0.5, &lon3, &lat3);

    // Take into account possible wrapping so that the pixel longitude/latitude
    // are close to the requested longitude/latitude
    if (lon - lon3 > M_PI)
      lon3 += 2 * M_PI;
    if (lon3 - lon > M_PI)
      lon3 -= 2 * M_PI;

    // Finally we can find the fourth pixel as before

    if (lon < lon3) { // Go to the left
      if (longitude_index == 0) {
        ring4 = ring3 + n_in_ring - 1;
      } else {
        ring4 = ring3 - 1;
      }
    } else { // Go to the right
      if (longitude_index == n_in_ring - 1) {
        ring4 = ring3 - n_in_ring + 1;
      } else {
        ring4 = ring3 + 1;
      }
    }

    xy_index = healpixl_ring_to_xy(ring4, Nside);

    healpixl_to_radec(xy_index, Nside, 0.5, 0.5, &lon4, &lat4);

    // Take into account possible wrapping so that the pixel longitude/latitude
    // are close to the requested longitude/latitude
    if (lon - lon4 > M_PI)
      lon4 += 2 * M_PI;
    if (lon4 - lon > M_PI)
      lon4 -= 2 * M_PI;

    // Determine the interpolation weights

    xfrac1 = (lon - lon1) / (lon2 - lon1);
    xfrac2 = (lon - lon3) / (lon4 - lon3);
    yfrac = (lat - lat1) / (lat3 - lat1);

    weights[0] = (1 - xfrac1) * (1 - yfrac);
    weights[1] = xfrac1 * (1 - yfrac);
    weights[2] = (1 - xfrac2) * yfrac;
    weights[3] = xfrac2 * yfrac;

  } else {

    // In the case where we are inside the four top/bottom-most
    // values, we effectively place a value at the pole that
    // is the average of the four values, and the interpolation
    // is the weighted average of this polar value and the
    // value interpolated along the ring.

    xfrac1 = (lon - lon1) / (lon2 - lon1);
    yfrac = (lat - lat1) / (0.5 * M_PI - lat1);

    if (ring_number == 0) {
      ring3 = (ring1 + 2) % 4;
      ring4 = (ring2 + 2) % 4;
      yfrac = (lat - lat1) / (0.5 * M_PI - lat1);
    } else {
      npix = 12 * (int64_t)Nside * (int64_t)Nside;
      ring3 = ((ring1 - (npix - 4)) + 2) % 4 + npix - 4;
      ring4 = ((ring2 - (npix - 4)) + 2) % 4 + npix - 4;
      yfrac = (lat - lat1) / (-0.5 * M_PI - lat1);
    }

    weights[0] = (1 - xfrac1) * (1 - yfrac) + 0.25 * yfrac;
    weights[1] = xfrac1 * (1 - yfrac) + 0.25 * yfrac;
    weights[2] = 0.25 * yfrac;
    weights[3] = 0.25 * yfrac;
  }

  ring_indices[0] = ring1;
  ring_indices[1] = ring2;
  ring_indices[2] = ring3;
  ring_indices[3] = ring4;
}
