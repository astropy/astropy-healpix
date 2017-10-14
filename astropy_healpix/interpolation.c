#include "interpolation.h"
#include "healpix.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void interpolate_weights(double lon, double lat, int64_t *ring_indices,
                         double *weights, int Nside) {

  int64_t xy_index, ring_index;

  int64_t ring1, ring2;
  double lon1, lat1, lon2, lat2;

  int64_t ring3, ring4;
  double lon3, lat3, lon4, lat4;

  double lon_tmp, lat_tmp;

  int ring_number, longitude_index;

  double xfrac1, xfrac2, yfrac;
  double x0, y0, dlon, frac;

  ring_neighbours(lon, lat, Nside, &ring1, &ring2, &lon1, &lat1, &lon2, &lat2);

  healpixl_decompose_ring(ring1, Nside, &ring_number, &longitude_index);

  if (lat > lat1) {
    ring_number -= 1;
  } else {
    ring_number += 1;
  }

  if (ring_number == 0) {
    xfrac1 = (lon - lon1) / (lon2 - lon1);
    yfrac = (lat - lat1) / (0.5 * M_PI - lat1);
    weights[0] = (1 - xfrac1) * (1 - yfrac) + 0.25 * yfrac;
    weights[1] = xfrac1 * (1 - yfrac) + 0.25 * yfrac;
    weights[2] = 0.25 * yfrac;
    weights[3] = 0.25 * yfrac;
    ring3 = (ring1 + 2) % 4;
    ring4 = (ring2 + 2) % 4;
  } else {
    ring_index = healpixl_compose_ring(ring_number, 1, Nside);

    xy_index = healpixl_ring_to_xy(ring_index, Nside);
    healpixl_to_radec(xy_index, Nside, 0.5, 0.5, &lon_tmp, &lat_tmp);

    dlon = 0.1 * (lon2 - lon1);

    ring_neighbours(lon + dlon, lat_tmp, Nside, &ring3, &ring4, &lon3, &lat3,
                    &lon4, &lat4);

    xfrac1 = (lon - lon1) / (lon2 - lon1);
    xfrac2 = (lon - lon3) / (lon4 - lon3);
    yfrac = (lat - lat1) / (lat3 - lat1);

    weights[0] = (1 - xfrac1) * (1 - yfrac);
    weights[1] = xfrac1 * (1 - yfrac);
    weights[2] = (1 - xfrac2) * yfrac;
    weights[3] = xfrac2 * yfrac;
  }

  ring_indices[0] = ring1;
  ring_indices[1] = ring2;
  ring_indices[2] = ring3;
  ring_indices[3] = ring4;
}

void ring_neighbours(double lon, double lat, int Nside, int64_t *ring1,
                     int64_t *ring2, double *lon1, double *lat1, double *lon2,
                     double *lat2) {

  int64_t xy_index;
  int ring_number, longitude_index, n_in_ring;
  double lon_tmp;
  int64_t ring_tmp;

  // Find the xy index of the pixel in which the coordinates fall
  xy_index = radec_to_healpixl(lon, lat, Nside);

  // Find the lon/lat of the center of that pixel
  healpixl_to_radec(xy_index, Nside, 0.5, 0.5, lon1, lat1);

  // Take into account possible wrapping so that the pixel longitude/latitude
  // are close to the requested longitude/latitude
  if (lon - *lon1 > M_PI)
    *lon1 += 2 * M_PI;
  if (*lon1 - lon > M_PI)
    *lon1 -= 2 * M_PI;

  // Also convert to a ring index and decompose into ring number and longitude
  // index
  *ring1 = healpixl_xy_to_ring(xy_index, Nside);
  healpixl_decompose_ring(*ring1, Nside, &ring_number, &longitude_index);

  // Figure out how many pixels are in the ring
  if (ring_number < Nside) {
    n_in_ring = 4 * ring_number;
  } else if (ring_number < 3 * Nside) {
    n_in_ring = 4 * Nside;
  } else {
    n_in_ring = 4 * (4 * Nside - ring_number);
  }

  // We now want to find the next index in the ring so that the point to
  // interpolate is between the two. First we check what direction to look in by
  // finding the longitude/latitude of the center of the HEALPix pixel.

  if (lon < *lon1) { // Go to the left
    if (longitude_index == 0) {
      *ring2 = *ring1 + n_in_ring - 1;
    } else {
      *ring2 = *ring1 - 1;
    }
  } else { // Go to the right
    if (longitude_index == n_in_ring - 1) {
      *ring2 = *ring1 - n_in_ring + 1;
    } else {
      *ring2 = *ring1 + 1;
    }
  }

  // Find the lon/lat of the new pixel
  xy_index = healpixl_ring_to_xy(*ring2, Nside);
  healpixl_to_radec(xy_index, Nside, 0.5, 0.5, lon2, lat2);

  // Take into account possible wrapping so that the pixel longitude/latitude
  // are close to the requested longitude/latitude
  if (lon - *lon2 > M_PI)
    *lon2 += 2 * M_PI;
  if (*lon2 - lon > M_PI)
    *lon2 -= 2 * M_PI;

  if (*lon1 > *lon2) {
    lon_tmp = *lon2;
    *lon2 = *lon1;
    *lon1 = lon_tmp;
    ring_tmp = *ring2;
    *ring2 = *ring1;
    *ring1 = ring_tmp;
  }
}
