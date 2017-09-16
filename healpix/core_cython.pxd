# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This file is needed in order to be able to cimport functions into other Cython
# files

from numpy cimport int64_t

cdef extern from "healpix.h":

    # Converts a healpix index from the XY scheme to the RING scheme.
    int64_t healpixl_xy_to_ring(int64_t hp, int Nside);

    # Converts a healpix index from the RING scheme to the XY scheme.
    int64_t healpixl_ring_to_xy(int64_t ring_index, int Nside);

    # Converts a healpix index from the XY scheme to the NESTED scheme.
    int64_t healpixl_xy_to_nested(int64_t hp, int Nside);

    # Converts a healpix index from the NESTED scheme to the XY scheme.
    int64_t healpixl_nested_to_xy(int64_t nested_index, int Nside);

    #  Converts a healpix index, plus fractional offsets (dx,dy), into (x,y,z)
    #  coordinates on the unit sphere.  (dx,dy) must be in [0, 1].  (0.5, 0.5)
    #  is the center of the healpix.  (0,0) is the southernmost corner, (1,1) is
    #  the northernmost corner, (1,0) is the easternmost, and (0,1) the
    #  westernmost.
    void healpixl_to_radec(int64_t hp, int Nside, double dx, double dy,
                           double* ra, double* dec);

    # Converts (RA, DEC) coordinates (in radians) to healpix XY index
    int64_t radec_to_healpixlf(double ra, double dec, int Nside, double* dx, double* dy);

    # Finds the healpixes neighbouring the given healpix, placing them in the
    # array "neighbour".  Returns the number of neighbours.  You must ensure
    # that "neighbour" has 8 elements.
    #
    # Healpixes in the interior of a large healpix will have eight neighbours;
    # pixels near the edges can have fewer.
    int healpixl_get_neighbours(int64_t hp, int64_t* neighbours, int Nside);

    int healpix_rangesearch_radec_simple(double ra, double dec, double radius, int Nside, int approx, int64_t** indices)
