# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This file is needed in order to be able to cimport functions into other Cython
# files

cdef extern from "healpix.h":

    # Converts a healpix index from the XY scheme to the RING scheme.
    int healpix_xy_to_ring(int hp, int Nside);

    # Converts a healpix index from the RING scheme to the XY scheme.
    int healpix_ring_to_xy(int ring_index, int Nside);

    # Converts a healpix index from the XY scheme to the NESTED scheme.
    int healpix_xy_to_nested(int hp, int Nside);

    # Converts a healpix index from the NESTED scheme to the XY scheme.
    int healpix_nested_to_xy(int nested_index, int Nside);

    #  Converts a healpix index, plus fractional offsets (dx,dy), into (x,y,z)
    #  coordinates on the unit sphere.  (dx,dy) must be in [0, 1].  (0.5, 0.5)
    #  is the center of the healpix.  (0,0) is the southernmost corner, (1,1) is
    #  the northernmost corner, (1,0) is the easternmost, and (0,1) the
    #  westernmost.
    void healpix_to_radec(int hp, int Nside, double dx, double dy,
                          double* ra, double* dec);

    # Converts (RA, DEC) coordinates (in radians) to healpix XY index
    int radectohealpixf(double ra, double dec, int Nside, double* dx, double* dy);
