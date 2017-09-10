# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This file is needed in order to be able to cimport functions into other Cython files
cdef extern from "healpix.h":
    int healpix_xy_to_ring(int hp, int Nside)
