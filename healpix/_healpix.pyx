import numpy as np
cimport numpy as np
import cython

ctypedef np.double_t DOUBLE_T

cdef extern from "healpix.h":
    int healpix_xy_to_ring(int hp, int Nside)


def heapix_xy_to_ring():
    return None
