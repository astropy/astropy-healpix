import numpy as np
cimport numpy as np
import cython

ctypedef np.double_t DOUBLE_T

from _healpix cimport healpix_xy_to_ring

def test():
    return healpix_xy_to_ring(1, 2)
