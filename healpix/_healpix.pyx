import numpy as np
cimport numpy as np
import cython

ctypedef np.double_t DOUBLE_T


def nested_to_lonlat(np.ndarray[int, ndim=1, mode="c"] nested_index, int n_side):

    cdef int n = nested_index.shape[0]
    cdef double dx, dy;
    cdef int xy_index
    cdef np.ndarray[double, ndim=1, mode="c"] lon = np.zeros(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1, mode="c"] lat = np.zeros(n, dtype=np.double)

    dx = 0
    dy = 0

    for i in range(n):
        xy_index = healpix_nested_to_xy(nested_index[i], n_side)
        healpix_to_radec(xy_index, n_side, dx, dy, &lon[i], &lat[i])
    return lon, lat


def nested_with_offset_to_lonlat(np.ndarray[int, ndim=1, mode="c"] nested_index,
                                 np.ndarray[double, ndim=1, mode="c"] dx,
                                 np.ndarray[double, ndim=1, mode="c"] dy,
                                 int n_side):

    cdef int n = nested_index.shape[0]
    cdef int xy_index
    cdef np.ndarray[double, ndim=1, mode="c"] lon = np.zeros(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1, mode="c"] lat = np.zeros(n, dtype=np.double)

    dx = 0
    dy = 0

    for i in range(n):
        xy_index = healpix_nested_to_xy(nested_index[i], n_side)
        healpix_to_radec(xy_index, n_side, dx[i], dy[i], &lon[i], &lat[i])

    return lon, lat


def lonlat_to_nested(np.ndarray[double, ndim=1, mode="c"] lon,
                     np.ndarray[double, ndim=1, mode="c"] lat,
                     int n_side):

    cdef int n = lon.shape[0]
    cdef int xy_index
    cdef double dx, dy;
    cdef np.ndarray[int, ndim=1, mode="c"] nested_index = np.zeros(n, dtype=np.int32)

    for i in range(n):
        xy_index = radectohealpixf(lon[i], lat[i], n_side, &dx, &dy)
        nested_index[i] = healpix_xy_to_nested(xy_index, n_side)

    return nested_index


def lonlat_to_nested_with_offset(np.ndarray[double, ndim=1, mode="c"] lon,
                                 np.ndarray[double, ndim=1, mode="c"] lat,
                                 int n_side):

    cdef int n = lon.shape[0]
    cdef int xy_index
    cdef np.ndarray[int, ndim=1, mode="c"] nested_index = np.zeros(n, dtype=np.int32)
    cdef np.ndarray[double, ndim=1, mode="c"] dx = np.zeros(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1, mode="c"] dy = np.zeros(n, dtype=np.double)

    for i in range(n):
        xy_index = radectohealpixf(lon[i], lat[i], n_side, &dx[i], &dy[i])
        nested_index[i] = healpix_xy_to_nested(xy_index, n_side)

    return nested_index, dx, dy
