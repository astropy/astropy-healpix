"""
This module contains vectorized Cython functions for common healpix operations.
Since they are written in Cython rather than Python, their input types are
strict and the functions will fail if the incorrect types are passed in.
"""

import numpy as np
cimport numpy as np
import cython

ctypedef np.double_t DOUBLE_T

cdef int ORDER_NESTED = 0
cdef int ORDER_RING = 1
cdef int ORDER_XY = 2


def nested_to_lonlat(np.ndarray[int, ndim=1, mode="c"] nested_index, int n_side):
    """
    Convert healpix indices to longitudes/latitudes using the 'nested' convention.

    This returns the longitudes/latitudes of the center of the healpix pixels.
    If you also want to provide relative offsets inside the pixels, see
    :func:`nested_with_offset_to_lonlat`.

    Parameters
    ----------
    nested_index : `~numpy.ndarray`
        1-D array of healpix indices using the 'nested' convention
    n_side : int
        Number of pixels along the side of each of the 12 top-level healpix tiles

    Returns
    -------
    lon, lat : `~numpy.ndarray`
        1-D arrays of longitude and latitude in radians
    """

    cdef int n = nested_index.shape[0]
    cdef double dx, dy;
    cdef int i, xy_index
    cdef np.ndarray[double, ndim=1, mode="c"] lon = np.zeros(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1, mode="c"] lat = np.zeros(n, dtype=np.double)

    dx = 0.5
    dy = 0.5

    for i in range(n):
        xy_index = healpix_nested_to_xy(nested_index[i], n_side)
        healpix_to_radec(xy_index, n_side, dx, dy, &lon[i], &lat[i])
    return lon, lat


def nested_with_offset_to_lonlat(np.ndarray[int, ndim=1, mode="c"] nested_index,
                                 np.ndarray[double, ndim=1, mode="c"] dx,
                                 np.ndarray[double, ndim=1, mode="c"] dy,
                                 int n_side):
    """
    Convert healpix indices to longitudes/latitudes using the 'nested' convention.

    This function takes relative offsets in x and y inside the healpix pixels.
    If you are only interested in the centers of the pixels, see
    `nested_to_lonlat`.

    Parameters
    ----------
    nested_index : `~numpy.ndarray`
        1-D array of healpix indices using the 'nested' convention
    dx, dy : `~numpy.ndarray`
        1-D arrays of offsets inside the healpix pixel, which should be in the
        range [0:1] (0.5 is the center of the healpix pixels)
    n_side : int
        Number of pixels along the side of each of the 12 top-level healpix tiles

    Returns
    -------
    lon, lat : `~numpy.ndarray`
        1-D arrays of longitude and latitude in radians
    """

    cdef int n = nested_index.shape[0]
    cdef int i, xy_index
    cdef np.ndarray[double, ndim=1, mode="c"] lon = np.zeros(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1, mode="c"] lat = np.zeros(n, dtype=np.double)

    for i in range(n):
        xy_index = healpix_nested_to_xy(nested_index[i], n_side)
        healpix_to_radec(xy_index, n_side, dx[i], dy[i], &lon[i], &lat[i])

    return lon, lat


def lonlat_to_nested(np.ndarray[double, ndim=1, mode="c"] lon,
                     np.ndarray[double, ndim=1, mode="c"] lat,
                     int n_side):
    """
    Convert longitudes/latitudes to healpix indices using the 'nested' convention.

    This returns only the healpix indices. If you also want to get relative
    offsets inside the pixels, see :func:`lonlat_to_nested_with_offset`.

    Parameters
    ----------
    lon, lat : `~numpy.ndarray`
        1-D arrays of longitude and latitude in radians
    n_side : int
        Number of pixels along the side of each of the 12 top-level healpix tiles

    Returns
    -------
    nested_index : `~numpy.ndarray`
        1-D array of healpix indices using the 'nested' convention
    """

    cdef int n = lon.shape[0]
    cdef int i, xy_index
    cdef double dx, dy;
    cdef np.ndarray[int, ndim=1, mode="c"] nested_index = np.zeros(n, dtype=np.int32)

    for i in range(n):
        xy_index = radectohealpixf(lon[i], lat[i], n_side, &dx, &dy)
        nested_index[i] = healpix_xy_to_nested(xy_index, n_side)

    return nested_index


def lonlat_to_nested_with_offset(np.ndarray[double, ndim=1, mode="c"] lon,
                                 np.ndarray[double, ndim=1, mode="c"] lat,
                                 int n_side):
    """
    Convert longitudes/latitudes to healpix indices using the 'nested' convention.

    This returns the healpix indices and relative offsets inside the pixels. If
    you want only the healpix indices, see :func:`lonlat_to_nested`.

    Parameters
    ----------
    lon, lat : `~numpy.ndarray`
        1-D arrays of longitude and latitude in radians
    n_side : int
        Number of pixels along the side of each of the 12 top-level healpix tiles

    Returns
    -------
    nested_index : `~numpy.ndarray`
        1-D array of healpix indices using the 'nested' convention
    dx, dy : `~numpy.ndarray`
        1-D arrays of offsets inside the healpix pixel in the range [0:1] (0.5
        is the center of the healpix pixels)
    """


    cdef int n = lon.shape[0]
    cdef int i, xy_index
    cdef np.ndarray[int, ndim=1, mode="c"] nested_index = np.zeros(n, dtype=np.int32)
    cdef np.ndarray[double, ndim=1, mode="c"] dx = np.zeros(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1, mode="c"] dy = np.zeros(n, dtype=np.double)

    for i in range(n):
        xy_index = radectohealpixf(lon[i], lat[i], n_side, &dx[i], &dy[i])
        nested_index[i] = healpix_xy_to_nested(xy_index, n_side)

    return nested_index, dx, dy


def ring_to_lonlat(np.ndarray[int, ndim=1, mode="c"] ring_index, int n_side):
    """
    Convert healpix indices to longitudes/latitudes using the 'ring' convention.

    This returns the longitudes/latitudes of the center of the healpix pixels.
    If you also want to provide relative offsets inside the pixels, see
    :func:`ring_with_offset_to_lonlat`.

    Parameters
    ----------
    ring_index : `~numpy.ndarray`
        1-D array of healpix indices using the 'ring' convention
    n_side : int
        Number of pixels along the side of each of the 12 top-level healpix tiles

    Returns
    -------
    lon, lat : `~numpy.ndarray`
        1-D arrays of longitude and latitude in radians
    """

    cdef int n = ring_index.shape[0]
    cdef double dx, dy;
    cdef int i, xy_index
    cdef np.ndarray[double, ndim=1, mode="c"] lon = np.zeros(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1, mode="c"] lat = np.zeros(n, dtype=np.double)

    dx = 0.5
    dy = 0.5

    for i in range(n):
        xy_index = healpix_ring_to_xy(ring_index[i], n_side)
        healpix_to_radec(xy_index, n_side, dx, dy, &lon[i], &lat[i])
    return lon, lat


def ring_with_offset_to_lonlat(np.ndarray[int, ndim=1, mode="c"] ring_index,
                                 np.ndarray[double, ndim=1, mode="c"] dx,
                                 np.ndarray[double, ndim=1, mode="c"] dy,
                                 int n_side):
    """
    Convert healpix indices to longitudes/latitudes using the 'ring' convention.

    This function takes relative offsets in x and y inside the healpix pixels.
    If you are only interested in the centers of the pixels, see
    `ring_to_lonlat`.

    Parameters
    ----------
    ring_index : `~numpy.ndarray`
        1-D array of healpix indices using the 'ring' convention
    dx, dy : `~numpy.ndarray`
        1-D arrays of offsets inside the healpix pixel, which should be in the
        range [0:1] (0.5 is the center of the healpix pixels)
    n_side : int
        Number of pixels along the side of each of the 12 top-level healpix tiles

    Returns
    -------
    lon, lat : `~numpy.ndarray`
        1-D arrays of longitude and latitude in radians
    """

    cdef int n = ring_index.shape[0]
    cdef int i, xy_index
    cdef np.ndarray[double, ndim=1, mode="c"] lon = np.zeros(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1, mode="c"] lat = np.zeros(n, dtype=np.double)

    for i in range(n):
        xy_index = healpix_ring_to_xy(ring_index[i], n_side)
        healpix_to_radec(xy_index, n_side, dx[i], dy[i], &lon[i], &lat[i])

    return lon, lat


def lonlat_to_ring(np.ndarray[double, ndim=1, mode="c"] lon,
                     np.ndarray[double, ndim=1, mode="c"] lat,
                     int n_side):
    """
    Convert longitudes/latitudes to healpix indices using the 'ring' convention.

    This returns only the healpix indices. If you also want to get relative
    offsets inside the pixels, see :func:`lonlat_to_ring_with_offset`.

    Parameters
    ----------
    lon, lat : `~numpy.ndarray`
        1-D arrays of longitude and latitude in radians
    n_side : int
        Number of pixels along the side of each of the 12 top-level healpix tiles

    Returns
    -------
    ring_index : `~numpy.ndarray`
        1-D array of healpix indices using the 'ring' convention
    """

    cdef int n = lon.shape[0]
    cdef int i, xy_index
    cdef double dx, dy;
    cdef np.ndarray[int, ndim=1, mode="c"] ring_index = np.zeros(n, dtype=np.int32)

    for i in range(n):
        xy_index = radectohealpixf(lon[i], lat[i], n_side, &dx, &dy)
        ring_index[i] = healpix_xy_to_ring(xy_index, n_side)

    return ring_index


def lonlat_to_ring_with_offset(np.ndarray[double, ndim=1, mode="c"] lon,
                                 np.ndarray[double, ndim=1, mode="c"] lat,
                                 int n_side):
    """
    Convert longitudes/latitudes to healpix indices using the 'ring' convention.

    This returns the healpix indices and relative offsets inside the pixels. If
    you want only the healpix indices, see :func:`lonlat_to_ring`.

    Parameters
    ----------
    lon, lat : `~numpy.ndarray`
        1-D arrays of longitude and latitude in radians
    n_side : int
        Number of pixels along the side of each of the 12 top-level healpix tiles

    Returns
    -------
    ring_index : `~numpy.ndarray`
        1-D array of healpix indices using the 'ring' convention
    dx, dy : `~numpy.ndarray`
        1-D arrays of offsets inside the healpix pixel in the range [0:1] (0.5
        is the center of the healpix pixels)
    """


    cdef int n = lon.shape[0]
    cdef int i, xy_index
    cdef np.ndarray[int, ndim=1, mode="c"] ring_index = np.zeros(n, dtype=np.int32)
    cdef np.ndarray[double, ndim=1, mode="c"] dx = np.zeros(n, dtype=np.double)
    cdef np.ndarray[double, ndim=1, mode="c"] dy = np.zeros(n, dtype=np.double)

    for i in range(n):
        xy_index = radectohealpixf(lon[i], lat[i], n_side, &dx[i], &dy[i])
        ring_index[i] = healpix_xy_to_ring(xy_index, n_side)

    return ring_index, dx, dy



def bilinear_interpolation_nested(np.ndarray[double, ndim=1, mode="c"] lon,
                                  np.ndarray[double, ndim=1, mode="c"] lat,
                                  np.ndarray[double, ndim=1, mode="c"] values,
                                  int n_side, int ordering):

    cdef int n = lon.shape[0]
    cdef int i, xy_index, ix, iy, i11, i12, i21, i22
    cdef double dx, dy, v11, v12, v21, v22, xfrac, yfrac
    cdef np.ndarray[double, ndim=1, mode="c"] result = np.zeros(n, dtype=np.double)
    cdef int neighbours[8]
    cdef double invalid = np.nan

    for i in range(n):

        xy_index = radectohealpixf(lon[i], lat[i], n_side, &dx, &dy)

        # We now need to identify the four pixels that surround the position
        # we've identified. The neighbours are ordered as follows:
        #
        #       3   2   1
        #       4   X   0
        #       5   6   7

        healpix_get_neighbours(xy_index, neighbours, n_side)

        if dx < 0.5:

            if dy < 0.5:
                i11 = neighbours[5]
                i12 = neighbours[4]
                i21 = neighbours[6]
                i22 = xy_index
                xfrac = 0.5 + dx
                yfrac = 0.5 + dy
            else:
                i11 = neighbours[4]
                i12 = neighbours[3]
                i21 = xy_index
                i22 = neighbours[2]
                xfrac = 0.5 + dx
                yfrac = dy - 0.5

        else:

            if dy < 0.5:
                i11 = neighbours[6]
                i12 = xy_index
                i21 = neighbours[7]
                i22 = neighbours[0]
                xfrac = dx - 0.5
                yfrac = 0.5 + dy
            else:
                i11 = xy_index
                i12 = neighbours[2]
                i21 = neighbours[0]
                i22 = neighbours[1]
                xfrac = dx - 0.5
                yfrac = dy - 0.5

        if i11 < 0 or i12 < 0 or i21 < 0 or i22 < 0:
            result[i] = invalid
            continue

        if ordering == ORDER_NESTED:
            i11 = healpix_xy_to_nested(i11, n_side)
            i12 = healpix_xy_to_nested(i12, n_side)
            i21 = healpix_xy_to_nested(i21, n_side)
            i22 = healpix_xy_to_nested(i22, n_side)
        elif ordering == ORDER_RING:
            i11 = healpix_xy_to_ring(i11, n_side)
            i12 = healpix_xy_to_ring(i12, n_side)
            i21 = healpix_xy_to_ring(i21, n_side)
            i22 = healpix_xy_to_ring(i22, n_side)

        v11 = values[i11]
        v12 = values[i12]
        v21 = values[i21]
        v22 = values[i22]

        result[i] = (v11 * (1 - xfrac) * (1 - yfrac) +
                     v12 * (1 - xfrac) * yfrac +
                     v21 * xfrac * (1 - yfrac) +
                     v22 * xfrac * yfrac)

    return result
