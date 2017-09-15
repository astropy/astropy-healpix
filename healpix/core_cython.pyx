"""
This module contains vectorized Cython functions for common HEALPix operations.
Since they are written in Cython rather than Python, their input types are
strict and the functions will fail if the incorrect types are passed in.
"""

import numpy as np
cimport numpy as np
import cython

ctypedef np.intp_t intp_t
ctypedef np.double_t double_t

npy_double = np.double
npy_intp = np.intp
npy_int64 = np.int64

cdef int ORDER_NESTED = 0
cdef int ORDER_RING = 1


def healpix_to_lonlat(np.ndarray[int64_t, ndim=1, mode="c"] healpix_index,
                      int nside, int order):
    """
    Convert HEALPix indices to longitudes/latitudes.

    This returns the longitudes/latitudes of the center of the HEALPix pixels.
    If you also want to provide relative offsets inside the pixels, see
    :func:`healpix_with_offset_to_lonlat`.

    Parameters
    ----------
    healpix_index : `~numpy.ndarray`
        1-D array of HEALPix indices
    nside : int
        Number of pixels along the side of each of the 12 top-level HEALPix tiles
    order : int
        Order of HEALPix pixels. Set this to 0 for nested order or 1 for ring
        order.

    Returns
    -------
    lon, lat : `~numpy.ndarray`
        1-D arrays of longitude and latitude in radians
    """

    cdef intp_t n = healpix_index.shape[0]
    cdef intp_t i
    cdef int64_t xy_index
    cdef double dx, dy;
    cdef np.ndarray[double_t, ndim=1, mode="c"] lon = np.zeros(n, dtype=npy_double)
    cdef np.ndarray[double_t, ndim=1, mode="c"] lat = np.zeros(n, dtype=npy_double)

    dx = 0.5
    dy = 0.5

    if order == ORDER_NESTED:
        for i in range(n):
            xy_index = healpixl_nested_to_xy(healpix_index[i], nside)
            if xy_index < 0:
                raise Exception("xy_index: " + str(xy_index))
            healpixl_to_radec(xy_index, nside, dx, dy, &lon[i], &lat[i])
    elif order == ORDER_RING:
        for i in range(n):
            xy_index = healpixl_ring_to_xy(healpix_index[i], nside)
            healpixl_to_radec(xy_index, nside, dx, dy, &lon[i], &lat[i])
    else:
        raise ValueError('order should be 0 or 1')

    return lon, lat


def healpix_with_offset_to_lonlat(np.ndarray[int64_t, ndim=1, mode="c"] healpix_index,
                                  np.ndarray[double_t, ndim=1, mode="c"] dx,
                                  np.ndarray[double_t, ndim=1, mode="c"] dy,
                                  int nside, int order):
    """
    Convert HEALPix indices to longitudes/latitudes

    This function takes relative offsets in x and y inside the HEALPix pixels.
    If you are only interested in the centers of the pixels, see
    `healpixl_to_lonlat`.

    Parameters
    ----------
    healpix_index : `~numpy.ndarray`
        1-D array of HEALPix indices
    dx, dy : `~numpy.ndarray`
        1-D arrays of offsets inside the HEALPix pixel, which should be in the
        range [0:1] (0.5 is the center of the HEALPix pixels)
    nside : int
        Number of pixels along the side of each of the 12 top-level HEALPix tiles
    order : int
        Order of HEALPix pixels. Set this to 0 for nested order or 1 for ring
        order.

    Returns
    -------
    lon, lat : `~numpy.ndarray`
        1-D arrays of longitude and latitude in radians
    """

    cdef intp_t n = healpix_index.shape[0]
    cdef intp_t i
    cdef int64_t xy_index
    cdef np.ndarray[double_t, ndim=1, mode="c"] lon = np.zeros(n, dtype=npy_double)
    cdef np.ndarray[double_t, ndim=1, mode="c"] lat = np.zeros(n, dtype=npy_double)

    if order == ORDER_NESTED:
        for i in range(n):
            xy_index = healpixl_nested_to_xy(healpix_index[i], nside)
            healpixl_to_radec(xy_index, nside, dx[i], dy[i], &lon[i], &lat[i])
    elif order == ORDER_RING:
        for i in range(n):
            xy_index = healpixl_ring_to_xy(healpix_index[i], nside)
            healpixl_to_radec(xy_index, nside, dx[i], dy[i], &lon[i], &lat[i])
    else:
        raise ValueError('order should be 0 or 1')

    return lon, lat


def lonlat_to_healpix(np.ndarray[double_t, ndim=1, mode="c"] lon,
                      np.ndarray[double_t, ndim=1, mode="c"] lat,
                      int nside, int order):
    """
    Convert longitudes/latitudes to HEALPix indices

    This returns only the HEALPix indices. If you also want to get relative
    offsets inside the pixels, see :func:`lonlat_to_healpix_with_offset`.

    Parameters
    ----------
    lon, lat : `~numpy.ndarray`
        1-D arrays of longitude and latitude in radians
    nside : int
        Number of pixels along the side of each of the 12 top-level HEALPix tiles
    order : int
        Order of HEALPix pixels. Set this to 0 for nested order or 1 for ring
        order.

    Returns
    -------
    healpix_index : `~numpy.ndarray`
        1-D array of HEALPix indices
    """

    cdef intp_t n = lon.shape[0]
    cdef intp_t i
    cdef int64_t xy_index
    cdef double dx, dy;
    cdef np.ndarray[int64_t, ndim=1, mode="c"] healpix_index = np.zeros(n, dtype=npy_int64)

    if order == ORDER_NESTED:
        for i in range(n):
            xy_index = radectohealpixlf(lon[i], lat[i], nside, &dx, &dy)
            healpix_index[i] = healpixl_xy_to_nested(xy_index, nside)
    elif order == ORDER_RING:
        for i in range(n):
            xy_index = radectohealpixlf(lon[i], lat[i], nside, &dx, &dy)
            healpix_index[i] = healpixl_xy_to_ring(xy_index, nside)
    else:
        raise ValueError('order should be 0 or 1')

    return healpix_index


def lonlat_to_healpix_with_offset(np.ndarray[double_t, ndim=1, mode="c"] lon,
                                  np.ndarray[double_t, ndim=1, mode="c"] lat,
                                  int nside, int order):
    """
    Convert longitudes/latitudes to healpix indices

    This returns the HEALPix indices and relative offsets inside the pixels. If
    you want only the HEALPix indices, see :func:`lonlat_to_healpix`.

    Parameters
    ----------
    lon, lat : `~numpy.ndarray`
        1-D arrays of longitude and latitude in radians
    nside : int
        Number of pixels along the side of each of the 12 top-level HEALPix tiles
    order : int
        Order of HEALPix pixels. Set this to 0 for nested order or 1 for ring
        order.

    Returns
    -------
    healpix_index : `~numpy.ndarray`
        1-D array of HEALPix indices
    dx, dy : `~numpy.ndarray`
        1-D arrays of offsets inside the HEALPix pixel in the range [0:1] (0.5
        is the center of the HEALPix pixels)
    """

    cdef intp_t n = lon.shape[0]
    cdef intp_t i
    cdef int64_t xy_index
    cdef np.ndarray[int64_t, ndim=1, mode="c"] healpix_index = np.zeros(n, dtype=npy_int64)
    cdef np.ndarray[double_t, ndim=1, mode="c"] dx = np.zeros(n, dtype=npy_double)
    cdef np.ndarray[double_t, ndim=1, mode="c"] dy = np.zeros(n, dtype=npy_double)

    if order == ORDER_NESTED:
        for i in range(n):
            xy_index = radectohealpixlf(lon[i], lat[i], nside, &dx[i], &dy[i])
            healpix_index[i] = healpixl_xy_to_nested(xy_index, nside)
    elif order == ORDER_RING:
        for i in range(n):
            xy_index = radectohealpixlf(lon[i], lat[i], nside, &dx[i], &dy[i])
            healpix_index[i] = healpixl_xy_to_ring(xy_index, nside)
    else:
        raise ValueError('order should be 0 or 1')

    return healpix_index, dx, dy


def nested_to_ring(np.ndarray[int64_t, ndim=1, mode="c"] nested_index, int nside):
    """
    Convert a HEALPix 'nested' index to a HEALPix 'ring' index

    Parameters
    ----------
    nested_index : `~numpy.ndarray`
        Healpix index using the 'nested' ordering
    nside : int
        Number of pixels along the side of each of the 12 top-level HEALPix tiles

    Returns
    -------
    ring_index : `~numpy.ndarray`
        Healpix index using the 'ring' ordering
    """

    cdef intp_t n = nested_index.shape[0]
    cdef intp_t i
    cdef np.ndarray[int64_t, ndim=1, mode="c"] ring_index = np.zeros(n, dtype=npy_int64)

    for i in range(n):
        ring_index[i] = healpixl_xy_to_ring(healpixl_nested_to_xy(nested_index[i], nside), nside)

    return ring_index


def ring_to_nested(np.ndarray[int64_t, ndim=1, mode="c"] ring_index, int nside):
    """
    Convert a HEALPix 'ring' index to a HEALPix 'nested' index

    Parameters
    ----------
    ring_index : `~numpy.ndarray`
        Healpix index using the 'ring' ordering
    nside : int
        Number of pixels along the side of each of the 12 top-level HEALPix tiles

    Returns
    -------
    nested_index : `~numpy.ndarray`
        Healpix index using the 'nested' ordering
    """

    cdef intp_t n = ring_index.shape[0]
    cdef intp_t i
    cdef np.ndarray[int64_t, ndim=1, mode="c"] nested_index = np.zeros(n, dtype=npy_int64)

    for i in range(n):
        nested_index[i] = healpixl_xy_to_nested(healpixl_ring_to_xy(ring_index[i], nside), nside)

    return nested_index


def interpolate_bilinear_lonlat(np.ndarray[double_t, ndim=1, mode="c"] lon,
                                np.ndarray[double_t, ndim=1, mode="c"] lat,
                                np.ndarray[double_t, ndim=1, mode="c"] values,
                                int order):
    """
    Interpolate values at specific longitudes/latitudes using bilinear interpolation

    Parameters
    ----------
    lon, lat : `~numpy.ndarray`
        1-D arrays of longitude and latitude in radians
    values : `~numpy.ndarray`
        1-D array with the values in each HEALPix pixel. This should have a
        length of the form 12 * nside ** 2 (and nside is determined
        automatically from this).
    order : int
        Order of HEALPix pixels. Set this to 0 for nested order or 1 for ring
        order.

    Returns
    -------
    result : `~numpy.ndarray`
        1-D array of interpolated values
    """

    cdef int nside
    cdef intp_t n = lon.shape[0]
    cdef intp_t i
    cdef int64_t xy_index, i11, i12, i21, i22
    cdef double dx, dy, xfrac, yfrac
    cdef double_t v11, v12, v21, v22
    cdef np.ndarray[double_t, ndim=1, mode="c"] result = np.zeros(n, dtype=npy_double)
    cdef int64_t neighbours[8]
    cdef double_t invalid = np.nan
    cdef intp_t npix = values.shape[0]
    cdef double square_root

    if npix % 12 != 0:
        raise ValueError('Number of pixels should be divisible by 12')

    square_root = (npix / 12.) ** 0.5

    if square_root ** 2 != npix / 12:
        raise ValueError('Number of pixels is not of the form 12 * nside ** 2')

    nside = int(square_root)

    for i in range(n):

        xy_index = radectohealpixlf(lon[i], lat[i], nside, &dx, &dy)

        # We now need to identify the four pixels that surround the position
        # we've identified. The neighbours are ordered as follows:
        #
        #       3   2   1
        #       4   X   0
        #       5   6   7

        healpix_get_neighboursl(xy_index, neighbours, nside)

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

        if order == ORDER_NESTED:
            i11 = healpixl_xy_to_nested(i11, nside)
            i12 = healpixl_xy_to_nested(i12, nside)
            i21 = healpixl_xy_to_nested(i21, nside)
            i22 = healpixl_xy_to_nested(i22, nside)
        elif order == ORDER_RING:
            i11 = healpixl_xy_to_ring(i11, nside)
            i12 = healpixl_xy_to_ring(i12, nside)
            i21 = healpixl_xy_to_ring(i21, nside)
            i22 = healpixl_xy_to_ring(i22, nside)

        v11 = values[i11]
        v12 = values[i12]
        v21 = values[i21]
        v22 = values[i22]

        result[i] = (v11 * (1 - xfrac) * (1 - yfrac) +
                     v12 * (1 - xfrac) * yfrac +
                     v21 * xfrac * (1 - yfrac) +
                     v22 * xfrac * yfrac)

    return result


def healpix_neighbors(np.ndarray[int64_t, ndim=1, mode="c"] healpix_index,
                       int nside, int order):
    """
    Find all the HEALPix pixels that are the neighbours of a HEALPix pixel

    Parameters
    ----------
    healpix_pixel : `~numpy.ndarray`
        1-D array of HEALPix pixels
    nside : int
        Number of pixels along the side of each of the 12 top-level HEALPix tiles
    order : int
        Order of HEALPix pixels. Set this to 0 for nested order or 1 for ring
        order.

    Returns
    -------
    neighbours : `~numpy.ndarray`
        2-D array with shape (8, N) giving the neighbours starting SW and
        rotating clockwise.
    """

    cdef intp_t n = healpix_index.shape[0]
    cdef intp_t i
    cdef int64_t xy_index
    cdef int j, k
    cdef np.ndarray[int64_t, ndim=2, mode="c"] neighbours = np.zeros((8, n), dtype=npy_int64)
    cdef int64_t neighbours_indiv[8]
    cdef int n_neighbours

    # The neighbours above are ordered as follows:
    #
    #       3   2   1
    #       4   X   0
    #       5   6   7
    #
    # but we want:
    #
    #       2   3   4
    #       1   X   5
    #       0   7   6
    #
    # so we reorder these on-the-fly

    if order == ORDER_NESTED:

        for i in range(n):

            xy_index = healpixl_nested_to_xy(healpix_index[i], nside)
            n_neighbours = healpix_get_neighboursl(xy_index, neighbours_indiv, nside)

            for j in range(8):
                k = 5 - j
                if k < 0:
                    k = k + 8
                if j >= n_neighbours or neighbours_indiv[k] < 0:
                    neighbours[j, i] = -1
                else:
                    neighbours[j, i] = healpixl_xy_to_nested(neighbours_indiv[k], nside)

    elif order == ORDER_RING:

        for i in range(n):

            xy_index = healpixl_ring_to_xy(healpix_index[i], nside)

            nn = healpix_get_neighboursl(xy_index, neighbours_indiv, nside)

            for j in range(8):
                k = 5 - j
                if k < 0:
                    k = k + 8
                if neighbours_indiv[k] < 0:
                    neighbours[j, i] = -1
                else:
                    neighbours[j, i] = healpixl_xy_to_ring(neighbours_indiv[k], nside)

    else:
        raise ValueError('order should be 0 or 1')


    return neighbours


def healpix_search(double ra, double dec, double radius, int nside):

    cdef intp_t i
    cdef int xy_index
    cdef int j, k
    cdef int *indices
    cdef int n_indices
    cdef double t
    cdef x=1, y=2

    print(1, ra, dec, radius, nside)
    # t = uniform_sample(x, y)
    # print(t)
    n_indices =  healpix_rangesearch_radec_simple(ra, dec, radius, nside, indices)
    print(2)

    cdef np.ndarray[int64_t, ndim=2, mode="c"] result = np.zeros(n_indices, dtype=npy_int64)

    print(3)
    # for i in range(n_indices):
    #     result[i] = indices[i]

    return result