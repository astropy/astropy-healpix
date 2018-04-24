# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains vectorized Cython functions for common HEALPix operations.
Since they are written in Cython rather than Python, their input types are
strict and the functions will fail if the incorrect types are passed in.
"""

import numpy as np
cimport numpy as np
import cython
from cython.parallel import parallel, prange
from libc.stdlib cimport abort, malloc, free
from libc.math cimport sin, cos, sqrt

ctypedef np.intp_t intp_t
ctypedef np.double_t double_t

npy_double = np.double
npy_intp = np.intp
npy_int64 = np.int64

np.import_array()
np.import_ufunc()


def _validate_order(str order):
    # We also support upper-case, to support directly the values
    # ORDERING = {'RING', 'NESTED'} in FITS headers
    # This is currently undocumented in the docstrings.
    if order == 'nested' or order == 'NESTED':
        return 'nested'
    elif order == 'ring' or order == 'RING':
        return 'ring'
    else:
        raise ValueError("order must be 'nested' or 'ring'")


@cython.boundscheck(False)
cdef void healpix_to_lonlat_ring_loop(char** args, intp_t* dimensions,
                                      intp_t* steps, void* data) nogil:
    cdef intp_t n = dimensions[0]
    cdef intp_t i
    cdef int64_t xy_index, healpix_index, nside
    cdef double dx, dy
    cdef double* lon
    cdef double* lat

    for i in prange(n, schedule='static'):
        healpix_index = (<int64_t*> &args[0][i * steps[0]])[0]
        nside = (<int64_t*> &args[1][i * steps[1]])[0]
        dx = (<double*> &args[2][i * steps[2]])[0]
        dy = (<double*> &args[3][i * steps[3]])[0]
        lon = <double*> &args[4][i * steps[5]]
        lat = <double*> &args[5][i * steps[5]]
        xy_index = healpixl_ring_to_xy(healpix_index, nside)
        healpixl_to_radec(xy_index, nside, dx, dy, lon, lat)


healpix_to_lonlat_ring = np.PyUFunc_FromFuncAndData(
    [healpix_to_lonlat_ring_loop], [NULL],
    [np.NPY_INT64, np.NPY_INT64,
     np.NPY_DOUBLE, np.NPY_DOUBLE, np.NPY_DOUBLE, np.NPY_DOUBLE],
    1, 4, 2, 0, "healpix_to_lonlat_ring", "", 0)


@cython.boundscheck(False)
cdef void healpix_to_lonlat_nested_loop(
        char** args, intp_t* dimensions, intp_t* steps, void* data) nogil:
    cdef intp_t n = dimensions[0]
    cdef intp_t i
    cdef int64_t xy_index, healpix_index, nside
    cdef double dx, dy
    cdef double* lon
    cdef double* lat

    for i in prange(n, schedule='static'):
        healpix_index = (<int64_t*> &args[0][i * steps[0]])[0]
        nside = (<int64_t*> &args[1][i * steps[1]])[0]
        dx = (<double*> &args[2][i * steps[2]])[0]
        dy = (<double*> &args[3][i * steps[3]])[0]
        lon = <double*> &args[4][i * steps[5]]
        lat = <double*> &args[5][i * steps[5]]
        xy_index = healpixl_nested_to_xy(healpix_index, nside)
        healpixl_to_radec(xy_index, nside, dx, dy, lon, lat)


healpix_to_lonlat_nested = np.PyUFunc_FromFuncAndData(
    [healpix_to_lonlat_nested_loop], [NULL],
    [np.NPY_INT64, np.NPY_INT64, np.NPY_DOUBLE, np.NPY_DOUBLE,
     np.NPY_DOUBLE, np.NPY_DOUBLE],
    1, 4, 2, 0, "healpix_to_lonlat_nested", "", 0)


@cython.boundscheck(False)
def lonlat_to_healpix(np.ndarray[double_t, ndim=1, mode="c"] lon,
                      np.ndarray[double_t, ndim=1, mode="c"] lat,
                      int nside, str order):
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
    order : { 'nested' | 'ring' }
        Order of HEALPix pixels


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

    order = _validate_order(order)

    if order == 'nested':
        for i in prange(n, nogil=True, schedule='static'):
            xy_index = radec_to_healpixlf(lon[i], lat[i], nside, &dx, &dy)
            healpix_index[i] = healpixl_xy_to_nested(xy_index, nside)
    elif order == 'ring':
        for i in prange(n, nogil=True, schedule='static'):
            xy_index = radec_to_healpixlf(lon[i], lat[i], nside, &dx, &dy)
            healpix_index[i] = healpixl_xy_to_ring(xy_index, nside)

    return healpix_index


@cython.boundscheck(False)
def lonlat_to_healpix_with_offset(np.ndarray[double_t, ndim=1, mode="c"] lon,
                                  np.ndarray[double_t, ndim=1, mode="c"] lat,
                                  int nside, str order):
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
    order : { 'nested' | 'ring' }
        Order of HEALPix pixels

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

    order = _validate_order(order)

    if order == 'nested':
        for i in prange(n, nogil=True, schedule='static'):
            xy_index = radec_to_healpixlf(lon[i], lat[i], nside, &dx[i], &dy[i])
            healpix_index[i] = healpixl_xy_to_nested(xy_index, nside)
    elif order == 'ring':
        for i in prange(n, nogil=True, schedule='static'):
            xy_index = radec_to_healpixlf(lon[i], lat[i], nside, &dx[i], &dy[i])
            healpix_index[i] = healpixl_xy_to_ring(xy_index, nside)

    return healpix_index, dx, dy


@cython.boundscheck(False)
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

    for i in prange(n, nogil=True, schedule='static'):
        ring_index[i] = healpixl_xy_to_ring(healpixl_nested_to_xy(nested_index[i], nside), nside)

    return ring_index


@cython.boundscheck(False)
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

    for i in prange(n, nogil=True, schedule='static'):
        nested_index[i] = healpixl_xy_to_nested(healpixl_ring_to_xy(ring_index[i], nside), nside)

    return nested_index



@cython.boundscheck(False)
def bilinear_interpolation_weights(np.ndarray[double_t, ndim=1, mode="c"] lon,
                                   np.ndarray[double_t, ndim=1, mode="c"] lat,
                                   int nside, str order):
    """
    Get the four neighbours for each (lon, lat) position and the weight
    associated with each one for bilinear interpolation.

    Parameters
    ----------
    lon, lat : `~numpy.ndarray`
        1-D arrays of longitude and latitude in radians
    nside : int
        Number of pixels along the side of each of the 12 top-level HEALPix tiles
    order : { 'nested' | 'ring' }
        Order of HEALPix pixels

    Returns
    -------
    indices : `~numpy.ndarray`
        2-D array with shape (4, N) giving the four indices to use for the
        interpolation
    weights : `~numpy.ndarray`
        2-D array with shape (4, N) giving the four weights to use for the
        interpolation
    """

    cdef intp_t n = lon.shape[0]
    cdef intp_t i, j
    cdef np.ndarray[int64_t, ndim=2, mode="c"] indices = np.zeros((4, n), dtype=npy_int64)
    cdef np.ndarray[double_t, ndim=2, mode="c"] weights = np.zeros((4, n), dtype=npy_double)
    cdef int order_int

    # Since we want to be able to use OpenMP in this function, we need to make
    # sure that any temporary buffers are allocated inside the parallel()
    # context. Note that we also need to do this for dx_buf and dy_buf otherwise
    # if we just passed &dx and &dy to radec_to_healpixlf, different threads
    # would be accessing the same location in memory, causing issues. We use
    # manual memory management with malloc as this appears to be the recommended
    # method at http://cython.readthedocs.io/en/latest/src/userguide/parallelism.html
    cdef double *weights_indiv
    cdef int64_t *indices_indiv

    order = _validate_order(order)

    if order == 'nested':
        order_int = 0
    elif order == 'ring':
        order_int = 1

    with nogil, parallel():

        indices_indiv = <int64_t *> malloc(sizeof(int64_t) * 4)
        if indices_indiv == NULL:
            abort()

        weights_indiv = <double *> malloc(sizeof(double) * 4)
        if weights_indiv == NULL:
            abort()

        for i in range(n):
              interpolate_weights(lon[i], lat[i], indices_indiv, weights_indiv, nside)
              for j in range(4):
                  if order_int == 0:
                      indices[j, i] = healpixl_xy_to_nested(healpixl_ring_to_xy(indices_indiv[j], nside), nside)
                  else:
                      indices[j, i] = indices_indiv[j]
                  weights[j, i] = weights_indiv[j]

    return indices, weights


@cython.boundscheck(False)
def neighbours(np.ndarray[int64_t, ndim=1, mode="c"] healpix_index,
               int nside, str order):
    """
    Find all the HEALPix pixels that are the neighbours of a HEALPix pixel

    Parameters
    ----------
    healpix_index : `~numpy.ndarray`
        1-D array of HEALPix indices
    nside : int
        Number of pixels along the side of each of the 12 top-level HEALPix tiles
    order : { 'nested' | 'ring' }
        Order of HEALPix pixels

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
    cdef int64_t * neighbours_indiv
    cdef int order_int

    order = _validate_order(order)

    if order == 'nested':
        order_int = 0
    elif order == 'ring':
        order_int = 1

    with nogil, parallel():

        neighbours_indiv = <int64_t *> malloc(sizeof(int64_t) * 8)
        if neighbours_indiv == NULL:
            abort()

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

        if order_int == 0:

            for i in prange(n, schedule='static'):

                xy_index = healpixl_nested_to_xy(healpix_index[i], nside)
                healpixl_get_neighbours(xy_index, neighbours_indiv, nside)

                for j in range(8):
                    k = 4 - j
                    if k < 0:
                        k = k + 8
                    if neighbours_indiv[k] < 0:
                        neighbours[j, i] = -1
                    else:
                        neighbours[j, i] = healpixl_xy_to_nested(neighbours_indiv[k], nside)

        elif order_int == 1:

            for i in prange(n, schedule='static'):

                xy_index = healpixl_ring_to_xy(healpix_index[i], nside)

                healpixl_get_neighbours(xy_index, neighbours_indiv, nside)

                for j in range(8):
                    k = 4 - j
                    if k < 0:
                        k = k + 8
                    if neighbours_indiv[k] < 0:
                        neighbours[j, i] = -1
                    else:
                        neighbours[j, i] = healpixl_xy_to_ring(neighbours_indiv[k], nside)

        free(neighbours_indiv)

    return neighbours


@cython.boundscheck(False)
def healpix_cone_search(double lon, double lat, double radius, int nside, str order):
    """
    Find all the HEALPix pixels within a given radius of a longitude/latitude.

    Note that this returns all pixels that overlap, including partially, with
    the search cone. This function can only be used for a single lon/lat pair at
    a time, since different calls to the function may result in a different
    number of matches.

    Parameters
    ----------
    lon, lat : float
        The longitude and latitude to search around, in degrees
    radius : float
        The search radius, in degrees
    nside : int
        Number of pixels along the side of each of the 12 top-level HEALPix tiles
    order : { 'nested' | 'ring' }
        Order of HEALPix pixels

    Returns
    -------
    healpix_index : `~numpy.ndarray`
        1-D array with all the matching HEALPix pixel indices.
    """
    cdef intp_t i
    cdef int64_t *indices
    cdef int64_t n_indices
    cdef int64_t index

    n_indices =  healpix_rangesearch_radec_simple(lon, lat, radius, nside, 0, &indices)

    cdef np.ndarray[int64_t, ndim=1, mode="c"] result = np.zeros(n_indices, dtype=npy_int64)

    order = _validate_order(order)

    if order == 'nested':
        for i in prange(n_indices, nogil=True, schedule='static'):
            index = indices[i]
            result[i] = healpixl_xy_to_nested(index, nside)
    elif order == 'ring':
        for i in prange(n_indices, nogil=True, schedule='static'):
            index = indices[i]
            result[i] = healpixl_xy_to_ring(index, nside)

    return result
