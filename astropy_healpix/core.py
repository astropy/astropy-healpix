# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, print_function, division

import math

import numpy as np

from astropy import units as u
from astropy.coordinates import Longitude, Latitude

from . import core_cython
from .core_cython import _validate_order

__all__ = [
    'nside_to_pixel_area',
    'nside_to_pixel_resolution',
    'nside_to_npix',
    'npix_to_nside',
    'lonlat_to_healpix',
    'healpix_to_lonlat',
    'bilinear_interpolation_weights',
    'interpolate_bilinear_lonlat',
    'neighbours',
]


def _restore_shape(*args, **kwargs):
    shape = kwargs['shape']
    if shape:
        if len(args) > 1:
            return [arg.reshape(shape) for arg in args]
        else:
            return args[0].reshape(shape)
    else:
        if len(args) > 1:
            return [np.asscalar(arg) for arg in args]
        else:
            return np.asscalar(args[0])


def _validate_healpix_index(label, healpix_index, nside):
    npix = nside_to_npix(nside)
    if np.any((healpix_index < 0) | (healpix_index > npix - 1)):
        raise ValueError('{0} must be in the range [0:{1}]'.format(label, npix))


def _validate_offset(label, offset):
    if np.any((offset < 0) | (offset > 1)):
        raise ValueError('d{0} must be in the range [0:1]'.format(label))


def _validate_level(level):
    if level < 0:
        raise ValueError('level must be positive')


def _validate_nside(nside):
    log_2_nside = np.round(np.log2(nside))
    if not np.all(2 ** log_2_nside == nside):
        raise ValueError('nside must be a power of two')


def level_to_nside(level):
    """
    Find the pixel dimensions of the top-level HEALPix tiles given the
    resolution level (this is given by 2**level).

    Parameters
    ----------
    level : int
        The resolution level

    Returns
    -------
    nside : int
        The number of pixels on the side of one of the 12 'top-level' HEALPix tiles.
    """
    _validate_level(level)
    return 2 ** level


def nside_to_pixel_area(nside):
    """
    Find the area of HEALPix pixels given the pixel dimensions of one of
    the 12 'top-level' HEALPix tiles.

    Parameters
    ----------
    nside : int
        The number of pixels on the side of one of the 12 'top-level' HEALPix tiles.

    Returns
    -------
    pixel_area : :class:`~astropy.units.Quantity`
        The area of the HEALPix pixels
    """
    nside = np.asanyarray(nside, dtype=np.int64)
    _validate_nside(nside)
    npix = 12 * nside * nside
    pixel_area = 4 * math.pi / npix * u.sr
    return pixel_area


def nside_to_pixel_resolution(nside):
    """
    Find the resolution of HEALPix pixels given the pixel dimensions of one of
    the 12 'top-level' HEALPix tiles.

    Parameters
    ----------
    nside : int
        The number of pixels on the side of one of the 12 'top-level' HEALPix tiles.

    Returns
    -------
    resolution : :class:`~astropy.units.Quantity`
        The resolution of the HEALPix pixels
    """
    nside = np.asanyarray(nside, dtype=np.int64)
    _validate_nside(nside)
    return (nside_to_pixel_area(nside) ** 0.5).to(u.arcmin)


def nside_to_npix(nside):
    """
    Find the number of pixels corresponding to a HEALPix resolution.

    Parameters
    ----------
    nside : int
        The number of pixels on the side of one of the 12 'top-level' HEALPix tiles.

    Returns
    -------
    npix : int
        The number of pixels in the HEALPix map.
    """
    nside = np.asanyarray(nside, dtype=np.int64)
    _validate_nside(nside)
    return 12 * nside ** 2


def npix_to_nside(npix):
    """
    Find the number of pixels on the side of one of the 12 'top-level' HEALPix
    tiles given a total number of pixels.

    Parameters
    ----------
    npix : int
        The number of pixels in the HEALPix map.

    Returns
    -------
    nside : int
        The number of pixels on the side of one of the 12 'top-level' HEALPix tiles.
    """

    npix = np.asanyarray(npix, dtype=np.int64)

    if not np.all(npix % 12 == 0):
        raise ValueError('Number of pixels must be divisible by 12')

    square_root = np.sqrt(npix / 12)
    if not np.all(square_root ** 2 == npix / 12):
        raise ValueError('Number of pixels is not of the form 12 * nside ** 2')

    return np.round(square_root).astype(int)


def healpix_to_lonlat(healpix_index, nside, dx=None, dy=None, order='ring'):
    """
    Convert HEALPix indices (optionally with offsets) to longitudes/latitudes.

    If no offsets (``dx`` and ``dy``) are provided, the coordinates will default
    to those at the center of the HEALPix pixels.

    Parameters
    ----------
    healpix_index : int or `~numpy.ndarray`
        HEALPix indices (as a scalar or array)
    nside : int
        Number of pixels along the side of each of the 12 top-level HEALPix tiles
    dx, dy : float or `~numpy.ndarray`, optional
        Offsets inside the HEALPix pixel, which must be in the range [0:1],
        where 0.5 is the center of the HEALPix pixels (as scalars or arrays)
    order : { 'nested' | 'ring' }, optional
        Order of HEALPix pixels

    Returns
    -------
    lon : :class:`~astropy.coordinates.Longitude`
        The longitude values
    lat : :class:`~astropy.coordinates.Latitude`
        The latitude values
    """

    if (dx is None and dy is not None) or (dx is not None and dy is None):
        raise ValueError('Either both or neither dx and dy must be specified')

    healpix_index = np.asarray(healpix_index, dtype=np.int64)

    if dx is None and dy is not None:
        dx = 0.5
    elif dx is not None and dy is None:
        dy = 0.5

    if dx is not None:
        dx = np.asarray(dx, dtype=np.float)
        dy = np.asarray(dy, dtype=np.float)
        _validate_offset('x', dx)
        _validate_offset('y', dy)
        healpix_index, dx, dy = np.broadcast_arrays(healpix_index, dx, dy)
        dx = dx.ravel()
        dy = dy.ravel()

    shape = healpix_index.shape
    healpix_index = healpix_index.ravel()
    nside = int(nside)

    _validate_healpix_index('healpix_index', healpix_index, nside)
    _validate_nside(nside)
    order = _validate_order(order)

    if dx is None:
        lon, lat = core_cython.healpix_to_lonlat(healpix_index, nside, order)
    else:
        lon, lat = core_cython.healpix_with_offset_to_lonlat(healpix_index, dx, dy, nside, order)

    lon = Longitude(lon, unit=u.rad, copy=False)
    lat = Latitude(lat, unit=u.rad, copy=False)

    return _restore_shape(lon, lat, shape=shape)


def lonlat_to_healpix(lon, lat, nside, return_offsets=False, order='ring'):
    """
    Convert longitudes/latitudes to HEALPix indices

    Parameters
    ----------
    lon, lat : :class:`~astropy.units.Quantity`
        The longitude and latitude values as :class:`~astropy.units.Quantity`
        instances with angle units.
    nside : int
        Number of pixels along the side of each of the 12 top-level HEALPix tiles
    order : { 'nested' | 'ring' }
        Order of HEALPix pixels
    return_offsets : bool, optional
        If `True`, the returned values are the HEALPix pixel indices as well as
        ``dx`` and ``dy``, the fractional positions inside the pixels. If
        `False` (the default), only the HEALPix pixel indices is returned.

    Returns
    -------
    healpix_index : int or `~numpy.ndarray`
        The HEALPix indices
    dx, dy : `~numpy.ndarray`
        Offsets inside the HEALPix pixel in the range [0:1], where 0.5 is the
        center of the HEALPix pixels
    """

    lon = lon.to(u.rad).value
    lat = lat.to(u.rad).value

    lon, lat = np.broadcast_arrays(lon, lat)

    shape = np.shape(lon)

    lon = lon.astype(float).ravel()
    lat = lat.astype(float).ravel()

    nside = int(nside)
    _validate_nside(nside)
    order = _validate_order(order)

    if return_offsets:
        healpix_index, dx, dy = core_cython.lonlat_to_healpix_with_offset(lon, lat, nside, order)
        return _restore_shape(healpix_index, dx, dy, shape=shape)
    else:
        healpix_index = core_cython.lonlat_to_healpix(lon, lat, nside, order)
        return _restore_shape(healpix_index, shape=shape)


def nested_to_ring(nested_index, nside):
    """
    Convert a HEALPix 'nested' index to a HEALPix 'ring' index

    Parameters
    ----------
    nested_index : int or `~numpy.ndarray`
        Healpix index using the 'nested' ordering
    nside : int
        Number of pixels along the side of each of the 12 top-level HEALPix tiles

    Returns
    -------
    ring_index : int or `~numpy.ndarray`
        Healpix index using the 'ring' ordering
    """

    nested_index = np.asarray(nested_index, dtype=np.int64)
    shape = nested_index.shape
    nested_index = nested_index.ravel()
    nside = int(nside)

    _validate_healpix_index('nested_index', nested_index, nside)
    _validate_nside(nside)

    ring_index = core_cython.nested_to_ring(nested_index, nside)
    return _restore_shape(ring_index, shape=shape)


def ring_to_nested(ring_index, nside):
    """
    Convert a HEALPix 'ring' index to a HEALPix 'nested' index

    Parameters
    ----------
    ring_index : int or `~numpy.ndarray`
        Healpix index using the 'ring' ordering
    nside : int
        Number of pixels along the side of each of the 12 top-level HEALPix tiles

    Returns
    -------
    nested_index : int or `~numpy.ndarray`
        Healpix index using the 'nested' ordering
    """

    ring_index = np.asarray(ring_index, dtype=np.int64)
    shape = ring_index.shape
    ring_index = ring_index.ravel()
    nside = int(nside)

    _validate_healpix_index('ring_index', ring_index, nside)
    _validate_nside(nside)

    nested_index = core_cython.ring_to_nested(ring_index, nside)
    return _restore_shape(nested_index, shape=shape)


def bilinear_interpolation_weights(lon, lat, nside, order='ring'):
    """
    Get the four neighbours for each (lon, lat) position and the weight
    associated with each one for bilinear interpolation.

    Parameters
    ----------
    lon, lat : :class:`~astropy.units.Quantity`
        The longitude and latitude values as :class:`~astropy.units.Quantity`
        instances with angle units.
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

    lon = lon.to(u.rad).value
    lat = lat.to(u.rad).value

    lon, lat = np.broadcast_arrays(lon, lat)

    shape = np.shape(lon)

    lon = lon.astype(float).ravel()
    lat = lat.astype(float).ravel()
    nside = int(nside)

    order = _validate_order(order)
    _validate_nside(nside)

    indices, weights = core_cython.bilinear_interpolation_weights(lon, lat, nside, order)
    indices = _restore_shape(indices, shape=(4,) + shape)
    weights = _restore_shape(weights, shape=(4,) + shape)

    return indices, weights


def interpolate_bilinear_lonlat(lon, lat, values, order='ring'):
    """
    Interpolate values at specific longitudes/latitudes using bilinear interpolation

    Parameters
    ----------
    lon, lat : :class:`~astropy.units.Quantity`
        The longitude and latitude values as :class:`~astropy.units.Quantity` instances
        with angle units.
    values : `~numpy.ndarray`
        Array with the values in each HEALPix pixel. The first dimension should
        have length 12 * nside ** 2 (and nside is determined automatically from
        this).
    order : { 'nested' | 'ring' }
        Order of HEALPix pixels

    Returns
    -------
    result : float `~numpy.ndarray`
        The interpolated values
    """
    nside = npix_to_nside(values.shape[0])
    indices, weights = bilinear_interpolation_weights(lon, lat, nside, order=order)
    values = values[indices]
    # At this point values has shape (N, M) where both N and M might be several
    # dimensions, and weights has shape (N,), so we need to transpose in order
    # to benefit from broadcasting, then transpose back so that the dimension
    # with length 4 is at the start again, ready for summing.
    result = (values.T * weights.T).T
    return result.sum(axis=0)


def neighbours(healpix_index, nside, order='ring'):
    """
    Find all the HEALPix pixels that are the neighbours of a HEALPix pixel

    Parameters
    ----------
    healpix_index : `~numpy.ndarray`
        Array of HEALPix pixels
    nside : int
        Number of pixels along the side of each of the 12 top-level HEALPix tiles
    order : { 'nested' | 'ring' }
        Order of HEALPix pixels

    Returns
    -------
    neigh : `~numpy.ndarray`
        Array giving the neighbours starting SW and rotating clockwise. This has
        one extra dimension compared to ``healpix_index`` - the first dimension -
        which is set to 8. For example if healpix_index has shape (2, 3),
        ``neigh`` has shape (8, 2, 3).
    """

    healpix_index = np.asarray(healpix_index, dtype=np.int64)
    nside = int(nside)

    shape = np.shape(healpix_index)

    healpix_index = healpix_index.ravel()

    _validate_healpix_index('healpix_index', healpix_index, nside)
    _validate_nside(nside)
    order = _validate_order(order)

    neigh = core_cython.neighbours(healpix_index, nside, order)
    return _restore_shape(neigh, shape=(8,) + shape)


def healpix_cone_search(lon, lat, radius, nside, order='ring'):
    """
    Find all the HEALPix pixels within a given radius of a longitude/latitude.

    Note that this returns all pixels that overlap, including partially, with
    the search cone. This function can only be used for a single lon/lat pair at
    a time, since different calls to the function may result in a different
    number of matches.

    Parameters
    ----------
    lon, lat : :class:`~astropy.units.Quantity`
        The longitude and latitude to search around
    radius : :class:`~astropy.units.Quantity`
        The search radius
    nside : int
        Number of pixels along the side of each of the 12 top-level HEALPix tiles
    order : { 'nested' | 'ring' }
        Order of HEALPix pixels

    Returns
    -------
    healpix_index : `~numpy.ndarray`
        1-D array with all the matching HEALPix pixel indices.
    """

    lon = float(lon.to(u.deg).value)
    lat = float(lat.to(u.deg).value)
    radius = float(radius.to(u.deg).value)
    nside = int(nside)

    _validate_nside(nside)
    order = _validate_order(order)

    return core_cython.healpix_cone_search(lon, lat, radius, nside, order)


def boundaries_lonlat(healpix_index, step, nside, order='ring'):
    """
    Return the longitude and latitude of the edges of HEALPix pixels

    This returns the longitude and latitude of points along the edge of each
    HEALPIX pixel. The number of points returned for each pixel is ``4 * step``,
    so setting ``step`` to 1 returns just the corners.

    Parameters
    ----------
    healpix_index : `~numpy.ndarray`
        1-D array of HEALPix pixels
    step : int
        The number of steps to take along each edge.
    nside : int
        Number of pixels along the side of each of the 12 top-level HEALPix tiles
    order : { 'nested' | 'ring' }
        Order of HEALPix pixels

    Returns
    -------
    lon, lat : :class:`~astropy.units.Quantity`
        The longitude and latitude, as 2-D arrays where the first dimension is
        the same as the ``healpix_index`` input, and the second dimension has
        size ``4 * step``.
    """

    healpix_index = np.asarray(healpix_index, dtype=np.int64)
    step = int(step)

    if step < 1:
        raise ValueError('step must be at least 1')

    # PERF: this could be optimized by writing a Cython routine to do this to
    # avoid allocating temporary arrays

    frac = np.linspace(0., 1., step + 1)[:-1]
    dx = np.hstack([1 - frac, np.repeat(0, step), frac, np.repeat(1, step)])
    dy = np.hstack([np.repeat(1, step), 1 - frac, np.repeat(0, step), frac])

    healpix_index, dx, dy = np.broadcast_arrays(healpix_index.reshape(-1, 1), dx, dy)

    lon, lat = healpix_to_lonlat(healpix_index.ravel(), nside, dx.ravel(), dy.ravel(), order=order)

    lon = lon.reshape(-1, 4 * step)
    lat = lat.reshape(-1, 4 * step)

    return lon, lat
