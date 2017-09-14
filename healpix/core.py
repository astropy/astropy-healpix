from __future__ import print_function, division

import math

import numpy as np

from astropy import units as u
from astropy.coordinates import Longitude, Latitude

from . import core_cython

__all__ = ['nside_to_pixel_area', 'nside_to_pixel_resolution', 'nside_to_npix',
           'npix_to_nside', 'lonlat_to_healpix', 'healpix_to_lonlat',
           'interpolate_bilinear', 'healpix_neighbors']

ORDER = {'nested': 0,
         'ring': 1}


def _validate_order(order):
    if order not in ORDER:
        raise ValueError("order should be 'nested' or 'ring'")


def _validate_healpix_index(label, healpix_index, nside):
    npix = nside_to_npix(nside)
    if np.any((healpix_index < 0) | (healpix_index > npix - 1)):
        raise ValueError('{0} should be in the range [0:{1}]'.format(label, npix))


def _validate_offset(label, offset):
    if np.any((offset < 0) | (offset > 1)):
        raise ValueError('d{0} should be in the range [0:1]'.format(label))


def _validate_nside(nside):
    log_2_nside = np.round(np.log2(nside))
    if not np.all(2 ** log_2_nside == nside):
        raise ValueError('nside should be a power of two')


def nside_to_pixel_area(nside):
    """
    Find the area of healpix pixels given the pixel dimensions of one of
    the 12 'top-level' healpix tiles.

    Parameters
    ----------
    nside : int
        The number of pixels on the side of one of the 12 'top-level' healpix tiles.

    Returns
    -------
    pixel_area : :class:`~astropy.units.Quantity`
        The area of the healpix pixels
    """
    nside = np.asanyarray(nside)
    _validate_nside(nside)
    npix = 12 * nside * nside
    pixel_area = 4 * math.pi / npix * u.sr
    return pixel_area


def nside_to_pixel_resolution(nside):
    """
    Find the resolution of healpix pixels given the pixel dimensions of one of
    the 12 'top-level' healpix tiles.

    Parameters
    ----------
    nside : int
        The number of pixels on the side of one of the 12 'top-level' healpix tiles.

    Returns
    -------
    resolution : :class:`~astropy.units.Quantity`
        The resolution of the healpix pixels
    """
    nside = np.asanyarray(nside)
    _validate_nside(nside)
    return (nside_to_pixel_area(nside) ** 0.5).to(u.arcmin)


def nside_to_npix(nside):
    """
    Find the number of pixels corresponding to a healpix resolution.

    Parameters
    ----------
    nside : int
        The number of pixels on the side of one of the 12 'top-level' healpix tiles.

    Returns
    -------
    npix : int
        The number of pixels in the healpix map.
    """
    nside = np.asanyarray(nside)
    _validate_nside(nside)
    return 12 * nside ** 2


def npix_to_nside(npix):
    """
    Find the number of pixels on the side of one of the 12 'top-level' healpix
    tiles given a total number of pixels.

    Parameters
    ----------
    npix : int
        The number of pixels in the healpix map.

    Returns
    -------
    nside : int
        The number of pixels on the side of one of the 12 'top-level' healpix tiles.
    """

    npix = np.asanyarray(npix)

    if not np.all(npix % 12 == 0):
        raise ValueError('Number of pixels should be divisible by 12')

    square_root = np.sqrt(npix / 12)
    if not np.all(square_root ** 2 == npix / 12):
        raise ValueError('Number of pixels is not of the form 12 * nside ** 2')

    return np.round(square_root).astype(int)


def healpix_to_lonlat(healpix_index, nside, dx=None, dy=None, order='nested'):
    """
    Convert healpix indices (optionally with offsets) to longitudes/latitudes.

    If no offsets (``dx`` and ``dy``) are provided, the coordinates will default
    to those at the center of the HEALPix pixels.

    Parameters
    ----------
    healpix_index : `~numpy.ndarray`
        1-D array of healpix indices
    nside : int
        Number of pixels along the side of each of the 12 top-level healpix tiles
    dx, dy : `~numpy.ndarray`, optional
        1-D arrays of offsets inside the healpix pixel, which should be in the
        range [0:1] (0.5 is the center of the healpix pixels)
    order : { 'nested' | 'ring' }, optional
        Order of healpix pixels

    Returns
    -------
    lon : :class:`~astropy.coordinates.Longitude`
        The longitude values
    lat : :class:`~astropy.coordinates.Latitude`
        The latitude values
    """

    if (dx is None and dy is not None) or (dx is not None and dy is None):
        raise ValueError('Either both or neither dx and dy should be specified')

    healpix_index = np.asarray(healpix_index, dtype=np.int64)
    nside = int(nside)

    _validate_healpix_index('healpix_index', healpix_index, nside)
    _validate_nside(nside)
    _validate_order(order)

    if dx is None:
        lon, lat = core_cython.healpix_to_lonlat(healpix_index, nside, ORDER[order])
    else:
        dx = np.asarray(dx, dtype=np.float)
        dy = np.asarray(dy, dtype=np.float)
        _validate_offset('x', dx)
        _validate_offset('y', dy)
        lon, lat = core_cython.healpix_with_offset_to_lonlat(healpix_index, dx, dy, nside, ORDER[order])

    lon = Longitude(lon, unit=u.rad, copy=False)
    lat = Latitude(lat, unit=u.rad, copy=False)

    return lon, lat


def lonlat_to_healpix(lon, lat, nside, return_offsets=False, order='nested'):
    """
    Convert longitudes/latitudes to healpix indices

    This returns only the healpix indices. If you also want to get relative
    offsets inside the pixels, see :func:`lonlat_to_healpix_with_offset`.

    Parameters
    ----------
    lon, lat : :class:`~astropy.units.Quantity`
        The longitude and latitude values as :class:`~astropy.units.Quantity` instances
        with angle units.
    nside : int
        Number of pixels along the side of each of the 12 top-level healpix tiles
    order : { 'nested' | 'ring' }
        Order of healpix pixels
    return_offsets : bool, optional
        If `True`, the returned values are the healpix pixel indices as well as
        ``dx`` and ``dy``, the fractional positions inside the pixels. If
        `False` (the default), only the HEALPix pixel indices is returned.

    Returns
    -------
    healpix_index : `~numpy.ndarray`
        1-D array of healpix indices
    dx, dy : `~numpy.ndarray`
        1-D arrays of offsets inside the healpix pixel in the range [0:1] (0.5
        is the center of the healpix pixels)
    """

    lon = lon.to(u.rad).value.astype(np.float)
    lat = lat.to(u.rad).value.astype(np.float)
    nside = int(nside)

    _validate_nside(nside)
    _validate_order(order)

    if return_offsets:
        return core_cython.lonlat_to_healpix_with_offset(lon, lat, nside, ORDER[order])
    else:
        return core_cython.lonlat_to_healpix(lon, lat, nside, ORDER[order])


def nested_to_ring(nested_index, nside):
    """
    Convert a healpix 'nested' index to a healpix 'ring' index

    Parameters
    ----------
    nested_index : `~numpy.ndarray`
        Healpix index using the 'nested' ordering
    nside : int
        Number of pixels along the side of each of the 12 top-level healpix tiles

    Returns
    -------
    ring_index : `~numpy.ndarray`
        Healpix index using the 'ring' ordering
    """

    nested_index = np.asarray(nested_index, dtype=np.int64)
    nside = int(nside)

    _validate_healpix_index('nested_index', nested_index, nside)
    _validate_nside(nside)

    return core_cython.nested_to_ring(nested_index, nside)


def ring_to_nested(ring_index, nside):
    """
    Convert a healpix 'ring' index to a healpix 'nested' index

    Parameters
    ----------
    ring_index : `~numpy.ndarray`
        Healpix index using the 'ring' ordering
    nside : int
        Number of pixels along the side of each of the 12 top-level healpix tiles

    Returns
    -------
    nested_index : `~numpy.ndarray`
        Healpix index using the 'nested' ordering
    """

    ring_index = np.asarray(ring_index, dtype=np.int64)
    nside = int(nside)

    _validate_healpix_index('ring_index', ring_index, nside)
    _validate_nside(nside)

    return core_cython.ring_to_nested(ring_index, nside)


def interpolate_bilinear(lon, lat, values, order='nested'):
    """
    Interpolate values at specific longitudes/latitudes using bilinear interpolation

    Parameters
    ----------
    lon, lat : :class:`~astropy.units.Quantity`
        The longitude and latitude values as :class:`~astropy.units.Quantity` instances
        with angle units.
    values : `~numpy.ndarray`
        1-D array with the values in each healpix pixel. This should have a
        length of the form 12 * nside ** 2 (and nside is determined
        automatically from this).
    order : { 'nested' | 'ring' }
        Order of healpix pixels

    Returns
    -------
    result : `~numpy.ndarray`
        1-D array of interpolated values
    """

    lon = lon.to(u.rad).value.astype(np.float)
    lat = lat.to(u.rad).value.astype(np.float)
    values = np.asarray(values, dtype=float)

    _validate_order(order)

    # TODO: in future we could potentially support higher-dimensional arrays
    if values.ndim != 1:
        raise ValueError("values should be a 1-dimensional array")

    return core_cython.interpolate_bilinear(lon, lat, values, ORDER[order])


def healpix_neighbors(healpix_index, nside, order='nested'):
    """
    Find all the healpix pixels that are the neighbours of a healpix pixel

    Parameters
    ----------
    healpix_pixel : `~numpy.ndarray`
        1-D array of healpix pixels
    nside : int
        Number of pixels along the side of each of the 12 top-level healpix tiles
    order : { 'nested' | 'ring' }
        Order of healpix pixels

    Returns
    -------
    neighbours : `~numpy.ndarray`
        2-D array with shape (8, N) giving the neighbours starting SW and
        rotating clockwise.
    """

    healpix_index = np.asarray(healpix_index, dtype=np.int64)
    nside = int(nside)

    _validate_healpix_index('healpix_index', healpix_index, nside)
    _validate_nside(nside)
    _validate_order(order)

    return core_cython.healpix_neighbors(healpix_index, nside, ORDER[order])
