from __future__ import print_function, division

import math

import numpy as np

from astropy import units as u
from astropy.coordinates import Longitude, Latitude

from . import _healpix

__all__ = ['n_side_to_pixel_area', 'n_side_to_resolution',
           'n_side_to_n_pix', 'n_pix_to_n_side']

ORDER = {'nested': 0,
         'ring': 1}


def _validate_order(order):
    if order not in ORDER:
        raise ValueError("order should be 'nested' or 'ring'")


def _validate_healpix_index(label, healpix_index, n_side):
    n_pix = n_side_to_n_pix(n_side)
    if np.any((healpix_index < 0) | (healpix_index > n_pix - 1)):
        raise ValueError('{0} should be in the range [0:{1}]'.format(label, n_pix))


def _validate_offset(label, offset):
    if np.any((offset < 0) | (offset > 1)):
        raise ValueError('d{0} should be in the range [0:1]'.format(label))


def _validate_n_side(n_side):
    log_2_n_side = np.round(np.log2(n_side))
    if not np.all(2 ** log_2_n_side == n_side):
        raise ValueError('n_side should be a power of two')


def n_side_to_pixel_area(n_side):
    """
    Find the area of healpix pixels given the pixel dimensions of one of
    the 12 'top-level' healpix tiles.

    Parameters
    ----------
    n_side : int
        The number of pixels on the side of one of the 12 'top-level' healpix tiles.

    Returns
    -------
    pixel_area : :class:`~astropy.units.Quantity`
        The area of the healpix pixels
    """
    n_side = np.asanyarray(n_side)
    _validate_n_side(n_side)
    n_pix = 12 * n_side * n_side
    pixel_area = 4 * math.pi / n_pix * u.sr
    return pixel_area


def n_side_to_resolution(n_side):
    """
    Find the resolution of healpix pixels given the pixel dimensions of one of
    the 12 'top-level' healpix tiles.

    Parameters
    ----------
    n_side : int
        The number of pixels on the side of one of the 12 'top-level' healpix tiles.

    Returns
    -------
    resolution : :class:`~astropy.units.Quantity`
        The resolution of the healpix pixels
    """
    n_side = np.asanyarray(n_side)
    _validate_n_side(n_side)
    return (n_side_to_pixel_area(n_side) ** 0.5).to(u.arcmin)


def n_side_to_n_pix(n_side):
    """
    Find the number of pixels corresponding to a healpix resolution.

    Parameters
    ----------
    n_side : int
        The number of pixels on the side of one of the 12 'top-level' healpix tiles.

    Returns
    -------
    n_pix : int
        The number of pixels in the healpix map.
    """
    n_side = np.asanyarray(n_side)
    _validate_n_side(n_side)
    return 12 * n_side ** 2


def n_pix_to_n_side(n_pix):
    """
    Find the number of pixels on the side of one of the 12 'top-level' healpix
    tiles given a total number of pixels.

    Parameters
    ----------
    n_pix : int
        The number of pixels in the healpix map.

    Returns
    -------
    n_side : int
        The number of pixels on the side of one of the 12 'top-level' healpix tiles.
    """

    n_pix = np.asanyarray(n_pix)

    if not np.all(n_pix % 12 == 0):
        raise ValueError('Number of pixels should be divisible by 12')

    square_root = np.sqrt(n_pix / 12)
    if not np.all(square_root ** 2 == n_pix / 12):
        raise ValueError('Number of pixels is not of the form 12 * n_side ** 2')

    return np.round(square_root).astype(int)


def healpix_to_lonlat(healpix_index, n_side, order='nested'):
    """
    Convert healpix indices to longitudes/latitudes.

    This returns the longitudes/latitudes of the center of the healpix pixels.
    If you also want to provide relative offsets inside the pixels, see
    :func:`healpix_with_offset_to_lonlat`.

    Parameters
    ----------
    healpix_index : `~numpy.ndarray`
        1-D array of healpix indices
    n_side : int
        Number of pixels along the side of each of the 12 top-level healpix tiles
    order : { 'nested' | 'ring' }
        Order of healpix pixels

    Returns
    -------
    lon : `~astropy.coordinates.Longitude`
        The longitude values
    lat : `~astropy.coordinates.Latitude`
        The latitude values
    """

    healpix_index = np.asarray(healpix_index, dtype=np.int64)
    n_side = int(n_side)

    _validate_healpix_index('healpix_index', healpix_index, n_side)
    _validate_n_side(n_side)
    _validate_order(order)

    lon, lat = _healpix.healpix_to_lonlat(healpix_index, n_side, ORDER[order])

    lon = Longitude(lon, unit=u.rad, copy=False)
    lat = Latitude(lat, unit=u.rad, copy=False)

    return lon, lat


def healpix_with_offset_to_lonlat(healpix_index, dx, dy, n_side, order='nested'):
    """
    Convert healpix indices to longitudes/latitudes

    This function takes relative offsets in x and y inside the healpix pixels.
    If you are only interested in the centers of the pixels, see
    `healpixl_to_lonlat`.

    Parameters
    ----------
    healpix_index : `~numpy.ndarray`
        1-D array of healpix indices
    dx, dy : `~numpy.ndarray`
        1-D arrays of offsets inside the healpix pixel, which should be in the
        range [0:1] (0.5 is the center of the healpix pixels)
    n_side : int
        Number of pixels along the side of each of the 12 top-level healpix tiles
    order : { 'nested' | 'ring' }
        Order of healpix pixels

    Returns
    -------
    lon : `~astropy.coordinates.Longitude`
        The longitude values
    lat : `~astropy.coordinates.Latitude`
        The latitude values
    """

    healpix_index = np.asarray(healpix_index, dtype=np.int64)
    dx = np.asarray(dx, dtype=np.float)
    dy = np.asarray(dy, dtype=np.float)
    n_side = int(n_side)

    _validate_healpix_index('healpix_index', healpix_index, n_side)
    _validate_offset('x', dx)
    _validate_offset('y', dy)
    _validate_n_side(n_side)
    _validate_order(order)

    lon, lat = _healpix.healpix_with_offset_to_lonlat(healpix_index, dx, dy, n_side, ORDER[order])

    lon = Longitude(lon, unit=u.rad, copy=False)
    lat = Latitude(lat, unit=u.rad, copy=False)

    return lon, lat


def lonlat_to_healpix(lon, lat, n_side, order='nested'):
    """
    Convert longitudes/latitudes to healpix indices

    This returns only the healpix indices. If you also want to get relative
    offsets inside the pixels, see :func:`lonlat_to_healpix_with_offset`.

    Parameters
    ----------
    lon, lat : `~astropy.units.Quantity`
        The longitude and latitude values as `~astropy.units.Quantity` instances
        with angle units.
    n_side : int
        Number of pixels along the side of each of the 12 top-level healpix tiles
    order : { 'nested' | 'ring' }
        Order of healpix pixels

    Returns
    -------
    healpix_index : `~numpy.ndarray`
        1-D array of healpix indices
    """

    lon = lon.to(u.rad).value.astype(np.float)
    lat = lat.to(u.rad).value.astype(np.float)
    n_side = int(n_side)

    _validate_n_side(n_side)
    _validate_order(order)

    return _healpix.lonlat_to_healpix(lon, lat, n_side, ORDER[order])


def lonlat_to_healpix_with_offset(lon, lat, n_side, order='nested'):
    """
    Convert longitudes/latitudes to healpix indices

    This returns the healpix indices and relative offsets inside the pixels. If
    you want only the healpix indices, see :func:`lonlat_to_healpix`.

    Parameters
    ----------
    lon, lat : `~astropy.units.Quantity`
        The longitude and latitude values as `~astropy.units.Quantity` instances
        with angle units.
    n_side : int
        Number of pixels along the side of each of the 12 top-level healpix tiles
    order : { 'nested' | 'ring' }
        Order of healpix pixels

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
    n_side = int(n_side)

    _validate_n_side(n_side)
    _validate_order(order)

    return _healpix.lonlat_to_healpix_with_offset(lon, lat, n_side, ORDER[order])


def nested_to_ring(nested_index, n_side):
    """
    Convert a healpix 'nested' index to a healpix 'ring' index

    Parameters
    ----------
    nested_index : `~numpy.ndarray`
        Healpix index using the 'nested' ordering
    n_side : int
        Number of pixels along the side of each of the 12 top-level healpix tiles

    Returns
    -------
    ring_index : `~numpy.ndarray`
        Healpix index using the 'ring' ordering
    """

    nested_index = np.asarray(nested_index, dtype=np.int64)
    n_side = int(n_side)

    _validate_healpix_index('nested_index', nested_index, n_side)
    _validate_n_side(n_side)

    return _healpix.nested_to_ring(nested_index, n_side)


def ring_to_nested(ring_index, n_side):
    """
    Convert a healpix 'ring' index to a healpix 'nested' index

    Parameters
    ----------
    ring_index : `~numpy.ndarray`
        Healpix index using the 'ring' ordering
    n_side : int
        Number of pixels along the side of each of the 12 top-level healpix tiles

    Returns
    -------
    nested_index : `~numpy.ndarray`
        Healpix index using the 'nested' ordering
    """

    ring_index = np.asarray(ring_index, dtype=np.int64)
    n_side = int(n_side)

    _validate_healpix_index('ring_index', ring_index, n_side)
    _validate_n_side(n_side)

    return _healpix.ring_to_nested(ring_index, n_side)


def bilinear_interpolation(lon, lat, values, order='nested'):
    """
    Interpolate values at specific longitudes/latitudes using bilinear interpolation

    Parameters
    ----------
    lon, lat : `~astropy.units.Quantity`
        The longitude and latitude values as `~astropy.units.Quantity` instances
        with angle units.
    values : `~numpy.ndarray`
        1-D array with the values in each healpix pixel. This should have a
        length of the form 12 * n_side ** 2 (and n_side is determined
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

    return _healpix.bilinear_interpolation(lon, lat, values, ORDER[order])


def healpix_neighbours(healpix_index, n_side, order='nested'):
    """
    Find all the healpix pixels that are the neighbours of a healpix pixel

    Parameters
    ----------
    healpix_pixel : `~numpy.ndarray`
        1-D array of healpix pixels
    n_side : int
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
    n_side = int(n_side)

    _validate_healpix_index('healpix_index', healpix_index, n_side)
    _validate_n_side(n_side)
    _validate_order(order)

    return _healpix.healpix_neighbors(healpix_index, n_side, ORDER[order])
