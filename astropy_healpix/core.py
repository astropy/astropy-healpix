# Licensed under a 3-clause BSD style license - see LICENSE.rst
import math

import numpy as np

from astropy import units as u
from astropy.coordinates import Longitude, Latitude

from . import _core

__all__ = [
    'nside_to_pixel_area',
    'nside_to_pixel_resolution',
    'pixel_resolution_to_nside',
    'nside_to_npix',
    'npix_to_nside',
    'level_to_nside',
    'nside_to_level',
    'level_ipix_to_uniq',
    'uniq_to_level_ipix',
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
            return [arg.item() for arg in args]
        else:
            return args[0].item()


def _validate_order(order):
    # We also support upper-case, to support directly the values
    # ORDERING = {'RING', 'NESTED'} in FITS headers
    # This is currently undocumented in the docstrings.
    if order == 'nested' or order == 'NESTED':
        return 'nested'
    elif order == 'ring' or order == 'RING':
        return 'ring'
    else:
        raise ValueError("order must be 'nested' or 'ring'")


def _validate_healpix_index(label, healpix_index, nside):
    npix = nside_to_npix(nside)
    if np.any((healpix_index < 0) | (healpix_index > npix - 1)):
        raise ValueError(f'{label} must be in the range [0:{npix}]')


def _validate_offset(label, offset):
    offset = np.asarray(offset)
    if np.any((offset < 0) | (offset > 1)):
        raise ValueError(f'd{label} must be in the range [0:1]')


def _validate_level(level):
    if np.any(level < 0):
        raise ValueError('level must be positive')


def _validate_nside(nside):
    log_2_nside = np.round(np.log2(nside))
    if not np.all(2 ** log_2_nside == nside):
        raise ValueError('nside must be a power of two')


def _validate_npix(level, ipix):
    if not np.all(ipix < (3 << 2*(level + 1))):
        raise ValueError('ipix for a specific level must be inferior to npix')


def level_to_nside(level):
    """
    Find the pixel dimensions of the top-level HEALPix tiles.

    This is given by ``nside = 2**level``.

    Parameters
    ----------
    level : int
        The resolution level

    Returns
    -------
    nside : int
        The number of pixels on the side of one of the 12 'top-level' HEALPix tiles.
    """
    level = np.asarray(level, dtype=np.int64)

    _validate_level(level)
    return 2 ** level


def nside_to_level(nside):
    """
    Find the HEALPix level for a given nside.

    This is given by ``level = log2(nside)``.

    This function is the inverse of `level_to_nside`.

    Parameters
    ----------
    nside : int
        The number of pixels on the side of one of the 12 'top-level' HEALPix tiles.
        Must be a power of two.

    Returns
    -------
    level : int
        The level of the HEALPix cells
    """
    nside = np.asarray(nside, dtype=np.int64)

    _validate_nside(nside)
    return np.log2(nside).astype(np.int64)


def uniq_to_level_ipix(uniq):
    """
    Convert a HEALPix cell uniq number to its (level, ipix) equivalent.

    A uniq number is a 64 bits integer equaling to : ipix + 4*(4**level). Please read
    this `paper <http://ivoa.net/documents/MOC/20140602/REC-MOC-1.0-20140602.pdf>`_
    for more details about uniq numbers.

    Parameters
    ----------
    uniq : int
        The uniq number of a HEALPix cell.

    Returns
    -------
    level, ipix: int, int
        The level and index of the HEALPix cell computed from ``uniq``.
    """
    uniq = np.asarray(uniq, dtype=np.int64)

    level = (np.log2(uniq//4)) // 2
    level = level.astype(np.int64)
    _validate_level(level)

    ipix = uniq - (1 << 2*(level + 1))
    _validate_npix(level, ipix)

    return level, ipix


def level_ipix_to_uniq(level, ipix):
    """
    Convert a level and HEALPix index into a uniq number representing the cell.

    This function is the inverse of `uniq_to_level_ipix`.

    Parameters
    ----------
    level : int
        The level of the HEALPix cell
    ipix : int
        The index of the HEALPix cell

    Returns
    -------
    uniq : int
        The uniq number representing the HEALPix cell.
    """
    level = np.asarray(level, dtype=np.int64)
    ipix = np.asarray(ipix, dtype=np.int64)

    _validate_level(level)
    _validate_npix(level, ipix)

    return ipix + (1 << 2*(level + 1))


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

    See also
    --------
    pixel_resolution_to_nside
    """
    nside = np.asanyarray(nside, dtype=np.int64)
    _validate_nside(nside)
    return (nside_to_pixel_area(nside) ** 0.5).to(u.arcmin)


def pixel_resolution_to_nside(resolution, round='nearest'):
    """Find closest HEALPix nside for a given angular resolution.

    This function is the inverse of `nside_to_pixel_resolution`,
    for the default rounding scheme of ``round='nearest'``.

    If you choose ``round='up'``, you'll get HEALPix pixels that
    have at least the requested resolution (usually a bit better
    due to rounding).

    Pixel resolution is defined as square root of pixel area.

    Parameters
    ----------
    resolution : `~astropy.units.Quantity`
        Angular resolution
    round : {'up', 'nearest', 'down'}
        Which way to round

    Returns
    -------
    nside : int
        The number of pixels on the side of one of the 12 'top-level' HEALPix tiles.
        Always a power of 2.

    Examples
    --------
    >>> from astropy import units as u
    >>> from astropy_healpix import pixel_resolution_to_nside
    >>> pixel_resolution_to_nside(13 * u.arcmin)
    256
    >>> pixel_resolution_to_nside(13 * u.arcmin, round='up')
    512
    """
    resolution = resolution.to(u.rad).value
    pixel_area = resolution * resolution
    npix = 4 * math.pi / pixel_area
    nside = np.sqrt(npix / 12)

    # Now we have to round to the closest ``nside``
    # Since ``nside`` must be a power of two,
    # we first compute the corresponding ``level = log2(nside)`
    # round the level and then go back to nside
    level = np.log2(nside)

    if round == 'up':
        level = np.ceil(level)
    elif round == 'nearest':
        level = np.round(level)
    elif round == 'down':
        level = np.floor(level)
    else:
        raise ValueError(f'Invalid value for round: {round!r}')

    # For very low requested resolution (i.e. large angle values), we
    # return ``level=0``, i.e. ``nside=1``, i.e. the lowest resolution
    # that exists with HEALPix
    level = np.clip(level.astype(int), 0, None)

    return level_to_nside(level)


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
    nside : int or `~numpy.ndarray`
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

    _validate_nside(nside)

    if _validate_order(order) == 'ring':
        func = _core.healpix_ring_to_lonlat
    else:  # _validate_order(order) == 'nested'
        func = _core.healpix_nested_to_lonlat

    if dx is None:
        dx = 0.5
    else:
        _validate_offset('x', dx)
    if dy is None:
        dy = 0.5
    else:
        _validate_offset('y', dy)

    nside = np.asarray(nside, dtype=np.intc)

    lon, lat = func(healpix_index, nside, dx, dy)

    lon = Longitude(lon, unit=u.rad, copy=False)
    lat = Latitude(lat, unit=u.rad, copy=False)

    return lon, lat


def lonlat_to_healpix(lon, lat, nside, return_offsets=False, order='ring'):
    """
    Convert longitudes/latitudes to HEALPix indices

    Parameters
    ----------
    lon, lat : :class:`~astropy.units.Quantity`
        The longitude and latitude values as :class:`~astropy.units.Quantity`
        instances with angle units.
    nside : int or `~numpy.ndarray`
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

    if _validate_order(order) == 'ring':
        func = _core.lonlat_to_healpix_ring
    else:  # _validate_order(order) == 'nested'
        func = _core.lonlat_to_healpix_nested

    nside = np.asarray(nside, dtype=np.intc)

    lon = lon.to_value(u.rad)
    lat = lat.to_value(u.rad)

    healpix_index, dx, dy = func(lon, lat, nside)

    if return_offsets:
        return healpix_index, dx, dy
    else:
        return healpix_index


def nested_to_ring(nested_index, nside):
    """
    Convert a HEALPix 'nested' index to a HEALPix 'ring' index

    Parameters
    ----------
    nested_index : int or `~numpy.ndarray`
        Healpix index using the 'nested' ordering
    nside : int or `~numpy.ndarray`
        Number of pixels along the side of each of the 12 top-level HEALPix tiles

    Returns
    -------
    ring_index : int or `~numpy.ndarray`
        Healpix index using the 'ring' ordering
    """

    nside = np.asarray(nside, dtype=np.intc)

    return _core.nested_to_ring(nested_index, nside)


def ring_to_nested(ring_index, nside):
    """
    Convert a HEALPix 'ring' index to a HEALPix 'nested' index

    Parameters
    ----------
    ring_index : int or `~numpy.ndarray`
        Healpix index using the 'ring' ordering
    nside : int or `~numpy.ndarray`
        Number of pixels along the side of each of the 12 top-level HEALPix tiles

    Returns
    -------
    nested_index : int or `~numpy.ndarray`
        Healpix index using the 'nested' ordering
    """

    nside = np.asarray(nside, dtype=np.intc)

    return _core.ring_to_nested(ring_index, nside)


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

    lon = lon.to_value(u.rad)
    lat = lat.to_value(u.rad)

    _validate_nside(nside)

    nside = np.asarray(nside, dtype=np.intc)

    result = _core.bilinear_interpolation_weights(lon, lat, nside)
    indices = np.stack(result[:4])
    weights = np.stack(result[4:])

    if _validate_order(order) == 'nested':
        indices = ring_to_nested(indices, nside)

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

    _validate_nside(nside)

    nside = np.asarray(nside, dtype=np.intc)

    if _validate_order(order) == 'ring':
        func = _core.neighbours_ring
    else:  # _validate_order(order) == 'nested'
        func = _core.neighbours_nested

    return np.stack(func(healpix_index, nside))


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

    lon = lon.to_value(u.deg)
    lat = lat.to_value(u.deg)
    radius = radius.to_value(u.deg)

    _validate_nside(nside)
    order = _validate_order(order)

    return _core.healpix_cone_search(lon, lat, radius, nside, order)


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
