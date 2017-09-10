from __future__ import print_function, division

import math

import numpy as np

from astropy import units as u

__all__ = ['n_side_to_pixel_area', 'n_side_to_resolution',
           'n_side_to_n_pix', 'n_pix_to_n_side']


def check_n_side_valid(n_side):
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
    check_n_side_valid(n_side)
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
    check_n_side_valid(n_side)
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
    check_n_side_valid(n_side)
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
