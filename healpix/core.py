from __future__ import print_function, division

import math
from astropy import units as u

__all__ = ['n_side_to_pixel_area', 'n_side_to_resolution']


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
    return (n_side_to_pixel_area(n_side) ** 0.5).to(u.arcmin)
