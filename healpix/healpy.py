# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
healpy compatible interface.
"""
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from astropy import units as u

from ._healpix import lonlat_to_healpix, healpix_to_lonlat
from .core import (n_side_to_resolution, n_side_to_pixel_area,
                   n_side_to_n_pix, n_pix_to_n_side)

__all__ = ['nside2resol', 'nside2pixarea', 'nside2npix', 'npix2nside',
           'pix2ang', 'ang2pix']


def nside2resol(nside, arcmin=False):
    resolution = n_side_to_resolution(nside)
    if arcmin:
        return resolution.to(u.arcmin).value
    else:
        return resolution.to(u.rad).value


def nside2pixarea(nside, degrees=False):
    area = n_side_to_pixel_area(nside)
    if degrees:
        return area.to(u.deg ** 2).value
    else:
        return area.to(u.sr).value


nside2npix = n_side_to_n_pix

npix2nside = n_pix_to_n_side


def pix2ang(nside, ipix, nest=False):
    ipix = np.atleast_1d(ipix).astype(np.int32)
    phi, theta = healpix_to_lonlat(ipix, nside, 1 - int(nest))
    return theta, phi


def ang2pix(nside, theta, phi, nest=False):
    theta = np.atleast_1d(theta)
    phi = np.atleast_1d(phi)
    return lonlat_to_healpix(phi, theta, nside, 1 - int(nest))
