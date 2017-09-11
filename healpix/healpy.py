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


def nside2npix(n_side):
    return n_side_to_n_pix(n_side)


def npix2nside(n_pix):
    return n_pix_to_n_side(n_pix)


def pix2ang(nside, ipix, nest=False):
    ipix = np.atleast_1d(ipix).astype(int)
    lon, lat = healpix_to_lonlat(ipix, nside, 1 - int(nest))
    return np.pi / 2 - lat, lon


def ang2pix(nside, theta, phi, nest=False):
    lat = np.pi / 2. - np.atleast_1d(theta)
    lon = np.atleast_1d(phi)
    return lonlat_to_healpix(lon, lat, nside, 1 - int(nest))
