# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This submodule provides a healpy-compatible interface.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from astropy import units as u

from .core_cython import lonlat_to_healpix, healpix_to_lonlat
from .core import (nside_to_pixel_resolution, nside_to_pixel_area,
                   nside_to_npix, npix_to_nside, level_to_nside)

__all__ = ['nside2resol', 'nside2pixarea', 'nside2npix', 'npix2nside',
           'pix2ang', 'ang2pix', 'order2nside']


def nside2resol(nside, arcmin=False):
    resolution = nside_to_pixel_resolution(nside)
    if arcmin:
        return resolution.to(u.arcmin).value
    else:
        return resolution.to(u.rad).value


def nside2pixarea(nside, degrees=False):
    area = nside_to_pixel_area(nside)
    if degrees:
        return area.to(u.deg ** 2).value
    else:
        return area.to(u.sr).value


def nside2npix(nside):
    return nside_to_npix(nside)


def npix2nside(npix):
    return npix_to_nside(npix)


def order2nside(order):
    return level_to_nside(order)


def pix2ang(nside, ipix, nest=False):
    ipix = np.atleast_1d(ipix).astype(np.int64, copy=False)
    lon, lat = healpix_to_lonlat(ipix, nside, 1 - int(nest))
    return np.pi / 2 - lat, lon


def ang2pix(nside, theta, phi, nest=False):
    lat = np.pi / 2. - np.atleast_1d(theta)
    lon = np.atleast_1d(phi)
    return lonlat_to_healpix(lon, lat, nside, 1 - int(nest))
