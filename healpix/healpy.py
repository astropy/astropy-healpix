# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This submodule provides a healpy-compatible interface.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from astropy import units as u

from .core_cython import lonlat_to_healpix, healpix_to_lonlat
from .core import (nside_to_pixel_resolution, nside_to_pixel_area,
                   nside_to_npix, npix_to_nside, nested_to_ring, ring_to_nested,
                   level_to_nside)

RAD2DEG = 180 / np.pi

__all__ = [
    'nside2resol',
    'nside2pixarea',
    'nside2npix',
    'npix2nside',
    'pix2ang',
    'ang2pix',
    'order2nside',
]


def nside2resol(nside, arcmin=False):
    """Drop-in replacement for healpy `~healpy.pixelfunc.nside2resol`."""
    resolution = nside_to_pixel_resolution(nside)
    if arcmin:
        return resolution.to(u.arcmin).value
    else:
        return resolution.to(u.rad).value


def nside2pixarea(nside, degrees=False):
    """Drop-in replacement for healpy `~healpy.pixelfunc.nside2pixarea`."""
    area = nside_to_pixel_area(nside)
    if degrees:
        return area.to(u.deg ** 2).value
    else:
        return area.to(u.sr).value


def nside2npix(nside):
    """Drop-in replacement for healpy `~healpy.pixelfunc.nside2npix`."""
    return nside_to_npix(nside)


def npix2nside(npix):
    """Drop-in replacement for healpy `~healpy.pixelfunc.npix2nside`."""
    return npix_to_nside(npix)


def order2nside(order):
    """Drop-in replacement for healpy `~healpy.pixelfunc.order2nside`."""
    return level_to_nside(order)


def pix2ang(nside, ipix, nest=False, lonlat=False):
    """Drop-in replacement for healpy `~healpy.pixelfunc.pix2ang`."""
    ipix = np.atleast_1d(ipix).astype(np.int64, copy=False)
    lon, lat = healpix_to_lonlat(ipix, nside, 1 - int(nest))
    # We use in-place operations below to avoid making temporary arrays - this
    # is safe because the lon/lat arrays returned from healpix_to_lonlat are
    # new and not used elsewhere.
    if lonlat:
        lon *= RAD2DEG
        lat *= RAD2DEG
        return lon, lat
    else:
        np.subtract(0.5 * np.pi, lat, out=lat)
        return lat, lon


def ang2pix(nside, theta, phi, nest=False, lonlat=False):
    """Drop-in replacement for healpy `~healpy.pixelfunc.ang2pix`."""
    # Unlike in pix2ang, we don't use in-place operations since we don't
    # want to modify theta and phi since the user may be using them elsewhere.
    if lonlat:
        lon, lat = np.atleast_1d(theta) / RAD2DEG, np.atleast_1d(phi) / RAD2DEG
    else:
        lat, lon = np.pi / 2. - np.atleast_1d(theta), np.atleast_1d(phi)
    return lonlat_to_healpix(lon, lat, nside, 1 - int(nest))


def nest2ring(nside, ipix):
    """Drop-in replacement for healpy `~healpy.pixelfunc.nest2ring`."""
    ipix = np.atleast_1d(ipix).astype(np.int64, copy=False)
    return nested_to_ring(ipix, nside)


def ring2nest(nside, ipix):
    """Drop-in replacement for healpy `~healpy.pixelfunc.ring2nest`."""
    ipix = np.atleast_1d(ipix).astype(np.int64, copy=False)
    return ring_to_nested(ipix, nside)
