# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
healpy compatible interface.
"""
from __future__ import absolute_import, division, print_function, unicode_literals

__all__ = [
    'nside2resol',
    'pix2ang',
    'ang2pix',
]


def nside2resol(nside, arcmin=False):
    import healpy as hp
    return hp.nside2resol(nside, arcmin)


def pix2ang(nside, ipix, nest=False):
    """Pixel to angles.

    Parameters
    ----------
    nside : int
        HEALPix ``nside`` parameter
    ipix : int or array-like
        HEALPix pixel indices
    nest : bool
        Use ``nested`` ordering if True, else ``ring`` ordering.

    Returns
    -------
    theta, phi : `~numpy.array`
        Colatitude and longitude in radians
    """
    import healpy as hp
    theta, phi = hp.pix2ang(nside, ipix, nest)
    return theta, phi


def ang2pix(nside, theta, phi, nest=False):
    import healpy as hp
    ipix = hp.ang2pix(nside, theta, phi, nest)
    return ipix
