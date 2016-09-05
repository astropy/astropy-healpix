# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np

__all__ = [
    'nside2npix',
    'npix2nside',
    # 'nside2resol',
    # 'pix2ang',
    # 'ang2pix',
]


def nside2npix(nside):
    """Convert ``nside`` to ``npix``.
    """
    nside = np.asanyarray(nside)

    nsidelog2 = np.round(np.log2(nside))
    if not np.all(2 ** nsidelog2 == nside):
        raise ValueError()

    return 12 * nside * nside


def npix2nside(npix):
    """Convert ``npix`` to ``nside``.
    """
    npix = np.asanyarray(npix)

    if not np.all(npix % 12 == 0):
        raise ValueError()

    square_root = np.sqrt(npix / 12)
    if not np.all(square_root * square_root == npix / 12):
        raise ValueError()

    return np.round(square_root)
