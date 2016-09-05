# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from numpy.testing import assert_equal, assert_allclose
from astropy.tests.helper import pytest
from ..pix import (
    nside2npix,
    npix2nside,
)


def test_nside2npix():
    npix = nside2npix(4)
    assert npix == 192

    npix = nside2npix([4, 4])
    assert_equal(npix, 192)

    with pytest.raises(ValueError):
        nside2npix(15)


def test_npix2nside():
    nside = npix2nside(192)
    assert nside == 4

    nside = npix2nside([192, 192])
    assert_equal(nside, 4)

    with pytest.raises(ValueError):
        npix2nside(7)

    with pytest.raises(ValueError):
        npix2nside(12 * 8 * 9)


