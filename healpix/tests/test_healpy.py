# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from numpy.testing import assert_equal, assert_allclose
from astropy.tests.helper import pytest
from .. import healpy as hp2

try:
    import healpy as hp

    HAS_HEALPY = True
except:
    HAS_HEALPY = False


@pytest.mark.skipif('not HAS_HEALPY')
def test_nside2resol():
    actual = hp2.nside2resol(nside=2)
    expected = hp.nside2resol(nside=2)
    assert_equal(actual, expected)


@pytest.mark.parametrize('nside,theta,phi,nest,ipix', [
    (256, 0.0000000000000000, 0.0000000000000000, True, 65535),
    (256, 0.0000000000000000, 1.2566370614359172, True, 65535),
    (256, 0.0000000000000000, 2.5132741228718345, True, 131071),
    (256, 0.0000000000000000, 3.7699111843077517, True, 196607),
    (256, 0.0000000000000000, 5.0265482457436690, True, 262143),
    (256, 0.0000000000000000, 6.2831853071795862, True, 65535),
])
def test_ang2pix(nside, theta, phi, nest, ipix):
    assert hp2.ang2pix(nside, theta, phi, nest) == ipix


@pytest.mark.parametrize('nside,ipix,nest,theta,phi', [
    (2, 0, True, 1.230959417340775, 0.785398163397448),
    (2, 1, True, 0.841068670567930, 1.178097245096172),
    (2, 2, True, 0.841068670567930, 0.392699081698724),
])
def test_pix2ang(nside, ipix, nest, theta, phi):
    assert_allclose(hp2.pix2ang(nside, ipix, nest), [theta, phi], rtol=1e-10)
