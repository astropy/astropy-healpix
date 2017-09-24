# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

from itertools import product

import pytest
import numpy as np

from numpy.testing import assert_equal, assert_allclose

from .. import healpy as hp_compat

hp = pytest.importorskip('healpy')

from hypothesis import given, settings
from hypothesis.strategies import integers, floats, booleans

NSIDE_VALUES = [2**n for n in range(1, 6)]


@pytest.mark.parametrize(('nside', 'degrees'), product(NSIDE_VALUES, (False, True)))
def test_nside2pixarea(nside, degrees):
    actual = hp_compat.nside2pixarea(nside=nside, degrees=degrees)
    expected = hp.nside2pixarea(nside=nside, degrees=degrees)
    assert_equal(actual, expected)


@pytest.mark.parametrize(('nside', 'arcmin'), product(NSIDE_VALUES, (False, True)))
def test_nside2resol(nside, arcmin):
    actual = hp_compat.nside2resol(nside=nside, arcmin=arcmin)
    expected = hp.nside2resol(nside=nside, arcmin=arcmin)
    assert_equal(actual, expected)


@pytest.mark.parametrize('nside', NSIDE_VALUES)
def test_nside2npix(nside):
    actual = hp_compat.nside2npix(nside)
    expected = hp.nside2npix(nside)
    assert_equal(actual, expected)


@pytest.mark.parametrize('level', [0, 3, 7])
def test_order2nside(level):
    actual = hp_compat.order2nside(level)
    expected = hp.order2nside(level)
    assert_equal(actual, expected)


@pytest.mark.parametrize('npix', [12 * 2**(2 * n) for n in range(1, 6)])
def test_npix2nside(npix):
    actual = hp_compat.npix2nside(npix)
    expected = hp.npix2nside(npix)
    assert_equal(actual, expected)


@given(nside_pow=integers(0, 16), nest=booleans(), lonlat=booleans(),
       lon=floats(0, 360, allow_nan=False, allow_infinity=False),
       lat=floats(-90, 90, allow_nan=False, allow_infinity=False))
@settings(max_examples=1000)
def test_ang2pix(nside_pow, lon, lat, nest, lonlat):
    nside = 2 ** nside_pow
    if lonlat:
        theta, phi = lon, lat
    else:
        theta, phi = np.pi / 2. - np.radians(lat), np.radians(lon)
    ipix1 = hp_compat.ang2pix(nside, theta, phi, nest=nest, lonlat=lonlat)
    ipix2 = hp.ang2pix(nside, theta, phi, nest=nest, lonlat=lonlat)
    assert ipix1 == ipix2


@given(nside_pow=integers(0, 16), nest=booleans(), lonlat=booleans(),
       frac=floats(0, 1, allow_nan=False, allow_infinity=False).filter(lambda x: x < 1))
@settings(max_examples=1000)
def test_pix2ang(nside_pow, frac, nest, lonlat):
    nside = 2 ** nside_pow
    ipix = int(frac * 12 * nside ** 2)
    theta1, phi1 = hp_compat.pix2ang(nside, ipix, nest=nest, lonlat=lonlat)
    theta2, phi2 = hp.pix2ang(nside, ipix, nest=nest, lonlat=lonlat)
    assert_allclose(phi1, phi2, atol=1e-10)
    assert_allclose(theta1, theta2, atol=1e-8)


@given(nside_pow=integers(0, 16),
       frac=floats(0, 1, allow_nan=False, allow_infinity=False).filter(lambda x: x < 1))
@settings(max_examples=1000)
def test_nest2ring(nside_pow, frac):
    nside = 2 ** nside_pow
    nest = int(frac * 12 * nside ** 2)
    ring1 = hp_compat.nest2ring(nside, nest)
    ring2 = hp.nest2ring(nside, nest)
    assert ring1 == ring2


@given(nside_pow=integers(0, 16),
       frac=floats(0, 1, allow_nan=False, allow_infinity=False).filter(lambda x: x < 1))
@settings(max_examples=1000)
def test_ring2nest(nside_pow, frac):
    nside = 2 ** nside_pow
    nest = int(frac * 12 * nside ** 2)
    nest1 = hp_compat.ring2nest(nside, nest)
    nest2 = hp.ring2nest(nside, nest)
    assert nest1 == nest2
