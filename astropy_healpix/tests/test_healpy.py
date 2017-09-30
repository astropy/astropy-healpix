# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, print_function, division

from itertools import product

from .six import integer_types

import pytest
import numpy as np

from numpy.testing import assert_equal, assert_allclose

from .. import healpy as hp_compat

# NOTE: If healpy is installed, we use it in these tests, but healpy is not a
# formal dependency of astropy-healpix.
hp = pytest.importorskip('healpy')

from hypothesis import given, settings
from hypothesis.strategies import integers, floats, booleans
from hypothesis.extra.numpy import arrays

NSIDE_VALUES = [2 ** n for n in range(1, 6)]


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


@pytest.mark.parametrize('npix', [12 * 2 ** (2 * n) for n in range(1, 6)])
def test_npix2nside(npix):
    actual = hp_compat.npix2nside(npix)
    expected = hp.npix2nside(npix)
    assert_equal(actual, expected)


# For the test below, we exclude latitudes that fall exactly on the pole or
# the equator since points that fall at exact boundaries are ambiguous.

@given(nside_pow=integers(0, 29), nest=booleans(), lonlat=booleans(),
       lon=floats(0, 360, allow_nan=False, allow_infinity=False).filter(lambda lon: abs(lon) > 1e-10),
       lat=floats(-90, 90, allow_nan=False, allow_infinity=False).filter(
           lambda lat: abs(lat) < 89.99 and abs(lat) > 1e-10))
@settings(max_examples=2000, derandomize=True)
def test_ang2pix(nside_pow, lon, lat, nest, lonlat):
    nside = 2 ** nside_pow
    if lonlat:
        theta, phi = lon, lat
    else:
        theta, phi = np.pi / 2. - np.radians(lat), np.radians(lon)
    ipix1 = hp_compat.ang2pix(nside, theta, phi, nest=nest, lonlat=lonlat)
    ipix2 = hp.ang2pix(nside, theta, phi, nest=nest, lonlat=lonlat)
    assert ipix1 == ipix2


def test_ang2pix_shape():
    ipix = hp_compat.ang2pix(8, 1., 2.)
    assert isinstance(ipix, integer_types)

    ipix = hp_compat.ang2pix(8, [[1., 2.], [3., 4.]], [[1., 2.], [3., 4.]])
    assert ipix.shape == (2, 2)


def test_pix2ang_shape():
    lon, lat = hp_compat.pix2ang(8, 1)
    assert isinstance(lon, float)
    assert isinstance(lat, float)

    lon, lat = hp_compat.pix2ang(8, [[1, 2, 3], [4, 5, 6]])
    assert lon.shape == (2, 3)
    assert lat.shape == (2, 3)


@given(nside_pow=integers(0, 29), nest=booleans(), lonlat=booleans(),
       frac=floats(0, 1, allow_nan=False, allow_infinity=False).filter(lambda x: x < 1))
@settings(max_examples=2000, derandomize=True)
def test_pix2ang(nside_pow, frac, nest, lonlat):
    nside = 2 ** nside_pow
    ipix = int(frac * 12 * nside ** 2)
    theta1, phi1 = hp_compat.pix2ang(nside, ipix, nest=nest, lonlat=lonlat)
    theta2, phi2 = hp.pix2ang(nside, ipix, nest=nest, lonlat=lonlat)
    if lonlat:
        assert_allclose(phi1, phi2, atol=1e-8)
        if abs(phi1) < 90:
            assert_allclose(theta1, theta2, atol=1e-10)
    else:
        assert_allclose(theta1, theta2, atol=1e-8)
        if theta1 > 0:
            assert_allclose(phi1, phi2, atol=1e-10)


@given(nside_pow=integers(0, 29),
       frac=floats(0, 1, allow_nan=False, allow_infinity=False).filter(lambda x: x < 1))
@settings(max_examples=2000, derandomize=True)
def test_nest2ring(nside_pow, frac):
    nside = 2 ** nside_pow
    nest = int(frac * 12 * nside ** 2)
    ring1 = hp_compat.nest2ring(nside, nest)
    ring2 = hp.nest2ring(nside, nest)
    assert ring1 == ring2


@given(nside_pow=integers(0, 29),
       frac=floats(0, 1, allow_nan=False, allow_infinity=False).filter(lambda x: x < 1))
@settings(max_examples=2000, derandomize=True)
def test_ring2nest(nside_pow, frac):
    nside = 2 ** nside_pow
    ring = int(frac * 12 * nside ** 2)
    nest1 = hp_compat.ring2nest(nside, ring)
    nest2 = hp.ring2nest(nside, ring)
    assert nest1 == nest2


@given(nside_pow=integers(0, 29), step=integers(1, 10), nest=booleans(),
       frac=floats(0, 1, allow_nan=False, allow_infinity=False).filter(lambda x: x < 1))
@settings(max_examples=500, derandomize=True)
def test_boundaries(nside_pow, frac, step, nest):
    nside = 2 ** nside_pow
    pix = int(frac * 12 * nside ** 2)
    b1 = hp_compat.boundaries(nside, pix, step=step, nest=nest)
    b2 = hp.boundaries(nside, pix, step=step, nest=nest)
    assert_allclose(b1, b2, atol=1e-8)


def test_boundaries_shape():
    pix = 1
    b1 = hp_compat.boundaries(8, pix, step=4)
    b2 = hp.boundaries(8, pix, step=4)
    assert b1.shape == b2.shape

    pix = [1, 2, 3, 4, 5]
    b1 = hp_compat.boundaries(8, pix, step=4)
    b2 = hp.boundaries(8, pix, step=4)
    assert b1.shape == b2.shape


@given(vectors=arrays(float, (3,), elements=floats(-1, 1)),
       lonlat=booleans(), ndim=integers(0, 4))
@settings(max_examples=500, derandomize=True)
def test_vec2ang(vectors, lonlat, ndim):
    vectors = np.broadcast_to(vectors, (2,) * ndim + (3,))
    theta1, phi1 = hp_compat.vec2ang(vectors, lonlat=lonlat)
    theta2, phi2 = hp.vec2ang(vectors, lonlat=lonlat)
    # Healpy sometimes returns NaNs for phi (somewhat incorrectly)
    phi2 = np.nan_to_num(phi2)
    assert_allclose(theta1, theta1, atol=1e-10)
    assert_allclose(phi1, phi2, atol=1e-10)
