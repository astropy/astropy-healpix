# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, print_function, division

from itertools import product

import pytest
import numpy as np

from numpy.testing import assert_equal, assert_allclose

from .. import healpy as hp_compat

# NOTE: If healpy is installed, we use it in these tests, but healpy is not a
# formal dependency of astropy-healpix.
hp = pytest.importorskip('healpy')

from hypothesis import given, settings, example
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
    assert np.can_cast(ipix, np.int64)

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
@example(nside_pow=29, frac=0.1666666694606345, nest=False, lonlat=False)
@example(nside_pow=27, frac=2./3., nest=True, lonlat=False)
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


@given(nside_pow=integers(0, 29), nest=booleans(),
       x=floats(-1, 1, allow_nan=False, allow_infinity=False).filter(lambda x: abs(x) > 1e-10),
       y=floats(-1, 1, allow_nan=False, allow_infinity=False).filter(lambda y: abs(y) > 1e-10),
       z=floats(-1, 1, allow_nan=False, allow_infinity=False).filter(lambda z: abs(z) > 1e-10))
@settings(max_examples=2000, derandomize=True)
def test_vec2pix(nside_pow, x, y, z, nest):
    nside = 2 ** nside_pow
    ipix1 = hp_compat.vec2pix(nside, x, y, z, nest=nest)
    ipix2 = hp.vec2pix(nside, x, y, z, nest=nest)
    assert ipix1 == ipix2


@given(nside_pow=integers(0, 29), nest=booleans(),
       frac=floats(0, 1, allow_nan=False, allow_infinity=False).filter(lambda x: x < 1))
@settings(max_examples=2000, derandomize=True)
@example(nside_pow=29, frac=0.1666666694606345, nest=False)
def test_pix2vec(nside_pow, frac, nest):
    nside = 2 ** nside_pow
    ipix = int(frac * 12 * nside ** 2)
    xyz1 = hp_compat.pix2vec(nside, ipix, nest=nest)
    xyz2 = hp.pix2vec(nside, ipix, nest=nest)
    assert_allclose(xyz1, xyz2, atol=1e-8)


def test_vec2pix_shape():
    ipix = hp_compat.vec2pix(8, 1., 2., 3.)
    assert np.can_cast(ipix, np.int64)

    ipix = hp_compat.vec2pix(8, [[1., 2.], [3., 4.]], [[5., 6.], [7., 8.]], [[9., 10.], [11., 12.]])
    assert ipix.shape == (2, 2)


def test_pix2vec_shape():
    x, y, z = hp_compat.pix2vec(8, 1)
    assert isinstance(x, float)
    assert isinstance(y, float)
    assert isinstance(z, float)

    x, y, z = hp_compat.pix2vec(8, [[1, 2, 3], [4, 5, 6]])
    assert x.shape == (2, 3)
    assert y.shape == (2, 3)
    assert z.shape == (2, 3)


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
@example(nside_pow=29, frac=0.16666666697710755)
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


def not_at_origin(vec):
    return np.linalg.norm(vec) > 0


@given(vectors=arrays(float, (3,), elements=floats(-1, 1)).filter(not_at_origin),
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


@given(lonlat=booleans(),
       lon=floats(0, 360, allow_nan=False, allow_infinity=False).filter(lambda lon: abs(lon) > 1e-10),
       lat=floats(-90, 90, allow_nan=False, allow_infinity=False).filter(
           lambda lat: abs(lat) < 89.99 and abs(lat) > 1e-10))
@settings(max_examples=2000, derandomize=True)
def test_ang2vec(lon, lat, lonlat):
    if lonlat:
        theta, phi = lon, lat
    else:
        theta, phi = np.pi / 2. - np.radians(lat), np.radians(lon)
    xyz1 = hp_compat.ang2vec(theta, phi, lonlat=lonlat)
    xyz2 = hp.ang2vec(theta, phi, lonlat=lonlat)
    assert_allclose(xyz1, xyz2, atol=1e-10)


# The following fails, need to investigate:
# @example(nside_pow=29, lon=1.0000000028043134e-05, lat=1.000000000805912e-05, nest=False, lonlat=False)
#

@given(nside_pow=integers(0, 28), nest=booleans(), lonlat=booleans(),
       lon=floats(0, 360, allow_nan=False, allow_infinity=False).filter(lambda lon: abs(lon) > 1e-5),
       lat=floats(-90, 90, allow_nan=False, allow_infinity=False).filter(
           lambda lat: abs(lat) < 89.99 and abs(lat) > 1e-5))
@settings(max_examples=500, derandomize=True)
@example(nside_pow=27, lon=1.0000000028043134e-05, lat=-41.81031451395941, nest=False, lonlat=False)
@example(nside_pow=6, lon=1.6345238095238293, lat=69.42254649458224, nest=False, lonlat=False)
@example(nside_pow=15, lon=1.0000000028043134e-05, lat=1.000000000805912e-05, nest=False, lonlat=False)
@example(nside_pow=0, lon=315.0000117809725, lat=1.000000000805912e-05, nest=False, lonlat=False)
@example(nside_pow=0, lon=1.0000000028043134e-05, lat=-41.81031489577861, nest=False, lonlat=False)
@example(nside_pow=0, lon=35.559942143736414, lat=-41.8103252622604, nest=False, lonlat=False)
@example(nside_pow=28, lon=359.9999922886491, lat=-41.81031470486902, nest=False, lonlat=False)
@example(nside_pow=0, lon=1.0000000028043134e-05, lat=-41.81031489577861, nest=False, lonlat=False)
@example(nside_pow=27, lon=1.0000000028043134e-05, lat=-41.81031451395941, nest=False, lonlat=False)
@example(nside_pow=26, lon=359.9999986588955, lat=41.81031489577861, nest=False, lonlat=False)
@example(nside_pow=27, lon=359.999997317791, lat=-41.81031451395943, nest=False, lonlat=False)
@example(nside_pow=27, lon=1.0000000028043134e-05, lat=89.80224636153702, nest=False, lonlat=False)
def test_interp_weights(nside_pow, lon, lat, nest, lonlat):
    nside = 2 ** nside_pow
    if lonlat:
        theta, phi = lon, lat
    else:
        theta, phi = np.pi / 2. - np.radians(lat), np.radians(lon)
    indices1, weights1 = hp_compat.get_interp_weights(nside, theta, phi, nest=nest, lonlat=lonlat)
    indices2, weights2 = hp.get_interp_weights(nside, theta, phi, nest=nest, lonlat=lonlat)

    # Ignore neighbours with weights < 1e-6 - we have to exclude these otherwise
    # in some corner cases there will be different low-probability neighbours.

    keep = weights1 > 1e-6
    indices1, weights1 = indices1[keep], weights1[keep]

    keep = weights2 > 1e-6
    indices2, weights2 = indices2[keep], weights2[keep]

    order1 = np.argsort(indices1)
    order2 = np.argsort(indices2)

    assert_equal(indices1[order1], indices2[order2])
    assert_allclose(weights1[order1], weights2[order2], atol=1e-5)


# Make an array that can be useful up to the highest nside tested below
NSIDE_POW_MAX = 8
VALUES = np.random.random(12 * NSIDE_POW_MAX ** 2)


@given(nside_pow=integers(0, NSIDE_POW_MAX), nest=booleans(), lonlat=booleans(),
       lon=floats(0, 360, allow_nan=False, allow_infinity=False).filter(lambda lon: abs(lon) > 1e-5),
       lat=floats(-90, 90, allow_nan=False, allow_infinity=False).filter(
           lambda lat: abs(lat) < 89.99 and abs(lat) > 1e-5))
@settings(max_examples=500, derandomize=True)
def test_interp_val(nside_pow, lon, lat, nest, lonlat):
    nside = 2 ** nside_pow
    if lonlat:
        theta, phi = lon, lat
    else:
        theta, phi = np.pi / 2. - np.radians(lat), np.radians(lon)
    m = VALUES[:12 * nside ** 2]
    value1 = hp_compat.get_interp_val(m, theta, phi, nest=nest, lonlat=lonlat)
    value2 = hp.get_interp_val(m, theta, phi, nest=nest, lonlat=lonlat)
    assert_allclose(value1, value2, rtol=0.1, atol=1.e-10)
