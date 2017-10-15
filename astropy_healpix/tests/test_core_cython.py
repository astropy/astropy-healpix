# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, print_function, division

from itertools import product

import pytest

import numpy as np
from numpy.testing import assert_equal, assert_allclose

from astropy import units as u
from astropy.coordinates.angle_utilities import angular_separation

from ..core import nside_to_pixel_resolution
from .. import core_cython

NSIDE_POWERS = range(0, 21)
ORDERS = ('nested', 'ring')


def get_test_indices(nside):
    # For large number of pixels, we only compute a random subset of points
    if nside > 2 ** 8:
        try:
            return np.random.randint(0, 12 * nside ** 2, 12 * 8 ** 2, dtype=np.int64)
        except TypeError:  # Numpy 1.9 and 1.10
            return (np.random.random(12 * 8 ** 2) * (12 * float(nside) ** 2)).astype(np.int64, copy=False)
    else:
        return np.arange(12 * nside ** 2, dtype=np.int64)


# NOTE: we use capfd in all tests here to make sure no errors/warnings are being
# raised by the C code.

@pytest.mark.parametrize(('order', 'nside_power'), product(ORDERS, NSIDE_POWERS))
def test_roundtrip_healpix_no_offsets(order, nside_power, capfd):
    nside = 2 ** nside_power
    index = get_test_indices(nside)
    lon, lat = core_cython.healpix_to_lonlat(index, nside, order)
    index_new = core_cython.lonlat_to_healpix(lon, lat, nside, order)
    assert_equal(index, index_new)
    out, err = capfd.readouterr()
    assert out == "" and err == ""


@pytest.mark.parametrize(('order', 'nside_power'), product(ORDERS, NSIDE_POWERS))
def test_roundtrip_healpix_with_offsets(order, nside_power, capfd):
    nside = 2 ** nside_power
    index = get_test_indices(nside)
    dx = np.random.random(index.shape)
    dy = np.random.random(index.shape)
    lon, lat = core_cython.healpix_with_offset_to_lonlat(index, dx, dy, nside, order)
    index_new, dx_new, dy_new = core_cython.lonlat_to_healpix_with_offset(lon, lat, nside, order)
    assert_equal(index, index_new)
    assert_allclose(dx, dx_new, atol=1e-8)
    assert_allclose(dy, dy_new, atol=1e-8)
    out, err = capfd.readouterr()
    assert out == "" and err == ""


@pytest.mark.parametrize('nside_power', NSIDE_POWERS)
def test_roundtrip_nested_ring(nside_power, capfd):
    nside = 2 ** nside_power
    nested_index = get_test_indices(nside)
    ring_index = core_cython.nested_to_ring(nested_index, nside)
    nested_index_new = core_cython.ring_to_nested(ring_index, nside)
    assert_equal(nested_index, nested_index_new)
    if nside == 1:
        assert np.all(nested_index == ring_index)
    else:
        assert not np.all(nested_index == ring_index)
    out, err = capfd.readouterr()
    assert out == "" and err == ""


@pytest.mark.parametrize(('order', 'nside_power'), product(ORDERS, NSIDE_POWERS))
def test_neighbours(order, nside_power, capfd):
    # This just makes sure things run, but doesn't check the validity of result
    nside = 2 ** nside_power
    index = get_test_indices(nside)
    neighbours = core_cython.neighbours(index, nside, order)
    assert np.all(neighbours >= -1) and np.all(neighbours < 12 * nside ** 2)
    out, err = capfd.readouterr()
    assert out == "" and err == ""


@pytest.mark.parametrize(('order', 'nside_power'), product(ORDERS, NSIDE_POWERS))
def test_healpix_cone_search(order, nside_power, capfd):
    nside = 2 ** nside_power
    lon0, lat0 = 12., 40.
    radius = nside_to_pixel_resolution(nside).to(u.degree).value * 10
    index_inside = core_cython.healpix_cone_search(lon0, lat0, radius, nside, order)
    n_inside = len(index_inside)
    dx = np.array([[0.0, 0.0, 1.0, 1.0]])
    dy = np.array([[0.0, 1.0, 1.0, 0.0]])
    dx = np.repeat(dx, n_inside, axis=0).ravel()
    dy = np.repeat(dy, n_inside, axis=0).ravel()
    index_inside = np.repeat(index_inside, 4).ravel()
    lon, lat = core_cython.healpix_with_offset_to_lonlat(index_inside, dx, dy, nside, order)
    lon, lat = lon.reshape((n_inside, 4)), lat.reshape((n_inside, 4))
    sep = angular_separation(lon0 * u.deg, lat0 * u.deg, lon * u.rad, lat * u.rad)
    sep = sep.min(axis=1)
    assert np.all(sep.to(u.degree).value < radius * 1.05)
    out, err = capfd.readouterr()
    assert out == "" and err == ""


# Regression tests for fixed bugs

def test_regression_healpix_to_lonlat_sqrt():

    # Regression test for a bug that caused the ring index decomposition to fail
    # and return a negative longitude.

    index = np.array([9007199120523263], dtype=np.int64)
    lon, lat = core_cython.healpix_to_lonlat(index, 2**26, order='ring')
    assert_allclose(lon, 6.283185295476241, rtol=1e-14)
    assert_allclose(lat, 0.729727669554970, rtol=1e-14)

    index = np.array([720575940916150240], dtype=np.int64)
    lon, lat = core_cython.healpix_to_lonlat(index, 2**28, order='ring')
    assert_allclose(lon, 6.283185122851909, rtol=1e-14)
    assert_allclose(lat, -0.729727656226966, rtol=1e-14)

    index = np.array([180143985363255292], dtype=np.int64)
    lon, lat = core_cython.healpix_to_lonlat(index, 2**27, order='ring')
    assert_allclose(lon, 6.283185266217880, rtol=1e-14)
    assert_allclose(lat, -0.729727656226966, rtol=1e-14)
