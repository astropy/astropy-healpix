from itertools import product

import pytest

import numpy as np
from numpy.testing import assert_equal, assert_allclose

from astropy import units as u
from astropy.coordinates.angle_utilities import angular_separation

from ..core import nside_to_pixel_resolution
from .. import core_cython


NSIDE_POWERS = range(0, 17)
ORDERS = (0, 1)


def get_test_indices(nside):
    # For large number of pixels, we only compute a random subset of points
    if nside > 2 ** 8:
        try:
            return np.random.randint(0, 12 * nside ** 2, 12 * 8 ** 2, dtype=np.int64)
        except TypeError:  # Numpy 1.9 and 1.10
            return (np.random.random(12 * 8 ** 2) * (12 * float(nside) ** 2)).astype(np.int64)
    else:
        return np.arange(12 * nside ** 2, dtype=np.int64)


def test_roundtrip_healpix_no_offsets(order=1, nside_power=14):
    nside = 2 ** nside_power
    index = get_test_indices(nside)
    lon, lat = core_cython.healpix_to_lonlat(index, nside, order)
    index_new = core_cython.lonlat_to_healpix(lon, lat, nside, order)
    assert_equal(index, index_new)


@pytest.mark.parametrize(('order', 'nside_power'), product(ORDERS, NSIDE_POWERS))
def test_roundtrip_healpix_with_offsets(order, nside_power):
    nside = 2 ** nside_power
    index = get_test_indices(nside)
    dx = np.random.random(index.shape)
    dy = np.random.random(index.shape)
    lon, lat = core_cython.healpix_with_offset_to_lonlat(index, dx, dy, nside, order)
    index_new, dx_new, dy_new = core_cython.lonlat_to_healpix_with_offset(lon, lat, nside, order)
    assert_equal(index, index_new)
    assert_allclose(dx, dx_new, atol=1e-10)
    assert_allclose(dy, dy_new, atol=1e-10)


@pytest.mark.parametrize('nside_power', NSIDE_POWERS)
def test_roundtrip_nested_ring(nside_power):
    nside = 2 ** nside_power
    nested_index = get_test_indices(nside)
    ring_index = core_cython.nested_to_ring(nested_index, nside)
    nested_index_new = core_cython.ring_to_nested(ring_index, nside)

    assert_equal(nested_index, nested_index_new)

    if nside == 1:
        assert np.all(nested_index == ring_index)
    else:
        assert not np.all(nested_index == ring_index)


@pytest.mark.parametrize(('order', 'nside_power'), product(ORDERS, NSIDE_POWERS))
def test_healpix_neighbors(order, nside_power):
    # This just makes sure things run, but doesn't check the validity of result
    nside = 2 ** nside_power
    index = get_test_indices(nside)
    neighbours = core_cython.healpix_neighbors(index, nside, order)
    assert np.all(neighbours >= -1) and np.all(neighbours < 12 * nside ** 2)


CASES = list()

# Also add a case with very high resolution to check things work properly
# with indices that would overflow a 32-bit int
CASES.append((0, 16, 0.1))

CASES = [(0, 6, 1)]


@pytest.mark.parametrize(('order', 'nside_power'), product(ORDERS, NSIDE_POWERS))
def test_healpix_cone_search(order, nside_power):
    nside = 2 ** nside_power
    lon0, lat0 = 12., 40.
    radius = nside_to_pixel_resolution(nside).to(u.degree).value * 10
    index_inside = core_cython.healpix_cone_search(lon0, lat0, radius, nside, order, 0)
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
    assert np.all(sep.to(u.degree).value < radius)
