from itertools import product

import pytest

import numpy as np
from numpy.testing import assert_equal, assert_allclose

from .. import core_cython


nside_POWERS = range(0, 16)
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


@pytest.mark.parametrize(('order', 'nside_power'), product(ORDERS, nside_POWERS))
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


@pytest.mark.parametrize('nside_power', nside_POWERS)
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


@pytest.mark.parametrize(('order', 'nside_power'), product(ORDERS, nside_POWERS))
def test_healpix_neighbors(order, nside_power):
    # This just makes sure things run, but doesn't check the validity of result
    nside = 2 ** nside_power
    index = get_test_indices(nside)
    neighbours = core_cython.healpix_neighbors(index, nside, order)
    assert np.all(neighbours >= -1) and np.all(neighbours < 12 * nside ** 2)
