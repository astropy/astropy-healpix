from itertools import product

import pytest

import numpy as np
from numpy.testing import assert_equal, assert_allclose

from .. import _healpix


N_SIDE_POWERS = range(0, 9)
ORDERS = (0, 1)


@pytest.mark.parametrize(('order', 'n_side_power'), product(ORDERS, N_SIDE_POWERS))
def test_roundtrip_healpix_no_offsets(order, n_side_power):
    n_side = 2 ** n_side_power
    index = np.arange(12 * n_side ** 2).astype(np.int32)
    lon, lat = _healpix.healpix_to_lonlat(index, n_side, order)
    index_new = _healpix.lonlat_to_healpix(lon, lat, n_side, order)
    assert_equal(index, index_new)


@pytest.mark.parametrize(('order', 'n_side_power'), product(ORDERS, N_SIDE_POWERS))
def test_roundtrip_healpix_with_offsets(order, n_side_power):
    n_side = 2 ** n_side_power
    index = np.arange(12 * n_side ** 2).astype(np.int32)
    dx = np.random.random(index.shape)
    dy = np.random.random(index.shape)
    lon, lat = _healpix.healpix_with_offset_to_lonlat(index, dx, dy, n_side, order)
    index_new, dx_new, dy_new = _healpix.lonlat_to_healpix_with_offset(lon, lat, n_side, order)
    assert_equal(index, index_new)
    assert_allclose(dx, dx_new, atol=1e-10)
    assert_allclose(dy, dy_new, atol=1e-10)


@pytest.mark.parametrize('n_side_power', N_SIDE_POWERS)
def test_roundtrip_nested_ring(n_side_power):
    n_side = 2 ** n_side_power
    nested_index = np.arange(12 * n_side ** 2).astype(np.int32)
    ring_index = _healpix.nested_to_ring(nested_index, n_side)
    nested_index_new = _healpix.ring_to_nested(ring_index, n_side)

    assert_equal(nested_index, nested_index_new)

    if n_side == 1:
        assert np.all(nested_index == ring_index)
    else:
        assert not np.all(nested_index == ring_index)


@pytest.mark.parametrize(('order', 'n_side_power'), product(ORDERS, N_SIDE_POWERS))
def test_neighbours_healpix(order, n_side_power):
    # This just makes sure things run, but doesn't check the validity of result
    n_side = 2 ** n_side_power
    index = np.arange(12 * n_side ** 2).astype(np.int32)
    neighbours = _healpix.neighbours_healpix(index, n_side, order)
    assert np.all(neighbours >= -1) and np.all(neighbours < 12 * n_side ** 2)
