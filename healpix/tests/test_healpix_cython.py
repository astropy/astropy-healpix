import pytest

import numpy as np
from numpy.testing import assert_equal, assert_allclose

from .. import _healpix


@pytest.mark.parametrize('n_side_power', range(0, 9))
def test_roundtrip_nested_no_offsets(n_side_power):
    n_side = 2**n_side_power
    index = np.arange(12 * n_side ** 2).astype(np.int32)
    lon, lat = _healpix.nested_to_lonlat(index, n_side)
    index_new = _healpix.lonlat_to_nested(lon, lat, n_side)
    assert_equal(index, index_new)


@pytest.mark.parametrize('n_side_power', range(0, 9))
def test_roundtrip_nested_with_offsets(n_side_power):
    n_side = 2**n_side_power
    index = np.arange(12 * n_side ** 2).astype(np.int32)
    dx = np.random.random(index.shape)
    dy = np.random.random(index.shape)
    lon, lat = _healpix.nested_with_offset_to_lonlat(index, dx, dy, n_side)
    index_new, dx_new, dy_new = _healpix.lonlat_to_nested_with_offset(lon, lat, n_side)
    assert_equal(index, index_new)
    assert_allclose(dx, dx_new, atol=1e-10)
    assert_allclose(dy, dy_new, atol=1e-10)


@pytest.mark.parametrize('n_side_power', range(0, 9))
def test_roundtrip_ring_no_offsets(n_side_power):
    n_side = 2**n_side_power
    index = np.arange(12 * n_side ** 2).astype(np.int32)
    lon, lat = _healpix.ring_to_lonlat(index, n_side)
    index_new = _healpix.lonlat_to_ring(lon, lat, n_side)
    assert_equal(index, index_new)


@pytest.mark.parametrize('n_side_power', range(0, 9))
def test_roundtrip_ring_with_offsets(n_side_power):
    n_side = 2**n_side_power
    index = np.arange(12 * n_side ** 2).astype(np.int32)
    dx = np.random.random(index.shape)
    dy = np.random.random(index.shape)
    lon, lat = _healpix.ring_with_offset_to_lonlat(index, dx, dy, n_side)
    index_new, dx_new, dy_new = _healpix.lonlat_to_ring_with_offset(lon, lat, n_side)
    assert_equal(index, index_new)
    assert_allclose(dx, dx_new, atol=1e-10)
    assert_allclose(dy, dy_new, atol=1e-10)
