# Licensed under a 3-clause BSD style license - see LICENSE.rst
from itertools import product

import pytest

import numpy as np
from numpy.testing import assert_allclose, assert_equal

from astropy import units as u
from astropy.coordinates import Longitude, Latitude

from ..core import (nside_to_pixel_area, nside_to_pixel_resolution, pixel_resolution_to_nside,
                    nside_to_npix, npix_to_nside, healpix_to_lonlat,
                    lonlat_to_healpix, interpolate_bilinear_lonlat,
                    neighbours, healpix_cone_search, boundaries_lonlat,
                    level_to_nside, nside_to_level,
                    nested_to_ring, ring_to_nested,
                    level_ipix_to_uniq, uniq_to_level_ipix,
                    bilinear_interpolation_weights)


def test_level_to_nside():
    assert level_to_nside(5) == 2 ** 5
    with pytest.raises(ValueError) as exc:
        level_to_nside(-1)
    assert exc.value.args[0] == 'level must be positive'


def test_nside_to_level():
    assert nside_to_level(1024) == 10
    with pytest.raises(ValueError) as exc:
        nside_to_level(511)
    assert exc.value.args[0] == 'nside must be a power of two'


def test_level_ipix_to_uniq():
    assert 11 + 4*4**0 == level_ipix_to_uniq(0, 11)
    assert 62540 + 4*4**15 == level_ipix_to_uniq(15, 62540)
    with pytest.raises(ValueError) as exc:
        level_ipix_to_uniq(1, 49)
    assert exc.value.args[0] == 'ipix for a specific level must be inferior to npix'


@pytest.mark.parametrize("level", [
    0, 5, 10, 15, 20, 22, 25, 26, 27, 28, 29
])
def test_uniq_to_level_ipix(level):
    npix = 3 << 2*(level + 1)
    # Take 10 pixel indices between 0 and npix - 1
    size = 10

    ipix = np.arange(size, dtype=np.int64) * (npix // size)
    level = np.ones(size) * level

    level_res, ipix_res = uniq_to_level_ipix(level_ipix_to_uniq(level, ipix))
    assert np.all(level_res == level) & np.all(ipix_res == ipix)


def test_nside_to_pixel_area():
    resolution = nside_to_pixel_area(256)
    assert_allclose(resolution.value, 1.5978966540475428e-05)
    assert resolution.unit == u.sr


def test_nside_to_pixel_resolution():
    resolution = nside_to_pixel_resolution(256)
    assert_allclose(resolution.value, 13.741945647269624)
    assert resolution.unit == u.arcmin


def test_pixel_resolution_to_nside():

    # Check the different rounding options
    nside = pixel_resolution_to_nside(13 * u.arcmin, round='nearest')
    assert nside == 256

    nside = pixel_resolution_to_nside(13 * u.arcmin, round='up')
    assert nside == 512

    nside = pixel_resolution_to_nside(13 * u.arcmin, round='down')
    assert nside == 256

    # Check that it works with arrays
    nside = pixel_resolution_to_nside([1e3, 10, 1e-3] * u.deg, round='nearest')
    assert_equal(nside, [1, 8, 65536])

    with pytest.raises(ValueError) as exc:
        pixel_resolution_to_nside(13 * u.arcmin, round='peaches')
    assert exc.value.args[0] == "Invalid value for round: 'peaches'"

    with pytest.raises(AttributeError) as exc:
        pixel_resolution_to_nside(13)
    assert exc.value.args[0] == "'int' object has no attribute 'to'"


def test_nside_to_npix():
    npix = nside_to_npix(4)
    assert npix == 192

    npix = nside_to_npix([4, 4])
    assert_equal(npix, 192)

    with pytest.raises(ValueError) as exc:
        nside_to_npix(15)
    assert exc.value.args[0] == 'nside must be a power of two'


def test_npix_to_nside():
    nside = npix_to_nside(192)
    assert nside == 4

    nside = npix_to_nside([192, 192])
    assert_equal(nside, 4)

    with pytest.raises(ValueError) as exc:
        npix_to_nside(7)
    assert exc.value.args[0] == 'Number of pixels must be divisible by 12'

    with pytest.raises(ValueError) as exc:
        npix_to_nside(12 * 8 * 9)
    assert exc.value.args[0] == 'Number of pixels is not of the form 12 * nside ** 2'


# For the following tests, the numerical accuracy of this function is already
# tested in test_cython_api.py, so we focus here on functionality specific to
# the Python functions.


@pytest.mark.parametrize('order', ['nested', 'ring'])
def test_healpix_to_lonlat(order):
    lon, lat = healpix_to_lonlat([1, 2, 3], 4, order=order)

    assert isinstance(lon, Longitude)
    assert isinstance(lat, Latitude)

    index = lonlat_to_healpix(lon, lat, 4, order=order)

    assert_equal(index, [1, 2, 3])

    lon, lat = healpix_to_lonlat([1, 2, 3], 4,
                                 dx=[0.1, 0.2, 0.3],
                                 dy=[0.5, 0.4, 0.7], order=order)

    assert isinstance(lon, Longitude)
    assert isinstance(lat, Latitude)

    index, dx, dy = lonlat_to_healpix(lon, lat, 4, order=order, return_offsets=True)

    assert_equal(index, [1, 2, 3])
    assert_allclose(dx, [0.1, 0.2, 0.3])
    assert_allclose(dy, [0.5, 0.4, 0.7])


def test_healpix_to_lonlat_invalid():
    dx = [0.1, 0.4, 0.9]
    dy = [0.4, 0.3, 0.2]

    with pytest.warns(RuntimeWarning, match='invalid value'):
        lon, lat = healpix_to_lonlat([-1, 2, 3], 4)

    with pytest.warns(RuntimeWarning, match='invalid value'):
        lon, lat = healpix_to_lonlat([192, 2, 3], 4)

    with pytest.raises(ValueError) as exc:
        lon, lat = healpix_to_lonlat([1, 2, 3], 5)
    assert exc.value.args[0] == 'nside must be a power of two'

    with pytest.raises(ValueError) as exc:
        lon, lat = healpix_to_lonlat([1, 2, 3], 4, order='banana')
    assert exc.value.args[0] == "order must be 'nested' or 'ring'"

    with pytest.raises(ValueError) as exc:
        lon, lat = healpix_to_lonlat([1, 2, 3], 4, dx=[-0.1, 0.4, 0.5], dy=dy)
    assert exc.value.args[0] == 'dx must be in the range [0:1]'

    with pytest.raises(ValueError) as exc:
        lon, lat = healpix_to_lonlat([1, 2, 3], 4, dx=dx, dy=[-0.1, 0.4, 0.5])
    assert exc.value.args[0] == 'dy must be in the range [0:1]'


def test_healpix_to_lonlat_shape():
    lon, lat = healpix_to_lonlat(2, 8)
    assert lon.isscalar and lat.isscalar

    lon, lat = healpix_to_lonlat([[1, 2, 3], [3, 4, 4]], 8)
    assert lon.shape == (2, 3) and lat.shape == (2, 3)

    lon, lat = healpix_to_lonlat([[1], [2], [3]], nside=8, dx=0.2, dy=[[0.1, 0.3]])
    assert lon.shape == (3, 2) and lat.shape == (3, 2)


def test_lonlat_to_healpix_shape():
    healpix_index = lonlat_to_healpix(2 * u.deg, 3 * u.deg, 8)
    assert np.can_cast(healpix_index, np.int64)

    lon, lat = np.ones((2, 4)) * u.deg, np.zeros((2, 4)) * u.deg
    healpix_index = lonlat_to_healpix(lon, lat, 8)
    assert healpix_index.shape == (2, 4)

    healpix_index, dx, dy = lonlat_to_healpix(2 * u.deg, 3 * u.deg, 8, return_offsets=True)
    assert np.can_cast(healpix_index, np.int64)
    assert isinstance(dx, float)
    assert isinstance(dy, float)

    lon, lat = np.ones((2, 4)) * u.deg, np.zeros((2, 4)) * u.deg
    healpix_index, dx, dy = lonlat_to_healpix(lon, lat, 8, return_offsets=True)
    assert healpix_index.shape == (2, 4)
    assert dx.shape == (2, 4)
    assert dy.shape == (2, 4)


def test_lonlat_to_healpix_invalid():
    """Check that if we pass NaN values for example, the index is set to -1"""
    with pytest.warns(RuntimeWarning, match='invalid value'):
        ipix = lonlat_to_healpix(np.nan * u.deg, np.nan * u.deg,
                                 nside=1, order='nested')
    assert ipix == -1


@pytest.mark.parametrize('function', [nested_to_ring, ring_to_nested])
def test_nested_ring_shape(function):
    index = function(1, 8)
    assert np.can_cast(index, np.int64)

    index = function([[1, 2, 3], [2, 3, 4]], 8)
    assert index.shape == (2, 3)


@pytest.mark.parametrize('order', ['nested', 'ring'])
def test_bilinear_interpolation_weights(order):

    indices, weights = bilinear_interpolation_weights(100 * u.deg, 10 * u.deg,
                                                      nside=4, order=order)
    if order == 'nested':
        indices = nested_to_ring(indices, nside=4)
    assert_equal(indices, [76, 77, 60, 59])
    assert_allclose(weights, [0.532723, 0.426179, 0.038815, 0.002283], atol=1e-6)


def test_bilinear_interpolation_weights_invalid():
    with pytest.raises(ValueError) as exc:
        bilinear_interpolation_weights(1 * u.deg, 2 * u.deg, nside=5)
    assert exc.value.args[0] == 'nside must be a power of two'

    with pytest.raises(ValueError) as exc:
        bilinear_interpolation_weights(3 * u.deg, 4 * u.deg,
                                       nside=4, order='banana')
    assert exc.value.args[0] == "order must be 'nested' or 'ring'"


def test_bilinear_interpolation_weights_shape():

    indices, weights = bilinear_interpolation_weights(3 * u.deg, 4 * u.deg, nside=8)
    assert indices.shape == (4,)
    assert weights.shape == (4,)

    indices, weights = bilinear_interpolation_weights([[1, 2, 3], [2, 3, 4]] * u.deg,
                                                      [[1, 2, 3], [2, 3, 4]] * u.deg, nside=8)
    assert indices.shape == (4, 2, 3)
    assert weights.shape == (4, 2, 3)


@pytest.mark.parametrize('order', ['nested', 'ring'])
def test_interpolate_bilinear_lonlat(order):
    values = np.ones(192) * 3
    result = interpolate_bilinear_lonlat([1, 3, 4] * u.deg, [3, 2, 6] * u.deg,
                                         values, order=order)
    assert_allclose(result, [3, 3, 3])


def test_interpolate_bilinear_invalid():
    values = np.ones(133)
    with pytest.raises(ValueError) as exc:
        interpolate_bilinear_lonlat([1, 3, 4] * u.deg, [3, 2, 6] * u.deg, values)
    assert exc.value.args[0] == 'Number of pixels must be divisible by 12'

    values = np.ones(192)
    with pytest.raises(ValueError) as exc:
        interpolate_bilinear_lonlat([1, 3, 4] * u.deg, [3, 2, 6] * u.deg,
                                    values, order='banana')
    assert exc.value.args[0] == "order must be 'nested' or 'ring'"

    with pytest.warns(RuntimeWarning, match='invalid value'):
        result = interpolate_bilinear_lonlat([0, np.nan] * u.deg,
                                             [0, np.nan] * u.deg, values,
                                             order='nested')
    assert result.shape == (2,)
    assert result[0] == 1
    assert np.isnan(result[1])


def test_interpolate_bilinear_lonlat_shape():

    values = np.ones(192) * 3

    result = interpolate_bilinear_lonlat(3 * u.deg, 4 * u.deg, values)
    assert isinstance(result, float)

    result = interpolate_bilinear_lonlat([[1, 2, 3], [2, 3, 4]] * u.deg,
                                         [[1, 2, 3], [2, 3, 4]] * u.deg, values)
    assert result.shape == (2, 3)

    values = np.ones((192, 50)) * 3

    lon = np.ones((3, 6, 5)) * u.deg
    lat = np.ones((3, 6, 5)) * u.deg

    result = interpolate_bilinear_lonlat(lon, lat, values)
    assert result.shape == (3, 6, 5, 50)


@pytest.mark.parametrize('order', ['nested', 'ring'])
def test_neighbours(order):
    neigh = neighbours([1, 2, 3], 4, order=order)

    if order == 'nested':
        expected = [[0, 71, 2],
                    [2, 77, 8],
                    [3, 8, 9],
                    [6, 9, 12],
                    [4, 3, 6],
                    [94, 1, 4],
                    [91, 0, 1],
                    [90, 69, 0]]
    else:

        expected = [[6, 8, 10],
                    [5, 7, 9],
                    [0, 1, 2],
                    [3, 0, 1],
                    [2, 3, 0],
                    [8, 10, 4],
                    [7, 9, 11],
                    [16, 19, 22]]

    assert_equal(neigh, expected)


def test_neighbours_invalid():
    with pytest.warns(RuntimeWarning, match='invalid value'):
        neighbours([-1, 2, 3], 4)

    with pytest.warns(RuntimeWarning, match='invalid value'):
        neighbours([192, 2, 3], 4)

    with pytest.raises(ValueError) as exc:
        neighbours([1, 2, 3], 5)
    assert exc.value.args[0] == 'nside must be a power of two'

    with pytest.raises(ValueError) as exc:
        neighbours([1, 2, 3], 4, order='banana')
    assert exc.value.args[0] == "order must be 'nested' or 'ring'"


def test_neighbours_shape():
    neigh = neighbours([[1, 2, 3], [2, 3, 4]], 4)
    assert neigh.shape == (8, 2, 3)


@pytest.mark.parametrize('order', ['nested', 'ring'])
def test_healpix_cone_search(order):
    indices = healpix_cone_search(10 * u.deg, 20 * u.deg, 1 * u.deg,
                                  nside=256, order=order)

    assert len(indices) == 80


@pytest.mark.parametrize(('step', 'order'), product([1, 4, 10], ['nested', 'ring']))
def test_boundaries_lonlat(step, order):
    lon, lat = boundaries_lonlat([10, 20, 30], step, 256, order=order)
    assert lon.shape == (3, 4 * step)
    assert lat.shape == (3, 4 * step)
