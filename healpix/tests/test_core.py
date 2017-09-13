from __future__ import print_function, division

import pytest

import numpy as np
from numpy.testing import assert_allclose, assert_equal

from astropy import units as u
from astropy.coordinates import Longitude, Latitude

from ..core import (n_side_to_pixel_area, n_side_to_resolution,
                    n_side_to_n_pix, n_pix_to_n_side, healpix_to_lonlat,
                    lonlat_to_healpix, healpix_with_offset_to_lonlat,
                    lonlat_to_healpix_with_offset, bilinear_interpolation,
                    healpix_neighbours)


def test_n_side_to_pixel_area():
    resolution = n_side_to_pixel_area(256)
    assert_allclose(resolution.value, 1.5978966540475428e-05)
    assert resolution.unit == u.sr


def test_n_side_to_resolution():
    resolution = n_side_to_resolution(256)
    assert_allclose(resolution.value, 13.741945647269624)
    assert resolution.unit == u.arcmin


def test_n_side_to_n_pix():
    npix = n_side_to_n_pix(4)
    assert npix == 192

    npix = n_side_to_n_pix([4, 4])
    assert_equal(npix, 192)

    with pytest.raises(ValueError) as exc:
        n_side_to_n_pix(15)
    assert exc.value.args[0] == 'n_side should be a power of two'


def test_npix2nside():
    nside = n_pix_to_n_side(192)
    assert nside == 4

    nside = n_pix_to_n_side([192, 192])
    assert_equal(nside, 4)

    with pytest.raises(ValueError) as exc:
        n_pix_to_n_side(7)
    assert exc.value.args[0] == 'Number of pixels should be divisible by 12'

    with pytest.raises(ValueError) as exc:
        n_pix_to_n_side(12 * 8 * 9)
    assert exc.value.args[0] == 'Number of pixels is not of the form 12 * n_side ** 2'


# For the following tests, the numerical accuracy of this function is already
# tested in test_cython_api.py, so we focus here on functionality specific to
# the high-level functions.


@pytest.mark.parametrize('order', ['nested', 'ring'])
def test_healpix_to_lonlat(order):

    lon, lat = healpix_to_lonlat([1, 2, 3], 4, order=order)

    assert isinstance(lon, Longitude)
    assert isinstance(lat, Latitude)

    index = lonlat_to_healpix(lon, lat, 4, order=order)

    assert_equal(index, [1, 2, 3])


@pytest.mark.parametrize('order', ['nested', 'ring'])
def test_healpix_with_to_lonlat(order):

    lon, lat = healpix_with_offset_to_lonlat([1, 2, 3],
                                             [0.1, 0.2, 0.3],
                                             [0.5, 0.4, 0.7], 4, order=order)

    assert isinstance(lon, Longitude)
    assert isinstance(lat, Latitude)

    index, dx, dy = lonlat_to_healpix_with_offset(lon, lat, 4, order=order)

    assert_equal(index, [1, 2, 3])
    assert_allclose(dx, [0.1, 0.2, 0.3])
    assert_allclose(dy, [0.5, 0.4, 0.7])


def test_healpix_to_lonlat_invalid():

    with pytest.raises(ValueError) as exc:
        lon, lat = healpix_to_lonlat([-1, 2, 3], 4)
    assert exc.value.args[0] == 'healpix_index should be in the range [0:192]'

    with pytest.raises(ValueError) as exc:
        lon, lat = healpix_to_lonlat([1, 2, 3], 5)
    assert exc.value.args[0] == 'n_side should be a power of two'

    with pytest.raises(ValueError) as exc:
        lon, lat = healpix_to_lonlat([1, 2, 3], 4, order='banana')
    assert exc.value.args[0] == "order should be 'nested' or 'ring'"


def test_healpix_with_offset_to_lonlat_invalid():

    dx, dy = [0.1, 0.2, 0.3], [0.5, 0.3, 0.2]

    with pytest.raises(ValueError) as exc:
        lon, lat = healpix_with_offset_to_lonlat([-1, 2, 3], dx, dy, 4)
    assert exc.value.args[0] == 'healpix_index should be in the range [0:192]'

    with pytest.raises(ValueError) as exc:
        lon, lat = healpix_with_offset_to_lonlat([1, 2, 3], dx, dy, 5)
    assert exc.value.args[0] == 'n_side should be a power of two'

    with pytest.raises(ValueError) as exc:
        lon, lat = healpix_with_offset_to_lonlat([1, 2, 3], dx, dy, 4, order='banana')
    assert exc.value.args[0] == "order should be 'nested' or 'ring'"

    with pytest.raises(ValueError) as exc:
        lon, lat = healpix_with_offset_to_lonlat([1, 2, 3], [-0.1, 0.4, 0.5], dy, 4)
    assert exc.value.args[0] == 'dx should be in the range [0:1]'

    with pytest.raises(ValueError) as exc:
        lon, lat = healpix_with_offset_to_lonlat([1, 2, 3], dx, [-0.1, 0.4, 0.5], 4)
    assert exc.value.args[0] == 'dy should be in the range [0:1]'


@pytest.mark.parametrize('order', ['nested', 'ring'])
def test_bilinear_interpolation(order):
    values = np.ones(192) * 3
    result = bilinear_interpolation([1, 3, 4] * u.deg, [3, 2, 6] * u.deg,
                                    values, order=order)
    assert_allclose(result, [3, 3, 3])


def test_bilinear_interpolation_invalid():

    values = np.ones(133)
    with pytest.raises(ValueError) as exc:
        bilinear_interpolation([1, 3, 4] * u.deg, [3, 2, 6] * u.deg, values)
    assert exc.value.args[0] == 'Number of pixels should be divisible by 12'

    values = np.ones(192)
    with pytest.raises(ValueError) as exc:
        bilinear_interpolation([1, 3, 4] * u.deg, [3, 2, 6] * u.deg,
                               values, order='banana')
    assert exc.value.args[0] == "order should be 'nested' or 'ring'"


@pytest.mark.parametrize('order', ['nested', 'ring'])
def test_healpix_neighbors(order):

    neighbours = healpix_neighbours([1, 2, 3], 4, order=order)

    if order == 'nested':
        expected = [[90, 69, 0],
                    [0, 71, 2],
                    [2, 77, 8],
                    [3, 8, 9],
                    [6, 9, 12],
                    [4, 3, 6],
                    [94, 1, 4],
                    [91, 0, 1]]
    else:
        expected = [[16, 19, 22],
                    [6, 8, 10],
                    [5, 7, 9],
                    [0, 1, 2],
                    [3, 0, 1],
                    [2, 3, 0],
                    [8, 10, 4],
                    [7, 9, 11]]

    assert_equal(neighbours, expected)


def test_healpix_neighbors_invalid():

    with pytest.raises(ValueError) as exc:
        healpix_neighbours([-1, 2, 3], 4)
    assert exc.value.args[0] == 'healpix_index should be in the range [0:192]'

    with pytest.raises(ValueError) as exc:
        healpix_neighbours([1, 2, 3], 5)
    assert exc.value.args[0] == 'n_side should be a power of two'

    with pytest.raises(ValueError) as exc:
        healpix_neighbours([1, 2, 3], 4, order='banana')
    assert exc.value.args[0] == "order should be 'nested' or 'ring'"
