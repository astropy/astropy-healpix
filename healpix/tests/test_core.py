from __future__ import print_function, division

import pytest

from numpy.testing import assert_allclose, assert_equal
from astropy import units as u

from ..core import (n_side_to_pixel_area, n_side_to_resolution,
                    n_side_to_n_pix, n_pix_to_n_side)


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
