from numpy.testing import assert_allclose
from astropy import units as u

from ..core import n_side_to_resolution


def test_n_side_to_resolution():
    resolution = n_side_to_resolution(256)
    assert_allclose(resolution.value, 13.741945647269624)
    assert resolution.unit == u.arcmin
