from __future__ import print_function, division

import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_equal

from astropy import units as u
from astropy.coordinates import Longitude, Latitude, Galactic, SkyCoord

from ..high_level import HEALPix, CelestialHEALPix


class TestHEALPix:

    def setup_class(self):
        self.pix = HEALPix(nside=256, order='nested')

    def test_pixel_area(self):
        pixel_area = self.pix.pixel_area
        assert_allclose(pixel_area.value, 1.5978966540475428e-05)
        assert pixel_area.unit == u.sr

    def test_pixel_resolution(self):
        pixel_resolution = self.pix.pixel_resolution
        assert_allclose(pixel_resolution.value, 13.741945647269624)
        assert pixel_resolution.unit == u.arcmin

    def test_npix(self):
        assert self.pix.npix == 12 * 256 ** 2

    # For the following tests, the numerical accuracy of this function is
    # already tested in test_cython_api.py, so we focus here on functionality
    # specific to the high-level functions.

    def test_healpix_to_lonlat(self):

        lon, lat = self.pix.healpix_to_lonlat([1, 2, 3])

        assert isinstance(lon, Longitude)
        assert isinstance(lat, Latitude)

        index = self.pix.lonlat_to_healpix(lon, lat)

        assert_equal(index, [1, 2, 3])

        lon, lat = self.pix.healpix_to_lonlat([1, 2, 3],
                                              dx=[0.1, 0.2, 0.3],
                                              dy=[0.5, 0.4, 0.7])

        assert isinstance(lon, Longitude)
        assert isinstance(lat, Latitude)

        index, dx, dy = self.pix.lonlat_to_healpix(lon, lat, return_offsets=True)

        assert_equal(index, [1, 2, 3])
        assert_allclose(dx, [0.1, 0.2, 0.3])
        assert_allclose(dy, [0.5, 0.4, 0.7])

    def test_nested_to_ring(self):
        nested_index_1 = [1, 3, 22]
        ring_index = self.pix.nested_to_ring(nested_index_1)
        nested_index_2 = self.pix.ring_to_nested(ring_index)
        assert_equal(nested_index_1, nested_index_2)

    def test_interpolate_bilinear_lonlat(self):
        values = np.ones(12 * 256 ** 2) * 3
        result = self.pix.interpolate_bilinear_lonlat([1, 3, 4] * u.deg,
                                                      [3, 2, 6] * u.deg, values)
        assert_allclose(result, [3, 3, 3])

    def test_interpolate_bilinear_lonlat_invalid(self):
        values = np.ones(222) * 3
        with pytest.raises(ValueError) as exc:
            self.pix.interpolate_bilinear_lonlat([1, 3, 4] * u.deg,
                                                 [3, 2, 6] * u.deg, values)
        assert exc.value.args[0] == 'values should be an array of length 786432 (got 222)'


class TestCelestialHEALPix:

    def setup_class(self):
        self.pix = CelestialHEALPix(nside=256, order='nested', frame=Galactic())

    def test_healpix_to_skycoord(self):

        coord = self.pix.healpix_to_skycoord([1, 2, 3])

        assert isinstance(coord, SkyCoord)
        assert isinstance(coord.frame, Galactic)

        # Make sure that the skycoord_to_healpix method converts coordinates
        # to the frame of the HEALPix
        coord = coord.transform_to('fk5')

        index = self.pix.skycoord_to_healpix(coord)

        assert_equal(index, [1, 2, 3])

        coord = self.pix.healpix_to_skycoord([1, 2, 3],
                                             dx=[0.1, 0.2, 0.3],
                                             dy=[0.5, 0.4, 0.7])

        assert isinstance(coord, SkyCoord)
        assert isinstance(coord.frame, Galactic)

        # Make sure that the skycoord_to_healpix method converts coordinates
        # to the frame of the HEALPix
        coord = coord.transform_to('fk5')

        index, dx, dy = self.pix.skycoord_to_healpix(coord, return_offsets=True)

        assert_equal(index, [1, 2, 3])
        assert_allclose(dx, [0.1, 0.2, 0.3])
        assert_allclose(dy, [0.5, 0.4, 0.7])

    def test_interpolate_bilinear_skycoord(self):
        values = np.ones(192) * 3
        coord = SkyCoord([1, 2, 3] * u.deg, [4, 3, 1] * u.deg, frame='fk4')
        result = self.pix.interpolate_bilinear_skycoord(coord, values)
        assert_allclose(result, [3, 3, 3])
