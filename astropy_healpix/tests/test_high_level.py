# Licensed under a 3-clause BSD style license - see LICENSE.rst
import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_equal

from astropy import units as u
from astropy.coordinates import Longitude, Latitude, Galactic, SkyCoord

from ..high_level import HEALPix


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

    def test_level(self):
        assert self.pix.level == 8

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

    def test_healpix_to_xyz(self):
        x, y, z = self.pix.healpix_to_xyz([1, 2, 3])

        assert isinstance(x, np.ndarray)
        assert isinstance(y, np.ndarray)
        assert isinstance(z, np.ndarray)

        index = self.pix.xyz_to_healpix(x, y, z)

        assert_equal(index, [1, 2, 3])

        x, y, z = self.pix.healpix_to_xyz([1, 2, 3],
                                          dx=[0.1, 0.2, 0.3],
                                          dy=[0.5, 0.4, 0.7])

        assert isinstance(x, np.ndarray)
        assert isinstance(y, np.ndarray)
        assert isinstance(z, np.ndarray)

        index, dx, dy = self.pix.xyz_to_healpix(x, y, z, return_offsets=True)

        assert_equal(index, [1, 2, 3])
        assert_allclose(dx, [0.1, 0.2, 0.3])
        assert_allclose(dy, [0.5, 0.4, 0.7])

    def test_nested_to_ring(self):
        nested_index_1 = [1, 3, 22]
        ring_index = self.pix.nested_to_ring(nested_index_1)
        nested_index_2 = self.pix.ring_to_nested(ring_index)
        assert_equal(nested_index_1, nested_index_2)

    def test_bilinear_interpolation_weights(self):
        indices, weights = self.pix.bilinear_interpolation_weights([1, 3, 4] * u.deg,
                                                                   [3, 2, 6] * u.deg)
        assert indices.shape == (4, 3)
        assert weights.shape == (4, 3)

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
        assert exc.value.args[0] == 'values must be an array of length 786432 (got 222)'

    def test_cone_search_lonlat(self):
        lon, lat = 1 * u.deg, 4 * u.deg
        result = self.pix.cone_search_lonlat(lon, lat, 1 * u.deg)
        assert len(result) == 77

    def test_cone_search_lonlat_invalid(self):
        lon, lat = [1, 2] * u.deg, [3, 4] * u.deg
        with pytest.raises(ValueError) as exc:
            self.pix.cone_search_lonlat(lon, lat, 1 * u.deg)
        assert exc.value.args[0] == ('The longitude, latitude and radius must '
                                     'be scalar Quantity objects')

    def test_boundaries_lonlat(self):
        lon, lat = self.pix.boundaries_lonlat([10, 20, 30], 4)
        assert lon.shape == (3, 16)
        assert lat.shape == (3, 16)

    def test_neighbours(self):
        neigh = self.pix.neighbours([10, 20, 30])
        assert neigh.shape == (8, 3)


class TestCelestialHEALPix:

    def setup_class(self):
        self.pix = HEALPix(nside=256, order='nested', frame=Galactic())

    def test_healpix_from_header(self):
        """Test instantiation from a FITS header.

        Notes
        -----
        We don't need to test all possible options, because
        :meth:`~astropy_healpix.HEALPix.from_header` is just a wrapper around
        :meth:`~astropy_healpix.utils.parse_input_healpix_data`, which is
        tested exhaustively in :mod:`~astropy_healpix.tests.test_utils`.
        """

        pix = HEALPix.from_header(
            (np.empty(self.pix.npix), 'G'),
            nested=self.pix.order == 'nested')

        assert pix.nside == self.pix.nside
        assert type(pix.frame) == type(self.pix.frame)  # noqa
        assert pix.order == self.pix.order

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
        values = np.ones(12 * 256 ** 2) * 3
        coord = SkyCoord([1, 2, 3] * u.deg, [4, 3, 1] * u.deg, frame='fk4')
        result = self.pix.interpolate_bilinear_skycoord(coord, values)
        assert_allclose(result, [3, 3, 3])

        # Make sure that coordinate system is correctly taken into account

        values = np.arange(12 * 256 ** 2) * 3
        coord = SkyCoord([1, 2, 3] * u.deg, [4, 3, 1] * u.deg, frame='fk4')

        result1 = self.pix.interpolate_bilinear_skycoord(coord, values)
        result2 = self.pix.interpolate_bilinear_skycoord(coord.icrs, values)

        assert_allclose(result1, result2)

    def test_cone_search_skycoord(self):
        coord = SkyCoord(1 * u.deg, 4 * u.deg, frame='galactic')
        result1 = self.pix.cone_search_skycoord(coord, 1 * u.deg)
        assert len(result1) == 77
        result2 = self.pix.cone_search_skycoord(coord.icrs, 1 * u.deg)
        assert_allclose(result1, result2)

    def test_boundaries_skycoord(self):
        coord = self.pix.boundaries_skycoord([10, 20, 30], 4)
        assert coord.shape == (3, 16)


class TestCelestialHEALPixFrameAsClass(TestCelestialHEALPix):

    def setup_class(self):
        self.pix = HEALPix(nside=256, order='nested', frame=Galactic)


class TestCelestialHEALPixFrameAsString(TestCelestialHEALPix):

    def setup_class(self):
        self.pix = HEALPix(nside=256, order='nested', frame='galactic')


def test_invalid_frame_name():
    with pytest.raises(ValueError, match='Coordinate frame name "foobar"'):
        HEALPix(nside=256, frame='foobar')


def test_invalid_frame_type():
    with pytest.raises(ValueError, match='Coordinate frame must be a'):
        HEALPix(nside=256, frame=('obviously', 'not', 'a', 'frame'))
