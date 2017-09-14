from __future__ import print_function, division

from astropy.coordinates import ICRS, SkyCoord
from astropy.coordinates.representation import UnitSphericalRepresentation

from .core import (nside_to_pixel_area, nside_to_pixel_resolution, nside_to_npix,
                   healpix_to_lonlat, lonlat_to_healpix, interpolate_bilinear,
                   ring_to_nested, nested_to_ring)


class HEALPix(object):
    """
    A HEALPix pixellization.

    Parameters
    ----------
    nside : int
        Number of pixels along the side of each of the 12 top-level healpix tiles
    order : { 'nested' | 'ring' }
        Order of healpix pixels
    """

    def __init__(self, nside=None, order='nested'):
        if nside is None:
            raise ValueError('nside has not been set')
        self.nside = nside
        self.order = order

    @property
    def pixel_area(self):
        """
        The area of a single HEALPix pixel.
        """
        return nside_to_pixel_area(self.nside)

    @property
    def pixel_resolution(self):
        """
        The resolution of a single HEALPix pixel.
        """
        return nside_to_pixel_resolution(self.nside)

    @property
    def npix(self):
        """
        The number of pixels in the pixellization of the sphere.
        """
        return nside_to_npix(self.nside)

    def healpix_to_lonlat(self, healpix_index, dx=None, dy=None):
        """
        Convert healpix indices (optionally with offsets) to longitudes/latitudes

        Parameters
        ----------
        healpix_index : `~numpy.ndarray`
            1-D array of healpix indices
        dx, dy : `~numpy.ndarray`, optional
            1-D arrays of offsets inside the healpix pixel, which should be in
            the range [0:1] (0.5 is the center of the healpix pixels). If not
            specified, the position at the center of the pixel is used.

        Returns
        -------
        lon : :class:`~astropy.coordinates.Longitude`
            The longitude values
        lat : :class:`~astropy.coordinates.Latitude`
            The latitude values
        """
        return healpix_to_lonlat(healpix_index, self.nside, dx=dx, dy=dy, order=self.order)

    def lonlat_to_healpix(self, lon, lat, return_offsets=False):
        """
        Convert longitudes/latitudes to healpix indices (optionally with offsets)

        Parameters
        ----------
        lon, lat : :class:`~astropy.units.Quantity`
            The longitude and latitude values as :class:`~astropy.units.Quantity` instances
            with angle units.
        return_offsets : bool
            If `True`, the returned values are the healpix pixel as well as
            ``dx`` and ``dy``, the fractional positions inside the pixel. If
            `False` (the default), only the HEALPix pixel is returned.

        Returns
        -------
        healpix_index : `~numpy.ndarray`
            1-D array of healpix indices
        dx, dy : `~numpy.ndarray`
            1-D arrays of offsets inside the healpix pixel in the range [0:1] (0.5
            is the center of the healpix pixels). This is returned if
            ``return_offsets`` is `True`.
        """
        return lonlat_to_healpix(lon, lat, self.nside,
                                 return_offsets=return_offsets, order=self.order)

    def nested_to_ring(self, nested_index):
        """
        Convert a healpix 'nested' index to a healpix 'ring' index

        Parameters
        ----------
        nested_index : `~numpy.ndarray`
            Healpix index using the 'nested' ordering

        Returns
        -------
        ring_index : `~numpy.ndarray`
            Healpix index using the 'ring' ordering
        """
        return nested_to_ring(nested_index, self.nside)

    def ring_to_nested(self, ring_index):
        """
        Convert a healpix 'ring' index to a healpix 'nested' index

        Parameters
        ----------
        ring_index : `~numpy.ndarray`
            Healpix index using the 'ring' ordering

        Returns
        -------
        nested_index : `~numpy.ndarray`
            Healpix index using the 'nested' ordering
        """
        return ring_to_nested(ring_index, self.nside)

    def interpolate_bilinear(self, lon, lat, values):
        """
        Interpolate values at specific longitudes/latitudes using bilinear interpolation

        Parameters
        ----------
        lon, lat : :class:`~astropy.units.Quantity`
            The longitude and latitude values as :class:`~astropy.units.Quantity` instances
            with angle units.
        values : `~numpy.ndarray`
            1-D array with the values in each healpix pixel. This should have a
            length of the form 12 * nside ** 2 (and nside is determined
            automatically from this).

        Returns
        -------
        result : `~numpy.ndarray`
            1-D array of interpolated values
        """
        if len(values) != self.npix:
            raise ValueError('values should be an array of length {0} (got {1})'.format(self.npix, len(values)))
        return interpolate_bilinear(lon, lat, values, order=self.order)


class CelestialHEALPix(HEALPix):
    """
    A HEALPix pixellization of the celestial sphere

    Parameters
    ----------
    nside : int
        Number of pixels along the side of each of the 12 top-level healpix tiles
    order : { 'nested' | 'ring' }
        Order of healpix pixels
    frame : :class:`~astropy.coordinates.BaseCoordinateFrame`, optional
        The coordinate frame of the pixellization, which defaults to ICRS
    """

    def __init__(self, nside=None, order='nested', frame=None):
        super(CelestialHEALPix, self).__init__(nside=nside, order=order)
        # Note that we can't do 'frame or ICRS() here since frames evaluate as False'
        self.frame = frame if frame is not None else ICRS()

    def healpix_to_skycoord(self, healpix_index, dx=None, dy=None):
        """
        Convert healpix indices (optionally with offsets) to celestial coordinates.

        Parameters
        ----------
        healpix_index : `~numpy.ndarray`
            1-D array of healpix indices
        dx, dy : `~numpy.ndarray`, optional
            1-D arrays of offsets inside the healpix pixel, which should be in
            the range [0:1] (0.5 is the center of the healpix pixels). If not
            specified, the position at the center of the pixel is used.

        Returns
        -------
        coord : :class:`~astropy.coordinates.SkyCoord`
            The resulting celestial coordinates
        """
        lon, lat = self.healpix_to_lonlat(healpix_index, dx=dx, dy=dy)
        representation = UnitSphericalRepresentation(lon, lat, copy=False)
        return SkyCoord(self.frame.realize_frame(representation))

    def skycoord_to_healpix(self, skycoord, return_offsets=False):
        """
        Convert celestial coordinates to healpix indices (optionally with offsets).

        Parameters
        ----------
        skycoord : :class:`~astropy.coordinates.SkyCoord`
            The celestial coordinates to convert
        return_offsets : bool
            If `True`, the returned values are the healpix pixel as well as
            ``dx`` and ``dy``, the fractional positions inside the pixel. If
            `False` (the default), only the HEALPix pixel is returned.

        Returns
        -------
        healpix_index : `~numpy.ndarray`
            1-D array of healpix indices
        dx, dy : `~numpy.ndarray`
            1-D arrays of offsets inside the healpix pixel in the range [0:1] (0.5
            is the center of the healpix pixels). This is returned if
            ``return_offsets`` is `True`.
        """
        skycoord = skycoord.transform_to(self.frame)
        representation = skycoord.represent_as(UnitSphericalRepresentation)
        lon, lat = representation.lon, representation.lat
        return self.lonlat_to_healpix(lon, lat, return_offsets=return_offsets)

    def interpolate_bilinear(self, skycoord, values):
        """
        Interpolate values at specific celestial coordinates using bilinear interpolation.

        Parameters
        ----------
        skycoord : :class:`~astropy.coordinates.SkyCoord`
            The celestial coordinates at which to interpolate
        values : `~numpy.ndarray`
            1-D array with the values in each healpix pixel. This should have a
            length of the form 12 * nside ** 2 (and nside is determined
            automatically from this).

        Returns
        -------
        result : `~numpy.ndarray`
            1-D array of interpolated values
        """
        representation = skycoord.represent_as(UnitSphericalRepresentation)
        lon, lat = representation.lon, representation.lat
        return interpolate_bilinear(lon, lat, values)
