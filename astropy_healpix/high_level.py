# Licensed under a 3-clause BSD style license - see LICENSE.rst
import inspect
import os

from astropy.coordinates import (BaseCoordinateFrame, frame_transform_graph,
                                 SkyCoord, UnitSphericalRepresentation)

from .core import (nside_to_pixel_area, nside_to_pixel_resolution,
                   nside_to_level, nside_to_npix, npix_to_nside,
                   healpix_to_lonlat, lonlat_to_healpix,
                   healpix_to_xyz, xyz_to_healpix,
                   bilinear_interpolation_weights, interpolate_bilinear_lonlat,
                   ring_to_nested, nested_to_ring, healpix_cone_search,
                   boundaries_lonlat, neighbours, _validate_order,
                   _NUMPY_COPY_IF_NEEDED)
from .utils import parse_input_healpix_data

__all__ = ['HEALPix']


NO_FRAME_MESSAGE = """
No frame was specified when initializing HEALPix, so SkyCoord objects cannot be
returned. Either specify a frame when initializing HEALPix or use the {0}
method.
""".replace(os.linesep, ' ').strip()


class NoFrameError(Exception):
    def __init__(self, alternative_method):
        super().__init__(NO_FRAME_MESSAGE.format(alternative_method))


def _get_frame(frame):
    """
    Get a frame instance or None, from the input `frame`, which could be a
    frame name string, frame instance, or frame class.

    Adapted from
    :meth:`astropy.coordinates.sky_coordinate_parsers._get_frame_class`.
    """

    if frame is None or isinstance(frame, BaseCoordinateFrame):
        return frame

    elif isinstance(frame, str):
        frame_cls = frame_transform_graph.lookup_name(frame)
        if frame_cls is None:
            frame_names = frame_transform_graph.get_names()
            raise ValueError('Coordinate frame name "{}" is not a known '
                             'coordinate frame ({})'
                             .format(frame, sorted(frame_names)))
        return frame_cls()

    elif inspect.isclass(frame) and issubclass(frame, BaseCoordinateFrame):
        return frame()

    else:
        raise ValueError("Coordinate frame must be a frame name, frame "
                         "instance, frame class, or None, not a '{}'"
                         .format(frame.__class__.__name__))


class HEALPix:
    """
    A HEALPix pixellization.

    Parameters
    ----------
    nside : int
        Number of pixels along the side of each of the 12 top-level HEALPix tiles
    order : { 'nested' | 'ring' }, optional
        Order of HEALPix pixels. Input string can be lower or upper case.
    frame : str or :class:`~astropy.coordinates.BaseCoordinateFrame`, optional
        The celestial coordinate frame of the pixellization. This can be
        ommitted, in which case the pixellization will not be attached to any
        particular celestial frame, and the methods ending in _skycoord will
        not work (but the _lonlat methods will still work and continue to
        return generic longitudes/latitudes). The frame may be passed as a
        string (such as ``galactic``), as a frame class, or as an instance of
        a frame class.

    Raises
    ------
    ValueError
        If 'order' is not one of the allowed options.
    """

    def __init__(self, nside=None, order='ring', frame=None):
        if nside is None:
            raise ValueError('nside has not been set')
        self.nside = nside
        self.order = _validate_order(order)
        self.frame = _get_frame(frame)

    @classmethod
    def from_header(cls, input_data, field=0, hdu_in=None, nested=None):
        """
        Parameters
        ----------
        input_data : str or `~astropy.io.fits.TableHDU` or `~astropy.io.fits.BinTableHDU` or tuple
            The input data to reproject. This can be:

                * The name of a HEALPIX FITS file
                * A `~astropy.io.fits.TableHDU` or `~astropy.io.fits.BinTableHDU`
                  instance
                * A tuple where the first element is a `~numpy.ndarray` and the
                  second element is a `~astropy.coordinates.BaseCoordinateFrame`
                  instance or a string alias for a coordinate frame.

        hdu_in : int or str, optional
            If ``input_data`` is a FITS file, specifies the HDU to use.
            (the default HDU for HEALPIX data is 1, unlike with image files where
            it is generally 0)
        nested : bool, optional
            The order of the healpix_data, either nested (True) or ring (False).
            If a FITS file is passed in, this is determined from the header.

        Returns
        -------
        healpix : `~astropy_healpix.HEALPix`
            A HEALPix pixellization corresponding to the input data.
        """
        array_in, frame, nested = parse_input_healpix_data(
            input_data, field=field, hdu_in=hdu_in, nested=nested)
        nside = npix_to_nside(len(array_in))
        order = 'nested' if nested else 'ring'
        return cls(nside=nside, order=order, frame=frame)

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

    @property
    def level(self):
        """
        The HEALPix level.
        """
        return nside_to_level(self.nside)

    def healpix_to_lonlat(self, healpix_index, dx=None, dy=None):
        """
        Convert HEALPix indices (optionally with offsets) to longitudes/latitudes

        Parameters
        ----------
        healpix_index : `~numpy.ndarray`
            1-D array of HEALPix indices
        dx, dy : `~numpy.ndarray`, optional
            1-D arrays of offsets inside the HEALPix pixel, which must be in
            the range [0:1] (0.5 is the center of the HEALPix pixels). If not
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
        Convert longitudes/latitudes to HEALPix indices (optionally with offsets)

        Parameters
        ----------
        lon, lat : :class:`~astropy.units.Quantity`
            The longitude and latitude values as :class:`~astropy.units.Quantity` instances
            with angle units.
        return_offsets : bool
            If `True`, the returned values are the HEALPix pixel as well as
            ``dx`` and ``dy``, the fractional positions inside the pixel. If
            `False` (the default), only the HEALPix pixel is returned.

        Returns
        -------
        healpix_index : `~numpy.ndarray`
            1-D array of HEALPix indices
        dx, dy : `~numpy.ndarray`
            1-D arrays of offsets inside the HEALPix pixel in the range [0:1] (0.5
            is the center of the HEALPix pixels). This is returned if
            ``return_offsets`` is `True`.
        """
        return lonlat_to_healpix(lon, lat, self.nside,
                                 return_offsets=return_offsets, order=self.order)

    def healpix_to_xyz(self, healpix_index, dx=None, dy=None):
        """
        Convert HEALPix indices (optionally with offsets) to Cartesian coordinates

        Parameters
        ----------
        healpix_index : `~numpy.ndarray`
            1-D array of HEALPix indices
        dx, dy : `~numpy.ndarray`, optional
            1-D arrays of offsets inside the HEALPix pixel, which must be in
            the range [0:1] (0.5 is the center of the HEALPix pixels). If not
            specified, the position at the center of the pixel is used.

        Returns
        -------
        x : :class:`~numpy.ndarray`
            The x coordinates
        y : :class:`~numpy.ndarray`
            The y coordinates
        z : :class:`~numpy.ndarray`
            The z coordinates
        """
        return healpix_to_xyz(healpix_index, self.nside, dx=dx, dy=dy, order=self.order)

    def xyz_to_healpix(self, x, y, z, return_offsets=False):
        """
        Convert Cartesian coordinates to HEALPix indices (optionally with offsets)

        Parameters
        ----------
        x : :class:`~numpy.ndarray`
            The x coordinates
        y : :class:`~numpy.ndarray`
            The y coordinates
        z : :class:`~numpy.ndarray`
            The z coordinates
        return_offsets : bool
            If `True`, the returned values are the HEALPix pixel as well as
            ``dx`` and ``dy``, the fractional positions inside the pixel. If
            `False` (the default), only the HEALPix pixel is returned.

        Returns
        -------
        healpix_index : `~numpy.ndarray`
            1-D array of HEALPix indices
        dx, dy : `~numpy.ndarray`
            1-D arrays of offsets inside the HEALPix pixel in the range [0:1] (0.5
            is the center of the HEALPix pixels). This is returned if
            ``return_offsets`` is `True`.
        """
        return xyz_to_healpix(x, y, z, self.nside,
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

    def bilinear_interpolation_weights(self, lon, lat):
        """
        Get the four neighbours for each (lon, lat) position and the weight
        associated with each one for bilinear interpolation.

        Parameters
        ----------
        lon, lat : :class:`~astropy.units.Quantity`
            The longitude and latitude values as
            :class:`~astropy.units.Quantity` instances with angle units.

        Returns
        -------
        indices : `~numpy.ndarray`
            2-D array with shape (4, N) giving the four indices to use for the
            interpolation
        weights : `~numpy.ndarray`
            2-D array with shape (4, N) giving the four weights to use for the
            interpolation
        """
        return bilinear_interpolation_weights(lon, lat, self.nside, order=self.order)

    def interpolate_bilinear_lonlat(self, lon, lat, values):
        """
        Interpolate values at specific longitudes/latitudes using bilinear interpolation

        If a position does not have four neighbours, this currently returns NaN.

        Parameters
        ----------
        lon, lat : :class:`~astropy.units.Quantity`
            The longitude and latitude values as :class:`~astropy.units.Quantity` instances
            with angle units.
        values : `~numpy.ndarray`
            1-D array with the values in each HEALPix pixel. This must have a
            length of the form 12 * nside ** 2 (and nside is determined
            automatically from this).

        Returns
        -------
        result : `~numpy.ndarray`
            1-D array of interpolated values
        """
        if len(values) != self.npix:
            raise ValueError('values must be an array of length {} (got {})'
                             .format(self.npix, len(values)))
        return interpolate_bilinear_lonlat(lon, lat, values, order=self.order)

    def cone_search_lonlat(self, lon, lat, radius):
        """
        Find all the HEALPix pixels within a given radius of a longitude/latitude.

        Note that this returns all pixels that overlap, including partially, with
        the search cone. This function can only be used for a single lon/lat pair at
        a time, since different calls to the function may result in a different
        number of matches.

        Parameters
        ----------
        lon, lat : :class:`~astropy.units.Quantity`
            The longitude and latitude to search around
        radius : :class:`~astropy.units.Quantity`
            The search radius

        Returns
        -------
        healpix_index : `~numpy.ndarray`
            1-D array with all the matching HEALPix pixel indices.
        """
        if not lon.isscalar or not lat.isscalar or not radius.isscalar:
            raise ValueError('The longitude, latitude and radius must be '
                             'scalar Quantity objects')
        return healpix_cone_search(lon, lat, radius, self.nside, order=self.order)

    def boundaries_lonlat(self, healpix_index, step):
        """
        Return the longitude and latitude of the edges of HEALPix pixels

        This returns the longitude and latitude of points along the edge of each
        HEALPIX pixel. The number of points returned for each pixel is ``4 * step``,
        so setting ``step`` to 1 returns just the corners.

        Parameters
        ----------
        healpix_index : `~numpy.ndarray`
            1-D array of HEALPix pixels
        step : int
            The number of steps to take along each edge.

        Returns
        -------
        lon, lat : :class:`~astropy.units.Quantity`
            The longitude and latitude, as 2-D arrays where the first dimension is
            the same as the ``healpix_index`` input, and the second dimension has
            size ``4 * step``.
        """
        return boundaries_lonlat(healpix_index, step, self.nside, order=self.order)

    def neighbours(self, healpix_index):
        """
        Find all the HEALPix pixels that are the neighbours of a HEALPix pixel

        Parameters
        ----------
        healpix_index : `~numpy.ndarray`
            Array of HEALPix pixels

        Returns
        -------
        neigh : `~numpy.ndarray`
            Array giving the neighbours starting SW and rotating clockwise. This has
            one extra dimension compared to ``healpix_index`` - the first dimension -
            which is set to 8. For example if healpix_index has shape (2, 3),
            ``neigh`` has shape (8, 2, 3).

        Notes
        -----
        Some HEALPix pixels do not have all 8 neighbours. In these cases, the
        corresponding entry in the returned array has the value of -1 and Numpy
        may print an invalid value warning. To suppress the warning, use
        :class:`numpy.errstate`.
        """
        return neighbours(healpix_index, self.nside, order=self.order)

    def healpix_to_skycoord(self, healpix_index, dx=None, dy=None):
        """
        Convert HEALPix indices (optionally with offsets) to celestial coordinates.

        Note that this method requires that a celestial frame was specified when
        initializing HEALPix. If you don't know or need the celestial frame, you
        can instead use :meth:`~astropy_healpix.HEALPix.healpix_to_lonlat`.

        Parameters
        ----------
        healpix_index : `~numpy.ndarray`
            1-D array of HEALPix indices
        dx, dy : `~numpy.ndarray`, optional
            1-D arrays of offsets inside the HEALPix pixel, which must be in
            the range [0:1] (0.5 is the center of the HEALPix pixels). If not
            specified, the position at the center of the pixel is used.

        Returns
        -------
        coord : :class:`~astropy.coordinates.SkyCoord`
            The resulting celestial coordinates
        """
        if self.frame is None:
            raise NoFrameError("healpix_to_skycoord")
        lon, lat = self.healpix_to_lonlat(healpix_index, dx=dx, dy=dy)
        representation = UnitSphericalRepresentation(lon, lat, copy=_NUMPY_COPY_IF_NEEDED)
        return SkyCoord(self.frame.realize_frame(representation))

    def skycoord_to_healpix(self, skycoord, return_offsets=False):
        """
        Convert celestial coordinates to HEALPix indices (optionally with offsets).

        Note that this method requires that a celestial frame was specified when
        initializing HEALPix. If you don't know or need the celestial frame, you
        can instead use :meth:`~astropy_healpix.HEALPix.lonlat_to_healpix`.

        Parameters
        ----------
        skycoord : :class:`~astropy.coordinates.SkyCoord`
            The celestial coordinates to convert
        return_offsets : bool
            If `True`, the returned values are the HEALPix pixel as well as
            ``dx`` and ``dy``, the fractional positions inside the pixel. If
            `False` (the default), only the HEALPix pixel is returned.

        Returns
        -------
        healpix_index : `~numpy.ndarray`
            1-D array of HEALPix indices
        dx, dy : `~numpy.ndarray`
            1-D arrays of offsets inside the HEALPix pixel in the range [0:1] (0.5
            is the center of the HEALPix pixels). This is returned if
            ``return_offsets`` is `True`.
        """
        if self.frame is None:
            raise NoFrameError("skycoord_to_healpix")
        skycoord = skycoord.transform_to(self.frame)
        representation = skycoord.represent_as(UnitSphericalRepresentation)
        lon, lat = representation.lon, representation.lat
        return self.lonlat_to_healpix(lon, lat, return_offsets=return_offsets)

    def interpolate_bilinear_skycoord(self, skycoord, values):
        """
        Interpolate values at specific celestial coordinates using bilinear interpolation.

        If a position does not have four neighbours, this currently returns NaN.

        Note that this method requires that a celestial frame was specified when
        initializing HEALPix. If you don't know or need the celestial frame, you
        can instead use :meth:`~astropy_healpix.HEALPix.interpolate_bilinear_lonlat`.

        Parameters
        ----------
        skycoord : :class:`~astropy.coordinates.SkyCoord`
            The celestial coordinates at which to interpolate
        values : `~numpy.ndarray`
            1-D array with the values in each HEALPix pixel. This must have a
            length of the form 12 * nside ** 2 (and nside is determined
            automatically from this).

        Returns
        -------
        result : `~numpy.ndarray`
            1-D array of interpolated values
        """
        if self.frame is None:
            raise NoFrameError("interpolate_bilinear_skycoord")
        skycoord = skycoord.transform_to(self.frame)
        representation = skycoord.represent_as(UnitSphericalRepresentation)
        lon, lat = representation.lon, representation.lat
        return self.interpolate_bilinear_lonlat(lon, lat, values)

    def cone_search_skycoord(self, skycoord, radius):
        """
        Find all the HEALPix pixels within a given radius of a celestial position.

        Note that this returns all pixels that overlap, including partially,
        with the search cone. This function can only be used for a single
        celestial position at a time, since different calls to the function may
        result in a different number of matches.

        This method requires that a celestial frame was specified when
        initializing HEALPix.  If you don't know or need the celestial frame,
        you can instead use :meth:`~astropy_healpix.HEALPix.cone_search_lonlat`.

        Parameters
        ----------
        skycoord : :class:`~astropy.coordinates.SkyCoord`
            The celestial coordinates to use for the cone search
        radius : :class:`~astropy.units.Quantity`
            The search radius

        Returns
        -------
        healpix_index : `~numpy.ndarray`
            1-D array with all the matching HEALPix pixel indices.
        """
        if self.frame is None:
            raise NoFrameError("cone_search_skycoord")
        skycoord = skycoord.transform_to(self.frame)
        representation = skycoord.represent_as(UnitSphericalRepresentation)
        lon, lat = representation.lon, representation.lat
        return self.cone_search_lonlat(lon, lat, radius)

    def boundaries_skycoord(self, healpix_index, step):
        """
        Return the celestial coordinates of the edges of HEALPix pixels

        This returns the celestial coordinates of points along the edge of each
        HEALPIX pixel. The number of points returned for each pixel is ``4 * step``,
        so setting ``step`` to 1 returns just the corners.

        This method requires that a celestial frame was specified when
        initializing HEALPix.  If you don't know or need the celestial frame,
        you can instead use :meth:`~astropy_healpix.HEALPix.boundaries_lonlat`.

        Parameters
        ----------
        healpix_index : `~numpy.ndarray`
            1-D array of HEALPix pixels
        step : int
            The number of steps to take along each edge.

        Returns
        -------
        skycoord : :class:`~astropy.coordinates.SkyCoord`
            The celestial coordinates of the HEALPix pixel boundaries
        """
        if self.frame is None:
            raise NoFrameError("boundaries_skycoord")
        lon, lat = self.boundaries_lonlat(healpix_index, step)
        representation = UnitSphericalRepresentation(lon, lat, copy=_NUMPY_COPY_IF_NEEDED)
        return SkyCoord(self.frame.realize_frame(representation))
