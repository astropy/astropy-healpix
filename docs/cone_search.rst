Seaching for pixels around a position (cone search)
===================================================

A common operation when using HEALPix maps is to try and find all pixels
that lie within a certain radius of a given longitude/latitude. One way to
do this would be to simply find the longitude/latitude of all pixels in the
HEALPix map then find the spherical distance to the requested longitude
and latitude, but in practice this would be very inefficient for high
resolution HEALPix maps where the number of pixels may become arbitrarily large.

Instead, the :meth:`~healpix.HEALPix.cone_search_lonlat` method can be used to
efficiently find all HEALpix pixels within a certain radius from a
longitude/latitude::

    >>> from astropy import units as u
    >>> from healpix import HEALPix
    >>> hp = HEALPix(nside=16, order='nested')
    >>> hp.cone_search_lonlat(10 * u.deg, 30 * u.deg, radius=10 * u.deg)
    array([1269,  160,  162, 1271, 1270, 1268, 1246, 1247,  138,  139,  161,
           1245,  136,  137,  140,  142,  130,  131, 1239, 1244, 1238, 1241,
           1243, 1265, 1267, 1276, 1273, 1277,  168,  169,  163,  166,  164])

Likewise, if using a :class:`~healpix.CelestialHEALPix` object, you can use the
:meth:`~healpix.CelestialHEALPix.cone_search_skycoord` method to query around
specific celestial coordinates::

    >>> from healpix import CelestialHEALPix
    >>> from astropy.coordinates import Galactic
    >>> hp = CelestialHEALPix(nside=16, order='nested', frame=Galactic())
    >>> from astropy.coordinates import SkyCoord
    >>> coord = SkyCoord('00h42m44.3503s +41d16m08.634s')
    >>> hp.cone_search_skycoord(coord, radius=5 * u.arcmin)
    array([2537, 2540, 2542, 2539, 2538, 2536, 2530, 2531, 2534, 2535, 2541,
           2529, 2532, 2507, 2528, 2506, 2485, 2487, 2486, 2493, 2492, 2495,
           1344, 1345, 2543])