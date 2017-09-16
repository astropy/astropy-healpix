.. include:: references.txt

.. doctest-skip-all

.. _using:

*************
Using healpix
*************

.. _using-intro:

What is HEALPix?
================

Description of [HEALPix from Wikipedia](https://en.wikipedia.org/wiki/HEALPix):

*HEALPix (sometimes written as Healpix), an acronym for Hierarchical Equal Area
isoLatitude Pixelisation of a 2-sphere, is an algorithm for pixelisation of the
2-sphere, and the associated class of map projections.*

More information about HEALPix can be found here:

* http://healpix.jpl.nasa.gov/
* http://adsabs.harvard.edu/abs/2005ApJ...622..759G
* http://adsabs.harvard.edu/abs/2007MNRAS.381..865C

Using the healpix package
=========================

High-level interface
--------------------

The cleanest way to use the functionality in this package is to make use of the
high-level :class:`~healpix.HEALPix` class (or the
:class:`~healpix.CelestialHEALPix` class if you are working with a HEALPix
pixellization of the sky). The :class:`~healpix.HEALPix` class should be
initialized with the ``nside`` parameter which controls the resolution of the
pixellization - it is the number of pixels on the side of each of the 12 top-level
HEALPix pixels::

    >>> from healpix import HEALPix
    >>> hp = HEALPix(nside=16)

As described in the references above, HEALPix pixel indices can follow two
different ordering conventions - the *nested* convention and the *ring*
convention. By default, the :class:`~healpix.HEALPix` class assumes the nested ordering
convention, but it is possible to explicitly specify the convention to use using
the ``order`` argument, for example::

    >>> hp = HEALPix(nside=16, order='ring')

or::

    >>> hp = HEALPix(nside=16, order='nested')

Once this class has been set up, you can access various properties and methods
related to the HEALPix pixellization. For example, you can calculate the
number of pixels as well as the pixel area or resolution::

    >>> hp.npix
    3072
    >>> hp.pixel_area
    <Quantity 0.0040906154343617095 sr>
    >>> hp.pixel_resolution
    <Quantity 219.87113035631398 arcmin>

As you can see, when appropriate the properties and the methods on the
:class:`~healpix.HEALPix` class return Astropy high-level classes such as
:class:`~astropy.units.Quantity`, :class:`~astropy.coordinates.SkyCoord`, and so
on.

The :meth:`~healpix.HEALPix.healpix_to_lonlat` method can be used to convert HEALPix indices to
:class:`~astropy.coordinates.Longitude` and
:class:`~astropy.coordinates.Latitude` objects:

    >>> lon, lat = hp.healpix_to_lonlat([1, 442, 2200])
    >>> lon
    <Longitude [ 0.83448555, 1.63624617, 0.4712389 ] rad>
    >>> lat
    <Latitude [ 0.08343009, 0.94842784,-0.78529135] rad>

The :class:`~astropy.coordinates.Longitude` and
:class:`~astropy.coordinates.Latitude` objects are fully-fledged
:class:`~astropy.units.Quantity` objects and also include shortcuts to get
the values in various units::

    >>> lon.hourangle
    array([ 3.1875,  6.25  ,  1.8   ])
    >>> lat.degree
    array([  4.78019185,  54.3409123 , -44.99388015])

Conversely, given longitudes and latitudes as :class:`~astropy.units.Quantity`
objects, it is possible to recover HEALPix pixel indices::

    >>> from astropy import units as u
    >>> hp.lonlat_to_healpix([1, 3, 4] * u.deg, [5, 6, 9] * u.deg)
    array([1217, 1217, 1222])

Note that :meth:`~healpix.HEALPix.healpix_to_lonlat` and
:meth:`~healpix.HEALPix.lonlat_to_healpix` can also take/give the fractional
position inside each HEALPix pixel, e.g.::

  >>> index, dx, dy = hp.lonlat_to_healpix([1, 3, 4] * u.deg, [5, 6, 9] * u.deg,
  ...                                      return_offsets=True)
  return_offsets=True)
  >>> index
  array([1217, 1217, 1222])
  >>> dx
  array([ 0.22364669,  0.78767489,  0.58832469])
  >>> dy
  array([ 0.86809114,  0.72100823,  0.16610247])

The :meth:`~healpix.HEALPix.interpolate_bilinear_lonlat` can be used to
interpolate a HEALPix map at given coordinates. A HEALPix map is typically given
as a 1-d array with as many values as pixels in the HEALPix map, and either in
nested or ring ordering. Assuming that we have an array of values in the correct
order, we can carry out bilinear interpolation at custom positions using::

    >>> import numpy as np
    >>> values = np.arange(3072)
    >>> hp.interpolate_bilinear_lonlat([1, 2, 3] * u.deg, [5, 8, 10] * u.deg, values)
    array([ 1217.45982896,  1220.20594161,  1222.41978026])

Finally, the :meth:`~healpix.HEALPix.cone_search_lonlat` method can be used to
find all HEALpix pixels within a certain radius from a longitude/latitude::

    >>> hp.cone_search_lonlat(10 * u.deg, 30 * u.deg, radius=10 * u.deg)
    array([1269,  160,  162, 1271, 1270, 1268, 1246, 1247,  138,  139,  161,
           1245,  136,  137,  140,  142,  130,  131, 1239, 1244, 1238, 1241,
           1243, 1265, 1267, 1276, 1273, 1277,  168,  169,  163,  166,  164])

Celestial HEALPix pixellization
-------------------------------

For cases where the HEALPix pixellization is of the celestial sphere, a
specialized class :class:`~healpix.CelestialHEALPix` is provided. This is a
sub-class of :class:`~healpix.HEALPix`, and in addition to the above
functionality, it is possible to convert HEALPix indices to celestial
coordinates (represented by :class:`~astropy.coordinates.SkyCoord`) and
vice-versa.

Initializing the :class:`~healpix.CelestialHEALPix` class is done as for
:class:`~healpix.HEALPix` but with an additional ``frame`` keyword argument
which specifies the frame in which the HEALPix pixellization is defined::

    >>> from healpix import CelestialHEALPix
    >>> from astropy.coordinates import Galactic
    >>> hp = CelestialHEALPix(nside=16, order='nested', frame=Galactic())

This can then be used to convert from HEALPix indices to celestial coordinates
using the :meth:`~healpix.CelestialHEALPix.healpix_to_skycoord` method::

    >>> hp.healpix_to_skycoord([144, 231])
    <SkyCoord (Galactic): (l, b) in deg
        [( 33.75      ,  32.7971683 ), ( 32.14285714,  69.42254649)]>

and from celestial coordinates to HEALPix indices using the
:meth:`~healpix.CelestialHEALPix.skycoord_to_healpix` method, e.g::

    >>> from astropy.coordinates import SkyCoord
    >>> coord = SkyCoord.from_name('m31')
    >>> hp.skycoord_to_healpix(coord)
    array([2537])

Finally, the :meth:`~healpix.CelestialHEALPix.interpolate_bilinear_skycoord` method can
be used for interpolation::

    >>> values = np.arange(3072)
    >>> hp.interpolate_bilinear_skycoord(coord, values)
    array([ 167.03780645])

and the :meth:`~healpix.CelestialHEALPix.cone_search_skycoord` method can be used for
cone searches::

    >>> hp.cone_search_skycoord(coord, values)
    array([ 167.03780645])

See the `High-level interface`_ section for more details on the interpolation
and the cone search.

Converting between ring and nested conventions
----------------------------------------------

The :class:`~healpix.HEALPix` class (and by extension the
:class:`~healpix.CelestialHEALPix` class) have methods that can be used to
convert HEALPix pixel indices between the ring and nested convention. These
are :meth:`~healpix.HEALPix.nested_to_ring`::

    >>> hp.nested_to_ring([30])
    array([873])

and :meth:`~healpix.HEALPix.ring_to_nested`::

    >>> hp.ring_to_nested([1, 2, 3])
    array([ 511,  767, 1023])

Low-level interface(s)
----------------------

If you would prefer to use a functional interface, you can use the functions
from `healpix.core`. These functions include sanity checking of the input, so if
performance is paramount and you want to access the Cython functions directly,
you can do so via the `healpix.core_cython` sub-package. Be sure to read the
documentation for the functions you want to use, since the Cython functions
require the data to be in specific numerical types in order to work properly.

Healpy-compatible interface
---------------------------

In addition to the above high- and low-level interfaces, we have provided
a `healpy <http://healpy.readthedocs.io>`_-compatible interface in
`healpix.healpy`. Note that this only includes a subset of the healpy functions.
This is not the recommended interface, and is only provided as a convenience
for packages that want to support both heapy and this package.
