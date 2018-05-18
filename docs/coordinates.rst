Coordinate conversions
======================

Converting between pixel indices and spherical coordinates
----------------------------------------------------------

As described in :doc:`getting_started`, coordinates in a HEALPix pixellization
can follow either the 'ring' or 'nested' convention. Let's start by setting up
an example pixellization::

    >>> from astropy_healpix import HEALPix
    >>> hp = HEALPix(nside=16, order='nested')

The :meth:`~astropy_healpix.HEALPix.healpix_to_lonlat` method can be used
to convert HEALPix indices to :class:`~astropy.coordinates.Longitude` and
:class:`~astropy.coordinates.Latitude` objects::

    >>> lon, lat = hp.healpix_to_lonlat([1, 442, 2200])
    >>> lon  # doctest: +FLOAT_CMP
    <Longitude [ 0.83448555, 1.63624617, 0.4712389 ] rad>
    >>> lat  # doctest: +FLOAT_CMP
    <Latitude [ 0.08343009, 0.94842784,-0.78529135] rad>

The :class:`~astropy.coordinates.Longitude` and
:class:`~astropy.coordinates.Latitude` objects are fully-fledged
:class:`~astropy.units.Quantity` objects and also include shortcuts to get
the values in various units::

    >>> lon.hourangle  # doctest: +FLOAT_CMP
    array([ 3.1875,  6.25  ,  1.8   ])
    >>> lat.degree  # doctest: +FLOAT_CMP
    array([  4.78019185,  54.3409123 , -44.99388015])

Conversely, given longitudes and latitudes as :class:`~astropy.units.Quantity`
objects, it is possible to recover HEALPix pixel indices::

    >>> from astropy import units as u
    >>> print(hp.lonlat_to_healpix([1, 3, 4] * u.deg, [5, 6, 9] * u.deg))
    [1217 1217 1222]

In these examples, what is being converted is the position of the center of each
pixel. In fact, the  :meth:`~astropy_healpix.HEALPix.lonlat_to_healpix` method can also
take or give the fractional position inside each HEALPix pixel, e.g.::

    >>> index, dx, dy = hp.lonlat_to_healpix([1, 3, 4] * u.deg, [5, 6, 9] * u.deg,
    ...                                      return_offsets=True)
    >>> print(index)
    [1217 1217 1222]
    >>> dx  # doctest: +FLOAT_CMP
    array([ 0.22364669,  0.78767489,  0.58832469])
    >>> dy  # doctest: +FLOAT_CMP
    array([ 0.86809114,  0.72100823,  0.16610247])

and the :meth:`~astropy_healpix.HEALPix.healpix_to_lonlat` method can take offset
positions - for example we can use this to find the position of the corners of
a given pixel::

    >>> dx = [0., 1., 1., 0.]
    >>> dy = [0., 0., 1., 1.]
    >>> lon, lat = hp.healpix_to_lonlat([133, 133, 133, 133], dx=dx, dy=dy)
    >>> lon  # doctest: +FLOAT_CMP
    <Longitude [ 0.53996124, 0.58904862, 0.53996124, 0.49087385] rad>
    >>> lat  # doctest: +FLOAT_CMP
    <Latitude [ 0.47611906, 0.52359878, 0.57241857, 0.52359878] rad>

.. _celestial:

Celestial coordinates
---------------------

For cases where the HEALPix pixellization is of the celestial sphere, a
``frame`` argument can be passed to :class:`~astropy_healpix.HEALPix`. This
argument should specify the celestial frame (using an `astropy.coordinates
<http://docs.astropy.org/en/stable/coordinates/index.html>`_ frame) in which the
HEALPix pixellization is defined::

    >>> from astropy_healpix import HEALPix
    >>> from astropy.coordinates import Galactic
    >>> hp = HEALPix(nside=16, order='nested', frame=Galactic())

Each method defined in :class:`~astropy_healpix.HEALPix` and ending in
``lonlat`` has an equivalent method ending in ``skycoord`` which can be used if
the frame is set. For example, to convert from HEALPix indices to celestial
coordinates, you can use the
:meth:`~astropy_healpix.HEALPix.healpix_to_skycoord` method::

    >>> hp.healpix_to_skycoord([144, 231])  # doctest: +FLOAT_CMP
    <SkyCoord (Galactic): (l, b) in deg
        [( 33.75      ,  32.7971683 ), ( 32.14285714,  69.42254649)]>

and to convert from celestial coordinates to HEALPix indices you can use the
:meth:`~astropy_healpix.HEALPix.skycoord_to_healpix` method, e.g::

    >>> from astropy.coordinates import SkyCoord
    >>> coord = SkyCoord('00h42m44.3503s +41d16m08.634s')
    >>> hp.skycoord_to_healpix(coord)
    2537

Converting between ring and nested conventions
----------------------------------------------

The :class:`~astropy_healpix.HEALPix` class has methods that can be used to
convert HEALPix pixel indices between the ring and nested convention. These are
:meth:`~astropy_healpix.HEALPix.nested_to_ring`::

    >>> print(hp.nested_to_ring([30]))
    [873]

and :meth:`~astropy_healpix.HEALPix.ring_to_nested`::

    >>> print(hp.ring_to_nested([1, 2, 3]))
    [ 511  767 1023]
