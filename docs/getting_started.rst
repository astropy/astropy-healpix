.. include:: references.txt

.. _using:

Getting started
===============

.. _using-intro:

The cleanest way to use the functionality in **healpix** is to make use of the
high-level :class:`~astropy_healpix.HEALPix` class. The
:class:`~astropy_healpix.HEALPix` class should be initialized with the ``nside``
parameter which controls the resolution of the pixellization - it is the number
of pixels on the side of each of the 12 top-level HEALPix pixels::

    >>> from astropy_healpix import HEALPix
    >>> hp = HEALPix(nside=16)

As described in the references above, HEALPix pixel indices can follow two
different ordering conventions - the *nested* convention and the *ring*
convention. By default, the :class:`~astropy_healpix.HEALPix` class assumes the ring
ordering convention, but it is possible to explicitly specify the convention to
use using the ``order`` argument, for example::

    >>> hp = HEALPix(nside=16, order='ring')

or::

    >>> hp = HEALPix(nside=16, order='nested')

Once this class has been set up, you can access various properties and methods
related to the HEALPix pixellization. For example, you can calculate the
number of pixels as well as the pixel area or resolution::

    >>> hp.npix
    3072
    >>> hp.pixel_area  # doctest: +FLOAT_CMP
    <Quantity 0.0040906154343617095 sr>
    >>> hp.pixel_resolution  # doctest: +FLOAT_CMP
    <Quantity 219.87113035631398 arcmin>

As you can see, when appropriate the properties and the methods on the
:class:`~astropy_healpix.HEALPix` class return Astropy high-level classes such as
:class:`~astropy.units.Quantity`, :class:`~astropy.coordinates.Longitude`, and
so on.

For example, the :meth:`~astropy_healpix.HEALPix.healpix_to_lonlat` method can be used
to convert HEALPix indices to :class:`~astropy.coordinates.Longitude` and
:class:`~astropy.coordinates.Latitude` objects::

    >>> lon, lat = hp.healpix_to_lonlat([1, 442, 2200])
    >>> lon  # doctest: +FLOAT_CMP
    <Longitude [ 0.83448555, 1.63624617, 0.4712389 ] rad>
    >>> lat  # doctest: +FLOAT_CMP
    <Latitude [ 0.08343009, 0.94842784,-0.78529135] rad>

The :class:`~astropy_healpix.HEALPix` class includes methods that take or return
:class:`~astropy.coordinates.SkyCoord` objects (we will take a look at this in
the :ref:`celestial` section).

In the subsequent sections of the documentation, we will take a closer look at
converting between coordinate systems, as well as more advanced features such as
interpolation and cone searches.
