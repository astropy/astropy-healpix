Healpy-compatible interface
===========================

In addition to the above high- and low-level interfaces, we have provided
a `healpy <http://healpy.readthedocs.io>`_-compatible interface in
:mod:`astropy_healpix.healpy`. Note that this only includes a subset of the healpy functions.
This is not the recommended interface, and is only provided as a convenience
for packages that want to support both healpy and this package.

Example
-------

As an example, the :func:`~astropy_healpix.healpy.pix2ang` function can be used to get
the longitude/latitude of a given HEALPix pixel (by default using the 'ring'
convention)::

  >>> from astropy_healpix.healpy import pix2ang
  >>> pix2ang(16, [100, 120])
  (array([ 0.35914432,  0.41113786]), array([ 3.70259134,  1.6689711 ]))

which agrees exactly with the healpy function::

.. doctest-requires:: healpy

  >>> from healpy import pix2ang
  >>> pix2ang(16, [100, 120])
  (array([ 0.35914432,  0.41113786]), array([ 3.70259134,  1.6689711 ]))

Migrate
-------

To migrate a script or package from using ``healpy`` to this ``healpix`` package,
to check if the required functionality is available by changing all::

    import healpy as hp

to::

    from astropy_healpix import healpy as hp

and see what's missing or breaks. Please file issues or feature requests!

As mentioned above, we then recommend that when you actually make the change,
you use the main API of this package instead of the ``healpy``-compatible interface.
