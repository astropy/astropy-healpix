Healpy-compatible interface
===========================

In addition to the above high- and low-level interfaces, we have provided
a `healpy <http://healpy.readthedocs.io>`_-compatible interface in
`healpix.healpy`. Note that this only includes a subset of the healpy functions.
This is not the recommended interface, and is only provided as a convenience
for packages that want to support both healpy and this package.

As an example, the :func:`~healpix.healpy.pix2ang` function can be used to get
the longitude/latitude of a given HEALPix pixel (by default using the 'ring'
convention)::

  >>> from healpix.healpy import pix2ang
  >>> pix2ang(16, [100, 120])
  (array([ 0.35914432,  0.41113786]), array([ 3.70259134,  1.6689711 ]))

which agrees exactly with the healpy function::

  >>> from healpy import pix2ang
  >>> pix2ang(16, [100, 120])
  (array([ 0.35914432,  0.41113786]), array([ 3.70259134,  1.6689711 ]))
