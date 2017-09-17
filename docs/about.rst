:orphan:

.. include:: references.txt

.. _about:

******************
About this package
******************

This is a BSD-licensed Python package for HEALPix, which is based on the C
HEALPix code written by Dustin Lang originally in `astrometry.net
<http://astrometry.net>`_, and was added here with a Cython wrapper and expanded
with a Python interface.

Why?
----

The `healpy <https://github.com/healpy/healpy>`__ package that is a wrapper
around the `HEALPIX <http://healpix.jpl.nasa.gov/>`__ C++ library has existed
for a long time.

So why this re-write?

The main motivation is that the original HEALPIX/healpy packages are
GPL-licensed, which is incompatible with the BSD license used by Astropy and
most Astropy-affiliated and scientific Python (Numpy, Scipy, ...) package, and
HEALPIX/healpy will not be relicensed (see `here
<https://sourceforge.net/p/healpix/mailman/message/34929929/>`__).

There are a few other reasons why the present package is useful:

* It doesn't have a big C++ package as a dependency, just a little C and
  Cython code like many other affiliated packages, i.e. it's easy to install
  everywhere (Linux, OS X, Windows).
* The rewrite offered the possibility to re-consider the implementation
  and API, e.g. we could use `astropy.coordinates.SkyCoord` in the API
  or support ``nside`` values that are not a power of 2 for the ring ordering
  scheme (see `here <https://github.com/healpy/healpy/issues/333>`__).

Plan
----

We plan to propose this package to be included in Astropy core,
so that it's available for example for MOC regions
(see https://github.com/astropy/regions/issues/62)
or Gammapy or other packages.

If this is accepted, this package will only be temporarily maintained as a
separate Python package, i.e. it's not recommended for production use
at this point.

There are no plans to implement all of the HEALPIX / healpy functionality.

- FITS I/O related functions -- Planned
- Other functions (e.g. spherical harmonics) -- Not planned.

Contributions welcome!

Note that code contributions can't be derived from the HEALPIX package
or healpy due to licensing reasons (see above).
