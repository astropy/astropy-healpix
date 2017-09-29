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

The `healpy <https://github.com/healpy/healpy>`__ package that is a wrapper
around the `HEALPix <http://healpix.jpl.nasa.gov/>`__ C++ library has existed
for a long time. So why this re-write?

The main motivation is that the original HEALPIX/healpy packages are
GPL-licensed, which is incompatible with the BSD license used by Astropy and
most Astropy-affiliated and scientific Python (Numpy, Scipy, ...) package, and
HEALPix/healpy will not be relicensed (see `here
<https://sourceforge.net/p/healpix/mailman/message/34929929/>`__).
In addition, the present package doesn't have a big C++ package as a dependency,
just a little C and Cython code, which makes it easy to install everywhere -- we
support Linux, MacOS X, and Windows.

However, this package is intended to be lightweight and we do not plan to fully
implement everything that healpy and the original HEALPix library support (such
as spherical harmonics).

Note that code contributions to this package can't be derived from the HEALPix
package or healpy due to licensing reasons (see above).
