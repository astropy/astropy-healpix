.. include:: references.txt

.. warning::
    This ``astropy_healpix`` package is in a very early stage of development.
    It is not feature complete or API stable!
    That said, please have a look and try to use it for your applications.
    Feedback and contributions welcome!

What is HEALPix?
================

**astropy_healpix** is a Python package that implements the `HEALPix
<https://en.wikipedia.org/wiki/HEALPix>`_ (Hierarchical Equal Area isoLatitude
Pixelisation) algorithm for pixellizing a sphere. This algorithm is sometimes
used in Astronomy to store data from all-sky surveys, but the general algorithm
can apply to any field that has to deal with representing data on a sphere.

More information about HEALPix can be found here:

* http://healpix.jpl.nasa.gov/
* http://adsabs.harvard.edu/abs/2005ApJ...622..759G
* http://adsabs.harvard.edu/abs/2007MNRAS.381..865C

About this package
==================

This package is a new BSD-licensed implementation that is separate from the
original GPL-licensed `HEALPix library <http://healpix.sourceforge.net>`_ and
associated `healpy <https://github.com/healpy/healpy>`__ Python wrapper. See
:ref:`about` for further information about the difference between this new
implementation and the original libraries.

The code can be found on `GitHub <https://dgithub.com/astropy/astropy_healpix>`__, along
with the list of `Contributors <https://github.com/astropy/astropy_healpix/graphs/contributors>`__.

User documentation
==================

.. toctree::
   :maxdepth: 1

   installation
   getting_started
   coordinates
   boundaries
   cone_search
   interpolation
   low_level_api
   healpy_compat
   api
