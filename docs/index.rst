.. include:: references.txt

What is HEALPix?
================

`HEALPix
<https://en.wikipedia.org/wiki/HEALPix>`_ (Hierarchical Equal Area isoLatitude
Pixelisation) is an algorithm for pixellizing a sphere that is sometimes
used in Astronomy to store data from all-sky surveys, but the general algorithm
can apply to any field that has to deal with representing data on a sphere.

More information about the HEALPix algorithm can be found here:

* https://healpix.jpl.nasa.gov/
* https://ui.adsabs.harvard.edu/abs/2005ApJ...622..759G
* https://ui.adsabs.harvard.edu/abs/2007MNRAS.381..865C

About this package
==================

**astropy-healpix** is a new BSD-licensed implementation that is separate from the
original GPL-licensed `HEALPix library <http://healpix.sourceforge.net>`_ and
associated `healpy <https://github.com/healpy/healpy>`__ Python wrapper. See
:ref:`about` for further information about the difference between this new
implementation and the original libraries.

The code can be found on `GitHub <https://github.com/astropy/astropy-healpix>`__, along
with the list of `Contributors <https://github.com/astropy/astropy-healpix/graphs/contributors>`__.

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
   performance
   healpy_compat
   api

Version history
===============

For a list of changes in each version, see the `CHANGES.rst
<https://github.com/astropy/astropy-healpix/blob/main/CHANGES.rst>`_ file.
