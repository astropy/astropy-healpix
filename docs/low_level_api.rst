Low-level functional interface(s)
=================================

The main interface available in **healpix** which has been described so far
is the object-oriented interface using the :class:`~healpix.HEALPix` or
:class:`~healpix.CelestialHEALPix` classes.

If you would prefer to use a functional interface, you can use the functions
from `healpix.core` (you can click on the module name to see the list of
available functions). These functions include checking of the input to make sure
it is sensible, so if performance is paramount and you want to access the Cython
functions directly, you can do so via the `healpix.core_cython` sub-package. Be
sure to read the documentation for the functions you want to use, since the
Cython functions require the data to be in specific numerical types in order to
work properly.
