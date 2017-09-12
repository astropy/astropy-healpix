.. include:: references.txt

.. doctest-skip-all

.. _gs:

***************
Getting started
***************

.. warning::
    This ``healpix`` package is in a very early stage of development.
    It is not feature complete or API stable!
    That said, please have a look and try to use it for your applications.
    Feedback and contributions welcome!

.. _gs-intro:

What is HEALPix?
================

Description of [HEALPix from Wikipedia](https://en.wikipedia.org/wiki/HEALPix):

> HEALPix (sometimes written as Healpix), an acronym for Hierarchical
> Equal Area isoLatitude Pixelisation of a 2-sphere, is an algorithm
> for pixelisation of the 2-sphere, and the associated class of map
> projections.

More information:

* http://healpix.jpl.nasa.gov/
* http://adsabs.harvard.edu/abs/2005ApJ...622..759G
* http://adsabs.harvard.edu/abs/2007MNRAS.381..865C

Using healpix
=============

The ``healpix`` API isn't stable yet.
We have to discuss what we want for Astropy and how to best
expose the functionality.

A major question is whether to use a healpy-compatible API and
whether to use `astropy.coordinates` in the API.

For now we're exposing the functionality in these namespaces:

- `healpix.healpy` -- `healpy`-compatible API
- `healpix.raw` -- Raw wrappers around C functionality

The top-level `healpix` namespace is still empty and we'll have to
discuss what API we want.

healpix.healpy
--------------

The API in `healpix.healpy` should be a drop-in replacement for `healpy <http://healpy.readthedocs.io/>`__ .

I.e. if your code currently uses ``healpy`` like this:

    >>> import healpy as hp

you should be able to replace that import statement with

    >>> import healpix.healpy as hp

and as long as we've implemented that functionality, the result should be the same.

Please file an issue or make a pull request if you find that:

- Some function from `healpix.healpy` doesn't give identical results as `healpy`.
- Some function you need is missing.
- Some function you need it too slow or uses too much memory.

Note: even if we don't want to support this API in the future, we probably should just make it private
(i.e. rename to `healpix._healpy`), because it's very useful for testing.

healpix.raw
-----------

The `healpix.raw` package exposes the functionality implemented in C via a Cython wrapper.


What next?
----------

Try out `healpix` for your applications.
Please give feedback or contribute!
