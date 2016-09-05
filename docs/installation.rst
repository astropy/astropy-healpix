.. include:: references.txt

.. _install:

************
Installation
************

* Python 2.7 and 3.4+ are supported.
* The only required dependency for ``healpix`` is Astropy (version 1.2 or later).

The ``healpix`` package works like most other Astropy affiliated packages.
Since it is planned to be merged into the Astropy core, we didn't put much
effort into writing up installation instructions for this separate package.

Stable version
==============

Installing the latest stable version is possible either using pip or conda.

Using pip
---------

To install regions with `pip <http://www.pip-installer.org/en/latest/>`_
from `PyPI <https://pypi.python.org/pypi/regions>`_
simply run::

    pip install --no-deps healpix

.. note::

    The ``--no-deps`` flag is optional, but highly recommended if you already
    have Numpy installed, since otherwise pip will sometimes try to "help" you
    by upgrading your Numpy installation, which may not always be desired.

Using conda
-----------

To install regions with `Anaconda <https://www.continuum.io/downloads>`_
from the `astropy channel on anaconda.org <https://anaconda.org/astropy/regions>`__
simply run::

    conda install -c astropy healpix


Testing installation
--------------------

To check if your install is OK, run the tests:

.. code-block:: bash

    python -c 'import healpix; healpix.test()'

Development version
===================

Install the latest development version from https://github.com/cdeil/healpix :

.. code-block:: bash

    git clone https://github.com/cdeil/healpix
    cd healpix
    python setup.py install
    python setup.py test
    python setup.py build_docs

Optional dependencies
=====================

The following packages are optional dependencies, install if needed:

* pytest and healpy for testing
* ...