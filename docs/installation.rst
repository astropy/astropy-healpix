.. include:: references.txt

.. _install:

************
Installation
************

Dependencies
============

Required dependencies
---------------------

The **healpix** package works with Python 2.7 or 3.5 and later, and requires the
following dependencies:

* `Numpy <http://www.numpy.org>`__ 1.10 or later
* `Astropy <http://www.astropy.org>`__ 1.2 or later

If you use :ref:`pip` or :ref:`conda`, these will be installed automatically.

Optional dependencies
---------------------

The following packages are optional dependencies, which can be installed if
needed:

* `pytest <http://www.pytest.org>`__ for testing
* `healpy <https://healpy.readthedocs.io>`__ for testing (but this is not required
  and the tests that require healpy will be skipped if healpy is not installed)

Stable version
==============

Installing the latest stable version is possible either using pip or conda.

.. _pip:

Using pip
---------

To install **healpix** with `pip <http://www.pip-installer.org/en/latest/>`__
from `PyPI <https://pypi.python.org/pypi/regions>`__
simply run::

    pip install --no-deps healpix

.. note::

    The ``--no-deps`` flag is optional, but highly recommended if you already
    have Numpy installed, since otherwise pip will sometimes try to "help" you
    by upgrading your Numpy installation, which may not always be desired.

.. _conda:

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
