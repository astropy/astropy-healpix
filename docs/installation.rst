.. include:: references.txt

.. doctest-skip-all

.. _install:

************
Installation
************

Dependencies
============

Required dependencies
---------------------

The **astropy-healpix** package works with Python 2.7 or 3.5 and later (on Linux, MacOS
X and Windows), and requires the following dependencies:

* `Numpy <http://www.numpy.org>`__ 1.10 or later
* `Astropy <http://www.astropy.org>`__ 1.2 or later

If you use :ref:`pip` or :ref:`conda`, these will be installed automatically.

Optional dependencies
---------------------

The following packages are optional dependencies, which can be installed if needed:

* `pytest <http://www.pytest.org>`__ for testing
* `healpy <https://healpy.readthedocs.io>`__ for testing (but this is not required
  and the tests that require healpy will be skipped if healpy is not installed)
* `hypothesis <hypothesis.readthedocs.io>`__ for the healpy-related tests.

Stable version
==============

Installing the latest stable version is possible either using pip or conda.

.. _pip:

Using pip
---------

To install **astropy-healpix** with `pip <http://www.pip-installer.org/en/latest/>`__
from `PyPI <https://pypi.python.org/pypi/astropy-healpix>`__
simply run::

    pip install --no-deps astropy-healpix

.. note::

    The ``--no-deps`` flag is optional, but highly recommended if you already
    have Numpy installed, since otherwise pip will sometimes try to "help" you
    by upgrading your Numpy installation, which may not always be desired.

.. _conda:

Using conda
-----------

To install healpix with `Anaconda <https://www.continuum.io/downloads>`_
from the `conda-forge channel on anaconda.org <https://anaconda.org/conda-forge/astropy-healpix>`__
simply run::

    conda install -c conda-forge astropy-healpix

Testing installation
--------------------

To check that you have this package installed and which version you're using,
start Python and execute the following code:

.. code-block:: bash

    $ python
    Python 3.6.2 |Continuum Analytics, Inc.| (default, Jul 20 2017, 13:14:59)
    [GCC 4.2.1 Compatible Apple LLVM 6.0 (clang-600.0.57)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import astropy_healpix
    >>> astropy_healpix.__version__
    0.1

To make sure that all functionality is working OK on your system, you can
run the automated tests of this package by executing the ``test`` function:

.. code-block:: bash

    python -c 'import astropy_healpix; astropy_healpix.test()'

Development version
===================

Install the latest development version from https://github.com/astropy/astropy-healpix :

.. code-block:: bash

    git clone https://github.com/astropy/astropy-healpix
    cd astropy-healpix
    pip install .

Contributing
============

This section contains some tips how to hack on **astropy-healpix**.

You can run the tests in a temp folder via::

    python setup.py test -V

Or build the C / Cython extensions in-place and run the tests from the source folder::

    python setup.py build_ext -i
    python -m pytest -v astropy_healpix

To build the docs::

    python setup.py build_docs
    open docs/_build/html/index.html

If you have any questions, just open an issue on Github and we'll help.
