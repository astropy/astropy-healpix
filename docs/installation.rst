.. include:: references.txt

.. _install:

************
Installation
************

.. important::

    This Python package has the name **healpix**.
    This name is used for importing from Python (``import healpix``),
    for installation from the Python package index PyPI (``pip install healpix``),
    and the git repository on Github is also called ``healpix``.

    However, the name was already taken on the Anaconda conda-forge channel,
    there **healpix** refers to the C++ HEALPix package.
    Therefore, since this package is built on Astropy and might be integrated into
    the Astropy core package at some point, we have chosen the name ``astropy-healpix``
    in the Anaconda conda-forge channel.

    We hope this makes it clear why in the commands below sometimes ``healpix``
    and sometimes ``astropy-healpix`` appears. For further information about this
    package and it's relation to other HEALPix packages, see :ref:`about`.


Dependencies
============

Required dependencies
---------------------

The **healpix** package works with Python 2.7 or 3.5 and later (on Linux, MacOS
X and Windows), and requires the following dependencies:

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
* `hypothesis <hypothesis.readthedocs.io>`__ for the healpy-related tests.

Stable version
==============

Installing the latest stable version is possible either using pip or conda.

.. _pip:

Using pip
---------

To install **healpix** with `pip <http://www.pip-installer.org/en/latest/>`__
from `PyPI <https://pypi.python.org/pypi/healpix>`__
simply run::

    pip install --no-deps healpix

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

To check that you have this package nstalled and which version you're using,
start Python and execute the following code:

.. code-block:: bash

    $ python
    >>> import healpix
    # The following line will print out the package install location
    >>> healpix
    # The following line will print the package version number
    >>> healpix.__version__


To make sure that all functionality is working OK on your system, you can
run the automated tests of this package by executing the ``test`` function:

.. code-block:: bash

    python -c 'import healpix; healpix.test()'

Development version
===================

Install the latest development version from https://github.com/cdeil/healpix :

.. code-block:: bash

    git clone https://github.com/cdeil/healpix
    cd healpix
    pip install .

Hacking on ``healpix``
======================

This section contains some tips how to hack on ``healpix``.

You can run the tests in a temp folder via::

    python setup.py test -V

Or build the C / Cython extensions in-place and run the tests from the source folder::

    python setup.py build_ext -i
    python -m pytest -v healpix

To build the docs::

    python setup.py build_docs
    open docs/_build/html/index.html

If you have any questions, just open an issue on Github and we'll help.
