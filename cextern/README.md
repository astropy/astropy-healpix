# astropy-healpix/cextern/

The `astropy-healpix` Python package is a wrapper around a C library.
See http://astropy-healpix.readthedocs.io/en/latest/about.html

This README gives some technical details on the C code here.

- The main file is `healpix.h` and `healpix.c`, start reading there first.
- For the Python `astropy-healpix` packge, the C code is built via `setup.py`
- However, to help work on the C code and test it directly, a `Makefile`
  is included here.
- For testing, a copy of `CuTest.h` and `CuTest.c` from here is bundled:
  https://github.com/asimjalis/cutest
