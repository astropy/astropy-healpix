# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

import numpy as np

HEALPIX_ROOT = os.path.relpath(os.path.dirname(__file__))

C_FILES = ['bl.c',
           'healpix-utils.c',
           'healpix.c',
           'mathutil.c',
           'permutedsort.c',
           'qsort_reentrant.c',
           'starutil.c']


C_DIR = os.path.join('cextern', 'astrometry.net')
C_DIRS = [np.get_include(), C_DIR, HEALPIX_ROOT,
          os.path.join('cextern', 'numpy')]


def get_extensions():
    from setuptools import Extension

    libraries = []

    sources = [os.path.join(C_DIR, filename) for filename in C_FILES]
    sources.append(os.path.join(HEALPIX_ROOT, 'interpolation.c'))
    sources.append(os.path.join(HEALPIX_ROOT, '_core.c'))

    extension = Extension(
        name="astropy_healpix._core",
        sources=sources,
        include_dirs=C_DIRS,
        libraries=libraries,
        language="c",
        extra_compile_args=['-O2'],
        py_limited_api=True,
        define_macros=[('Py_LIMITED_API', 0x030A0000),
                       ('NPY_TARGET_VERSION', 'NPY_1_19_API_VERSION'),
                       ('NPY_NO_DEPRECATED_API', 'NPY_1_19_API_VERSION')])

    return [extension]
