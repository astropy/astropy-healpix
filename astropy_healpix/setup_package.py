# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

from distutils.core import Extension

HEALPIX_ROOT = os.path.relpath(os.path.dirname(__file__))

C_FILES = ['bl.c',
           'healpix-utils.c',
           'healpix.c',
           'mathutil.c',
           'permutedsort.c',
           'qsort_reentrant.c',
           'starutil.c']


C_DIR = os.path.join('cextern', 'astrometry.net')
C_DIRS = ['numpy', C_DIR, HEALPIX_ROOT,
          os.path.join('cextern', 'lalsuite'),
          os.path.join('cextern', 'numpy')]


def get_extensions():

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
        extra_compile_args=['-O2'])

    return [extension]
