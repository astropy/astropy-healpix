# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, print_function, division

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


def get_extensions():

    libraries = []

    sources = [os.path.join(C_DIR, filename) for filename in C_FILES]
    sources.append(os.path.join(HEALPIX_ROOT, 'core_cython.pyx'))

    include_dirs = ['numpy', C_DIR]

    extension = Extension(
        name="astropy_healpix.core_cython",
        sources=sources,
        include_dirs=include_dirs,
        libraries=libraries,
        language="c",
        extra_compile_args=['-O2'])

    return [extension]
