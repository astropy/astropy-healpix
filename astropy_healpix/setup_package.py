# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

from distutils.core import Extension
from astropy_helpers.openmp_helpers import add_openmp_flags_if_available

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
    sources.append(os.path.join(HEALPIX_ROOT, 'interpolation.c'))
    sources.append(os.path.join(HEALPIX_ROOT, 'core_cython.pyx'))

    include_dirs = ['numpy', C_DIR, HEALPIX_ROOT]

    extension = Extension(
        name="astropy_healpix.core_cython",
        sources=sources,
        include_dirs=include_dirs,
        libraries=libraries,
        language="c",
        extra_compile_args=['-O2'])

    add_openmp_flags_if_available(extension)

    return [extension]
