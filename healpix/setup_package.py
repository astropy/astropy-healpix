import os
from distutils.core import Extension

HEALPIX_ROOT = os.path.relpath(os.path.dirname(__file__))

C_FILES = ['bl.c',
           'errors.c',
           'healpix-utils.c',
           'healpix.c',
           'ioutils.c',
           'log.c',
           'mathutil.c',
           'permutedsort.c',
           'qsort_reentrant.c',
           'starutil.c',
           'tic.c']


C_DIR = os.path.join('astrometry.net')


def get_extensions():

    libraries = []

    sources = [os.path.join(C_DIR, filename) for filename in C_FILES]
    sources.append(os.path.join(HEALPIX_ROOT, '_healpix.pyx'))

    include_dirs = ['numpy', C_DIR]

    extension = Extension(
        name="healpix._healpix",
        sources=sources,
        include_dirs=include_dirs,
        libraries=libraries,
        language="c",
        extra_compile_args=['-O2'])

    return [extension]
