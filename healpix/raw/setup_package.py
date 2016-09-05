import os
from distutils.core import Extension

CWD = os.path.relpath(os.path.dirname(__file__))


def get_extensions():
    libraries = []

    sources = []
    sources.append(os.path.join(CWD, 'wrap.c'))
    sources.append(os.path.join(CWD, 'example.c'))

    include_dirs = ['numpy']
    include_dirs.append(CWD)

    extension = Extension(
        name='healpix.raw.wrap',
        sources=sources,
        include_dirs=include_dirs,
        libraries=libraries,
        language='c',
        # extra_compile_args=['-O2'],
    )

    return [extension]


def get_package_data():
    header_files = ['example.h']

    return {'healpix.raw': header_files}
