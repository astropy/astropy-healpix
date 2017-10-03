# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This file incldues code adapted from astroscrappy, originally released under
# the following license:
#
# Copyright (c) 2015, Curtis McCully
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice, this
#   list of conditions and the following disclaimer in the documentation and/or
#   other materials provided with the distribution.
# * Neither the name of the Astropy Team nor the names of its contributors may be
#   used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from __future__ import absolute_import, print_function, division

import os
import tempfile

from distutils import log
from distutils.core import Extension
from distutils.ccompiler import new_compiler
from distutils.sysconfig import customize_compiler


HEALPIX_ROOT = os.path.relpath(os.path.dirname(__file__))

C_FILES = ['bl.c',
           'healpix-utils.c',
           'healpix.c',
           'mathutil.c',
           'permutedsort.c',
           'qsort_reentrant.c',
           'starutil.c']


C_DIR = os.path.join('cextern', 'astrometry.net')

CCODE = """
#include <omp.h>
#include <stdio.h>
int main() {
#pragma omp parallel
printf("nthreads=%d\\n", omp_get_num_threads());
}
"""


def get_openmp_flags():

    ccompiler = new_compiler()
    customize_compiler(ccompiler)

    filename = tempfile.mktemp(suffix='.c')

    with open(filename, 'w') as f:
        f.write(CCODE)

    for flag in ['-fopenmp', '-openmp']:

        try:
            ccompiler.compile([filename], extra_postargs=[flag])
            return [flag]
        except:
            pass

    return None


def get_extensions():

    libraries = []

    sources = [os.path.join(C_DIR, filename) for filename in C_FILES]
    sources.append(os.path.join(HEALPIX_ROOT, 'core_cython.pyx'))

    include_dirs = ['numpy', C_DIR]

    compile_args = ['-O2']
    link_args = []

    openmp_flags = get_openmp_flags()

    if openmp_flags is None:
        log.warn("Cannot compile Cython extension with OpenMP, reverting to non-parallel code")
    else:
        log.info("Compiling Cython extension with OpenMP support")
        compile_args.extend(openmp_flags)
        link_args.extend(openmp_flags)

    extension = Extension(
        name="astropy_healpix.core_cython",
        sources=sources,
        include_dirs=include_dirs,
        libraries=libraries,
        language="c",
        extra_compile_args=compile_args,
        extra_link_args=link_args)

    return [extension]
