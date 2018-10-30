/* -*- c -*- */
/*
 * vim:syntax=c
 *
 * Low-level routines related to IEEE-754 format
 *
 * Adapted from https://github.com/numpy/numpy/blob/master/numpy/core/src/npymath/ieee754.c.src,
 * removing all functions except for npy_set_floatstatus_invalid, which is
 * made static and has is renamed to _npy_set_floatstatus_invalid.
 *
 *
 * Copyright (c) 2005-2017, NumPy Developers.
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *        copyright notice, this list of conditions and the following
 *        disclaimer in the documentation and/or other materials provided
 *        with the distribution.
 * 
 *     * Neither the name of the NumPy Developers nor the names of any
 *        contributors may be used to endorse or promote products derived
 *        from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include "numpy/utils.h"


#if (defined(__unix__) || defined(unix)) && !defined(USG)
#include <sys/param.h>
#endif


/* Solaris --------------------------------------------------------*/
/* --------ignoring SunOS ieee_flags approach, someone else can
**         deal with that! */
#if defined(sun) || defined(__BSD__) || defined(__OpenBSD__) || \
    (defined(__FreeBSD__) && (__FreeBSD_version < 502114)) || \
    defined(__NetBSD__)
#include <ieeefp.h>

static void _npy_set_floatstatus_invalid(void)
{
    fpsetsticky(FP_X_INV);
}

#elif defined(_AIX)
#include <float.h>
#include <fpxcp.h>

static void _npy_set_floatstatus_invalid(void)
{
    fp_raise_xcp(FP_INVALID);
}

#elif defined(_MSC_VER) || (defined(__osf__) && defined(__alpha))

/*
 * By using a volatile floating point value,
 * the compiler is forced to actually do the requested
 * operations because of potential concurrency.
 *
 * We shouldn't write multiple values to a single
 * global here, because that would cause
 * a race condition.
 */
static volatile double _npy_floatstatus_x, _npy_floatstatus_inf;

static void _npy_set_floatstatus_invalid(void)
{
    _npy_floatstatus_inf = NPY_INFINITY;
    _npy_floatstatus_x = _npy_floatstatus_inf - NPY_INFINITY;
}

#else
/* General GCC code, should work on most platforms */
#  include <fenv.h>

static void _npy_set_floatstatus_invalid(void)
{
    feraiseexcept(FE_INVALID);
}

#endif
