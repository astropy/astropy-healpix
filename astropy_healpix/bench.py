# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, print_function, division

"""Benchmarks for this package.

To run all benchmarks and print a report to the console::

    python -m healpix.bench

You can also run the benchmarks first, save the results dict
to disk as a JSON file (or share it with others) and then
print the results later, or compare them with other results.

We should now that this is not very comprehensive / flexible.

If your application depends on performance of HEALPix computations,
you should write benchmarks with cases relevant for that application
and check if HEALPix computations are really the bottleneck and if
this package is fast enough for you or not.
"""

import timeit

from astropy.table import Table

# NOTE: If healpy is installed, we use it in the benchmarks, but healpy is not
# a formal dependency of astropy-healpix.
try:
    import healpy as hp  # noqa
except ImportError:
    HEALPY_INSTALLED = False
else:
    HEALPY_INSTALLED = True


# Copied from https://github.com/kwgoodman/bottleneck/blob/master/bottleneck/benchmark/autotimeit.py
def autotimeit(stmt, setup='pass', repeat=3, mintime=0.2):
    timer = timeit.Timer(stmt, setup)
    number, time1 = autoscaler(timer, mintime)
    time2 = timer.repeat(repeat=repeat - 1, number=number)
    return min(time2 + [time1]) / number


# Copied from https://github.com/kwgoodman/bottleneck/blob/master/bottleneck/benchmark/autotimeit.py
def autoscaler(timer, mintime):
    number = 1
    for i in range(12):
        time = timer.timeit(number)
        if time > mintime:
            return number, time
        number *= 10
    raise RuntimeError('function is too fast to test')


def get_import(package, fct):
    if package == 'astropy_healpix':
        return 'from astropy_healpix.healpy import {}'.format(fct)
    else:
        return 'from healpy import {}'.format(fct)


def bench_pix2ang(size, nside, nest, package, fast=False):
    shape = (int(size), )

    setup = '\n'.join([
        get_import(package, 'pix2ang'),
        'import numpy as np',
        'nside={}'.format(int(nside)),
        'ipix=np.zeros({}, dtype=np.int64)'.format(shape),
        'nest={}'.format(nest)])

    stmt = 'pix2ang(nside, ipix, nest)'

    return autotimeit(stmt=stmt, setup=setup, repeat=1, mintime=0 if fast else 0.1)


def bench_run(fast=False):
    """Run all benchmarks. Return results as a dict."""
    results = []

    for nest in [True, False]:
        for size in [10, 1e3, 1e6]:
            for nside in [1, 128]:

                time_self = bench_pix2ang(size=size, nside=nside,
                                          nest=nest, package='astropy_healpix',
                                          fast=fast)

                results_single = dict(fct='pix2ang', size=int(size),
                                      nside=nside, time_self=time_self)

                if HEALPY_INSTALLED:
                    time_healpy = bench_pix2ang(size=size, nside=nside,
                                                nest=nest, package='healpy',
                                                fast=fast)
                    results_single['time_healpy'] = time_healpy

                results.append(results_single)

    return results


def bench_report(results):
    """Print a report for given benchmark results to the console."""
    table = Table(rows=results)

    table['time_self'].format = '10.7f'

    if HEALPY_INSTALLED:
        table['ratio'] = table['time_self'] / table['time_healpy']
        table['time_healpy'].format = '10.7f'
        table['ratio'].format = '7.2f'

    table.pprint(max_lines=-1)


def main(fast=False):
    """Run all benchmarks and print report to the console."""
    print('Running benchmarks...\n')
    results = bench_run(fast=fast)
    bench_report(results)


if __name__ == '__main__':
    main()
