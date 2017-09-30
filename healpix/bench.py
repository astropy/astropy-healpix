# Licensed under a 3-clause BSD style license - see LICENSE.rst
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
from textwrap import dedent
import timeit
from astropy.table import Table
import healpy as hp
from . import healpy as hp_compat


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
    if package == 'self':
        return 'from healpix.healpy import {}'.format(fct)
    else:
        return 'from healpy import {}'.format(fct)


def bench_pix2ang(size, nside, nest, package):
    shape = (int(size), )

    setup = '\n'.join([
        get_import(package, 'pix2ang'),
        'import numpy as np',
        'nside={}'.format(int(nside)),
        'ipix=np.zeros({}, dtype=np.int64)'.format(shape),
        'nest={}'.format(nest),
    ])

    stmt = 'pix2ang(nside, ipix, nest)'

    return autotimeit(stmt=stmt, setup=setup, repeat=1, mintime=0.1)


def bench_run():
    """Run all benchmarks. Return results as a dict."""
    results = []

    for nest in [True, False]:
        for size in [10, 1e3, 1e6]:
            for nside in [1, 128]:
                time_self = bench_pix2ang(size=size, nside=nside, nest=nest, package='self')
                time_healpy = bench_pix2ang(size=size, nside=nside, nest=nest, package='healpy')

                results.append(dict(
                    fct='pix2ang', size=size, nside=nside,
                    time_self=time_self, time_healpy=time_healpy,
                ))

    return results


def bench_report(results):
    """Print a report for given benchmark results to the console."""
    table = Table(rows=results)
    table['ratio'] = table['time_self'] / table['time_healpy']
    table.pprint(max_lines=-1)


def main():
    """Run all benchmarks and print report to the console."""
    print('Running benchmarks ...\n')
    results = bench_run()
    print('Printint benchmark report ...\n')
    bench_report(results)


if __name__ == '__main__':
    main()
