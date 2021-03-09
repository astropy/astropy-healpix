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


def get_import(package, function):
    if package == 'astropy_healpix':
        return f'from astropy_healpix.healpy import {function}'
    else:
        return f'from healpy import {function}'


def bench_pix2ang(size=None, nside=None, nest=None, package=None, fast=False):
    setup = '\n'.join([
        get_import(package, 'pix2ang'),
        'import numpy as np',
        f'nside={nside}',
        f'ipix=(np.random.random({size}) * 12 * nside ** 2).astype(np.int64)',
        f'nest={nest}'])

    stmt = 'pix2ang(nside, ipix, nest)'

    return autotimeit(stmt=stmt, setup=setup, repeat=1, mintime=0 if fast else 0.1)


def bench_ang2pix(size=None, nside=None, nest=None, package=None, fast=False):
    setup = '\n'.join([
        get_import(package, 'ang2pix'),
        'import numpy as np',
        f'nside={nside}',
        f'lon=360 * np.random.random({size})',
        f'lat=180 * np.random.random({size}) - 90',
        f'nest={nest}'])

    stmt = 'ang2pix(nside, lon, lat, nest, lonlat=True)'

    return autotimeit(stmt=stmt, setup=setup, repeat=1, mintime=0 if fast else 0.1)


def bench_nest2ring(size=None, nside=None, package=None, fast=False):
    setup = '\n'.join([
        get_import(package, 'nest2ring'),
        'import numpy as np',
        f'nside={nside}',
        f'ipix=(np.random.random({size}) * 12 * nside ** 2).astype(np.int64)'])

    stmt = 'nest2ring(nside, ipix)'

    return autotimeit(stmt=stmt, setup=setup, repeat=1, mintime=0 if fast else 0.1)


def bench_ring2nest(size=None, nside=None, package=None, fast=False):
    setup = '\n'.join([
        get_import(package, 'ring2nest'),
        'import numpy as np',
        f'nside={nside}',
        f'ipix=(np.random.random({size}) * 12 * nside ** 2).astype(np.int64)'])

    stmt = 'ring2nest(nside, ipix)'

    return autotimeit(stmt=stmt, setup=setup, repeat=1, mintime=0 if fast else 0.1)


def bench_get_interp_weights(size=None, nside=None, nest=None, package=None, fast=False):
    setup = '\n'.join([
        get_import(package, 'get_interp_weights'),
        'import numpy as np',
        f'nside={nside}',
        f'lon=360 * np.random.random({size})',
        f'lat=180 * np.random.random({size}) - 90',
        f'nest={nest}'])

    stmt = 'get_interp_weights(nside, lon, lat, nest=nest, lonlat=True)'

    return autotimeit(stmt=stmt, setup=setup, repeat=1, mintime=0 if fast else 0.1)


def run_single(name, benchmark, fast=False, **kwargs):

    time_self = benchmark(package='astropy_healpix', fast=fast, **kwargs)
    results_single = dict(function=name, time_self=time_self, **kwargs)

    if HEALPY_INSTALLED:
        time_healpy = bench_ang2pix(package='healpy', fast=fast, **kwargs)
        results_single['time_healpy'] = time_healpy

    return results_single


def bench_run(fast=False):
    """Run all benchmarks. Return results as a dict."""
    results = []

    if fast:
        SIZES = [10, 1_000, 100_000]
    else:
        SIZES = [10, 1_000, 1_000_000]

    for nest in [True, False]:
        for size in SIZES:
            for nside in [1, 128]:
                results.append(run_single('pix2ang', bench_pix2ang, fast=fast,
                                          size=size, nside=nside, nest=nest))

    for nest in [True, False]:
        for size in SIZES:
            for nside in [1, 128]:
                results.append(run_single('ang2pix', bench_ang2pix, fast=fast,
                                          size=size, nside=nside, nest=nest))

    for size in SIZES:
        for nside in [1, 128]:
            results.append(run_single('nest2ring', bench_nest2ring, fast=fast,
                                      size=size, nside=nside))

    for size in SIZES:
        for nside in [1, 128]:
            results.append(run_single('ring2nest', bench_ring2nest, fast=fast,
                                      size=size, nside=nside))

    for nest in [True, False]:
        for size in SIZES:
            for nside in [1, 128]:
                results.append(run_single('get_interp_weights', bench_get_interp_weights,
                                          fast=fast, size=size,
                                          nside=nside, nest=nest))

    return results


def bench_report(results):
    """Print a report for given benchmark results to the console."""

    table = Table(names=['function', 'nest', 'nside', 'size',
                         'time_healpy', 'time_self', 'ratio'],
                  dtype=['S20', bool, int, int, float, float, float], masked=True)
    for row in results:
        table.add_row(row)

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
