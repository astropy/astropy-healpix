Performance
===========

At this time, we have focused mostly on implementing functionality into the
**astropy-healpix** package, performance is not as good in most
cases as the `healpy <https://github.com/healpy/healpy>`__ library. Once
the API is stable, we will focus on improving performance.

Benchmark
---------

To get an idea of
how the performance of the two packages compare, we have included some simple
benchmarks that compare the healpy-compatible interface of **astropy-healpix**
with healpy itself. These benchmarks are run with:

.. code-block:: bash

    $ python -m astropy_healpix.bench
    Running benchmarks...

      fct    nest nside   size  time_healpy time_self   ratio
    ------- ----- ----- ------- ----------- ---------- -------
    pix2ang  True     1      10   0.0000081  0.0003575   43.91
    pix2ang  True   128      10   0.0000082  0.0003471   42.52
    pix2ang  True     1    1000   0.0000399  0.0004751   11.92
    pix2ang  True   128    1000   0.0000345  0.0004575   13.28
    pix2ang  True     1 1000000   0.0434032  0.1589150    3.66
    pix2ang  True   128 1000000   0.0364285  0.1383810    3.80
    pix2ang False     1      10   0.0000080  0.0004040   50.30
    pix2ang False   128      10   0.0000082  0.0003322   40.63
    pix2ang False     1    1000   0.0000400  0.0005005   12.50
    pix2ang False   128    1000   0.0000548  0.0005045    9.21
    pix2ang False     1 1000000   0.0342841  0.1429310    4.17
    pix2ang False   128 1000000   0.0478645  0.1405270    2.94

For small arrays, ``pix2ang`` in **astropy-healpix** performs worse, but in both
caes the times are less than a millisecond, and such differences may therefore
not matter. For larger arrays, the difference is a factor of a few at most.
We will add more benchmarks over time to provide a more complete picture.
