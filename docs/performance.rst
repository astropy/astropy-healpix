Performance and OpenMP
======================

Measuring performance
---------------------

At this time, we have focused mostly on implementing functionality into the
**astropy-healpix** package, and performance may therefore not be as good in
some cases as the `healpy <https://github.com/healpy/healpy>`__ library. Once
the API is stable, we will focus on improving performance. To get an idea of
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
not matter. For larger arrays, the difference is a factor of a few at most
(without multi-threading). We will add more benchmarks over time to provide a
more complete picture.

Controlling the number of threads
---------------------------------

By default, OpenMP is used if available to parallelize many of the functions in
**astropy-healpix** using multi-threading (the results above do). For example, when calling
:meth:`~astropy_healpix.HEALPix.healpix_to_skycoord` with an array of HEALPix indices, the array
is split up into multiple chunks which are distributed over threads. This should
speed up the calculation by a factor depending on the number of available cores.

However, in some cases (such as in a shared computing environment) you may want
to restrict this behavior to only ever use one core. You can control the number
of threads used by setting the ``OMP_NUM_THREADS`` environment variable (set
this to 1 to disable multi-threading).

Handling compilers that don't support OpenMP
--------------------------------------------

If you are compiling the **astropy-healpix** package yourself, note that the
default C compiler on MacOS X (clang) does not support OpenMP. You will need to
explicitly set a different compiler - for example if you have GCC 6 installed
via MacPorts, you can do:

.. code-block:: bash

    $ CC=gcc-mp-6 python setup.py install

You can check the installation log for the following message which indicates
that OpenMP is being used::

    Compiling Cython extension with OpenMP support

If your compiler does not support OpenMP, you will instead see the following
message::

    Cannot compile Cython extension with OpenMP, reverting to non-parallel code
