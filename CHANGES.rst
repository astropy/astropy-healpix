0.2 (2017-10-15)
================

- Expand benchmarks to include ang2pix, nest2ring and ring2nest. [#62]

- Use OpenMP to parallelize the Cython wrappers. [#59]

- Renamed the ``healpix_neighbours`` function to ``neighbours`` and added
  a wrapper to the high-level class. [#61]

- Fix bilinear interpolation which was being done incorrectly, and added
  a new ``bilinear_interpolation_weights`` function to get the interpolation
  weights. [#63]

0.1 (2017-10-01)
================

- Initial release

