.. _changes:

*******
Changes
*******

0.4 (2018-12-18)
================

- Healpix rangesearch cleanup [#113]
- Update astropy-helpers to v2.0.8 [#112]
- Rewrite core module in C to make ``healpix_to_lonlat`` and
  ``lonlat_to_healpix`` broadcastable over both pixel index and nside. [#110]

0.3.1 (2018-10-24)
==================

- Ensure .c files are included in tar file.

0.3 (2018-10-24)
================

- Remove OpenMP from astropy-healpix [#108]
- Fix bilinear interpolation of invalid values [#106]
- Add uniq to (level, ipix) and inverse function [#105]
- compute z more stably; improve on z2dec [#101]
- use more stable cos(Dec) term [#94]
- Fix get_interp_weights for phi=None case [#89]
- Add pix2vec, vec2pix, ang2vec [#73]
- Add ``pixel_resolution_to_nside`` function. [#31]

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

