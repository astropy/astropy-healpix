.. _changes:

*******
Changes
*******

1.0.3 (2024-04-05)
==================

- Fixed compatibility with upcoming astropy 6.1.0 and numpy 2.0.0 releases. [#214, #215]

1.0.2 (2023-12-12)
==================

- lonlat_to_healpix now correctly returns -1 if the longitude or latitude is
  NaN or infinite. [#208]

1.0.1 (2023-11-28)
==================

- Allow building using any version of Numpy between 1.25 and 2. [#201]

- Build wheels for PyPI. [#200]

1.0.0 (2023-08-21)
==================

- Drop support for Python 3.7 and 3.8, which are not supported by the latest
  minor release of Numpy (1.25).

- Build binary wheels using the Python limited API.

- Remove warning about API stability. The API is now considered stable.

0.7 (2022-09-15)
================

- Added new methods ``healpix_to_xyz`` and ``xyz_to_healpix`` to
  the high level interface. [#153]
- The ``frame`` keyword argument for the high-level ``HEALPix`` class may now
  be a frame name, frame instance, or frame class. [#156]
- On instantiation, the ``HEALPix`` class checks the ``order`` argument. [#162]
- Drop support for Python 3.6, which has passed end-of-life. [#166]

0.6 (2021-03-10)
================

- Update package infrastructure to follow APE17 guidelines. [#142]
- Added new functions ``healpix_to_xyz`` and ``xyz_to_healpix`` to
  convert to/from cartesian coordinates. [#141]
- Add ``HEALPix.level`` property to the high-level interface.
  This is a shortcut for the ``nside_to_level`` function. [#147]

0.5 (2019-11-25)
================

- Update package infrastructure to use ``setup.cfg``. [#134]
- Make sure that Numpy is declared as a build-time dependency. [#134]
- Update astropy-helpers to v3.2.2. [#134]
- Update minimum required Python version to 3.6. [#125]
- Add ``HEALPix.from_header``. [#127]
- Clean up C code to avoid compilation warnings. [#118, #119, #120, #121, #122, #123]
- Fix unit tests on 32-bit architectures. [#117]
- Fix compatibility with Numpy 1.16 and later. [#116]

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
