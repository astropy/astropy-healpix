# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
BSD-licensed HEALPix for Astropy.

We plan to propose this package to be included in Astropy core.
See the about page for in the docs for further information.

If this is accepted, this package will only be temporarily maintained as a
separate Python package, i.e. it's not recommended for production use
at this point.

* Code : https://github.com/cdeil/healpix
* Docs : http://healpix.readthedocs.io
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
from ._astropy_init import *

# For egg_info test builds to pass, put package imports here.
if not _ASTROPY_SETUP_:
    from .high_level import *  # noqa
