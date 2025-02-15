# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
BSD-licensed HEALPix for Astropy.
"""

from .high_level import *  # noqa
from .core import *  # noqa

try:
    from .version import version as __version__
except ImportError:
    __version__ = ''
