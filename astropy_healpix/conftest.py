# This file is used to configure the behavior of pytest when using the Astropy
# test infrastructure. It needs to live inside the package in order for it to
# get picked up when running the tests inside an interpreter using
# astropy_helpers.test

import os

import numpy as np

from astropy.version import version as astropy_version

# For Astropy 3.0 and later, we can use the standalone pytest plugin
if astropy_version < '3.0':
    from astropy.tests.pytest_plugins import *  # noqa
    del pytest_report_header
    ASTROPY_HEADER = True
else:
    try:
        from pytest_astropy_header.display import PYTEST_HEADER_MODULES, TESTED_VERSIONS
        ASTROPY_HEADER = True
    except ImportError:
        ASTROPY_HEADER = False


def pytest_configure(config):

    if ASTROPY_HEADER:

        config.option.astropy_header = True

        PYTEST_HEADER_MODULES.pop('h5py', None)
        PYTEST_HEADER_MODULES.pop('Pandas', None)
        PYTEST_HEADER_MODULES['Astropy'] = 'astropy'
        PYTEST_HEADER_MODULES['healpy'] = 'healpy'

        from . import __version__
        packagename = os.path.basename(os.path.dirname(__file__))
        TESTED_VERSIONS[packagename] = __version__

    # Set the Numpy print style to a fixed version to make doctest outputs
    # reproducible.
    try:
        np.set_printoptions(legacy='1.13')
    except TypeError:
        # On older versions of Numpy, the unrecognized 'legacy' option will
        # raise a TypeError.
        pass

from astropy.tests.helper import enable_deprecations_as_exceptions  # noqa
enable_deprecations_as_exceptions()
