# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

import numpy as np

from astropy.version import version as astropy_version
if astropy_version < '3.0':
    from astropy.tests.pytest_plugins import *
    del pytest_report_header
else:
    from pytest_astropy_header.display import PYTEST_HEADER_MODULES, TESTED_VERSIONS

from astropy.tests.helper import enable_deprecations_as_exceptions

enable_deprecations_as_exceptions()


def pytest_configure(config):

    config.option.astropy_header = True

    PYTEST_HEADER_MODULES.pop('h5py', None)
    PYTEST_HEADER_MODULES.pop('Pandas', None)
    PYTEST_HEADER_MODULES['Astropy'] = 'astropy'
    PYTEST_HEADER_MODULES['healpy'] = 'healpy'

    from .version import version, astropy_helpers_version
    packagename = os.path.basename(os.path.dirname(__file__))
    TESTED_VERSIONS[packagename] = version
    TESTED_VERSIONS['astropy_helpers'] = astropy_helpers_version

    # Set the Numpy print style to a fixed version to make doctest outputs
    # reproducible.
    try:
        np.set_printoptions(legacy='1.13')
    except TypeError:
        # On older versions of Numpy, the unrecognized 'legacy' option will
        # raise a TypeError.
        pass
