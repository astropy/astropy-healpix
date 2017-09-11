# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import pytest
import numpy as np

from numpy.testing import assert_equal, assert_allclose

from .. import healpy as hp_compat

hp = pytest.importorskip('healpy')


def test_nside2resol():
    actual = hp_compat.nside2resol(nside=2)
    expected = hp.nside2resol(nside=2)
    assert_equal(actual, expected)


@pytest.mark.parametrize('nside,theta,phi,nest', [
    (256, 0.0000000000000000, 0.0000000000000000, True),
    (256, 0.0000000000000000, 1.2566370614359172, True),
    (256, 0.0000000000000000, 2.5132741228718345, True),
    (256, 0.0000000000000000, 3.7699111843077517, True),
    (256, 0.0000000000000000, 5.0265482457436690, True),
    (256, 0.0000000000000000, 6.2831853071795862, True)])
def test_ang2pix(nside, theta, phi, nest):
    ipix1 = hp_compat.ang2pix(nside, theta, phi, nest=nest)
    ipix2 = hp.ang2pix(nside, theta, phi, nest=nest)
    assert ipix1 == ipix2


@pytest.mark.parametrize('nside,ipix,nest', [
    (2, 0, True),
    (2, 1, True),
    (2, 2, True)])
def test_pix2ang(nside, ipix, nest):
    theta1, phi1 = hp_compat.pix2ang(nside, ipix, nest=nest)
    theta2, phi2 = hp.pix2ang(nside, ipix, nest=nest)
    assert_allclose(phi1, phi2, rtol=1e-10)
    assert_allclose(theta1, theta2, rtol=1e-10)
