# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from numpy.testing import assert_equal, assert_allclose
from .. import wrap


def test_fourty_two():
    actual = wrap.fourty_two()
    expected = 42
    assert_equal(actual, expected)
