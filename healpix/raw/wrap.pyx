# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

__all__ = [
    'fourty_two',
]

cdef extern from "example.h":
    int fourty_two_c_function(int)


def fourty_two(int x=1):
    """Answer to the Ultimate Question of Life, the Universe, and Everything.

    Let's see if this works (tests run, docstring shows up in the docs.

    https://en.wikipedia.org/wiki/42_(number)#The_Hitchhiker.27s_Guide_to_the_Galaxy
    """
    return fourty_two_c_function(x)
