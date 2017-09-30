"""Python 2 / 3 compatibility helpers.

Usually we would import and depend on
https://github.com/benjaminp/six/blob/master/six.py
but since in this package we need so little,
we just copied the few lines over from there
"""
import sys
PY2 = sys.version_info[0] == 2
PY3 = sys.version_info[0] == 3


if PY3:
    integer_types = int,
else:
    integer_types = (int, long)
