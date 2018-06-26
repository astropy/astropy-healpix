# Playground for alternative implementations

I'd like to experiment a bit with alternative implementations
of HEALPix functionality here.

The goals are for me to get a better understanding how HEALPix
works, and to explore the option of a rewrite in Numpy,
Cython or C, or to play around with Numba a bit.

For now there's `healpix_py.py`, which is a pure Python port
from the Typescript lib https://github.com/michitaro/healpix.

Translating to C should be easy.
How does it compare to the current implementation?

I think even a translation to Numpy that works for arrays
directly should be possible?
It probalby wouldn't have good performance though, because
the expression would make more array copies and pass over
the array values more than once, i.e. use more memory and
be slower than a C/Cython/Numba implementation.
