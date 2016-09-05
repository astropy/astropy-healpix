# Makefile with some convenient quick ways to do common things

PROJECT = healpix
CYTHON ?= cython
PYTHON_VERSION=`python -c 'import sys; print("%i" % (sys.hexversion<0x03000000))'`

help:
	@echo 'See Makefile for available targets'

clean:
	rm -rf build docs/_build docs/api htmlcov MANIFEST $(PROJECT).egg-info .coverage
	find . -name "*.pyc" -exec rm {} \;
	find . -name "*.so" -exec rm {} \;
	find . -name __pycache__ | xargs rm -fr

clean-repo:
	@git clean -f -x -d

cython:
	find $(PROJECT) -name "*.pyx" -exec $(CYTHON) {} \;

trailing-spaces:
	find $(PROJECT) examples docs -name "*.py" -exec perl -pi -e 's/[ \t]*$$//' {} \;

code-analysis: flake8 pylint

flake8:
	flake8 $(PROJECT) | grep -v __init__ | grep -v external

# TODO: once the errors are fixed, remove the -E option and tackle the warnings
pylint:
	pylint -E $(PROJECT)/ -d E1103,E0611,E1101 \
	       --ignore=_astropy_init.py -f colorized \
	       --msg-template='{C}: {path}:{line}:{column}: {msg} ({symbol})'
