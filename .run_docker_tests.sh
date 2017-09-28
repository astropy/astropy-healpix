#!/bin/bash

docker info

cat << EOF | docker run -i \
                        -v ${PWD}:/regions_src \
                        -a stdin -a stdout -a stderr \
                        astropy/affiliated-32bit-test-env:1.7 \
                        bash || exit $?

cd /regions_src

echo "Output of uname -m:"
uname -m

echo "Output of sys.maxsize in Python:"
python -c 'import sys; print(sys.maxsize)'

pip install six

# We only test the docs output on 64-bit Linux so as not to have to 
# degrade the output with ellipses to work on all platforms.
python setup.py test -V -a "-s" --skip-docs

EOF
