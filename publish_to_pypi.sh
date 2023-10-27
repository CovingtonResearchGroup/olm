#!/bin/bash

#
rm -r ./dist
rm -r ./build
#python -m build
python setup.py sdist bdist_wheel
twine upload dist/*


