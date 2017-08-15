#!/bin/sh

#backup files we want to manually maintain
mkdir ./tmp
cp ./docs/conf.py ./tmp/
cp ./docs/index.rst ./tmp/
cp ./docs/olm.rst ./tmp/
cp ./docs/olm.USGS.rst ./tmp/
cp ./docs/olm.general.solution.rst ./tmp/
cp ./docs/read-the-docs-requirements.txt ./tmp/

#remove previous docs
rm -r ./docs/*

#generate new documents using sphinx
sphinx-apidoc -d 2 -F -o ./docs/ olm

#copy manually maintained files back into docs
cp ./tmp/* ./docs/

#make new docs
cd ./docs/
make html
cd ..

#remove temporary directory
rm -r ./tmp/
