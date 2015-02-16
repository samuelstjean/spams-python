#!/bin/bash
case $1 in -x) set -x;shift;;esac

VC="/c/Program Files (x86)/Microsoft Visual Studio 9.0/VC/bin/amd64"
PATH="$VC:/c/Program Files (x86)/Microsoft Visual Studio 9.0/Common7/IDE":$PATH
/c/Python27/python setup.py build

/c/Python27/python setup.py bdist_wininst
