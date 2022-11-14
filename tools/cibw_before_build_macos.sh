# https://github.com/scipy/scipy/blob/main/tools/wheels/cibw_before_build_macos.sh

set -xe

# Install Openblas
basedir=$(python tools/openblas_support.py)
cp -r $basedir/lib/* /usr/local/lib
cp $basedir/include/* /usr/local/include
