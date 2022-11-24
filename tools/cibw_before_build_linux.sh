# https://github.com/scipy/scipy/blob/main/tools/wheels/cibw_before_build_linux.sh

# set -xe

# # Install Openblas
# basedir=$(python tools/openblas_support.py)
# cp -r $basedir/lib/* /usr/local/lib
# cp $basedir/include/* /usr/local/include

# conda activate openblas
# conda install openblas libopenblas

yum install openblas-devel lapack-devel -y