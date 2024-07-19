# From https://github.com/numpy/numpy/blob/main/tools/wheels/cibw_before_build.sh
set -xe

if [[ $RUNNER_OS == "Windows" ]]; then 
    outdir=C:\openblas
    pip install delvewheel
fi 

if [[ $RUNNER_OS == "Linux" ]]; then
    outdir=/usr/local/lib
    sudo chmod 777 /usr/local/lib
fi

if [[ $RUNNER_OS == "macOS" ]]; then 
    outdir=/usr/local/lib
    sudo mkdir -p /usr/local/lib
    sudo chmod 777 /usr/local/lib
fi

pip install scipy-openblas64

python <<EOF
import os, scipy_openblas64, shutil
srcdir = os.path.join(os.path.dirname(scipy_openblas64.__file__), "lib")
shutil.copytree(srcdir, "$outdir", dirs_exist_ok=True)
srcdir = os.path.join(os.path.dirname(scipy_openblas64.__file__), ".dylibs")
if os.path.exists(srcdir):  # macosx delocate
    shutil.copytree(srcdir, "$outdir", dirs_exist_ok=True)
EOF

ls $outdir
echo $outdir
