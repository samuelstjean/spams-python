# From https://github.com/numpy/numpy/blob/main/tools/wheels/cibw_before_build.sh
set -xe

PROJECT_DIR="${1:-$PWD}"

if [[ $RUNNER_OS == "Windows" ]]; then 
    outdir = C:\openblas
fi 

if [[ $RUNNER_OS == "Linux" ]]; then
    outdir = /usr/local/lib
fi

if [[ $RUNNER_OS == "MACOS" ]]; then 
    outdir = /usr/local/lib
fi

python <<EOF
import os, scipy_openblas64, shutil
srcdir = os.path.join(os.path.dirname(scipy_openblas64.__file__), "lib")
shutil.copytree(srcdir, os.path.join("$outdir", "lib"))
srcdir = os.path.join(os.path.dirname(scipy_openblas64.__file__), ".dylibs")
if os.path.exists(srcdir):  # macosx delocate
    shutil.copytree(srcdir, os.path.join("$outdir", ".dylibs"))
EOF

