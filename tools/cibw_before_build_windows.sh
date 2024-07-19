# From https://github.com/numpy/numpy/blob/main/tools/wheels/cibw_before_build.sh
set -xe

set outdir=C:\openblas
pip install delvewheel
pip install scipy-openblas64

ls $outdir
echo $outdir

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
