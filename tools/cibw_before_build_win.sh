# https://github.com/scipy/scipy/blob/main/tools/wheels/cibw_before_build_win.sh

set -xe

# Install Openblas
# PYTHONPATH=tools python -c "import openblas_support; openblas_support.make_init('spams')"
mkdir -p /c/opt/openblas/if_32/64/lib/pkgconfig

# delvewheel is the equivalent of delocate/auditwheel for windows.
python -m pip install delvewheel

# make the DLL available for tools/wheels/repair_windows.sh. If you change
# this location you need to alter that script.
mkdir -p /c/opt/openblas/openblas_dll
which strip

# 64-bit openBLAS
curl -L https://github.com/scipy/scipy-ci-artifacts/raw/main/openblas_32_if.zip -o openblas_32_if.zip
unzip openblas_32_if.zip -d /c
cp /c/opt/openblas/if_32/64/bin/*.dll /c/opt/openblas/openblas_dll
cp -r /c/opt/openblas/if_32/64/include /c/opt/openblas/
cp -r /c/opt/openblas/if_32/64/lib /c/opt/openblas/
