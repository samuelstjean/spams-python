# https://github.com/scipy/scipy/blob/main/tools/wheels/cibw_before_build_win.sh

set -xe

# Install Openblas
# PYTHONPATH=tools python -c "import openblas_support; openblas_support.make_init('spams')"
# mkdir -p C:/opt/openblas/if_32/64/lib/pkgconfig

# delvewheel is the equivalent of delocate/auditwheel for windows.
python -m pip install delvewheel

# make the DLL available for tools/wheels/repair_windows.sh. If you change
# this location you need to alter that script.
mkdir -p C:/opt/openblas/openblas_dll
# which strip

# 64-bit openBLAS
curl -L https://github.com/xianyi/OpenBLAS/releases/download/v0.3.21/OpenBLAS-0.3.21-x64.zip -o openblas_32_if.zip
unzip -f openblas_32_if.zip -d C:/opt/openblas
cp C:/opt/openblas/bin/*.dll C:/opt/openblas/openblas_dll
# cp -r C:/opt/openblas/if_32/64/include C:/opt/openblas/
# cp -r C:/opt/openblas/if_32/64/lib C:/opt/openblas/
# cp C:/opt/openblas/if_32/64/bin/*.dll C:/opt/openblas/lib
