# From https://github.com/scipy/scipy/blob/main/tools/wheels/cibw_before_build_win.sh
set -xe

# Install Openblas
mkdir -p /c/opt/64/lib/pkgconfig
mkdir -p /c/opt/openblas/openblas_dll

# delvewheel is the equivalent of delocate/auditwheel for windows.
python -m pip install delvewheel

target=$(python -c "import tools.openblas_support as obs; plat=obs.get_plat(); ilp64=obs.get_ilp64(); target=f'openblas_{plat}.zip'; obs.download_openblas(target, plat, ilp64);print(target)")

# 64-bit openBLAS
unzip $target -d /c/opt/
cp /c/opt/64/bin/*.dll /c/opt/openblas/openblas_dll
