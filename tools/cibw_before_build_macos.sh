# https://github.com/scipy/scipy/blob/main/tools/wheels/cibw_before_build_macos.sh

set -xe

# Openblas
basedir=$(python tools/openblas_support.py)

# copy over the OpenBLAS library stuff first
cp -r $basedir/lib/* /usr/local/lib
cp $basedir/include/* /usr/local/include

#GFORTRAN=$(type -p gfortran-9)
#sudo ln -s $GFORTRAN /usr/local/bin/gfortran
# same version of gfortran as the openblas-libs and scipy-wheel builds
curl -L https://github.com/MacPython/gfortran-install/raw/master/archives/gfortran-4.9.0-Mavericks.dmg -o gfortran.dmg
GFORTRAN_SHA256=$(shasum -a 256 gfortran.dmg)
KNOWN_SHA256="d2d5ca5ba8332d63bbe23a07201c4a0a5d7e09ee56f0298a96775f928c3c4b30  gfortran.dmg"
if [ "$GFORTRAN_SHA256" != "$KNOWN_SHA256" ]; then
    echo sha256 mismatch
    exit 1
fi

hdiutil attach -mountpoint /Volumes/gfortran gfortran.dmg
sudo installer -pkg /Volumes/gfortran/gfortran.pkg -target /
otool -L /usr/local/gfortran/lib/libgfortran.3.dylib
