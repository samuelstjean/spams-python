#!/bin/bash
case $1 in -x) set -x;shift;;esac


# modify the following variable if you have a different version of R
Rdir="/c/Program Files/R/R-2.15.1"
###################
# NB: this script needs GnuWin zip and unzip

die () {
    echo "$*"
    exit 1
}

ldir=spams/inst/libs
SVPATH=$PATH
##
[  -d "$Rdir" ] || die "I don't find the R directory $Rdir!"

[ -d spams ] || die "missing ./spams directory!"
cdir=`/bin/pwd`

rm -f *.zip
arch=i386
archname=win32
rm -f spams/src/spams.o
rm -rf spams/inst/libs

ZPATH="/c/Program Files (x86)/GnuWin32/bin:/c/Program Files/GnuWin32/bin"
#PATH=$SVPATH:"$Rdir/bin":/c/cygwin/bin
PATH=$SVPATH:"$Rdir/bin":"$ZPATH"
dstd=$ldir/$arch;mkdir -p $dstd
srcd=/c/MinGW/bin
cp -p $srcd/libgomp-1.dll $srcd/"libstdc++-6.dll" $srcd/libgcc_s_*.dll $srcd/pthread*.dll $dstd
R --arch i386 CMD INSTALL --html --no-multiarch --build spams
zipname=`echo spams_*.zip`
rm -rf $arch;mkdir $arch;mv $zipname $arch
if [  -d "$Rdir/bin/x64" ];then
    arch=x64
    archname=win32+amd64
#PATH=/c/mingw64/bin:$SVPATH:"/c/Program Files/R/R-2.15.1/bin/x64":/c/cygwin/bin:/c/mingw64/x86_64-w64-mingw32/bin
    PATH=/c/mingw64/bin:$SVPATH:"/c/Program Files/R/R-2.15.1/bin/x64":"$ZPATH":/c/mingw64/x86_64-w64-mingw32/bin
    rm -f spams/src/spams.o
    rm -rf spams/inst/libs
    dstd=$ldir/$arch;mkdir -p $dstd
    srcd=/c/mingw64/bin
    cp -p $srcd/libgomp-1.dll $srcd/"libstdc++-6.dll" $srcd/libgcc_s_*.dll $dstd
    cp -p /c/mingw64/x86_64-w64-mingw32/bin/pthread*.dll $dstd
    R --arch $arch CMD INSTALL --no-multiarch --build spams
    rm -rf $arch;mkdir $arch;mv $zipname $arch
    
# merge 32 and 64 bits packages 
    cd i386;unzip -x $zipname
    cd ../x64;unzip -x $zipname
    mv ../i386/spams/libs/i386 spams/libs
    cd spams
    rm -f MD5
    md5sum `find . -type f` | sed 's:\./:*:' >../MD5
    mv ../MD5 .
    cd ..
    find spams | zip -@ ../$zipname
else
    cp -p $arch/$zipname .
fi
cd $cdir
x="`ls -d  /c/Program\ Files*/Inno\ Setup*`"
if [ -z "$x" ]; then
    echo "InnoSetup is not installed. I cannot create spams-R.exe"
    exit 0
fi
iss="$x/Compil32"
dst=/c/WINDOWS/Temp/spams-R
if [ -r Release-name ]; then
    vers=`cat Release-name`
    mkdir $dst
    cp -p $zipname spams/LICENSE INSTALL-windows $dst
    "$iss" //cc spams-R.iss
    [ -r Output/spams-R.exe ] || {
	echo "Error"
	exit 1
    }
    mv Output/spams-R.exe spams-R-$vers.$archname.exe
fi
exit 0

