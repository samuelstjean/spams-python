set PATH=%PATH%;C:\MinGW\bin;C:\MinGW\msys\1.0\bin
set srcd=C:\MinGW\bin
set dstd=build\lib.win32-2.7
C:\Python27\python setup.py build -c mingw32
copy %srcd%\libgomp-1.dll %dstd%
copy %srcd%\"libstdc++-6.dll" %dstd%
copy %srcd%\libgcc_s_dw2-1.dll %dstd%
copy %srcd%\pthreadGC2.dll %dstd%
set srcd="C:\Program Files\R\R-2.15.1\bin\i386"
copy %srcd%\*.dll %dstd%
C:\Python27\python setup.py bdist_wininst
