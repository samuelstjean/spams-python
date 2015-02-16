set PATH=%PATH%;C:\MinGW\bin;C:\MinGW\msys\1.0\bin;"C:\Program Files\R\R-2.15.1\bin";C:\cygwin\bin
mkdir spams\inst\libs\i386
set srcd=C:\MinGW\bin
set dstd=spams\inst\libs\i386
copy %srcd%\libgomp-1.dll %dstd%
copy %srcd%\"libstdc++-6.dll" %dstd%
copy %srcd%\libgcc_s_dw2-1.dll %dstd%
copy %srcd%\pthreadGC2.dll %dstd%
R CMD INSTALL --build spams

