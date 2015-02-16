rm -rf spams-matlab
svn export ./ spams-matlab
rm -rf spams-matlab/swig/
rm spams-matlab/make_matlab_package.sh
rm spams-matlab/TODO
da=$(``date +%F)
echo $da
tar -czf spams-matlab-v2.5-svn$da.tar.gz spams-matlab
rm -rf spams-matlab

