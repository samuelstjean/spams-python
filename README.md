spams-python
============

The source code for this fork is also available at https://github.com/samuelstjean/spams-python/

A swig-regenerated rehost of the python version of SPArse Modeling Software (SPAMS) 2.6, available at http://spams-devel.gforge.inria.fr/downloads.html

Disclaimer : I am not the author of the package, I just host a patched version which builds on recent gcc and python 3 for (mostly my) convenience. Most of the patches were taken from https://aur.archlinux.org/packages/python-spams-svn/

As such, I can probably help out with small stuff, but for technical and theoretical details please contact the original authors at http://spams-devel.gforge.inria.fr/contacts.html

You can find the original swig wrapper to re-generate these files on the branch swig_generator at https://github.com/samuelstjean/spams-python/tree/swig_generator

This (unofficial) version includes pre-built wheels for python 3 on windows, mac (with openmp support) and linux and can be installed with ``pip install spams-bin``.

# Improvements from the original version

+ Automatic build system with meson
+ Pytest for testing
+ Newer swig wrappers
+ Fixed internal naming clash for newer openmp namespace and the orthogonal matching pursuit (omp) function
+ Fixed C++20 deprecation issues with templates
+ Fixed undefined behavior with openmp and passing 0 as the number of cores (it now defaults to 1)
+ Fixed a positivity bug with spams.lassoWeighted
+ Probably other stuff
