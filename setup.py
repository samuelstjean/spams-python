import os
import platform

from setuptools import setup, Extension, find_packages
from distutils.sysconfig import get_python_inc
from openmp_helpers import add_openmp_flags_if_available

import numpy as np
from numpy.distutils.system_info import blas_info

import distro


# for np >= 1.22
try:
    blasinfo = np.distutils.__config__.blas_ilp64_opt_info
    lapackinfo = np.distutils.__config__.lapack_ilp64_opt_info
except Exception:
    blasinfo = np.__config__.blas_opt_info
    lapackinfo = np.__config__.lapack_opt_info

def get_config():

    incs = ['spams_wrap']
    for x in ['linalg', 'prox', 'decomp', 'dictLearn']:
        incs.append(os.path.join('spams_wrap', 'spams', x))
    incs.append(np.get_include())
    incs.append(get_python_inc())
    incs.extend(blas_info().get_include_dirs())

    cc_flags = ['-fPIC', '-m64']

    for _ in blasinfo.get('extra_compile_args', []):
        if _ not in cc_flags:
            cc_flags.append(_)
    for _ in lapackinfo.get('extra_compile_args', []):
        if _ not in cc_flags:
            cc_flags.append(_)

    link_flags = []
    for _ in blasinfo.get('extra_link_args', []):
        if _ not in link_flags:
            link_flags.append(_)
    for _ in lapackinfo.get('extra_link_args', []):
        if _ not in link_flags:
            link_flags.append(_)

    if platform.system() == 'Windows':
        libs = []
        is_mkl = False
    else:
        libs = ['stdc++']

        is_mkl = False
        for lib in blasinfo.get('libraries', []):
            if 'mkl' in lib:
                is_mkl = True
                break

    libdirs = blas_info().get_lib_dirs()
    if is_mkl:
        for _ in blasinfo.get('include_dirs', []):
            if _ not in incs:
                incs.append(_)
        for _ in blasinfo.get('library_dirs', []):
            if _ not in libdirs:
                libdirs.append(_)
        libs.extend(['mkl_rt'])
    else:
        if 'centos' in distro.id():
            libs.extend(['openblaso', 'lapack'])  # for openmp support in openblas
        else:
            libs.extend(['openblas'])

    # Check for openmp flag, mac is done later
    if platform.system() != 'Darwin':
        if platform.system() == 'Windows':
            cc_flags.append('-openmp')
            # link_flags.append('-openmp')
        else:
            cc_flags.append('-fopenmp')
            link_flags.append('-fopenmp')

    if platform.system() == 'Darwin':
        cc_flags.append('-I/usr/local/opt/openblas/include')
        link_flags.append('-L/usr/local/opt/openblas/lib')

    if platform.system() == 'Windows':
        # dir_path = os.path.dirname(os.path.realpath(__file__))
        # Look for local intel mkl
        # libpath = os.path.join(dir_path, 'lib', 'native', 'win-x64')
        libs.append('openblas')
        libpath = os.path.join('C:/Miniconda/envs/openblas/Library/lib')
        libdirs.append(libpath)
        incs.append('C:/Miniconda/envs/openblas/Library/include')
        incs.append('C:/Miniconda/envs/openblas/Library/include/openblas')

    return incs, libs, libdirs, cc_flags, link_flags


incs, libs, libdirs, cc_flags, link_flags = get_config()
source = ['spams_wrap/spams_wrap.cpp']

spams_wrap = Extension(
    'spams_wrap._spams_wrap',
    sources=source,
    include_dirs=incs,
    extra_compile_args=['-DNDEBUG', '-DUSE_BLAS_LIB'] + cc_flags,
    library_dirs=libdirs,
    libraries=libs,
    extra_link_args=link_flags,
    language='c++',
    depends=['spams_wrap/spams.h'],
)

if platform.system() == 'Darwin':
    add_openmp_flags_if_available(spams_wrap)

long_description = """Python interface for SPArse Modeling Software (SPAMS),
an optimization toolbox for solving various sparse estimation problems.
This (unofficial) version includes pre-built wheels for python 3 on windows, mac (with openmp support) and linux.

In addition, building from source explicitly requires openblas on mac/linux and intel mkl on windows,
unlike the official version which can use any blas implementation.

The source code for this fork is also available at https://github.com/samuelstjean/spams-python/"""

setup(name='spams-bin',
      version='2.6.5',
      description='Python interface for SPAMS - binary wheels with openblas (mac, linux, windows)',
      long_description=long_description,
      author='Julien Mairal',
      author_email='spams.dev@inria.fr',
      url='http://spams-devel.gforge.inria.fr/',
      ext_modules=[spams_wrap],
      packages=find_packages(),
      install_requires=['numpy>=1.12',
                        'scipy>=0.19'],
      extras_require={"test": ['Pillow>=6.0']})
