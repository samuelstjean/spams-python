import os
import platform
import subprocess

from setuptools import setup, Extension, find_packages
from sysconfig import get_paths
from openmp_helpers import add_openmp_flags_if_available
from tempfile import gettempdir

import numpy as np
from numpy.distutils.system_info import blas_info


def get_config():

    incs = ['spams_wrap']
    for x in ['linalg', 'prox', 'decomp', 'dictLearn']:
        incs.append(os.path.join('spams_wrap', x))
    incs.append(np.get_include())
    incs.append(get_paths()['include'])
    incs.extend(blas_info().get_include_dirs())

    cc_flags = ['-fPIC', '-m64']

    for _ in np.__config__.blas_opt_info.get('extra_compile_args', []):
        if _ not in cc_flags:
            cc_flags.append(_)
    for _ in np.__config__.lapack_opt_info.get('extra_compile_args', []):
        if _ not in cc_flags:
            cc_flags.append(_)

    link_flags = []
    for _ in np.__config__.blas_opt_info.get('extra_link_args', []):
        if _ not in link_flags:
            link_flags.append(_)
    for _ in np.__config__.lapack_opt_info.get('extra_link_args', []):
        if _ not in link_flags:
            link_flags.append(_)

    libs = ['stdc++']

    is_mkl = False
    for lib in np.__config__.blas_opt_info.get('libraries', []):
        if 'mkl' in lib:
            is_mkl = True
            break

    if not is_mkl:
        # Grab a fresh openblas for the current platform
        cmd = 'python', 'openblas_support.py'
        subprocess.run(cmd)
        openblasdir = os.path.join(gettempdir(), 'openblas')
        incs.append(os.path.join(openblasdir, 'include'))

    libdirs = blas_info().get_lib_dirs()
    if is_mkl:
        for _ in np.__config__.blas_opt_info.get('include_dirs', []):
            if _ not in incs:
                incs.append(_)
        for _ in np.__config__.blas_opt_info.get('library_dirs', []):
            if _ not in libdirs:
                libdirs.append(_)
        libs.extend(['mkl_rt'])
    else:
        libs.extend(['openblas'])
        libdirs.append(os.path.join(openblasdir, 'lib'))

    # Check for openmp flag, mac is done later
    if platform.system() != 'Darwin':
        if platform.system() == 'Windows':
            cc_flags.append('-openmp')
            # link_flags.append('-openmp')
        else:
            cc_flags.append('-fopenmp')
            link_flags.append('-fopenmp')

    return incs, libs, libdirs, cc_flags, link_flags


incs, libs, libdirs, cc_flags, link_flags = get_config()

if platform.system() == 'Windows':
    source = ['spams_wrap/spams_wrap-windows.cpp']
else:
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
      version='2.6.4',
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
