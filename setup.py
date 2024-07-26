import os
import platform
import distro

from setuptools import setup, Extension, find_packages
from distutils.sysconfig import get_python_inc
from extension_helpers import add_openmp_flags_if_available, pkg_config

import numpy as np

ismac = platform.system() == 'Darwin'
iswindows = platform.system() == 'Windows'

def get_config():

    blasinfo = pkg_config(['openblas'], ['openblas64'])
    incs = ['spams_wrap']
    for x in ['linalg', 'prox', 'decomp', 'dictLearn']:
        incs.append(os.path.join('spams_wrap', 'spams', x))
    incs.append(np.get_include())
    incs.append(get_python_inc())
    incs.extend(blasinfo.get('include_dirs', []))

    cc_flags = ['-fPIC', '-m64']

    for _ in blasinfo.get('extra_compile_args', []):
        if _ not in cc_flags:
            cc_flags.append(_)
    # for _ in lapackinfo.get('extra_compile_args', []):
    #     if _ not in cc_flags:
    #         cc_flags.append(_)

    link_flags = []
    for _ in blasinfo.get('extra_link_args', []):
        if _ not in link_flags:
            link_flags.append(_)
    # for _ in lapackinfo.get('extra_link_args', []):
    #     if _ not in link_flags:
    #         link_flags.append(_)

    if iswindows:
        libs = []
        is_mkl = False
    else:
        libs = ['stdc++']

        is_mkl = False
        for lib in blasinfo.get('libraries', []):
            if 'mkl' in lib:
                is_mkl = True
                break

    libdirs = blasinfo.get('library_dirs', [])
    if is_mkl:
        for _ in blasinfo.get('include_dirs', []):
            if _ not in incs:
                incs.append(_)
        for _ in blasinfo.get('library_dirs', []):
            if _ not in libdirs:
                libdirs.append(_)
        libs.extend(['mkl_rt'])
    else:
        if 'centos' in distro.id() or 'alma' in distro.id():
            libs.extend(['openblaso'])  # for openmp support in openblas under redhat
        else:
            libs.extend(['openblas'])

    # Check for openmp flag
    if iswindows:
        cc_flags.append('/openmp')
    elif ismac:
        cc_flags.extend(['-Xpreprocessor', '-fopenmp'])
        link_flags.append('-lomp')
    else:
        cc_flags.append('-fopenmp')
        link_flags.append('-fopenmp')

    if ismac:
        # homebrew path x64
        cc_flags.append('-I/usr/local/opt/openblas/include')
        link_flags.append('-L/usr/local/opt/openblas/lib')

        cc_flags.append('-I/usr/local/opt/libomp/include')
        link_flags.append('-L/usr/local/opt/libomp/lib')

        # homebrew openblas path arm64
        cc_flags.append('-I/opt/homebrew/opt/openblas/include')
        link_flags.append('-L/opt/homebrew/opt/openblas/lib')

        cc_flags.append('-I/opt/homebrew/opt/libomp/include')
        link_flags.append('-L/opt/homebrew/opt/libomp/lib')

        # # use accelerate
        # link_flags.append('-framework accelerate')
        # cc_flags.append('-I/usr/local/opt/openblas/include')

    if iswindows:
        # dir_path = os.path.dirname(os.path.realpath(__file__))
        # Look for local intel mkl
        # libpath = os.path.join(dir_path, 'lib', 'native', 'win-x64')
        libs.append('openblas')
        libpath = os.path.join('C:/Miniconda/envs/openblas/Library/lib')
        # libpath = os.path.join('C:/opt/openblas/openblas_dll')
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

long_description = """Python interface for SPArse Modeling Software (SPAMS),
an optimization toolbox for solving various sparse estimation problems.
This (unofficial) version includes pre-built wheels for python 3 on windows, mac (with openmp support) and linux.
The source code for this fork is also available at https://github.com/samuelstjean/spams-python/"""

setup(name='spams-bin',
      python_requires='>=3.9.0',
      version='2.6.10',
      description='Python interface for SPAMS - binary wheels with openblas (mac, linux, windows)',
      long_description=long_description,
      author='Julien Mairal',
      author_email='spams.dev@inria.fr',
      url='http://spams-devel.gforge.inria.fr/',
      ext_modules=[spams_wrap],
      packages=find_packages(),
      install_requires=['numpy>=1.21.3',
                        'scipy>=1.5'],
      extras_require={"test": ['Pillow>=6.0']})
