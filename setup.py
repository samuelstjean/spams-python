import os
import platform
import sys

from setuptools import setup, Extension
# from setuptools.command.build_ext import build_ext
from distutils.sysconfig import get_python_inc


def get_config():
    # Import numpy here, only when headers are needed
    import numpy as np
    from numpy.distutils.system_info import blas_info

    incs = ['spams']
    for x in ['linalg', 'prox', 'decomp', 'dictLearn']:
        incs.append(os.path.join('spams', x))
    incs.append(np.get_include())
    incs.append(get_python_inc())
    incs.extend(blas_info().get_include_dirs())

    cc_flags = ['-fPIC']
    if sys.maxsize > 2**32:
        cc_flags.append('-m64')
    else:
        cc_flags.append('-m32')

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

    if platform.system() == 'Windows':
        libs = []
    else:
        libs = ['stdc++']

    is_mkl = False
    for lib in np.__config__.blas_opt_info.get('libraries', []):
        if 'mkl' in lib:
            is_mkl = True
            break

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
        libs.extend(['blas', 'lapack'])

    if platform.system() != 'Darwin':
        cc_flags.append('-fopenmp')
        link_flags.append('-fopenmp')

    return incs, libs, libdirs, cc_flags, link_flags


incs, libs, libdirs, cc_flags, link_flags = get_config()

if platform.system() == 'Windows':
    source = ['spams_wrap-windows.cpp']
else:
    source = ['spams_wrap.cpp']

spams_wrap = Extension(
    '_spams_wrap',
    sources=source,
    include_dirs=incs,
    extra_compile_args=['-DNDEBUG', '-DUSE_BLAS_LIB'] + cc_flags,
    library_dirs=libdirs,
    libraries=libs,
    # strip the .so
    extra_link_args=link_flags,
    language='c++',
    depends=['spams.h'],
)


def mkhtml(d=None, base='sphinx'):
    if d is None:
        d = base
    else:
        d = os.path.join(base, d)
    if not os.path.isdir(base):
        return []
    hdir = d

    l1 = os.listdir(hdir)
    l = []
    for s in l1:
        s = os.path.join(d, s)
        if not os.path.isdir(s):
            l.append(s)
    return l


setup(name='spams',
      version='2.6',
      description='Python interface for SPAMS',
      author='Julien Mairal',
      # author_email='nomail',
      # url='http://',
      ext_modules=[spams_wrap, ],
      py_modules=['spams', 'spams_wrap', 'myscipy_rand'],
      #       scripts = ['test_spams.py'],
      data_files=[('test', ['test_spams.py', 'test_decomp.py', 'test_dictLearn.py',
                            'test_linalg.py', 'test_prox.py', 'test_utils.py'])],
      )
