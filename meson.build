project('_spams_wrap',
        'cython',
        'cpp',
        version: '2.6.11',
        meson_version: '>= 1.2.3',
        default_options : [
                            'warning_level=1',
                            'buildtype=release'
                          ])

cpp_compiler = meson.get_compiler('cpp')

# Find some macos libs
m_dep = cpp_compiler.find_library('m', required : false)
dl_dep = cpp_compiler.find_library('dl', required : false)

extra_args = ['-fPIC',
              '-DNDEBUG']

python_sources = [
    'spams/__init__.py',
    'spams/spams.py',
]

python_wrap_sources = [
    'spams_wrap/__init__.py',
    'spams_wrap/spams_wrap.py',
]

test_sources = [
  'test/test_decomp.py',
  'test/test_dictLearn.py',
  'test/test_linalg.py',
  'test/test_prox.py',
  'test/boat.png',
  'test/lena.png',
]

spams_include = [
                'spams_wrap',
                'spams_wrap/spams/decomp',
                'spams_wrap/spams/dictLearn',
                'spams_wrap/spams/linalg',
                'spams_wrap/spams/prox'
                ]

py = import('python').find_installation(pure: false)
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

np_dep = declare_dependency(include_directories: incdir_numpy)

openmp = dependency('openmp')
# blas = dependency('blas', required: false)
# lapack = dependency('lapack', required: false)

# # centos has no pkg-config file apparently, so look for openblas if it was not found previously
# if not blas.found()
#   blas = cpp_compiler.find_library('openblas', required: true)
#   if not blas.found()
#     blas = dependency('openblas', required: true)
#   endif
# endif

# if not lapack.found()
#   lapack = cpp_compiler.find_library('openblas', required: true)
#   if not lapack.found()
#     lapack = dependency('openblas', required: true)
#   endif
# endif

cython_cpp_lib = static_library(
    'blas_lapack',
    [
      'spams_wrap/spams/linalg/cython_blas_lapack.h',
      'spams_wrap/spams/linalg/cython_blas_lapack.pyx'
    ],
    override_options : ['cython_language=cpp'],
    dependencies : [
                  py_dep,
                  np_dep
    ]
)

py.extension_module('_spams_wrap',
                     [
                      'spams_wrap/spams/linalg/cython_blas_lapack.h',
                      # 'spams_wrap/spams/linalg/cython_blas_lapack.c',
                      'spams_wrap/spams/linalg/cython_blas_lapack.pyx',
                      'spams_wrap/spams_wrap.cpp'
                     ],
                     include_directories: spams_include,
                     override_options: ['cython_language=cpp'],
                     link_with : [cython_cpp_lib],
                     dependencies : [
                     py_dep,
                     np_dep,
                     m_dep,
                     dl_dep,
                     openmp,
                    #  blas,
                    #  lapack,
                     ],
                     c_args: extra_args,
                     install : true,
                     subdir: 'spams_wrap')

py.install_sources(
  python_sources,
  pure: false,             # Will be installed next to binaries
  subdir: 'spams'  # Folder relative to site-packages to install to
)

py.install_sources(
  python_wrap_sources,
  pure: false,             # Will be installed next to binaries
  subdir: 'spams_wrap'  # Folder relative to site-packages to install to
)

py.install_sources(
  test_sources,
  pure: false,             # Will be installed next to binaries
  subdir: 'test'  # Folder relative to site-packages to install to
)
