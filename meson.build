project('_spams_wrap',
        ['cpp'],
        version: '2.6.12',
        meson_version: '>= 1.2.3',
        default_options : [
                            'warning_level=1',
                            'buildtype=release'
                          ])

numpy_nodepr_api = '-DNPY_NO_DEPRECATED_API=NPY_1_23_API_VERSION'
extra_args = ['-fPIC',
              '-DNDEBUG',
              '-m64']

python_sources = [
    'spams/__init__.py',
    'spams/spams.py',
    'spams_wrap/spams_wrap.py',
]

test_sources = [
  'test/test_decomp.py',
  'test/test_dictLearn.py',
  'test/test_linalg.py',
  'test/test_prox.py',
  'test/test_utils.py',
  'test/boat.png',
  'test/lena.png',
]

spams_include = [
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

np_dep = declare_dependency(include_directories: incdir_numpy,
                            compile_args: numpy_nodepr_api)
omp = dependency('openmp')
blas = dependency('blas')
lapack = dependency('lapack')

# headers =
                    #  ['''spams_wrap/spams/decomp/decomp.h''',
                    #  '''spams_wrap/spams/decomp/lsqsplx.h''',
                    #  '''spams_wrap/spams/decomp/projsplx.h''',
                    #  '''spams_wrap/spams/dictLearn/arch.h''',
                    #  '''spams_wrap/spams/dictLearn/dicts.h''',
                    #  '''spams_wrap/spams/linalg/cblas_alt_template.h''',
                    #  '''spams_wrap/spams/linalg/cblas_defvar.h''',
                    #  '''spams_wrap/spams/linalg/cblas_template.h''',
                    #  '''spams_wrap/spams/linalg/linalg.h''',
                    #  '''spams_wrap/spams/linalg/list.h''',
                    #  '''spams_wrap/spams/linalg/mexutils.h''',
                    #  '''spams_wrap/spams/linalg/misc.h''',
                    #  '''spams_wrap/spams/linalg/utils.h''',
                    #  '''spams_wrap/spams/prox/fista.h''',
                    #  '''spams_wrap/spams/prox/groups-graph.h''',
                    #  '''spams_wrap/spams/prox/mexgrouputils.h''',
                    #  '''spams_wrap/spams/prox/project.h''',
                    #  '''spams_wrap/spams/prox/surrogate.h''',
                    #  '''spams_wrap/spams/prox/svm.h''']

py.extension_module('spams_wrap',
                     [
                      'spams_wrap/spams_wrap.cpp',
                     ],
                     include_directories: spams_include,
                     dependencies : [
                     py_dep,
                     np_dep,
                     omp,
                     blas,
                     lapack,
                     'spams_wrap/spams.h'
                     ],
                     extra_args: extra_args,
                     install : true,
                     subdir: 'spams_wrap')

py.install_sources(
  python_sources,
  pure: false,             # Will be installed next to binaries
  subdir: 'spams_wrap'  # Folder relative to site-packages to install to
)

py.install_sources(
  test_sources,
  pure: false,             # Will be installed next to binaries
  preserve_path: true
)
