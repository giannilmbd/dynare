project('dynare',
        'cpp', 'fortran', 'c',
        version : '8-unstable',
        license : 'GPL-3.0-or-later',
        license_files : [ 'license.txt', 'COPYING' ],
        # NB: update C++ standard in .clang-format whenever the following is modified
        default_options : [ 'cpp_std=gnu++20', 'fortran_std=f2018',
                            'c_std=gnu17', 'warning_level=2' ],
        meson_version : '>=1.3.0')

add_global_arguments('-DPACKAGE_VERSION="' + meson.project_version() + '"', language : 'cpp')

if get_option('warning_level').to_int() >= 2
  add_global_arguments('-Wold-style-cast', language : 'cpp')
  add_global_arguments('-Wimplicit-interface', '-Wno-compare-reals', language : 'fortran')
endif

cpp_compiler = meson.get_compiler('cpp')
fortran_compiler = meson.get_compiler('fortran')
c_compiler = meson.get_compiler('c')

# If this is GCC, enforce minimal version
# NB: when modifiying this value, also update the gcc-*.ini native file under scripts/
# Other compilers are unsupported, but are nevertheless accepted by the present meson file
gcc_minimal_version = 13
if ((cpp_compiler.get_id() == 'gcc'
     and cpp_compiler.version().version_compare('<' + gcc_minimal_version.to_string()))
    or (fortran_compiler.get_id() == 'gcc'
        and fortran_compiler.version().version_compare('<' + gcc_minimal_version.to_string()))
    or (c_compiler.get_id() == 'gcc'
        and c_compiler.version().version_compare('<' + gcc_minimal_version.to_string())))
  error('GCC is too old. Version ' + gcc_minimal_version.to_string() + ' or higher is required.')
endif

# Put trampolines on the heap, for Fortran MEX files, on platforms where this is supported
# (as of GCC 14 and 15: x86_64 Darwin, x86_64 and aarch64 Linux)
# Seems only needed under Linux (likely the default under Darwin; not needed under Windows)
# Seems also only needed under Octave (MATLAB likely requests an executable stack to the OS)
# For more context, see https://git.dynare.org/Dynare/dynare/-/merge_requests/2428#note_25854
# Requires GCC ⩾ 14, which is why it is only conditionally enabled here after compiler checks
if (fortran_compiler.get_id() == 'gcc' and fortran_compiler.version().version_compare('>=14')
    and ((host_machine.system() == 'linux' and (host_machine.cpu_family() == 'x86_64' or
                                                host_machine.cpu_family() == 'aarch64'))
         or (host_machine.system() == 'darwin' and host_machine.cpu_family() == 'x86_64')))
  add_global_arguments('-ftrampoline-impl=heap', language : 'fortran')
elif get_option('build_for') == 'octave' and host_machine.system() == 'linux'
  warning('Your Fortran compiler does not support the -ftrampoline-impl=heap flag. As a consequence, you may need to set the environment variable GLIBC_TUNABLES=glibc.rtld.execstack=2 before running Octave, otherwise some MEX files may not load.')
endif


### Preprocessor

subdir('preprocessor/src')

# NB: the following two symlinks may be improved if this wishlist item is done:
#  https://github.com/mesonbuild/meson/issues/11519

install_symlink('dynare-preprocessor', install_dir : 'lib/dynare/preprocessor/',
                pointing_to : '../../../bin/dynare-preprocessor'
                              + (host_machine.system() == 'windows' ? '.exe' : ''))

# Compatibility symlink
install_symlink('dynare_m' + (host_machine.system() == 'windows' ? '.exe' : ''),
                install_dir : 'lib/dynare/matlab/preprocessor64',
                pointing_to : '../../../../bin/dynare-preprocessor'
                              + (host_machine.system() == 'windows' ? '.exe' : ''))

### M files

install_subdir('matlab', install_dir : 'lib/dynare',
               exclude_files : [ 'utilities/tests/.git' ,
                                 'utilities/tests/.gitignore',
                                 'dseries/.git',
                                 'dseries/.gitignore',
                                 'dseries/src/modules/matlab-fame-io/.git',
                                 'dseries/src/modules/matlab-fame-io/.gitignore',
                                 'dseries/.gitmodules',
                                 'dseries/.gitlab-ci.yml' ])

sed_exe = find_program('sed')
custom_target(output : 'dynare_version.m', input : 'matlab/dynare_version.m.in',
              command : [ sed_exe, 's/@PACKAGE_VERSION@/' + meson.project_version() + '/', '@INPUT@' ],
              capture : true,
              install : true,
              install_dir : 'lib/dynare/matlab')

install_subdir('contrib/ms-sbvar/TZcode/MatlabFiles',
               install_dir : 'lib/dynare/contrib/ms-sbvar/TZcode')

### MEX files

mex_incdir = include_directories('mex/sources')

## Various dependencies

if host_machine.system() != 'windows'
  dl_dep = dependency('dl')
else
  # Under Windows, we don’t use dlopen but rely on LoadLibrary
  dl_dep = []
endif

openmp_dep = dependency('openmp')

pthread_t_sizeof = c_compiler.sizeof('pthread_t', prefix : '#include <pthread.h>')
pthread_fortran_dep = declare_dependency(compile_args : '-DSIZEOF_PTHREAD_T=' + pthread_t_sizeof.to_string(), sources : 'mex/sources/pthread.F08')

## Determine MEX compilation options

if get_option('build_for') == 'matlab'
  matlab_path = get_option('matlab_path')
  if matlab_path == ''
    error('The “matlab_path” option must be provided when doing a MATLAB build')
  endif
  matlab_version = run_command('scripts/get-matlab-version', matlab_path, check : true).stdout().strip()
  matlab_minimal_version = [ '9.8', 'R2020a' ]

  if matlab_version.version_compare('<' + matlab_minimal_version[0])
    error('MATLAB is too old (version ' + matlab_version + '), please upgrade to version ' + matlab_minimal_version[0] + ' (' + matlab_minimal_version[1] + ') at least.')
  endif

  matlab_version_hex = run_command('scripts/get-matlab-version', '--hex', matlab_path, check : true).stdout().strip()
  matlab_exe = find_program(matlab_path / 'bin' / 'matlab', required : not meson.is_cross_build(), disabler : true)

  if host_machine.system() == 'linux' and host_machine.cpu_family() == 'x86_64'
    mexext = 'mexa64'
    matlab_arch = 'glnxa64'
    export_file = matlab_path / 'extern/lib/glnxa64/mexFunction.map'
    export_link_arg = '-Wl,--version-script,' + export_file
  elif host_machine.system() == 'windows' and host_machine.cpu_family() == 'x86_64'
    mexext = 'mexw64'
    matlab_arch = 'win64'
    export_file = matlab_path / 'extern/lib/win64/mingw64/mexFunction.def'
    export_link_arg = export_file
    if get_option('build_for') == 'matlab'
      arch_fortran_args = [ '-fno-underscoring' ]
    endif
  elif host_machine.system() == 'darwin'
    if host_machine.cpu_family() == 'x86_64'
      mexext = 'mexmaci64'
      matlab_arch = 'maci64'
    elif host_machine.cpu_family() == 'aarch64'
      mexext = 'mexmaca64'
      matlab_arch = 'maca64'
    else
      error('Unsupported platform')
    endif
    export_file = matlab_path / 'extern/lib' / matlab_arch / 'mexFunction.map'
    export_link_arg = '-Wl,-exported_symbols_list,' + export_file
  else
    error('Unsupported platform')
  endif

  matlab_defs = [ '-DMATLAB_MEX_FILE', '-DMATLAB_VERSION=' + matlab_version_hex, '-DMEXEXT=".' + mexext + '"' ]
  matlab_incdir = include_directories(matlab_path / 'extern/include')

  mex_kwargs = { 'name_prefix' : '',
                 'name_suffix' : mexext,
                 'include_directories' : [ mex_incdir, matlab_incdir ],
                 'cpp_args' : matlab_defs,
                 'fortran_args' : matlab_defs + [ '-fexceptions' ] + get_variable('arch_fortran_args', []),
                 'c_args' : matlab_defs + [ '-fexceptions' ],
                 'link_args' : [ export_link_arg, '-L' + (matlab_path / 'bin' / matlab_arch), '-lmx', '-lmex', '-lmat' ],
                 'link_depends' : export_file,
                 'install' : true,
                 'install_dir' : 'lib/dynare/mex/matlab',
                 'sources' : [] }

  # For unit tests
  exe_rpath = matlab_path / 'bin' / matlab_arch
  exe_link_args = [ '-L' + exe_rpath, '-lmx', '-lmex', '-lmat' ]

  # No need to use find_library() for the following libraries, since they are always shipped with MATLAB
  blas_dep = declare_dependency(link_args : '-lmwblas')
  lapack_dep = declare_dependency(link_args : '-lmwlapack', dependencies : blas_dep)
  umfpack_dep = declare_dependency(link_args : '-lmwumfpack', dependencies : blas_dep)
  ut_dep = declare_dependency(link_args : '-lut')
  extra_octave_dep = []

  if meson.version().version_compare('<1.7.0') and get_option('prefer_static') and host_machine.system() == 'linux'
    # Workaround for Meson bug https://github.com/mesonbuild/meson/issues/12757 (fixed in 1.7.0)
    # Use the C compiler as a fallback for detecting SLICOT under Linux with
    # prefer_static=true (but still try the Fortran compiler to honour the -B
    # option in fortran_args, as documented). Needed for building the MATLAB
    # Online package.
    slicot_dep_tmp = fortran_compiler.find_library('slicot64_pic', required : false)
    if not slicot_dep_tmp.found()
      slicot_dep_tmp = c_compiler.find_library('slicot64_pic')
    endif
    slicot_dep = declare_dependency(dependencies : [ slicot_dep_tmp, blas_dep, lapack_dep ])
  else
    slicot_dep = declare_dependency(dependencies : [ fortran_compiler.find_library('slicot64_pic'), blas_dep, lapack_dep ])
  endif
else # Octave build
  octave_exe = find_program('octave', required : not meson.is_cross_build(), disabler : true)
  mkoctfile_exe = find_program('mkoctfile')
  octave_minimal_version = '8.4.0'
  octave_version = run_command(mkoctfile_exe, '-v', check : true).stdout().replace('mkoctfile, version ', '').strip()

  if octave_version.version_compare('<' + octave_minimal_version)
    error('Octave is too old (version ' + octave_version + '), please upgrade to version ' + octave_minimal_version + ' at least.')
  endif

  octave_incflags = run_command(mkoctfile_exe, '-p', 'INCFLAGS', check : true).stdout().split()
  octlibdir = run_command(mkoctfile_exe, '-p', 'OCTLIBDIR', check : true).stdout().strip()
  if octave_version.version_compare('<10')
    octave_libs = run_command(mkoctfile_exe, '-p', 'OCTAVE_LIBS', check : true).stdout().split()
  else
    octave_libs = run_command(mkoctfile_exe, '-p', 'LIBOCTMEX', check : true).stdout().split()
  endif

  octave_link_args = []

  # Under Linux, ensure that MEX files use the symbols defined within them,
  # preferably to symbols exported by Octave or by other MEX.
  # This matters for the two MS-SBVAR MEX: they both define a main() symbol, and
  # without the following linker option, the first of these two MEX that is loaded
  # will override the main() symbol of the other one.
  # With MATLAB, the same effect is achieved by restricting the exported symbols.
  # This is not needed under Windows and macOS, which apparently have the correct
  # behaviour by default.
  # See also https://savannah.gnu.org/bugs/?55363 which is related
  if host_machine.system() == 'linux'
    octave_link_args += '-Wl,-Bsymbolic'
  endif

  # Determine whether to link MEX files against the Octave libraries. Octave <
  # 10 is not doing this by default but in practice it is needed for Windows
  # and macOS (under Windows, --enable-link-all-dependencies is hardcoded in
  # src/mkoctfile.cc.in; under macOS, the Homebrew formula passes
  # --enable-link-all-dependencies to the configure script).
  if octave_version.version_compare('>=10') or (host_machine.system() == 'windows' or host_machine.system() == 'darwin')
    octave_link_args += [ '-L' + octlibdir ] + octave_libs
  endif

  # For unit tests
  exe_rpath = octlibdir
  exe_link_args = [ '-L' + octlibdir ] + octave_libs

  octave_defs = [ '-DOCTAVE_MEX_FILE', '-DMEXEXT=".mex"', '-DMX_HAS_INTERLEAVED_COMPLEX=1' ]
  if octave_version.version_compare('>=10')
    octave_defs += '-DOCTAVE_MEX_SOVERSION=1'
  endif

  mex_kwargs = { 'name_prefix' : '',
                 'name_suffix' : 'mex',
                 'include_directories' : [ mex_incdir ],
                 'cpp_args' : octave_incflags + octave_defs,
                 'fortran_args' : octave_incflags + octave_defs,
                 'c_args' : octave_incflags + octave_defs,
                 'link_args' : octave_link_args,
                 'install' : true,
                 'install_dir' : 'lib/dynare/mex/octave',
                 'sources' : [ 'mex/sources/octave_interleaved.c' ]}

  if octave_version.version_compare('>=10')
    mex_kwargs += { 'sources' : mex_kwargs['sources'] + [ 'mex/sources/octave_mex_soversion.c' ] }
  endif

  # The -L argument is useful when cross-compiling.
  blas_dep = declare_dependency(link_args : [ '-L' + (octlibdir / '../..') ] + run_command(mkoctfile_exe, '-p', 'BLAS_LIBS', check : true).stdout().split())
  lapack_dep = declare_dependency(link_args : run_command(mkoctfile_exe, '-p', 'LAPACK_LIBS', check : true).stdout().split(),
                                  dependencies : blas_dep)

  # Do not enforce static linking even if prefer_static is true, since that library is shipped
  # with Octave.
  umfpack_dep = dependency('UMFPACK', static : false)

  # This library does not exist under Octave
  ut_dep = []

  # Extra Octave libraries (needed in particular when using the OCTAVE_QUIT macro)
  if octave_version.version_compare('>=10')
    extra_octave_dep = declare_dependency(link_args : [ '-L' + octlibdir ] + run_command(mkoctfile_exe, '-p', 'OCTAVE_LIBS', check : true).stdout().split())
  else
    extra_octave_dep = []
  endif

  # First look for libslicot, then if needed fallback on libslicot_pic
  slicot_dep_tmp = fortran_compiler.find_library('slicot', required : false)
  if not slicot_dep_tmp.found()
    slicot_dep_tmp = fortran_compiler.find_library('slicot_pic')
  endif
  slicot_dep = declare_dependency(dependencies : [ slicot_dep_tmp, blas_dep, lapack_dep ])
endif

# Fortran-specific variables, that include interface modules
mex_fortran_dep = declare_dependency(sources : [ 'mex/sources/matlab_mex.F08' ])
blas_fortran_dep = declare_dependency(sources : [ 'mex/sources/blas.F08' ], dependencies : blas_dep)
lapack_fortran_dep = declare_dependency(sources : [ 'mex/sources/blas.F08', 'mex/sources/lapack.F08' ], dependencies : lapack_dep) # blas.F08 must be included for the integer types

# When static linking is preferred, try to statically link the compiler libraries
if get_option('prefer_static')
  static_link_args = [ '-static-libgcc', '-static-libstdc++' ]

  # Flags to statically link libgfortran and libquadmath
  # Under Debian 13, libgfortran.a and libquadmath.a are not compiled with -fPIC, so can’t be
  # linked into a MEX, and we do not try to statically link those.
  if host_machine.system() == 'windows' or host_machine.system() == 'darwin'
    static_link_args += [ '-static-libgfortran', '-static-libquadmath' ]
  endif

  # Hacks for statically linking libgfortran and libquadmath in MEX files that are a mix of
  # Fortran and another language (C++ or C). Those are needed because for
  # such MEX, meson uses gcc or g++ as linker and manually adds -lgfortran.
  mixed_static_link_args = []
  if host_machine.system() == 'windows'
    # Under Windows, the -static-libgfortran and -static-libquadmath options do not have any
    # effect on g++, hence enforce static linking of them.
    mixed_static_link_args += [ '-Wl,-Bstatic', '-lgfortran', '-lquadmath', '-Wl,-Bdynamic' ]
  endif
  if host_machine.system() == 'darwin'
    # Under macOS, the -static-libgfortran and -static-libquadmath options work as expected when
    # g++ is used as a linker.
    # But meson does not pass -lquadmath to g++, hence leading to undefined symbols (only happens
    # when statically linking, because in dynamic context libquadmath is loaded recursively after
    # libgfortran).
    # Fixes https://forum.dynare.org/t/invalid-mex-file-whe-running-sensitivity-analysis-m1-mac-6-1/25829
    mixed_static_link_args += [ '-lquadmath' ]
  endif

  if get_option('build_for') == 'octave'
    # Under Octave, we add a .c source file in all MEX files, so any MEX with some Fortran is
    # mixed language.
    # Ideally we would add these flags only for MEX with some Fortran, but doing it for all MEX
    # does no harm (as long as b_asneeded remains to its default value, which is true).
    static_link_args += mixed_static_link_args
  endif

  if host_machine.system() == 'windows'
    # Hack to avoid dynamically linking against libwinpthread DLL (which is
    # pulled in by libstdc++, even without using threads, since we are using
    # the POSIX threads version of MinGW).
    static_link_args += [ '-Wl,-Bstatic,--whole-archive', '-lwinpthread', '-Wl,-Bdynamic,--no-whole-archive' ]
  endif

  mex_kwargs += { 'link_args' : static_link_args + mex_kwargs.get('link_args', []) }

  # NB: constructing a dependency object with link_args : ['-Wl,-Bstatic', '-lgomp', '-Wl,-Bdynamic'] does not work,
  # because it reorders the three arguments and puts -lgomp at the end
  if host_machine.system() != 'linux'
    # Under Debian 13, trying to link (static) libgomp.a in a MEX fails.
    openmp_dep_tmp = cpp_compiler.find_library('gomp', static : true)
    openmp_dep = declare_dependency(dependencies : [ openmp_dep, openmp_dep_tmp ])
  endif
endif

# For use when creating intermediate static libraries to be incorporated in MEX files
static_library_kwargs = mex_kwargs + { 'name_prefix' : [], 'name_suffix' : [], 'link_args' : [], 'pic' : true, 'install' : false }


## Various core MEX

shared_module('mjdgges', [ 'mex/sources/mjdgges/mjdgges.f08' ], kwargs : mex_kwargs, dependencies : [ mex_fortran_dep, lapack_fortran_dep ])
shared_module('num_procs', 'mex/sources/num_procs/num_procs.cc', kwargs : mex_kwargs)

perfect_foresight_problem_src = [ 'mex/sources/perfect_foresight_problem/perfect_foresight_problem.cc',
                                  'mex/sources/perfect_foresight_problem/DynamicModelCaller.cc' ]
shared_module('perfect_foresight_problem', perfect_foresight_problem_src, kwargs : mex_kwargs, dependencies : [ openmp_dep, dl_dep ])

perfect_foresight_block_problem_src = [ 'mex/sources/perfect_foresight_problem/perfect_foresight_block_problem.cc',
                                        'mex/sources/perfect_foresight_problem/DynamicModelCaller.cc' ]
shared_module('perfect_foresight_block_problem', perfect_foresight_block_problem_src, kwargs : mex_kwargs, dependencies : [ openmp_dep, dl_dep ])

block_trust_region_src = [ 'mex/sources/block_trust_region/dulmage_mendelsohn.f08',
                           'mex/sources/block_trust_region/matlab_fcn_closure.f08',
                           'mex/sources/block_trust_region/trust_region.f08',
                           'mex/sources/block_trust_region/mexFunction.f08' ]
shared_module('block_trust_region', block_trust_region_src, kwargs : mex_kwargs, dependencies : [ mex_fortran_dep, lapack_fortran_dep ])

trust_region_test_exe = executable('trust_region_test',
                                   [ 'mex/sources/block_trust_region/trust_region.f08',
                                     'mex/sources/block_trust_region/test/trust_region_test.f08' ],
                                   dependencies : lapack_fortran_dep,
                                   link_args : exe_link_args,
                                   build_rpath : exe_rpath,
                                   fortran_args : mex_kwargs['fortran_args'])
test('trust_region', trust_region_test_exe, timeout : 0)

bytecode_src = [ 'mex/sources/bytecode/bytecode.cc',
                 'mex/sources/bytecode/Interpreter.cc',
                 'mex/sources/bytecode/Mem_Mngr.cc',
                 'mex/sources/bytecode/Evaluate.cc',
                 'mex/sources/bytecode/BasicSymbolTable.cc' ]
preprocessor_headers_dep = declare_dependency(include_directories : include_directories('preprocessor/src'))
shared_module('bytecode', bytecode_src, kwargs : mex_kwargs, dependencies : [ umfpack_dep, ut_dep, extra_octave_dep, preprocessor_headers_dep ])

shared_module('sparse_hessian_times_B_kronecker_C', 'mex/sources/kronecker/sparse_hessian_times_B_kronecker_C.cc',
              kwargs : mex_kwargs, dependencies : openmp_dep)

shared_module('A_times_B_kronecker_C', [ 'mex/sources/kronecker/A_times_B_kronecker_C.f08' ],
              kwargs : mex_kwargs, dependencies : [ mex_fortran_dep, blas_fortran_dep ])

shared_module('cycle_reduction', [ 'mex/sources/cycle_reduction/mexFunction.f08' ],
              kwargs : mex_kwargs, dependencies : [ mex_fortran_dep, blas_fortran_dep, lapack_fortran_dep ])

shared_module('logarithmic_reduction', [ 'mex/sources/logarithmic_reduction/mexFunction.f08' ],
              kwargs : mex_kwargs, dependencies : [ mex_fortran_dep, blas_fortran_dep, lapack_fortran_dep ])

shared_module('disclyap_fast', [ 'mex/sources/disclyap_fast/disclyap_fast.f08' ],
              kwargs : mex_kwargs, dependencies : [ mex_fortran_dep, blas_fortran_dep, lapack_fortran_dep ])

qmc_sequence_src = [ 'mex/sources/sobol/qmc_sequence.cc',
                     'mex/sources/sobol/sobol.f08' ]
qmc_sequence_mex_kwargs = mex_kwargs
if get_option('prefer_static')
  # This MEX has both Fortran and C++, so use the additional static linking args for mixed languages
  qmc_sequence_mex_kwargs += { 'link_args' : qmc_sequence_mex_kwargs.get('link_args') + mixed_static_link_args }
endif
# Workaround for cross-compiler bug in Debian 13
# With LTO, one gets a compilation failure with messages such as:
#  /usr/bin/x86_64-w64-mingw32-ld: cannot export _Z10icdfmSigmaIdEviiPT_PKd: symbol wrong type (4 vs 3)
if get_option('build_for') == 'octave' and host_machine.system() == 'windows'
  qmc_sequence_mex_kwargs += { 'override_options' : [ 'b_lto=false' ] }
endif
shared_module('qmc_sequence', qmc_sequence_src, kwargs : qmc_sequence_mex_kwargs, dependencies : [ blas_dep, openmp_dep ])

shared_module('kalman_steady_state', 'mex/sources/kalman_steady_state/kalman_steady_state.cc', kwargs : mex_kwargs, dependencies : slicot_dep)

shared_module('kalman_filter_mex', [ 'mex/sources/kalman_filter/kalman_filter.f08' ],
              kwargs : mex_kwargs, dependencies : [ mex_fortran_dep, blas_fortran_dep, lapack_fortran_dep ])

## k-order simulation stuff

kordersim_src = [ 'mex/sources/libkordersim/pascal.f08',
                  'mex/sources/libkordersim/sort.f08',
                  'mex/sources/libkordersim/partitions.f08',
                  'mex/sources/libkordersim/tensors.f08',
                  'mex/sources/libkordersim/simulation.f08',
                  'mex/sources/libkordersim/struct.f08' ]
kordersim_lib = static_library('kordersim', kordersim_src, kwargs : static_library_kwargs, dependencies : [ mex_fortran_dep, blas_fortran_dep, pthread_fortran_dep ])

shared_module('folded_to_unfolded_dr', 'mex/sources/folded_to_unfolded_dr/mexFunction.f08', kwargs : mex_kwargs, link_with : kordersim_lib)
shared_module('k_order_mean', 'mex/sources/k_order_mean/mexFunction.f08', kwargs : mex_kwargs, link_with : kordersim_lib)
shared_module('k_order_simul', 'mex/sources/k_order_simul/mexFunction.f08', kwargs : mex_kwargs, link_with : kordersim_lib)

shared_module('local_state_space_iteration_2', 'mex/sources/local_state_space_iterations/local_state_space_iteration_2.cc', kwargs : mex_kwargs, dependencies : openmp_dep)
shared_module('local_state_space_iteration_3', 'mex/sources/local_state_space_iterations/local_state_space_iteration_3.f08', kwargs : mex_kwargs, link_with : kordersim_lib)
shared_module('local_state_space_iteration_k', 'mex/sources/local_state_space_iterations/local_state_space_iteration_k.f08', kwargs : mex_kwargs, link_with : kordersim_lib)


## k-order resolution stuff

korder_src = [ 'mex/sources/libkorder/kord/approximation.cc',
               'mex/sources/libkorder/kord/decision_rule.cc',
               'mex/sources/libkorder/kord/dynamic_model.cc',
               'mex/sources/libkorder/kord/faa_di_bruno.cc',
               'mex/sources/libkorder/kord/first_order.cc',
               'mex/sources/libkorder/kord/korder.cc',
               'mex/sources/libkorder/kord/korder_stoch.cc',
               'mex/sources/libkorder/kord/journal.cc',
               'mex/sources/libkorder/sylv/BlockDiagonal.cc',
               'mex/sources/libkorder/sylv/GeneralMatrix.cc',
               'mex/sources/libkorder/sylv/GeneralSylvester.cc',
               'mex/sources/libkorder/sylv/IterativeSylvester.cc',
               'mex/sources/libkorder/sylv/KronUtils.cc',
               'mex/sources/libkorder/sylv/KronVector.cc',
               'mex/sources/libkorder/sylv/QuasiTriangular.cc',
               'mex/sources/libkorder/sylv/QuasiTriangularZero.cc',
               'mex/sources/libkorder/sylv/SchurDecomp.cc',
               'mex/sources/libkorder/sylv/SchurDecompEig.cc',
               'mex/sources/libkorder/sylv/SimilarityDecomp.cc',
               'mex/sources/libkorder/sylv/SylvException.cc',
               'mex/sources/libkorder/sylv/SylvMatrix.cc',
               'mex/sources/libkorder/sylv/SylvParams.cc',
               'mex/sources/libkorder/sylv/SymSchurDecomp.cc',
               'mex/sources/libkorder/sylv/TriangularSylvester.cc',
               'mex/sources/libkorder/sylv/Vector.cc',
               'mex/sources/libkorder/tl/equivalence.cc',
               'mex/sources/libkorder/tl/fine_container.cc',
               'mex/sources/libkorder/tl/fs_tensor.cc',
               'mex/sources/libkorder/tl/gs_tensor.cc',
               'mex/sources/libkorder/tl/int_sequence.cc',
               'mex/sources/libkorder/tl/kron_prod.cc',
               'mex/sources/libkorder/tl/normal_moments.cc',
               'mex/sources/libkorder/tl/permutation.cc',
               'mex/sources/libkorder/tl/ps_tensor.cc',
               'mex/sources/libkorder/tl/pyramid_prod.cc',
               'mex/sources/libkorder/tl/pyramid_prod2.cc',
               'mex/sources/libkorder/tl/rfs_tensor.cc',
               'mex/sources/libkorder/tl/sparse_tensor.cc',
               'mex/sources/libkorder/tl/stack_container.cc',
               'mex/sources/libkorder/tl/symmetry.cc',
               'mex/sources/libkorder/tl/t_container.cc',
               'mex/sources/libkorder/tl/t_polynomial.cc',
               'mex/sources/libkorder/tl/tensor.cc',
               'mex/sources/libkorder/tl/tl_static.cc',
               'mex/sources/libkorder/tl/twod_matrix.cc',
               'mex/sources/libkorder/utils/pascal_triangle.cc',
               'mex/sources/libkorder/utils/int_power.cc',
               'mex/sources/libkorder/utils/thread_pool.cc',
               'mex/sources/libkorder/k_ord_dynare.cc',
               'mex/sources/libkorder/dynamic_dll.cc',
               'mex/sources/libkorder/dynamic_m.cc' ]
korder_incdir = include_directories('mex/sources/libkorder', 'mex/sources/libkorder/tl', 'mex/sources/libkorder/sylv',
                                    'mex/sources/libkorder/kord', 'mex/sources/libkorder/utils')
korder_lib = static_library('korder', korder_src,
                            kwargs : static_library_kwargs + { 'include_directories' : static_library_kwargs.get('include_directories') + korder_incdir},
                            dependencies : [ blas_dep, lapack_dep, dl_dep ])

korder_mex_kwargs = mex_kwargs + { 'include_directories' : mex_kwargs.get('include_directories') + korder_incdir}
shared_module('gensylv', 'mex/sources/gensylv/gensylv.cc', kwargs : korder_mex_kwargs, link_with : korder_lib)
shared_module('k_order_perturbation', 'mex/sources/k_order_perturbation/k_order_perturbation.cc', kwargs : korder_mex_kwargs, link_with : korder_lib)
k_order_welfare_src = [ 'mex/sources/k_order_welfare/k_order_welfare.cc',
                        'mex/sources/k_order_welfare/approximation_welfare.cc',
                        'mex/sources/k_order_welfare/k_ord_objective.cc',
                        'mex/sources/k_order_welfare/objective_m.cc' ]
shared_module('k_order_welfare', k_order_welfare_src, kwargs : korder_mex_kwargs, link_with : korder_lib)

# Unit tests

k_order_test_kwargs = { 'include_directories' : mex_kwargs['include_directories'] + korder_incdir,
                        'cpp_args' : mex_kwargs['cpp_args'],
                        'link_args' : exe_link_args,
                        'build_rpath' : exe_rpath,
                        'link_with' : korder_lib }

k_order_sylv_test_exe = executable('k_order_sylv_test', [ 'mex/sources/libkorder/sylv/tests/MMMatrix.cc',
                                                          'mex/sources/libkorder/sylv/tests/tests.cc'],
                                   kwargs : k_order_test_kwargs)
test('k_order_sylv', k_order_sylv_test_exe, workdir : meson.current_source_dir() / 'mex/sources/libkorder/sylv/tests', timeout : 0, suite : 'k_order')

k_order_tl_test_exe = executable('k_order_tl_test', [ 'mex/sources/libkorder/tl/tests/factory.cc',
                                                      'mex/sources/libkorder/tl/tests/monoms.cc',
                                                      'mex/sources/libkorder/tl/tests/tests.cc'],
                                 kwargs : k_order_test_kwargs)
test('k_order_tl', k_order_tl_test_exe, timeout : 0, suite : 'k_order')

k_order_kord_test_exe = executable('k_order_kord_test', 'mex/sources/libkorder/kord/tests/tests.cc', kwargs : k_order_test_kwargs)
test('k_order_kord', k_order_kord_test_exe, timeout : 0, suite : 'k_order')


## MS-SBVAR stuff

ms_sbvar_src = [ 'contrib/ms-sbvar/utilities_dw/arrays/dw_array.c',
                 'contrib/ms-sbvar/utilities_dw/arrays/dw_matrix_array.c',
                 'contrib/ms-sbvar/utilities_dw/ascii/dw_ascii.c',
                 'contrib/ms-sbvar/utilities_dw/ascii/dw_parse_cmd.c',
                 'contrib/ms-sbvar/utilities_dw/elliptical/dw_elliptical.c',
                 'contrib/ms-sbvar/utilities_dw/error/dw_error.c',
                 'contrib/ms-sbvar/utilities_dw/histogram/dw_histogram.c',
                 'contrib/ms-sbvar/utilities_dw/math/dw_math.c',
                 'contrib/ms-sbvar/utilities_dw/matrix/dw_matrix.c',
                 'contrib/ms-sbvar/utilities_dw/matrix/bmatrix.c',
                 'contrib/ms-sbvar/utilities_dw/sort/dw_matrix_sort.c',
                 'contrib/ms-sbvar/utilities_dw/stat/dw_rand.cc',
                 'contrib/ms-sbvar/utilities_dw/stat/dw_matrix_rand.c',
                 'contrib/ms-sbvar/switch_dw/switching/dw_switch.c',
                 'contrib/ms-sbvar/switch_dw/switching/dw_switchio.c',
                 'contrib/ms-sbvar/switch_dw/switching/dw_dirichlet_restrictions.c',
                 'contrib/ms-sbvar/switch_dw/switching/dw_metropolis_theta.c',
                 'contrib/ms-sbvar/switch_dw/state_space/sbvar/VARbase.c',
                 'contrib/ms-sbvar/switch_dw/state_space/sbvar/VARio.c',
                 'contrib/ms-sbvar/switch_dw/state_space/sbvar/write-mat.cc',
                 'mex/sources/ms-sbvar/modify_for_mex.cc' ]
ms_sbvar_defs = [ '-DSTRUCTURED_COLUMN_MAJOR' ]
ms_sbvar_incdir = include_directories('contrib/ms-sbvar/utilities_dw/include', 'contrib/ms-sbvar/switch_dw/switching', 'mex/sources/ms-sbvar')
ms_sbvar_lib = static_library('ms_sbvar', ms_sbvar_src,
                              kwargs : static_library_kwargs + { 'c_args' : static_library_kwargs.get('c_args') + ms_sbvar_defs,
                                                                 'cpp_args' : static_library_kwargs.get('cpp_args') + ms_sbvar_defs,
                                                                 'include_directories' : static_library_kwargs.get('include_directories') + ms_sbvar_incdir },
                              dependencies : [ blas_dep, lapack_dep, boost_dep, ut_dep, extra_octave_dep ])

# NB: in the following, mex_top_level.cc is put as extra source and not in ms_sbvar_lib,
# to ensure that mexFunction() is present in the MEX
mex_ms_sbvar_kwargs = mex_kwargs + { 'c_args' : mex_kwargs.get('c_args') + ms_sbvar_defs,
                                     'cpp_args' : mex_kwargs.get('cpp_args') + ms_sbvar_defs,
                                     'include_directories' : mex_kwargs.get('include_directories') + ms_sbvar_incdir,
                                     'link_with' : ms_sbvar_lib,
                                     'sources' : mex_kwargs['sources'] + [ 'mex/sources/ms-sbvar/mex_top_level.cc' ] }
ms_sbvar_create_init_file_src = [ 'contrib/ms-sbvar/switch_dw/state_space/sbvar/create_init_file.c',
                                  'contrib/ms-sbvar/switch_dw/state_space/sbvar/VARio_matlab.c' ]
shared_module('ms_sbvar_create_init_file', ms_sbvar_create_init_file_src, kwargs : mex_ms_sbvar_kwargs)
ms_sbvar_command_line_src = [ 'contrib/ms-sbvar/switch_dw/switching/dw_switch_opt.c',
                              'contrib/ms-sbvar/switch_dw/switching/dw_mdd_switch.c',
                              'contrib/ms-sbvar/switch_dw/state_space/sbvar/dw_sbvar_command_line.c',
                              'contrib/ms-sbvar/switch_dw/state_space/sbvar/sbvar_estimate.c',
                              'contrib/ms-sbvar/switch_dw/state_space/sbvar/sbvar_simulate.c',
                              'contrib/ms-sbvar/switch_dw/state_space/sbvar/sbvar_probabilities.c',
                              'contrib/ms-sbvar/switch_dw/state_space/sbvar/sbvar_mdd.c',
                              'contrib/ms-sbvar/switch_dw/state_space/sbvar/sbvar_forecast.c',
                              'contrib/ms-sbvar/switch_dw/state_space/sbvar/sbvar_variance_decomposition.c',
                              'contrib/ms-sbvar/switch_dw/state_space/sbvar/sbvar_impulse_responses.c',
                              'contrib/ms-sbvar/switch_dw/state_space/sbvar/dw_csminwel.c' ]
shared_module('ms_sbvar_command_line', ms_sbvar_command_line_src, kwargs : mex_ms_sbvar_kwargs)

heterogeneity_src = [
   'mex/sources/block_trust_region/trust_region.f08',
   'mex/sources/heterogeneity/utils/gray_code.f08',
   'mex/sources/heterogeneity/utils/interpolation.f08',
   'mex/sources/heterogeneity/utils/markov.f08',
   'mex/sources/heterogeneity/utils/expectations.f08',
   'mex/sources/heterogeneity/utils/curlyDs.f08',
   'mex/sources/heterogeneity/utils/aggregate_residuals_tensor.f08',
   'mex/sources/heterogeneity/utils/fischer_burmeister.f08',
   'mex/sources/heterogeneity/utils/broyden_solver.f08',
   'mex/sources/heterogeneity/utils/calibration_types.f08',
   'mex/sources/heterogeneity/utils/calibration.f08',
   'mex/sources/heterogeneity/utils/pack_output.f08',
   'mex/sources/heterogeneity/utils/time_iteration_tensor.f08'
]
heterogeneity_lib = static_library('heterogeneity', heterogeneity_src, kwargs : static_library_kwargs, dependencies : [ mex_fortran_dep, blas_fortran_dep, lapack_fortran_dep ])
shared_module('find_bracket_linear_weight',  'mex/sources/heterogeneity/find_bracket_linear_weight.f08', kwargs : mex_kwargs, link_with : heterogeneity_lib)
shared_module('rouwenhorst',  'mex/sources/heterogeneity/rouwenhorst.f08', kwargs : mex_kwargs, link_with : heterogeneity_lib)
shared_module('compute_expected_dx',  'mex/sources/heterogeneity/compute_expected_dx.f08', kwargs : mex_kwargs, link_with : heterogeneity_lib)
shared_module('compute_expected_y',  'mex/sources/heterogeneity/compute_expected_y.f08', kwargs : mex_kwargs, link_with : heterogeneity_lib)
shared_module('compute_curlyDs',  'mex/sources/heterogeneity/compute_curlyDs.f08', kwargs : mex_kwargs, link_with : heterogeneity_lib)
shared_module('compute_transition_matrices',  'mex/sources/heterogeneity/compute_transition_matrices.f08', kwargs : mex_kwargs, link_with : heterogeneity_lib)
shared_module('compute_steady_state_tensor',  'mex/sources/heterogeneity/compute_steady_state_tensor.f08', kwargs : mex_kwargs, link_with : heterogeneity_lib)
shared_module('compute_Phi_tilde_e',  'mex/sources/heterogeneity/compute_Phi_tilde_e.f08', kwargs : mex_kwargs, link_with : heterogeneity_lib)


## Optionally build PARDISO MEX files

pardiso_dep = fortran_compiler.find_library('pardiso', required : get_option('pardiso'), disabler : true)
shared_module('pardiso_init', [ 'mex/sources/pardiso/pardiso_init.f08', 'mex/sources/pardiso/pardiso.f08' ], kwargs : mex_kwargs, dependencies : [ mex_fortran_dep, pardiso_dep ])
shared_module('pardiso_free', [ 'mex/sources/pardiso/pardiso_free.f08', 'mex/sources/pardiso/pardiso.f08' ], kwargs : mex_kwargs, dependencies : [ mex_fortran_dep, pardiso_dep ])
shared_module('pardiso_solve', [ 'mex/sources/pardiso/pardiso_solve.f08', 'mex/sources/pardiso/pardiso.f08' ], kwargs : mex_kwargs, dependencies : [ mex_fortran_dep, pardiso_dep ])


## Optionally build ParU MEX
# Build instructions inspired from ParU/MATLAB/paru_make.m in the SuiteSparse distribution

if get_option('suitesparse_src_path') != ''
  paru_inc = [
    'include',
    'ParU/Include',
    'SuiteSparse_config',
    'AMD/Source',
    'AMD/Include',
    'CAMD/Source',
    'CAMD/Include',
    'COLAMD/Source',
    'COLAMD/Include',
    'CCOLAMD/Source',
    'CCOLAMD/Include',
    'UMFPACK/Source',
    'UMFPACK/Include',
    'CHOLMOD/Include',
    'CHOLMOD/Check',
    'CHOLMOD/Cholesky',
    'CHOLMOD/Utility',
    'CHOLMOD/MATLAB',
    'CHOLMOD/MatrixOps',
    'CHOLMOD/Partition',
    'CHOLMOD',
    'CHOLMOD/SuiteSparse_metis/GKlib',
    'CHOLMOD/SuiteSparse_metis/libmetis',
    'CHOLMOD/SuiteSparse_metis/include',
    'CHOLMOD/Supernodal'
  ]

  paru_src = [
    'SuiteSparse_config/SuiteSparse_config.c',
    'AMD/Source/amd_l1.c',
    'AMD/Source/amd_l2.c',
    'AMD/Source/amd_l_aat.c',
    'AMD/Source/amd_l_control.c',
    'AMD/Source/amd_l_defaults.c',
    'AMD/Source/amd_l_dump.c',
    'AMD/Source/amd_l_info.c',
    'AMD/Source/amd_l_order.c',
    'AMD/Source/amd_l_postorder.c',
    'AMD/Source/amd_l_post_tree.c',
    'AMD/Source/amd_l_preprocess.c',
    'AMD/Source/amd_l_valid.c',
    'CAMD/Source/camd_l1.c',
    'CAMD/Source/camd_l2.c',
    'CAMD/Source/camd_l_aat.c',
    'CAMD/Source/camd_l_control.c',
    'CAMD/Source/camd_l_defaults.c',
    'CAMD/Source/camd_l_dump.c',
    'CAMD/Source/camd_l_info.c',
    'CAMD/Source/camd_l_order.c',
    'CAMD/Source/camd_l_postorder.c',
    'CAMD/Source/camd_l_preprocess.c',
    'CAMD/Source/camd_l_valid.c',
    'COLAMD/Source/colamd_l.c',
    'CCOLAMD/Source/ccolamd_l.c',
    'CHOLMOD/MATLAB/sputil2.c',
    'CHOLMOD/Utility/cholmod_l_aat.c',
    'CHOLMOD/Utility/cholmod_l_add.c',
    'CHOLMOD/Utility/cholmod_l_add_size_t.c',
    'CHOLMOD/Utility/cholmod_l_allocate_dense.c',
    'CHOLMOD/Utility/cholmod_l_allocate_factor.c',
    'CHOLMOD/Utility/cholmod_l_allocate_sparse.c',
    'CHOLMOD/Utility/cholmod_l_allocate_triplet.c',
    'CHOLMOD/Utility/cholmod_l_allocate_work.c',
    'CHOLMOD/Utility/cholmod_l_alloc_factor.c',
    'CHOLMOD/Utility/cholmod_l_alloc_work.c',
    'CHOLMOD/Utility/cholmod_l_band.c',
    'CHOLMOD/Utility/cholmod_l_band_nnz.c',
    'CHOLMOD/Utility/cholmod_l_calloc.c',
    'CHOLMOD/Utility/cholmod_l_change_factor.c',
    'CHOLMOD/Utility/cholmod_l_clear_flag.c',
    'CHOLMOD/Utility/cholmod_l_copy.c',
    'CHOLMOD/Utility/cholmod_l_copy_dense2.c',
    'CHOLMOD/Utility/cholmod_l_copy_dense.c',
    'CHOLMOD/Utility/cholmod_l_copy_factor.c',
    'CHOLMOD/Utility/cholmod_l_copy_sparse.c',
    'CHOLMOD/Utility/cholmod_l_copy_triplet.c',
    'CHOLMOD/Utility/cholmod_l_cumsum.c',
    'CHOLMOD/Utility/cholmod_l_dbound.c',
    'CHOLMOD/Utility/cholmod_l_defaults.c',
    'CHOLMOD/Utility/cholmod_l_dense_nnz.c',
    'CHOLMOD/Utility/cholmod_l_dense_to_sparse.c',
    'CHOLMOD/Utility/cholmod_l_divcomplex.c',
    'CHOLMOD/Utility/cholmod_l_ensure_dense.c',
    'CHOLMOD/Utility/cholmod_l_error.c',
    'CHOLMOD/Utility/cholmod_l_eye.c',
    'CHOLMOD/Utility/cholmod_l_factor_to_sparse.c',
    'CHOLMOD/Utility/cholmod_l_finish.c',
    'CHOLMOD/Utility/cholmod_l_free.c',
    'CHOLMOD/Utility/cholmod_l_free_dense.c',
    'CHOLMOD/Utility/cholmod_l_free_factor.c',
    'CHOLMOD/Utility/cholmod_l_free_sparse.c',
    'CHOLMOD/Utility/cholmod_l_free_triplet.c',
    'CHOLMOD/Utility/cholmod_l_free_work.c',
    'CHOLMOD/Utility/cholmod_l_hypot.c',
    'CHOLMOD/Utility/cholmod_l_malloc.c',
    'CHOLMOD/Utility/cholmod_l_maxrank.c',
    'CHOLMOD/Utility/cholmod_l_mult_size_t.c',
    'CHOLMOD/Utility/cholmod_l_nnz.c',
    'CHOLMOD/Utility/cholmod_l_ones.c',
    'CHOLMOD/Utility/cholmod_l_pack_factor.c',
    'CHOLMOD/Utility/cholmod_l_ptranspose.c',
    'CHOLMOD/Utility/cholmod_l_reallocate_column.c',
    'CHOLMOD/Utility/cholmod_l_reallocate_factor.c',
    'CHOLMOD/Utility/cholmod_l_reallocate_sparse.c',
    'CHOLMOD/Utility/cholmod_l_reallocate_triplet.c',
    'CHOLMOD/Utility/cholmod_l_realloc.c',
    'CHOLMOD/Utility/cholmod_l_realloc_multiple.c',
    'CHOLMOD/Utility/cholmod_l_sbound.c',
    'CHOLMOD/Utility/cholmod_l_score_comp.c',
    'CHOLMOD/Utility/cholmod_l_set_empty.c',
    'CHOLMOD/Utility/cholmod_l_sort.c',
    'CHOLMOD/Utility/cholmod_l_sparse_to_dense.c',
    'CHOLMOD/Utility/cholmod_l_sparse_to_triplet.c',
    'CHOLMOD/Utility/cholmod_l_speye.c',
    'CHOLMOD/Utility/cholmod_l_spzeros.c',
    'CHOLMOD/Utility/cholmod_l_start.c',
    'CHOLMOD/Utility/cholmod_l_transpose.c',
    'CHOLMOD/Utility/cholmod_l_transpose_sym.c',
    'CHOLMOD/Utility/cholmod_l_transpose_unsym.c',
    'CHOLMOD/Utility/cholmod_l_triplet_to_sparse.c',
    'CHOLMOD/Utility/cholmod_l_version.c',
    'CHOLMOD/Utility/cholmod_l_xtype.c',
    'CHOLMOD/Utility/cholmod_l_zeros.c',
    'CHOLMOD/Utility/cholmod_mult_uint64_t.c',
    'CHOLMOD/Utility/cholmod_memdebug.c',
    'CHOLMOD/Check/cholmod_l_check.c',
    'CHOLMOD/Check/cholmod_l_read.c',
    'CHOLMOD/Check/cholmod_l_write.c',
    'CHOLMOD/Cholesky/cholmod_l_amd.c',
    'CHOLMOD/Cholesky/cholmod_l_analyze.c',
    'CHOLMOD/Cholesky/cholmod_l_colamd.c',
    'CHOLMOD/Cholesky/cholmod_l_etree.c',
    'CHOLMOD/Cholesky/cholmod_l_factorize.c',
    'CHOLMOD/Cholesky/cholmod_l_postorder.c',
    'CHOLMOD/Cholesky/cholmod_l_rcond.c',
    'CHOLMOD/Cholesky/cholmod_l_resymbol.c',
    'CHOLMOD/Cholesky/cholmod_l_rowcolcounts.c',
    'CHOLMOD/Cholesky/cholmod_l_rowfac.c',
    'CHOLMOD/Cholesky/cholmod_l_solve.c',
    'CHOLMOD/Cholesky/cholmod_l_spsolve.c',
    'CHOLMOD/MatrixOps/cholmod_l_drop.c',
    'CHOLMOD/MatrixOps/cholmod_l_horzcat.c',
    'CHOLMOD/MatrixOps/cholmod_l_norm.c',
    'CHOLMOD/MatrixOps/cholmod_l_scale.c',
    'CHOLMOD/MatrixOps/cholmod_l_sdmult.c',
    'CHOLMOD/MatrixOps/cholmod_l_ssmult.c',
    'CHOLMOD/MatrixOps/cholmod_l_submatrix.c',
    'CHOLMOD/MatrixOps/cholmod_l_vertcat.c',
    'CHOLMOD/MatrixOps/cholmod_l_symmetry.c',
    'CHOLMOD/Supernodal/cholmod_l_super_numeric.c',
    'CHOLMOD/Supernodal/cholmod_l_super_solve.c',
    'CHOLMOD/Supernodal/cholmod_l_super_symbolic.c',
    'CHOLMOD/Partition/cholmod_metis_wrapper.c',
    'CHOLMOD/Partition/cholmod_l_ccolamd.c',
    'CHOLMOD/Partition/cholmod_l_csymamd.c',
    'CHOLMOD/Partition/cholmod_l_camd.c',
    'CHOLMOD/Partition/cholmod_l_metis.c',
    'CHOLMOD/Partition/cholmod_l_nesdis.c',
    'UMFPACK/Source2/umf_l_analyze.c',
    'UMFPACK/Source2/umf_l_apply_order.c',
    'UMFPACK/Source2/umf_l_cholmod.c',
    'UMFPACK/Source2/umf_l_colamd.c',
    'UMFPACK/Source2/umf_l_free.c',
    'UMFPACK/Source2/umf_l_fsize.c',
    'UMFPACK/Source2/umf_l_is_permutation.c',
    'UMFPACK/Source2/umf_l_malloc.c',
    'UMFPACK/Source2/umf_l_realloc.c',
    'UMFPACK/Source2/umf_l_report_perm.c',
    'UMFPACK/Source2/umf_l_singletons.c',
    'UMFPACK/Source2/umf_dl_assemble.c',
    'UMFPACK/Source2/umf_dl_assemble_fixq.c',
    'UMFPACK/Source2/umf_dl_blas3_update.c',
    'UMFPACK/Source2/umf_dl_build_tuples.c',
    'UMFPACK/Source2/umf_dl_create_element.c',
    'UMFPACK/Source2/umf_dl_dump.c',
    'UMFPACK/Source2/umf_dl_extend_front.c',
    'UMFPACK/Source2/umf_dl_garbage_collection.c',
    'UMFPACK/Source2/umf_dl_get_memory.c',
    'UMFPACK/Source2/umf_dl_grow_front.c',
    'UMFPACK/Source2/umf_dl_init_front.c',
    'UMFPACK/Source2/umf_dl_kernel.c',
    'UMFPACK/Source2/umf_dl_kernel_init.c',
    'UMFPACK/Source2/umf_dl_kernel_wrapup.c',
    'UMFPACK/Source2/umf_dl_lhsolve.c',
    'UMFPACK/Source2/umf_dl_local_search.c',
    'UMFPACK/Source2/umf_dl_lsolve.c',
    'UMFPACK/Source2/umf_dl_ltsolve.c',
    'UMFPACK/Source2/umf_dl_mem_alloc_element.c',
    'UMFPACK/Source2/umf_dl_mem_alloc_head_block.c',
    'UMFPACK/Source2/umf_dl_mem_alloc_tail_block.c',
    'UMFPACK/Source2/umf_dl_mem_free_tail_block.c',
    'UMFPACK/Source2/umf_dl_mem_init_memoryspace.c',
    'UMFPACK/Source2/umf_dl_report_vector.c',
    'UMFPACK/Source2/umf_dl_row_search.c',
    'UMFPACK/Source2/umf_dl_scale.c',
    'UMFPACK/Source2/umf_dl_scale_column.c',
    'UMFPACK/Source2/umf_dl_set_stats.c',
    'UMFPACK/Source2/umf_dl_solve.c',
    'UMFPACK/Source2/umf_dl_start_front.c',
    'UMFPACK/Source2/umf_dl_store_lu.c',
    'UMFPACK/Source2/umf_dl_store_lu_drop.c',
    'UMFPACK/Source2/umf_dl_symbolic_usage.c',
    'UMFPACK/Source2/umf_dl_transpose.c',
    'UMFPACK/Source2/umf_dl_triplet_map_nox.c',
    'UMFPACK/Source2/umf_dl_triplet_map_x.c',
    'UMFPACK/Source2/umf_dl_triplet_nomap_nox.c',
    'UMFPACK/Source2/umf_dl_triplet_nomap_x.c',
    'UMFPACK/Source2/umf_dl_tuple_lengths.c',
    'UMFPACK/Source2/umf_dl_uhsolve.c',
    'UMFPACK/Source2/umf_dl_usolve.c',
    'UMFPACK/Source2/umf_dl_utsolve.c',
    'UMFPACK/Source2/umf_dl_valid_numeric.c',
    'UMFPACK/Source2/umf_dl_valid_symbolic.c',
    'UMFPACK/Source2/umfpack_gn_tictoc.c',
    'UMFPACK/Source2/umfpack_gn_timer.c',
    'UMFPACK/Source2/umfpack_dl_col_to_triplet.c',
    'UMFPACK/Source2/umfpack_dl_defaults.c',
    'UMFPACK/Source2/umfpack_dl_free_numeric.c',
    'UMFPACK/Source2/umfpack_dl_free_symbolic.c',
    'UMFPACK/Source2/umfpack_dl_get_determinant.c',
    'UMFPACK/Source2/umfpack_dl_get_lunz.c',
    'UMFPACK/Source2/umfpack_dl_get_numeric.c',
    'UMFPACK/Source2/umfpack_dl_get_symbolic.c',
    'UMFPACK/Source2/umfpack_dl_load_numeric.c',
    'UMFPACK/Source2/umfpack_dl_load_symbolic.c',
    'UMFPACK/Source2/umfpack_dl_numeric.c',
    'UMFPACK/Source2/umfpack_dl_qsymbolic.c',
    'UMFPACK/Source2/umfpack_dl_report_control.c',
    'UMFPACK/Source2/umfpack_dl_report_info.c',
    'UMFPACK/Source2/umfpack_dl_report_matrix.c',
    'UMFPACK/Source2/umfpack_dl_report_numeric.c',
    'UMFPACK/Source2/umfpack_dl_report_perm.c',
    'UMFPACK/Source2/umfpack_dl_report_status.c',
    'UMFPACK/Source2/umfpack_dl_report_symbolic.c',
    'UMFPACK/Source2/umfpack_dl_report_triplet.c',
    'UMFPACK/Source2/umfpack_dl_report_vector.c',
    'UMFPACK/Source2/umfpack_dl_save_numeric.c',
    'UMFPACK/Source2/umfpack_dl_save_symbolic.c',
    'UMFPACK/Source2/umfpack_dl_scale.c',
    'UMFPACK/Source2/umfpack_dl_solve.c',
    'UMFPACK/Source2/umfpack_dl_symbolic.c',
    'UMFPACK/Source2/umfpack_dl_transpose.c',
    'UMFPACK/Source2/umfpack_dl_triplet_to_col.c',
    'UMFPACK/Source2/umfpack_dl_wsolve.c',
    'ParU/Source/ParU_Factorize.cpp',
    'ParU/Source/ParU_Analyze.cpp',
    'ParU/Source/paru_assemble.cpp',
    'ParU/Source/paru_assemble_row2U.cpp',
    'ParU/Source/paru_bin_search.cpp',
    'ParU/Source/ParU_C.cpp',
    'ParU/Source/paru_create_element.cpp',
    'ParU/Source/paru_cumsum.cpp',
    'ParU/Source/paru_dgemm.cpp',
    'ParU/Source/paru_diag_update.cpp',
    'ParU/Source/paru_exec_tasks.cpp',
    'ParU/Source/paru_finalize_perm.cpp',
    'ParU/Source/ParU_FreeNumeric.cpp',
    'ParU/Source/ParU_FreeSymbolic.cpp',
    'ParU/Source/ParU_FreeControl.cpp',
    'ParU/Source/ParU_InitControl.cpp',
    'ParU/Source/paru_front.cpp',
    'ParU/Source/paru_fs_factorize.cpp',
    'ParU/Source/paru_full_summed.cpp',
    'ParU/Source/paru_gaxpy.cpp',
    'ParU/Source/ParU_Get.cpp',
    'ParU/Source/paru_hash.cpp',
    'ParU/Source/paru_heap.cpp',
    'ParU/Source/paru_init_rel.cpp',
    'ParU/Source/paru_init_rowFronts.cpp',
    'ParU/Source/paru_intersection.cpp',
    'ParU/Source/ParU_InvPerm.cpp',
    'ParU/Source/ParU_LSolve.cpp',
    'ParU/Source/paru_mem.cpp',
    'ParU/Source/paru_free_work.cpp',
    'ParU/Source/paru_memcpy.cpp',
    'ParU/Source/paru_memset.cpp',
    'ParU/Source/paru_norms.cpp',
    'ParU/Source/paru_nthreads.cpp',
    'ParU/Source/ParU_Perm.cpp',
    'ParU/Source/paru_pivotal.cpp',
    'ParU/Source/paru_prior_assemble.cpp',
    'ParU/Source/ParU_Residual.cpp',
    'ParU/Source/ParU_Set.cpp',
    'ParU/Source/ParU_Solve.cpp',
    'ParU/Source/paru_tasked_dgemm.cpp',
    'ParU/Source/paru_tasked_dtrsm.cpp',
    'ParU/Source/paru_dtrsm.cpp',
    'ParU/Source/paru_tuples.cpp',
    'ParU/Source/paru_umfpack_info.cpp',
    'ParU/Source/paru_update_rel_ind.cpp',
    'ParU/Source/paru_update_rowDeg.cpp',
    'ParU/Source/ParU_USolve.cpp',
    'ParU/Source/ParU_Version.cpp',
    'ParU/MATLAB/paru.c' ]

  paru_abs_src = []
  foreach p : paru_src
    paru_abs_src += get_option('suitesparse_src_path') / p
  endforeach

  paru_abs_inc = []
  foreach p : paru_inc
    paru_abs_inc += get_option('suitesparse_src_path') / p
  endforeach

  paru_deps = [ blas_dep, lapack_dep ]

  if get_option('build_for') == 'matlab'
    paru_compile_args = [ '-DBLAS64' ]
    if host_machine.system() == 'linux' or (host_machine.system() == 'darwin' and host_machine.cpu_family() == 'x86_64')
      # Under Windows, there is no mkl_serv_set_num_threads symbol in mkl.dll
      # (and the ParU manual says that Windows+MATLAB is not supported with OpenMP)
      # Under macOS/ARM, there is no MKL at all.
      paru_compile_args += [ '-DBLAS_Intel10_64ilp', '-DBLAS_Intel10' ]
      paru_deps += [ openmp_dep ]
      if host_machine.system() == 'linux'
        paru_deps += declare_dependency(link_args : '-l:mkl.so')
      else
        paru_deps += declare_dependency(link_args : matlab_path / 'bin' / matlab_arch / 'mkl.dylib')
      endif
    endif
    if host_machine.system() == 'windows'
      paru_compile_args += [ '-DBLAS_NO_UNDERSCORE' ]
    endif
  else # Octave
    paru_compile_args = [ '-DBLAS_OpenBLAS' ]
    paru_deps += [ openmp_dep ]
    # Do not try static linking since that OpenBLAS is shipped with Octave
    paru_deps += dependency('openblas', static : false)
  endif

  paru_mex_kwargs = mex_kwargs
  paru_mex_kwargs += { 'include_directories' : paru_mex_kwargs.get('include_directories') + paru_abs_inc }
  paru_mex_kwargs += { 'c_args' : paru_mex_kwargs.get('c_args') + paru_compile_args + [ '-Wno-sign-compare', '-Wno-unused', '-Wno-array-parameter', '-Wno-implicit-fallthrough', '-Wno-parentheses', '-Wno-empty-body', '-Wno-type-limits', '-Wno-uninitialized', '-Wno-unknown-pragmas', '-Wno-clobbered' ]}
  paru_mex_kwargs += { 'cpp_args' : paru_mex_kwargs.get('cpp_args') + paru_compile_args + [ '-Wno-old-style-cast', '-Wno-unused', '-Wno-sized-deallocation', '-Wno-unknown-pragmas' ] }

  shared_module('paru', paru_abs_src, kwargs : paru_mex_kwargs, dependencies : paru_deps)
endif


### Documentation

# Documentation is not built by default. A special-purpose target named “doc” triggers it.
# Executables that are needed for building documentation are not marked as required, to
# facilitate compilation of binaries in environments that cannot build the documentation.

## Reference manual

sphinx_build_exe = find_program('sphinx-build', required : false, disabler : true)

sphinx_src = [ 'doc/manual/source/conf.py',
               'doc/manual/source/time-series.rst',
               'doc/manual/source/parallel-execution.rst',
               'doc/manual/source/index.rst',
               'doc/manual/source/dynare-misc-commands.rst',
               'doc/manual/source/bibliography.rst',
               'doc/manual/source/running-dynare.rst',
               'doc/manual/source/the-model-file.rst',
               'doc/manual/source/introduction.rst',
               'doc/manual/source/reporting.rst',
               'doc/manual/source/examples.rst',
               'doc/manual/source/installation-and-configuration.rst',
               'doc/manual/source/_templates/about.html',
               'doc/manual/source/_templates/navigation.html',
               'doc/manual/source/_static/dlogo.svg',
               'doc/manual/source/_static/custom.css',
               'doc/manual/source/_static/report.png',
               'doc/manual/utils/dynare_dom.py',
               'doc/manual/utils/dynare_lex.py' ]

sphinx_defines = [ '-D', 'release=' + meson.project_version(),
                   '-D', 'version=' + meson.project_version() ]

if get_option('mathjax_path') != ''
  sphinx_defines += [ '-D', 'mathjax_path=' + get_option('mathjax_path') ]
endif

dynare_manual_html = custom_target('dynare-manual.html',
                                   output : 'dynare-manual.html', input : sphinx_src,
                                   command : [ sphinx_build_exe, '-n', '-b', 'html', sphinx_defines,
                                               '-d', '@PRIVATE_DIR@',
                                               meson.current_source_dir() / 'doc/manual/source',
                                               '@OUTPUT@' ],
                                   build_by_default : false,
                                   install : true, install_dir : 'share/doc/dynare')

# For the PDF, we are forced to use a wrapper around sphinx-build, because
# the PDF is created in a subdirectory, and meson does not provide an easy way to
# do this. See https://github.com/mesonbuild/meson/discussions/10488
dynare_manual_pdf = custom_target('dynare-manual.pdf',
                                  output : 'dynare-manual.pdf', input : sphinx_src,
                                  command : [ 'scripts/sphinx-build-pdf-wrapper',
                                              sphinx_build_exe.full_path(),
                                              meson.current_source_dir() / 'doc/manual/source',
                                              '@OUTPUT@', '@PRIVATE_DIR@', sphinx_defines ],
                                  build_by_default : false,
                                  install : true, install_dir : 'share/doc/dynare')

## Various PDFs

# The following defines the latexmk_command, macroprocessor_pdf and preprocessor_pdf variables
subdir('preprocessor/doc')

guide_pdf = custom_target('guide.pdf', output : 'guide.pdf', input : 'doc/guide.tex',
                          command : latexmk_command,
                          build_by_default : false,
                          install : true, install_dir : 'share/doc/dynare')

bvar_a_la_sims_pdf = custom_target('bvar-a-la-sims.pdf',
                                   output : 'bvar-a-la-sims.pdf', input : 'doc/bvar-a-la-sims.tex',
                                   command : latexmk_command,
                                   build_by_default : false,
                                   install : true, install_dir : 'share/doc/dynare')

dr_pdf = custom_target('dr.pdf', output : 'dr.pdf', input : 'doc/dr.tex',
                       command : latexmk_command,
                       env : { 'BIBINPUTS': meson.current_source_dir() / 'doc' },
                       depend_files : 'doc/dynare.bib',
                       build_by_default : false,
                       install : true, install_dir : 'share/doc/dynare')

sylvester_pdf = custom_target('sylvester.pdf', output : 'sylvester.pdf', input : 'doc/sylvester.tex',
                              command : latexmk_command,
                              build_by_default : false,
                              install : true, install_dir : 'share/doc/dynare')

tl_pdf = custom_target('tl.pdf', output : 'tl.pdf', input : 'doc/tl.tex',
                       command : latexmk_command,
                       build_by_default : false,
                       install : true, install_dir : 'share/doc/dynare')

parallel_pdf = custom_target('parallel.pdf',
                             output : 'parallel.pdf', input : 'doc/parallel/parallel.tex',
                             command : latexmk_command,
                             env : { 'TEXINPUTS': meson.current_source_dir() + '/doc/parallel:',
                                     'BIBINPUTS': meson.current_source_dir() / 'doc/parallel' },
                             depend_files : [ 'doc/parallel/marco.bib',
                                              'doc/parallel/AvenueParadigm.pdf',
                                              'doc/parallel/iVaNo_gain.pdf',
                                              'doc/parallel/iVaNo_time_comp.pdf',
	                                      'doc/parallel/netbook_complete_comp.pdf',
                                              'doc/parallel/netbook_complete_openclose.pdf',
	                                      'doc/parallel/netbook_partial_comp.pdf',
                                              'doc/parallel/netbook_partial_openclose.pdf',
                                              'doc/parallel/quest_complete_comp.pdf',
                                              'doc/parallel/quest_complete_openclose.pdf',
                                              'doc/parallel/quest_partial_comp.pdf',
                                              'doc/parallel/quest_partial_openclose.pdf',
                                              'doc/parallel/RWMH_quest1_PriorsAndPosteriors1Comp.pdf',
	                                      'doc/parallel/RWMH_quest1_PriorsAndPosteriors2Comp.pdf',
	                                      'doc/parallel/RWMH_quest1_PriorsAndPosteriors3Comp.pdf',
	                                      'doc/parallel/RWMH_quest1_PriorsAndPosteriors4Comp.pdf',
	                                      'doc/parallel/RWMH_quest1_PriorsAndPosteriors5Comp.pdf',
	                                      'doc/parallel/RWMH_quest1_PriorsAndPosteriors6Comp.pdf',
	                                      'doc/parallel/RWMH_quest1_PriorsAndPosteriors7Comp.pdf',
                                              'doc/parallel/waitbars1.pdf',
                                              'doc/parallel/waitbars2.pdf',
	                                      'doc/parallel/waitbarsP.pdf' ],
                             build_by_default : false,
                             install : true, install_dir : 'share/doc/dynare')

gsa_pdf = custom_target('gsa.pdf', output : 'gsa.pdf', input : 'doc/gsa/gsa.tex',
                        command : latexmk_command,
                        env : { 'BIBINPUTS': meson.current_source_dir() / 'doc/gsa' },
                        depend_files : 'doc/gsa/marco.bib',
                        build_by_default : false,
                        install : true, install_dir : 'share/doc/dynare')

dseries_reporting_pdf = custom_target('dseriesReporting.pdf', output : 'dseriesReporting.pdf',
                                      input : 'doc/dseries-and-reporting/dseriesReporting.tex',
                                      command : latexmk_command,
                                      build_by_default : false,
                                      install : true, install_dir : 'share/doc/dynare')

## Special-purpose target for building documentation

alias_target('doc',
             dynare_manual_html,
             dynare_manual_pdf,
             macroprocessor_pdf,
             preprocessor_pdf,
             guide_pdf,
             bvar_a_la_sims_pdf,
             dr_pdf,
             sylvester_pdf,
             tl_pdf,
             parallel_pdf,
             gsa_pdf,
             dseries_reporting_pdf)


### Tests implemented as .mod or .m (both integration and unit tests)

mod_and_m_tests = [
  # Particle files (moved to beginning as they take the most time in parallel execution)
  { 'test' : [ 'particle/dsge_base2.mod' ],
    'extra' : [ 'particle/mysample.m',
                'particle/dsge_base.inc' ] },
  { 'test' : [ 'particle/dsge_base_DSMH.mod' ],
    'extra' : [ 'particle/mysample.m',
                'particle/dsge_base.inc' ] },
  { 'test' : [ 'particle/dsge_base_HSSMC.mod' ],
    'extra' : [ 'particle/mysample.m',
                'particle/dsge_base.inc' ] },
  { 'test' : [ 'particle/dsge_base_Liu_West.mod' ],
    'extra' : [ 'particle/mysample.m',
                'particle/dsge_base.inc' ] },
  { 'test' : [ 'particle/dsge_base_apf.mod' ],
    'extra' : [ 'particle/mysample.m',
                'particle/dsge_base.inc' ] },
  { 'test' : [ 'particle/dsge_base_cpf.mod' ],
    'extra' : [ 'particle/mysample.m',
                'particle/dsge_base.inc' ] },
  { 'test' : [ 'particle/dsge_base_gf.mod' ],
    'extra' : [ 'particle/mysample.m',
                'particle/dsge_base.inc' ] },
  { 'test' : [ 'particle/dsge_base_gmf.mod' ],
    'extra' : [ 'particle/mysample.m',
                'particle/dsge_base.inc' ] },
  { 'test' : [ 'particle/dsge_base_nlkf.mod' ],
    'extra' : [ 'particle/mysample.m',
                'particle/dsge_base.inc' ] },
  { 'test' : [ 'particle/dsge_base_sis.mod' ],
    'extra' : [ 'particle/mysample.m',
                'particle/dsge_base.inc' ] },
  { 'test' : [ 'particle/first_spec_hssmc.mod' ],
    'extra' : [ 'particle/first_spec_common.inc' ] },
  { 'test' : [ 'particle/first_spec.mod',
               'particle/local_state_space_iteration_k_test.mod',
               'particle/local_state_space_iteration_3_test.mod' ],
    'extra' : [ 'particle/first_spec_common.inc' ] },
  { 'test' : [ 'particle/first_spec_MCMC.mod' ],
    'extra' : [ 'particle/first_spec_common.inc' ] },
  { 'test' : [ 'occbin/ppf/swpkfppf.mod'],
    'extra' : [ 'occbin/ppf/swpkfest_mh_mode.mat',
                'occbin/ppf/usdata_pkf.mat',
                'occbin/ppf/swpkfest/metropolis/swpkfest_mh1_blck1.mat',
                'occbin/ppf/swpkfest/metropolis/swpkfest_mh1_blck2.mat',
                'occbin/ppf/swpkfest/metropolis/swpkfest_mh1_blck3.mat',
                'occbin/ppf/swpkfest/metropolis/swpkfest_mh1_blck4.mat',
                'occbin/ppf/swpkfest/metropolis/swpkfest_mh_history_0.mat',
                'occbin/ppf/swpkfest/metropolis/swpkfest_mh_history_1.mat',
                'occbin/ppf/swpkfest/metropolis/swpkfest_mh_history_2.mat'] },
  { 'test' : [ 'decision_rules/first_order/walsh.mod' ],
    'extra' : [ 'decision_rules/first_order/+matlab/+namespace/y_k.m'] },
  { 'test' : [ 'occbin/model_irrcap_twoconstraints/dynrbc.mod',
               'occbin/model_irrcap_twoconstraints/dynrbc_0_std_shocks.mod' ],
    'extra' : [ 'occbin/model_irrcap_twoconstraints/dynrbc_common.inc' ] },
  { 'test' : [ 'occbin/model_borrcon/borrcon.mod',
               'occbin/model_borrcon/borrcon_0_std_shocks.mod' ],
    'extra' : [ 'occbin/model_borrcon/borrcon_common.inc' ] },
  { 'test' : [ 'occbin/filter/NKM.mod' ],
    'extra' : [ 'occbin/filter/dataobsfile.mat',
                'occbin/filter/NKM_mh_mode_saved.mat' ] },
  { 'test' : [ 'occbin/filter/NKM_calib_smoother.mod' ],
    'extra' : [ 'occbin/filter/dataobsfile.mat' ] },
  { 'test' : [ 'occbin/filter/NKM_0_std_shocks.mod' ],
    'extra' : [ 'occbin/filter/dataobsfile.mat',
                'occbin/filter/NKM_mh_mode_saved.mat' ] },
  { 'test' : [ 'moments/example1_hp_test.mod' ],
    'extra' : [ 'moments/example1_common.inc' ] },
  { 'test' : [ 'moments/fs2000_post_moments.mod' ],
    'extra' : [ 'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'moments/example1_order2_pruning.mod' ] },
  { 'test' : [ 'lmmcp/rbcii.mod' ] },
  { 'test' : [ 'lmmcp/purely_backward.mod' ] },
  { 'test' : [ 'lmmcp/purely_forward.mod' ] },
  { 'test' : [ 'lmmcp/MCP_ramsey.mod' ] },
  { 'test' : [ 'lmmcp/SMOPEC.mod' ] },
  { 'test' : [ 'estimation/TaRB/fs2000_tarb.mod' ],
    'extra' : [ 'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/pct/runpcttests.m' ],
    'extra' : [ 'estimation/pct/SW_2007_fmincon_parallel.mod',
                'estimation/pct/fs2000.mod',
                'estimation/pct/fs2000_data.m',
                'estimation/pct/usmodel_data.mat',
                'estimation/pct/tParallel.m' ] },
  { 'test' : [ 'parallel/ls2003.mod',
               'parallel/serial_mhblocks_1.mod',
               'parallel/serial_mhblocks_2.mod',
               'parallel/pct_mhblocks_1.mod',
               'parallel/pct_mhblocks_2.mod',
               'parallel/serial_vs_pct.mod' ],
    'extra' : [ 'parallel/_estimation_commands_with_different_seeds.inc',
                'parallel/_ls2003.inc',
                'parallel/data_ca1.m',
                'parallel/folder/endovars.mod' ] },
  { 'test' : [ 'observation_trends_and_prefiltering/MCMC/Trend_loglinear_no_prefilter_MC.mod',
	       'observation_trends_and_prefiltering/MCMC/Trend_no_prefilter_first_obs_MC.mod',
	       'observation_trends_and_prefiltering/MCMC/Trend_no_prefilter_MC.mod',
	       'observation_trends_and_prefiltering/MCMC/Trend_loglinear_prefilter_MC.mod',
	       'observation_trends_and_prefiltering/MCMC/Trend_loglin_prefilt_first_obs_MC.mod',
               'observation_trends_and_prefiltering/MCMC/Trend_prefilter_first_obs_MC.mod',
	       'observation_trends_and_prefiltering/MCMC/Trend_prefilter_MC.mod',
               'observation_trends_and_prefiltering/MCMC/Trend_loglin_no_prefilt_first_obs_MC.mod' ],
    'extra' : [ 'observation_trends_and_prefiltering/Trend_exp_model_no_prefilter_common.inc',
                'observation_trends_and_prefiltering/Trend_exp_model_prefilter_common.inc',
                'observation_trends_and_prefiltering/Trend_model_no_prefilter_common.inc',
                'observation_trends_and_prefiltering/Trend_model_prefilter_common.inc',
                'observation_trends_and_prefiltering/Trend_load_data_common.inc',
                'observation_trends_and_prefiltering/Trend_diagnostics_MCMC_common.inc',
                'observation_trends_and_prefiltering/generate_trend_stationary_AR1.m' ] },
  { 'test' : [ 'dates/dseries_interact.mod' ],
    'extra' : [ 'dates/data_uav.xlsx' ] },
  { 'test' : [ 'dates/fs2000.mod' ],
    'extra' : [ 'dates/fsdat_simul.m' ] },
  { 'test' : [ 'optimizers/fs2000_1.mod' ],
    'extra' : [ 'optimizers/fs2000.common.inc',
                'optimizers/optimizer_function_count.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'optimizers/fs2000_2.mod' ],
    'extra' : [ 'optimizers/fs2000.common.inc',
                'optimizers/optimizer_function_count.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'optimizers/fs2000_3.mod' ],
    'extra' : [ 'optimizers/fs2000.common.inc',
                'optimizers/optimizer_function_count.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'optimizers/fs2000_4.mod' ],
    'extra' : [ 'optimizers/fs2000.common.inc',
                'optimizers/optimizer_function_count.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'optimizers/fs2000_4_with_optim.mod' ],
    'extra' : [ 'optimizers/fs2000.common.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'optimizers/fs2000_5.mod' ],
    'extra' : [ 'optimizers/fs2000.common.inc',
                'optimizers/optimizer_function_count.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'optimizers/fs2000_6.mod' ],
    'extra' : [ 'optimizers/fs2000.common.inc',
                'optimizers/optimizer_function_count.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'optimizers/fs2000_7.mod' ],
    'extra' : [ 'optimizers/fs2000.common.inc',
                'optimizers/optimizer_function_count.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'optimizers/fs2000_8.mod' ],
    'extra' : [ 'optimizers/fs2000.common.inc',
                'optimizers/optimizer_function_count.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'optimizers/fs2000_8_with_optim.mod' ],
    'extra' : [ 'optimizers/fs2000.common.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'optimizers/fs2000_9.mod' ],
    'extra' : [ 'optimizers/fs2000.common.inc',
                'optimizers/optimizer_function_count.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'optimizers/fs2000_10.mod' ],
    'extra' : [ 'optimizers/fs2000.common.inc',
                'optimizers/optimizer_function_count.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'optimizers/fs2000_101.mod' ],
    'extra' : [ 'optimizers/fs2000.common.inc',
                'optimizers/optimizer_function_count.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'optimizers/fs2000_102.mod' ],
    'extra' : [ 'optimizers/fs2000.common.inc',
                'optimizers/optimizer_function_count.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'optimizers/fs2000_w.mod' ],
    'extra' : [ 'optimizers/fs2000.common.inc',
                'optimizers/optimizer_function_count.inc',
                'optimizers/optimizer_function_wrapper.m',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/slice/fs2000_slice.mod'],
    'extra' : [ 'estimation/slice/fs2000.inc',
                'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/slice/fs2000_slice_use_prior.mod'],
    'extra' : [ 'estimation/slice/fs2000.inc' ,
                'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/slice/bf.mod' ],
    'extra' : [ 'estimation/slice/datafile.mat' ] },
  { 'test' : [ 'estimation/slice/fs2000_estimate_initial_states_log.mod' ],
    'extra' : [ 'estimation/slice/fs2000.inc',
                'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/slice/fs2000_estimate_initial_states.mod' ],
    'extra' : [ 'estimation/slice/fs2000.inc'] },
  { 'test' : [ 'analytic_derivatives/fs2000_analytic_derivation.mod' ],
    'extra' : [ 'analytic_derivatives/fs2000.inc' ] },
  { 'test' : [ 'analytic_derivatives/fs2000_analytic_derivation_diag_H.mod' ],
    'extra' : [ 'analytic_derivatives/fs2000.inc' ] },
  { 'test' : [ 'analytic_derivatives/fs2000_analytic_derivation_newrat.mod' ] },
  { 'test' : [ 'analytic_derivatives/BrockMirman_PertParamsDerivs.mod' ],
    'extra' : [ 'analytic_derivatives/nBrockMirmanSYM.mat' ] },
  { 'test' : [ 'analytic_derivatives/burnside_3_order_PertParamsDerivs.mod' ] },
  { 'test' : [ 'stochastic_simulations/pruning/AnSchorfheide_pruned_state_space.mod',
               'stochastic_simulations/pruning/AS_pruned_state_space_red_shock.mod' ],
    'extra' : [ 'stochastic_simulations/pruning/Andreasen_et_al_2018_Dynare44Pruning_v2.mat'] },
  { 'test' : [ 'stochastic_simulations/pruning/fs2000_pruning.mod' ] },
  { 'test' : [ 'measurement_errors/fs2000_corr_me_ml_mcmc/fs2000_corr_ME.mod' ] },
  { 'test' : [ 'TeX/fs2000_corr_ME.mod' ],
    'extra' : [ 'fs2000/fsdat_simul.m'] },
  { 'test' : [ 'estimation/MH_recover/fs2000_recover.mod',
               'estimation/MH_recover/fs2000_recover_tarb.mod' ],
    'extra' : [ 'estimation/MH_recover/fs2000.common.inc',
                'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/MH_recover/fs2000_recover_2.mod' ],
    'extra' : [ 'estimation/MH_recover/fs2000.common.inc',
                'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/MH_recover/fs2000_recover_3.mod' ],
    'extra' : [ 'estimation/MH_recover/fs2000.common.inc',
                'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/MH_recover/fs2000_slice_recover.mod' ],
    'extra' : [ 'estimation/MH_recover/fs2000_slice_recover_mh_tmp_blck1.mat',
                'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/MH_recover/fs2000_recover_tmp.mod' ],
    'extra' : [ 'estimation/MH_recover/fs2000_recover_tmp_mh_tmp_blck1.mat',
                'estimation/MH_recover/fs2000_recover_tmp_mode.mat',
                'estimation/MH_recover/fs2000.common.inc',
                'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/fs2000.mod' ],
    'extra' : [ 'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/dime_mcmc/fs2000.mod' ],
    'extra' : [ 'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/hssmc/fs2000.mod' ],
    'extra' : [ 'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'gsa/ls2003a.mod' ],
    'extra' : [ 'gsa/data_ca1.m',
                'gsa/ls2003_model.inc' ] },
  { 'test' : [ 'gsa/ls2003.mod' ],
    'extra' : [ 'gsa/data_ca1.m',
                'gsa/ls2003_model.inc' ] },
  { 'test' : [ 'gsa/ls2003scr.mod' ],
    'extra' : [ 'gsa/data_ca1.m',
                'gsa/ls2003_model.inc' ] },
  { 'test' : [ 'gsa/ls2003ide.mod' ],
    'extra' : [ 'gsa/data_ca1.m',
                'gsa/ls2003_model.inc' ] },
  { 'test' : [ 'gsa/ls2003prior.mod' ],
    'extra' : [ 'gsa/data_ca1.m',
                'gsa/ls2003_model.inc' ] },
  { 'test' : [ 'gsa/cod_ML_morris/cod_ML_morris.mod' ] },
  { 'test' : [ 'gsa/morris/nk_est.mod' ],
    'extra' : [ 'gsa/morris/nk_est_data.m' ] },
  { 'test' : [ 'estimation/fs2000_MCMC_jumping_covariance.mod' ],
    'extra' : [ 'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'ms-sbvar/test_ms_variances_repeated_runs.mod' ],
    'extra' : [ 'ms-sbvar/msdata.m' ] },
  { 'test' : [ 'fs2000/fs2000.mod' ],
    'extra' : [ 'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'ls2003/ls2003_hessian_zero.mod' ] },
  { 'test' : [ 'exogenous-observed-variables/preprocessor.mod' ] },
  { 'test' : [ 'estimation/fs2000_with_weibull_prior.mod' ],
    'extra' : [ 'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/fs2000.mod',
               'estimation/fs2000_initialize_from_calib.mod',
               'estimation/fs2000_calibrated_covariance.mod',
               'estimation/fs2000_model_comparison.mod' ],
    'extra' : [ 'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/fs2000_estimated_params_init.mod' ] },
  { 'test' : [ 'estimation/fs2000_estimated_params_remove.mod' ],
    'extra' : [ 'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/fs2000_fast.mod' ],
    'extra' : [ 'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/fs2000.mod',
               'estimation/fs2000_init_from_previous.mod' ],
    'extra' : [ 'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/ls2003_endog_prior_restrict_estimation.mod' ],
    'extra' : [ 'gsa/data_ca1.m' ] },
  { 'test' : [ 'estimation/independent_mh/fs2000_independent_mh.mod' ],
    'extra' : [ 'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/heteroskedastic_shocks/fs2000_het.mod' ],
    'extra' : [ 'estimation/heteroskedastic_shocks/fs2000_het_model.inc',
                'estimation/heteroskedastic_shocks/fs2000_het_check.inc',
                'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/heteroskedastic_shocks/fs2000_het_corr.mod' ],
    'extra' : [ 'estimation/heteroskedastic_shocks/fs2000_het_model.inc',
                'estimation/heteroskedastic_shocks/fs2000_het_check.inc',
                'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/heteroskedastic_shocks/fs2000_het_sample_restriction.mod' ],
    'extra' : [ 'estimation/heteroskedastic_shocks/fs2000_het_model.inc',
                'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/heteroskedastic_shocks/fs2000_het_dates.mod' ],
    'extra' : [ 'estimation/heteroskedastic_shocks/fs2000_het_model.inc',
                'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/t_proposal/fs2000_student.mod' ],
    'extra' : [ 'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/tune_mh_jscale/fs2000.mod' ],
    'extra' : [ 'estimation/tune_mh_jscale/fs2000.inc',
                'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/method_of_moments/AnScho/AnScho_matched_moments.mod' ],
    'extra' : [ 'estimation/method_of_moments/AnScho/AnScho_MoM_common.inc',
                'estimation/method_of_moments/AnScho/AnScho_MoM_data_2.mat' ] },
  { 'test' : [ 'estimation/method_of_moments/AnScho/AnScho_GMM_estimParams0.mod' ],
    'extra' : [ 'estimation/method_of_moments/AnScho/AnScho_MoM_common.inc',
                'estimation/method_of_moments/AnScho/AnScho_MoM_data_2.mat' ] },
  { 'test' : [ 'estimation/method_of_moments/AnScho/AnScho_GMM_estimParams1.mod' ],
    'extra' : [ 'estimation/method_of_moments/AnScho/AnScho_MoM_common.inc',
                'estimation/method_of_moments/AnScho/AnScho_MoM_data_2.mat' ] },
  { 'test' : [ 'estimation/method_of_moments/AnScho/AnScho_GMM_estimParams2.mod' ],
    'extra' : [ 'estimation/method_of_moments/AnScho/AnScho_MoM_common.inc',
                'estimation/method_of_moments/AnScho/AnScho_MoM_data_2.mat' ] },
  { 'test' : [ 'estimation/method_of_moments/AnScho/AnScho_SMM_estimParams0.mod' ],
    'extra' : [ 'estimation/method_of_moments/AnScho/AnScho_MoM_common.inc',
                'estimation/method_of_moments/AnScho/AnScho_MoM_data_2.mat' ] },
  { 'test' : [ 'estimation/method_of_moments/AnScho/AnScho_SMM_estimParams1.mod' ],
    'extra' : [ 'estimation/method_of_moments/AnScho/AnScho_MoM_common.inc',
                'estimation/method_of_moments/AnScho/AnScho_MoM_data_2.mat' ] },
  { 'test' : [ 'estimation/method_of_moments/AnScho/AnScho_SMM_estimParams2.mod' ],
    'extra' : [ 'estimation/method_of_moments/AnScho/AnScho_MoM_common.inc',
                'estimation/method_of_moments/AnScho/AnScho_MoM_data_2.mat' ] },
  { 'test' : [ 'estimation/method_of_moments/RBC/RBC_MoM_Andreasen.mod' ],
    'extra' : [ 'estimation/method_of_moments/RBC/RBC_MoM_common.inc',
                'estimation/method_of_moments/RBC/RBC_Andreasen_Data_2.mat',
                'estimation/method_of_moments/RBC/RBC_MoM_steady_helper.m' ] },
  { 'test' : [ 'estimation/method_of_moments/RBC/RBC_MoM_SMM_ME.mod' ],
    'extra' : [ 'estimation/method_of_moments/RBC/RBC_MoM_common.inc',
                'estimation/method_of_moments/RBC/RBC_MoM_steady_helper.m' ] },
  { 'test' : [ 'estimation/method_of_moments/RBC/RBC_MoM_prefilter.mod' ],
    'extra' : [ 'estimation/method_of_moments/RBC/RBC_MoM_common.inc',
                'estimation/method_of_moments/RBC/RBC_MoM_steady_helper.m' ] },
  { 'test' : [ 'estimation/method_of_moments/RBC/RBC_MoM_SMM_ME_filter.mod' ],
    'extra' : [ 'estimation/method_of_moments/RBC/RBC_MoM_common.inc',
                'estimation/method_of_moments/RBC/RBC_MoM_steady_helper.m' ] },
  { 'test' : [ 'estimation/method_of_moments/RBC/RBC_MoM_optimizer.mod' ],
    'extra' : [ 'estimation/method_of_moments/RBC/RBC_MoM_common.inc',
                'estimation/method_of_moments/RBC/RBC_Andreasen_Data_2.mat',
                'estimation/method_of_moments/RBC/RBC_MoM_steady_helper.m' ] },
  { 'test' : [ 'estimation/method_of_moments/RBC/RBC_MoM_GMM_gradient_optim.mod' ],
    'extra' : [ 'estimation/method_of_moments/RBC/RBC_MoM_common.inc',
                'estimation/method_of_moments/RBC/RBC_Andreasen_Data_2.mat',
                'estimation/method_of_moments/RBC/RBC_MoM_steady_helper.m' ] },
  { 'test' : [ 'estimation/method_of_moments/AFVRR/AFVRR_M0.mod' ],
    'extra' : [ 'estimation/method_of_moments/AFVRR/AFVRR_common.inc',
                'estimation/method_of_moments/AFVRR/AFVRR_data.mat',
                'estimation/method_of_moments/AFVRR/AFVRR_steady_helper.m' ] },
  { 'test' : [ 'estimation/method_of_moments/AFVRR/AFVRR_MFB.mod' ],
    'extra' : [ 'estimation/method_of_moments/AFVRR/AFVRR_common.inc',
                'estimation/method_of_moments/AFVRR/AFVRR_data.mat',
                'estimation/method_of_moments/AFVRR/AFVRR_steady_helper.m' ] },
  { 'test' : [ 'estimation/method_of_moments/AFVRR/AFVRR_MFB_RRA.mod' ],
    'extra' : [ 'estimation/method_of_moments/AFVRR/AFVRR_common.inc',
                'estimation/method_of_moments/AFVRR/AFVRR_data.mat',
                'estimation/method_of_moments/AFVRR/AFVRR_steady_helper.m' ] },
  { 'test' : [ 'estimation/method_of_moments/CET/cet_interface.mod' ],
    'extra' : [ 'estimation/method_of_moments/CET/cet_data.mat',
                'estimation/method_of_moments/CET/cet_irf_matching_file.m',
                'estimation/method_of_moments/CET/cet_model.inc',
                'estimation/method_of_moments/CET/cet_original_mode.mat',
                'estimation/method_of_moments/CET/cet_steady_helper.m' ] },
  { 'test' : [ 'estimation/method_of_moments/CET/cet_imh.mod' ],
    'extra' : [ 'estimation/method_of_moments/CET/cet_data.mat',
                'estimation/method_of_moments/CET/cet_irf_matching_file.m',
                'estimation/method_of_moments/CET/cet_model.inc',
                'estimation/method_of_moments/CET/cet_original_mode.mat',
                'estimation/method_of_moments/CET/cet_steady_helper.m' ] },
  { 'test' : [ 'estimation/method_of_moments/CET/cet_mle.mod' ],
    'extra' : [ 'estimation/method_of_moments/CET/cet_data.mat',
                'estimation/method_of_moments/CET/cet_irf_matching_file.m',
                'estimation/method_of_moments/CET/cet_model.inc',
                'estimation/method_of_moments/CET/cet_original_mode.mat',
                'estimation/method_of_moments/CET/cet_steady_helper.m' ] },
  { 'test' : [ 'estimation/method_of_moments/CET/cet_rwmh.mod' ],
    'extra' : [ 'estimation/method_of_moments/CET/cet_data.mat',
                'estimation/method_of_moments/CET/cet_irf_matching_file.m',
                'estimation/method_of_moments/CET/cet_model.inc',
                'estimation/method_of_moments/CET/cet_original_mode.mat',
                'estimation/method_of_moments/CET/cet_steady_helper.m' ] },
  { 'test' : [ 'estimation/method_of_moments/CET/cet_mh_posterior_mode.mod' ],
    'extra' : [ 'estimation/method_of_moments/CET/cet_data.mat',
                'estimation/method_of_moments/CET/cet_irf_matching_file.m',
                'estimation/method_of_moments/CET/cet_model.inc',
                'estimation/method_of_moments/CET/cet_original_mode.mat',
                'estimation/method_of_moments/CET/cet_steady_helper.m' ] },
  { 'test' : [ 'estimation/method_of_moments/CET/cet_slice.mod' ],
    'extra' : [ 'estimation/method_of_moments/CET/cet_data.mat',
                'estimation/method_of_moments/CET/cet_irf_matching_file.m',
                'estimation/method_of_moments/CET/cet_model.inc',
                'estimation/method_of_moments/CET/cet_original_mode.mat',
                'estimation/method_of_moments/CET/cet_steady_helper.m' ] },
  { 'test' : [ 'estimation/method_of_moments/CET/cet_tarb.mod' ],
    'extra' : [ 'estimation/method_of_moments/CET/cet_data.mat',
                'estimation/method_of_moments/CET/cet_irf_matching_file.m',
                'estimation/method_of_moments/CET/cet_model.inc',
                'estimation/method_of_moments/CET/cet_original_mode.mat',
                'estimation/method_of_moments/CET/cet_steady_helper.m' ] },
  { 'test' : [ 'estimation/conditional-likelihood/1/fs2000_estimation_exact.mod',
               'estimation/conditional-likelihood/1/fs2000_estimation_conditional.mod' ],
    'extra' : [ 'estimation/fsdat_simul.m' ] },
  { 'test' : [ 'estimation/system_prior_restriction/Gali_2015.mod' ],
    'extra' : [ 'estimation/system_prior_restriction/Gali_2015_PC_slope.m',
                'estimation/system_prior_restriction/Gali_2015_prior_restrictions.m' ] },
  { 'test' : [ 'estimation/no_init_estimation_check_first_obs/fs2000_init_check.mod' ],
    'extra' : [ 'estimation/no_init_estimation_check_first_obs/fsdat_mat.m',
                'estimation/no_init_estimation_check_first_obs/fsdat_simul.m'] },
  { 'test' : [ 'estimation/univariate/nls/staticmodel.mod' ] },
  { 'test' : [ 'estimation/univariate/nls/dynamicmodel1.mod' ] },
  { 'test' : [ 'estimation/univariate/nls/dynamicmodel2.mod' ] },
  { 'test' : [ 'model_info/NKM_3_back.mod' ] },
  { 'test' : [ 'model_info/NKM_3_forward.mod' ] },
  { 'test' : [ 'model_info/bkk.mod' ] },
  { 'test' : [ 'moments/example1_var_decomp.mod' ],
    'extra' : [ 'moments/example1_common.inc' ] },
  { 'test' : [ 'moments/example1_bp_test.mod' ],
    'extra' : [ 'moments/example1_common.inc' ] },
  { 'test' : [ 'moments/example1_theoretical_simulated.mod' ],
    'extra' : [ 'moments/example1_common.inc' ] },
  { 'test' : [ 'moments/test_AR1_spectral_density.mod' ] },
  { 'test' : [ 'moments/example1_one_sided_hp_test.mod' ],
    'extra' : [ 'moments/example1_common.inc' ] },
  { 'test' : [ 'on-the-fly/ex1.mod' ] },
  { 'test' : [ 'on-the-fly/ex2.mod' ] },
  { 'test' : [ 'on-the-fly/ex3.mod' ] },
  { 'test' : [ 'on-the-fly/ex4.mod' ] },
  { 'test' : [ 'on-the-fly/ex5.mod' ] },
  { 'test' : [ 'on-the-fly/ex6.mod' ] },
  { 'test' : [ 'equation_tags/example1_with_tags.mod' ] },
  { 'test' : [ 'equation_tags/ramst_static_tag.mod' ] },
  { 'test' : [ 'equation_tags/ramst_static_tag_block.mod' ] },
  { 'test' : [ 'stochastic_simulations/example1.mod' ] },
  { 'test' : [ 'stochastic_simulations/example1_irf_shocks.mod' ] },
  { 'test' : [ 'stochastic_simulations/example1long.mod',
               'stochastic_simulations/example2long.mod',
               'stochastic_simulations/example2long_use_dll.mod' ] },
  { 'test' : [ 'stochastic_simulations/example2_simul_replic.mod' ] },
  { 'test' : [ 'stochastic_simulations/histval_sto.mod' ] },
  { 'test' : [ 'stochastic_simulations/histval_predetermined.mod' ] },
  { 'test' : [ 'macro_processor/test_ifndef.mod' ] },
  { 'test' : [ 'macro_processor/example1_macro.mod' ] },
  { 'test' : [ 'irfs/example1_unit_std.mod' ] },
  { 'test' : [ 'optimal_policy/OSR/osr_example.mod' ] },
  { 'test' : [ 'optimal_policy/OSR/osr_example_objective_correctness.mod' ] },
  { 'test' : [ 'optimal_policy/OSR/osr_objective_correctness_anal_deriv.mod' ] },
  { 'test' : [ 'optimal_policy/OSR/osr_example_obj_corr_non_stat_vars.mod' ] },
  { 'test' : [ 'optimal_policy/OSR/osr_example_param_bounds.mod' ] },
  { 'test' : [ 'optimal_policy/OSR/osr_obj_corr_algo_1.mod' ] },
  { 'test' : [ 'optimal_policy/OSR/osr_obj_corr_algo_3.mod' ] },
  { 'test' : [ 'optimal_policy/OSR/osr_obj_corr_algo_4.mod' ] },
  { 'test' : [ 'optimal_policy/OSR/osr_obj_corr_algo_7.mod' ] },
  { 'test' : [ 'optimal_policy/OSR/osr_obj_corr_algo_8.mod' ] },
  { 'test' : [ 'optimal_policy/OSR/osr_obj_corr_algo_9.mod' ] },
  { 'test' : [ 'optimal_policy/OSR/OSR_planner_objective.mod' ] },
  { 'test' : [ 'optimal_policy/ramsey_.mod' ] },
  { 'test' : [ 'optimal_policy/nk_ramsey.mod' ] },
  { 'test' : [ 'optimal_policy/nk_ramsey_model.mod' ] },
  { 'test' : [ 'optimal_policy/nk_ramsey_det.mod' ] },
  { 'test' : [ 'optimal_policy/neo_ramsey_det_bytecode.mod' ] },
  { 'test' : [ 'optimal_policy/nk_ramsey_expectation.mod',
               'optimal_policy/nk_ramsey_expectation_a.mod' ] },
  { 'test' : [ 'optimal_policy/mult_elimination_test.mod' ] },
  { 'test' : [ 'optimal_policy/neo_growth.mod',
               'optimal_policy/neo_growth_ramsey.mod' ],
    'extra' : [ 'optimal_policy/neo_growth_common.inc',
                'optimal_policy/neo_growth_ramsey_common.inc' ] },
  { 'test' : [ 'optimal_policy/neo_growth_k_order.mod',
               'optimal_policy/neo_growth_ramsey_k_order.mod' ],
    'extra' : [ 'optimal_policy/neo_growth_common.inc',
                'optimal_policy/neo_growth_ramsey_common.inc' ] },
  { 'test' : [ 'optimal_policy/Ramsey/ramsey_ex_initval.mod',
               'optimal_policy/Ramsey/ramsey_ex.mod' ],
    'extra' : [ 'optimal_policy/Ramsey/find_c.m' ] },
  { 'test' : [ 'optimal_policy/Ramsey/ramsey_ex_initval_AR2.mod' ] },
  { 'test' : [ 'optimal_policy/Ramsey/ramsey_ex_aux.mod' ],
    'extra' : [ 'optimal_policy/Ramsey/find_c.m' ] },
  { 'test' : [ 'optimal_policy/Ramsey/ramsey_histval.mod' ],
    'extra' : [ 'optimal_policy/Ramsey/find_c.m' ] },
  { 'test' : [ 'optimal_policy/Ramsey/Gali_commitment.mod' ] },
  { 'test' : [ 'optimal_policy/RamseyConstraints/test1.mod' ] },
  { 'test' : [ 'optimal_policy/Ramsey/Ramsey_Example_estimation.mod' ] },
  { 'test' : [ 'discretionary_policy/dennis_1.mod',
               'discretionary_policy/dennis_1_estim.mod',
               'discretionary_policy/dennis_1_estim_MoM.mod' ] },
  { 'test' : [ 'discretionary_policy/Gali_discretion.mod' ] },
  { 'test' : [ 'discretionary_policy/NK_discretion_forward.mod' ] },
  { 'test' : [ 'discretionary_policy/Gali_2015_chapter_3.mod',
               'discretionary_policy/Gali_2015_chapter_3_nonlinear.mod' ] },
  { 'test' : [ 'histval_initval_file/ramst_data_generate.mod',
               'histval_initval_file/ramst_initval_file_data.m',
               'histval_initval_file/ramst_initval_file.mod' ],
    'extra' : [ 'histval_initval_file/ramst_data.xlsx' ] },
  { 'test' : [ 'histval_initval_file/ramst_initval_file_with_dseries.mod' ] },
  { 'test' : [ 'histval_initval_file/ramst_initval_file_data.m',
               'histval_initval_file/ramst_datafile.mod' ] },
  { 'test' : [ 'histval_initval_file/sim_exo_lead_lag.mod',
               'histval_initval_file/sim_exo_lead_lag_initvalf.mod' ] },
  { 'test' : [ 'stat_functions/ramst_normcdf_and_friends_mfile.mod' ],
    'extra' : [ 'stat_functions/ramst_normcdf_and_friends.inc' ] },
  { 'test' : [ 'stat_functions/ramst_normcdf_and_friends_use_dll.mod' ],
    'extra' : [ 'stat_functions/ramst_normcdf_and_friends.inc' ] },
  { 'test' : [ 'stat_functions/ramst_normcdf_and_friends_bytecode.mod' ],
    'extra' : [ 'stat_functions/ramst_normcdf_and_friends.inc' ] },
  { 'test' : [ 'stat_functions/example1_abs_sign.mod' ] },
  { 'test' : [ 'preprocessor_checks/comments.mod' ] },
  { 'test' : [ 'preprocessor_checks/example1_mlv.mod' ] },
  { 'test' : [ 'deterministic_simulations/ramst.mod',
               'preprocessor_checks/ramst_model_edit.mod' ] },
  { 'test' : [ 'auxiliary_variables/test1.mod' ] },
  { 'test' : [ 'expectations/expectation.mod' ] },
  { 'test' : [ 'expectations/expectation_ss.mod' ] },
  { 'test' : [ 'expectations/expectation_ss_old.mod' ],
    'extra' : [ 'expectations/expectation_ss_old_steadystate.m' ] },
  { 'test' : [ 'expectations/expectation_nested.mod' ] },
  { 'test' : [ 'steady_state/walsh1_initval.mod' ] },
  { 'test' : [ 'steady_state/walsh1_old_ss.mod' ],
    'extra' : [ 'steady_state/walsh1_old_ss_steadystate.m' ] },
  { 'test' : [ 'steady_state/walsh1_ssm.mod' ] },
  { 'test' : [ 'steady_state/walsh1_ssm_block.mod' ] },
  { 'test' : [ 'steady_state/multi_leads.mod' ] },
  { 'test' : [ 'steady_state/example1_trust_region.mod' ] },
  { 'test' : [ 'steady_state/example1_block_trust_region.mod' ] },
  { 'test' : [ 'steady_state/Gali_2015_chapter_6_4.mod' ] },
  { 'test' : [ 'steady_state/ramst_solve_algo_0.mod' ] },
  { 'test' : [ 'steady_state/fs2000_ssfile.mod' ],
    'extra' : [ 'steady_state/fs2000_ssfile_aux.m' ] },
  { 'test' : [ 'steady_state_operator/standard.mod' ] },
  { 'test' : [ 'steady_state_operator/use_dll.mod' ] },
  { 'test' : [ 'steady_state_operator/block.mod' ] },
  { 'test' : [ 'steady_state_operator/bytecode_test.mod' ] },
  { 'test' : [ 'block_bytecode/ireland.mod' ] },
  { 'test' : [ 'block_bytecode/lola_solve_one_boundary.mod',
               'block_bytecode/lola_solve_one_boundary_mfs1.mod',
               'block_bytecode/lola_solve_one_boundary_mfs2.mod',
               'block_bytecode/lola_solve_one_boundary_mfs3.mod' ],
    'extra' : [ 'block_bytecode/lola_common.inc',
                'block_bytecode/lola_data.mat' ] },
  { 'test' : [ 'decision_rules/k_order/fs2000k2a.mod',
               'decision_rules/k_order/fs2000k2_use_dll.mod',
               'decision_rules/k_order/fs2000k_1_use_dll.mod',
               'decision_rules/k_order/fs2000k3_use_dll.mod',
               'decision_rules/k_order/fs2000k2_m.mod',
               'decision_rules/k_order/fs2000k_1_m.mod',
               'decision_rules/k_order/fs2000k3_m.mod' ] },
  { 'test' : [ 'stochastic_simulations/pruning/fs2000k3_p.mod' ] },
  { 'test' : [ 'decision_rules/k_order/burnside_k_order.mod' ] },
  { 'test' : [ 'stochastic_simulations/fs2000_k_order_simul.mod' ] },
  { 'test' : [ 'partial_information/PItest3aHc0PCLsimModPiYrVarobsAll.mod' ] },
  { 'test' : [ 'partial_information/PItest3aHc0PCLsimModPiYrVarobsCNR.mod' ] },
  { 'test' : [ 'partial_information/Partial_ramsey.mod' ] },
  { 'test' : [ 'partial_information/Partial_steady_state_model.mod' ] },
  { 'test' : [ 'arima/mod1.mod',
               'arima/mod1a.mod',
               'arima/mod1b.mod',
               'arima/mod1c.mod' ] },
  { 'test' : [ 'arima/mod2.mod',
               'arima/mod2a.mod',
               'arima/mod2b.mod',
               'arima/mod2c.mod' ] },
  { 'test' : [ 'data/mod1a.mod' ],
    'extra' : [ 'data/test.xlsx' ] },
  { 'test' : [ 'fs2000/fs2000a.mod' ],
    'extra' : [ 'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'fs2000/fs2000c.mod' ] },
  { 'test' : [ 'fs2000/fs2000_calib.mod' ],
    'extra' : [ 'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'fs2000/fs2000_calib_dseries.mod' ],
    'extra' : [ 'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'fs2000/fs2000_missing_data.mod' ],
    'extra' : [ 'fs2000/fsdat_simul_missing_obs.m' ] },
  { 'test' : [ 'fs2000/fs2000_dseries_a.mod' ],
    'extra' : [ 'fs2000/fsdat_simul_dseries.m' ] },
  { 'test' : [ 'fs2000/fs2000_dseries_b.mod' ],
    'extra' : [ 'fs2000/fsdat_simul_dseries.m' ] },
  { 'test' : [ 'steady_state/homotopy/homotopy1_test.mod' ],
    'extra' : [ 'steady_state/homotopy/common.mod' ] },
  { 'test' : [ 'steady_state/homotopy/homotopy2_test.mod' ],
    'extra' : [ 'steady_state/homotopy/common.mod' ] },
  { 'test' : [ 'steady_state/homotopy/homotopy3_test.mod' ],
    'extra' : [ 'steady_state/homotopy/common.mod' ] },
  { 'test' : [ 'steady_state/homotopy/homotopy_from_initval_to_endval.mod' ],
    'extra' : [ 'steady_state/homotopy/common.mod' ] },
  { 'test' : [ 'bvar_a_la_sims/bvar_standalone.mod' ],
    'extra' : [ 'bvar_a_la_sims/bvar_sample.m' ] },
  { 'test' : [ 'bvar_a_la_sims/bvar_and_dsge.mod' ],
    'extra' : [ 'bvar_a_la_sims/bvar_sample.m' ] },
  { 'test' : [ 'AIM/fs2000x10L9_L.mod',
               'AIM/fs2000x10L9_L_AIM.mod' ] },
  { 'test' : [ 'AIM/fs2000x10_L9_L.mod',
               'AIM/fs2000x10_L9_L_AIM.mod' ] },
  { 'test' : [ 'AIM/fs2000_b1L1L.mod',
               'AIM/fs2000_b1L1L_AIM.mod' ] },
  { 'test' : [ 'AIM/ls2003_2L0L.mod',
               'AIM/ls2003_2L0L_AIM.mod' ] },
  { 'test' : [ 'AIM/ls2003_2L2L.mod',
               'AIM/ls2003_2L2L_AIM.mod' ] },
  { 'test' : [ 'conditional_variance_decomposition/example1.mod' ] },
  { 'test' : [ 'conditional_variance_decomposition/example3.mod' ],
    'extra' : [ 'conditional_variance_decomposition/example3_steady_state_helper.m' ] },
  { 'test' : [ 'dsge-var/simul_hybrid.mod',
               'dsge-var/dsgevar_forward_calibrated_lambda.mod',
               'dsge-var/dsgevar_forward_estimated_lambda.mod' ] },
  { 'test' : [ 'external_function/benchmark.mod',
               'external_function/first_deriv_given.mod',
               'external_function/first_deriv_given_dll.mod',
               'external_function/first_and_2nd_deriv_given.mod',
               'external_function/first_and_2nd_deriv_given_namespace.mod',
               'external_function/first_and_2nd_deriv_given_dll.mod',
               'external_function/first_and_2nd_deriv_given_dll_namespace.mod',
               'external_function/no_deriv_given.mod',
               'external_function/no_deriv_given_dll.mod',
               'external_function/extfun_in_mlv.mod' ],
    'extra' : [ 'external_function/extFunDeriv.m',
                'external_function/extFunNoDerivs.m',
                'external_function/extFunWithFirstAndSecondDerivs.m',
                'external_function/+matlab/+namespace/extFunWithFirstAndSecondDerivsNamespace.m' ] },
  { 'test' : [ 'minimal_state_space_system/as2007_minimal.mod' ] },
  { 'test' : [ 'minimal_state_space_system/sw_minimal.mod' ] },
  { 'test' : [ 'identification/kim/kim2.mod' ] },
  { 'test' : [ 'identification/as2007/as2007.mod' ] },
  { 'test' : [ 'identification/as2007/as2007_kronflags.mod' ] },
  { 'test' : [ 'identification/as2007/as2007_QT.mod' ],
    'extra' : [ 'identification/as2007/G_QT.mat' ] },
  { 'test' : [ 'identification/as2007/as2007_QT_equal_autocorr.mod' ] },
  { 'test' : [ 'identification/as2007/as2007_order_1_2_3.mod' ] },
  { 'test' : [ 'identification/as2007/as2007_order3_no_stoch_simul.mod' ] },
  { 'test' : [ 'identification/BrockMirman/BrockMirman.mod' ] },
  { 'test' : [ 'identification/cgg/cgg_criteria_differ.mod' ] },
  { 'test' : [ 'identification/ident_unit_root/ident_unit_root.mod' ] },
  { 'test' : [ 'identification/rbc_ident/rbc_ident_std_as_structural_par.mod',
               'identification/rbc_ident/rbc_ident_varexo_only.mod' ] },
  { 'test' : [ 'identification/correlated_errors/fs2000_corr.mod' ] },
  { 'test' : [ 'identification/forward_looking/forward_looking_empty_ghx.mod' ],
    'extra' : [ 'identification/forward_looking/forward_looking_common.inc' ] },
  { 'test' : [ 'identification/forward_looking/forward_looking_varobs_x.mod' ],
    'extra' : [ 'identification/forward_looking/forward_looking_common.inc' ] },
  { 'test' : [ 'model_diagnostics/example2_block_model_diag.mod' ] },
  { 'test' : [ 'deterministic_simulations/initval_endval_blocks/initval_endval_test.mod' ] },
  { 'test' : [ 'deterministic_simulations/initval_endval_blocks/histval_initval_test.mod' ] },
  { 'test' : [ 'conditional_forecasts/1/fs2000_cal.mod' ] },
  { 'test' : [ 'conditional_forecasts/2/fs2000_est.mod' ],
    'extra' : [ 'conditional_forecasts/2/fsdat_simul.m' ] },
  { 'test' : [ 'conditional_forecasts/3/fs2000_conditional_forecast_initval.mod' ] },
  { 'test' : [ 'conditional_forecasts/4/fs2000_conditional_forecast_histval.mod' ] },
  { 'test' : [ 'conditional_forecasts/6/fs2000_cal.mod' ] },
  { 'test' : [ 'recursive/ls2003.mod' ],
    'extra' : [ 'recursive/data_ca1.m' ] },
  { 'test' : [ 'recursive/ls2003_bayesian.mod' ],
    'extra' : [ 'recursive/data_ca1.m' ] },
  { 'test' : [ 'recursive/ls2003_bayesian_xls.mod' ],
    'extra' : [ 'recursive/data_ca1_xls.xlsx' ] },
  { 'test' : [ 'recursive/ls2003_bayesian_csv.mod' ],
    'extra' : [ 'recursive/data_ca1_csv.csv' ] },
  { 'test' : [ 'ms-sbvar/test_exclusions.mod' ],
    'extra' : [ 'ms-sbvar/msdata.m' ] },
  { 'test' : [ 'ms-sbvar/test_exclusions_nc.mod' ],
    'extra' : [ 'ms-sbvar/msdata.m' ] },
  { 'test' : [ 'ms-sbvar/test_lower_cholesky.mod' ],
    'extra' : [ 'ms-sbvar/msdata.m' ] },
  { 'test' : [ 'ms-sbvar/test_lower_cholesky_a.mod' ],
    'extra' : [ 'ms-sbvar/msdata.m' ] },
  { 'test' : [ 'ms-sbvar/test_lower_cholesky_nc.mod' ],
    'extra' : [ 'ms-sbvar/msdata.m' ] },
  { 'test' : [ 'ms-sbvar/test_upper_cholesky.mod' ],
    'extra' : [ 'ms-sbvar/msdata.m' ] },
  { 'test' : [ 'ms-sbvar/test_upper_cholesky_nc.mod' ],
    'extra' : [ 'ms-sbvar/msdata.m' ] },
  { 'test' : [ 'ms-sbvar/test_ms_variances.mod' ],
    'extra' : [ 'ms-sbvar/msdata.m' ] },
  { 'test' : [ 'ms-dsge/test_ms_dsge.mod' ] },
  { 'test' : [ 'kalman/lyapunov/fs2000_lyap.mod' ],
    'extra' : [ 'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'kalman/lik_init/fs2000_ns_lik_init_2.mod' ],
    'extra' : [ 'kalman/lik_init/fs2000_ns_common.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'kalman/lik_init/fs2000_ns_lik_init_3.mod' ],
    'extra' : [ 'kalman/lik_init/fs2000_ns_common.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'kalman/lik_init/fs2000_ns_lik_init_5.mod' ],
    'extra' : [ 'kalman/lik_init/fs2000_ns_common.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'kalman/lik_init/fs2000_lik_init_1.mod' ],
    'extra' : [ 'kalman/lik_init/fs2000_common.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'kalman/lik_init/fs2000_lik_init_2.mod' ],
    'extra' : [ 'kalman/lik_init/fs2000_common.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'kalman/lik_init/fs2000_lik_init_3.mod' ],
    'extra' : [ 'kalman/lik_init/fs2000_common.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'kalman/lik_init/fs2000_lik_init_4.mod' ],
    'extra' : [ 'kalman/lik_init/fs2000_common.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'kalman/lik_init/fs2000_lik_init_5.mod' ],
    'extra' : [ 'kalman/lik_init/fs2000_common.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'kalman/block/fs2000.mod' ],
    'extra' : [ 'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'kalman/block/fs2000_missing_data.mod' ],
    'extra' : [ 'fs2000/fsdat_simul_missing_obs.m' ] },
  { 'test' : [ 'kalman_initial_state/fs2000_smoother_only_initial_state.mod' ],
    'extra' : [ 'kalman_initial_state/fs2000_common.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'kalman_initial_state/fs2000_ns_smoother_only_initial_state.mod' ],
    'extra' : [ 'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'kalman_initial_state/fs2000_kalman_initial.mod' ],
    'extra' : [ 'kalman_initial_state/fs2000_common.inc',
                'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'kalman_initial_state/fs2000_kalman_initial_2_lag.mod' ],
    'extra' : [ 'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'kalman_filter_smoother/gen_data.mod',
               'kalman_filter_smoother/algo1.mod',
               'kalman_filter_smoother/algo2.mod',
               'kalman_filter_smoother/algo3.mod',
               'kalman_filter_smoother/algo4.mod',
               'kalman_filter_smoother/algo4a.mod',
               'kalman_filter_smoother/algo4b.mod',
               'kalman_filter_smoother/algoH1.mod',
               'kalman_filter_smoother/algoH2.mod',
               'kalman_filter_smoother/algoH3.mod',
               'kalman_filter_smoother/fs2000.mod',
               'kalman_filter_smoother/fs2000_1.mod',
               'kalman_filter_smoother/fs2000_2.mod',
               'kalman_filter_smoother/fs2000a.mod' ],
    'extra' : [ 'kalman_filter_smoother/fsdat_simul.m' ] },
  { 'test' : [ 'kalman_filter_smoother/fs2000_smoother_only.mod' ],
    'extra' : [ 'kalman_filter_smoother/fsdat_simul.m' ] },
  { 'test' : [ 'kalman_filter_smoother/fs2000_smoother_only_ns.mod' ],
    'extra' : [ 'kalman_filter_smoother/fsdat_simul.m' ] },
  { 'test' : [ 'kalman_filter_smoother/fs2000_smoother_redux.mod' ],
    'extra' : [ 'kalman_filter_smoother/fsdat_simul.m' ] },
  { 'test' : [ 'kalman_filter_smoother/test_compute_Pinf_Pstar_data.mod',
               'kalman_filter_smoother/test_compute_Pinf_Pstar.mod' ] },
  { 'test' : [ 'kalman_filter_smoother/SOE.mod' ],
    'extra' : [ 'kalman_filter_smoother/SOE_data_file.m' ] },
  { 'test' : [ 'kalman_filter_smoother/check_variable_dimensions/fs2000.mod' ],
    'extra' : [ 'kalman_filter_smoother/fsdat_simul.m' ] },
  { 'test' : [ 'kalman_filter_smoother/check_variable_dimensions/fs2000_ML.mod' ],
    'extra' : [ 'kalman_filter_smoother/fsdat_simul.m' ] },
  { 'test' : [ 'kalman_filter_smoother/rbc_non_stationary_obs.mod' ],
    'extra' : [ 'kalman_filter_smoother/data_rbc_non_stationary_obs.m' ] },
  { 'test' : [ 'kalman/likelihood_from_dynare/fs2000_uncorr_ME.mod',
               'kalman/likelihood_from_dynare/fs2000_corr_ME.mod',
               'kalman/likelihood_from_dynare/fs2000_corr_ME_missing.mod',
               'kalman/likelihood_from_dynare/fs2000_uncorr_ME_missing.mod' ],
    'extra' : [ 'kalman/likelihood_from_dynare/fs2000_model.inc',
                'kalman/likelihood_from_dynare/fs2000_estimation_check.inc' ] },
  { 'test' : [ 'kalman/likelihood_from_dynare/fs2000ns_uncorr_ME.mod',
               'kalman/likelihood_from_dynare/fs2000ns_corr_ME.mod',
               'kalman/likelihood_from_dynare/fs2000ns_corr_ME_missing.mod',
               'kalman/likelihood_from_dynare/fs2000ns_uncorr_ME_missing.mod' ],
    'extra' : [ 'kalman/likelihood_from_dynare/fs2000ns_model.inc',
                'kalman/likelihood_from_dynare/fs2000ns_estimation_check.inc' ] },
  { 'test' : [ 'decision_rules/second_order/burnside_1.mod' ] },
  { 'test' : [ 'kalman_filter_smoother/compare_results_simulation/fs2000_ML.mod' ],
    'extra' : [ 'kalman_filter_smoother/compare_results_simulation/fsdat_simul_logged.m' ] },
  { 'test' : [ 'kalman_filter_smoother/compare_results_simulation/fs2000_ML_loglinear.mod' ],
    'extra' : [ 'kalman_filter_smoother/fsdat_simul.m',
                'kalman_filter_smoother/compare_results_simulation/fsdat_simul_logged.m' ] },
  { 'test' : [ 'kalman_filter_smoother/compare_results_simulation/fs2000.mod' ],
    'extra' : [ 'kalman_filter_smoother/compare_results_simulation/fsdat_simul_logged.m' ] },
  { 'test' : [ 'kalman_filter_smoother/compare_results_simulation/fs2000_loglinear.mod' ],
    'extra' : [ 'kalman_filter_smoother/fsdat_simul.m',
                'kalman_filter_smoother/compare_results_simulation/fsdat_simul_logged.m' ] },
  { 'test' : [ 'kalman_skewed/skew_in_shocks_block.mod',
               'kalman_skewed/skew_in_estimated_params.mod',
               'kalman_skewed/ireland2004_stoch_simul.mod',
               'kalman_skewed/ireland2004_ml.mod',
               'kalman_skewed/ireland2004_bayes.mod' ],
    'extra' : [ 'kalman_skewed/ireland2004_model.inc',
                'kalman_skewed/ireland2004_post1980.csv' ] },
  { 'test' : [ 'decision_rules/second_order/burnside_1.mod' ] },
  { 'test' : [ 'decision_rules/second_order/ds1.mod',
               'decision_rules/second_order/ds2.mod' ] },
  { 'test' : [ 'ep/rbc2.mod' ] },
  { 'test' : [ 'ep/rbc_mc.mod' ] },
  { 'test' : [ 'ep/ar_backward_xfail.mod' ],
    'should_fail' : true },
  { 'test' : [ 'ep/pfm_xfail.mod' ],
    'should_fail' : true },
  { 'test' : [ 'ep/ar_backward.mod' ]},
  { 'test' : [ 'ep/rbcii.mod',
               'ep/rbcii_MCP.mod' ],
    'extra' : [ 'ep/rbcii-calibration.inc',
                'ep/rbcii_steadystate.m' ] },
  { 'test' : [ 'ep/zeroshocks/rstrue.mod',
               'ep/zeroshocks/rs.mod' ] },
  { 'test' : [ 'ep/zeroshocks/rstrue.mod',
               'ep/zeroshocks/rs.mod' ] },
  { 'test' : [ 'ep/linearmodel0.mod' ] },
  { 'test' : [ 'ep/linearmodel1.mod' ] },
  { 'test' : [ 'ep/rbc_bytecode.mod' ] },
  { 'test' : [ 'ep/burnside.mod' ] },
  { 'test' : [ 'ep/hybrid/burnside_xfail.mod' ],
    'should_fail' : true },
  { 'test' : [ 'ep/hybrid/rs.mod' ] },
  { 'test' : [ 'ep/integration/burnside.mod' ] },
  { 'test' : [ 'ep/integration/burnside2.mod' ] },

  { 'test' : [ 'stochastic_simulations/example1_noprint.mod' ] },
  { 'test' : [ 'stochastic_simulations/example_seeds.mod' ] },
  { 'test' : [ 'stochastic-backward-models/solow_cd.mod' ] },
  { 'test' : [ 'stochastic-backward-models/solow_ces.mod' ] },
  { 'test' : [ 'stochastic-backward-models/solow_cd_with_steadystate.mod' ] },
  { 'test' : [ 'deterministic_simulations/ramst.mod' ] },
  { 'test' : [ 'deterministic_simulations/ramst_dseries.mod' ] },
  { 'test' : [ 'deterministic_simulations/ramst_a.mod' ] },
  { 'test' : [ 'deterministic_simulations/ramst_mshocks.mod' ] },
  { 'test' : [ 'deterministic_simulations/ramst_mshocks_vec.mod' ] },
  { 'test' : [ 'deterministic_simulations/ramst_mshocks_relative_to_initval.mod' ] },
  { 'test' : [ 'deterministic_simulations/predetermined_variables.mod' ] },
  { 'test' : [ 'deterministic_simulations/histval_det.mod' ] },
  { 'test' : [ 'deterministic_simulations/ramst_vec.mod' ] },
  { 'test' : [ 'deterministic_simulations/Solow_no_varexo.mod' ] },
  { 'test' : [ 'deterministic_simulations/simul_ZLB_purely_forward.mod' ] },
  { 'test' : [ 'deterministic_simulations/simul_ZLB_purely_forward_xfail.mod' ],
    'should_fail' : true  },
  { 'test' : [ 'deterministic_simulations/Irreversible_investment.mod',
               'deterministic_simulations/Irreversible_investment_leadlag.mod' ] },
  { 'test' : [ 'deterministic_simulations/linear_state_space_arma.mod' ] },
  { 'test' : [ 'deterministic_simulations/linear_state_space_arma_condpaths.mod' ] },
  { 'test' : [ 'deterministic_simulations/purely_forward/ar1.mod' ] },
  { 'test' : [ 'deterministic_simulations/purely_forward/nk.mod' ] },
  { 'test' : [ 'deterministic_simulations/purely_backward/ar1.mod' ] },
  { 'test' : [ 'deterministic_simulations/purely_static/toto.mod' ] },
  { 'test' : [ 'deterministic_simulations/rbc_det1.mod' ] },
  { 'test' : [ 'deterministic_simulations/rbc_det2.mod' ] },
  { 'test' : [ 'deterministic_simulations/rbc_det3.mod' ] },
  { 'test' : [ 'deterministic_simulations/rbc_det4.mod' ] },
  { 'test' : [ 'deterministic_simulations/rbc_det4_endval_initval_xfail.mod' ],
    'should_fail' : true },
  { 'test' : [ 'deterministic_simulations/rbc_det5.mod' ] },
  { 'test' : [ 'deterministic_simulations/rbc_det6.mod' ] },
  { 'test' : [ 'deterministic_simulations/homotopy.mod' ] },
  { 'test' : [ 'deterministic_simulations/homotopy_histval.mod' ] },
  { 'test' : [
    'deterministic_simulations/homotopy_endval_steady_linearization.mod' ],
    'extra' : [ 'deterministic_simulations/homotopy_linearization.inc' ] },
  { 'test' : [ 'deterministic_simulations/homotopy_marginal_linearization.mod' ],
    'extra' : [ 'deterministic_simulations/homotopy_linearization.inc' ] },
  { 'test' : [ 'deterministic_simulations/homotopy_linearization_condpaths.mod' ],
    'extra' : [ 'deterministic_simulations/homotopy_linearization.inc' ] },
  { 'test' : [ 'deterministic_simulations/homotopy_marglinearization_condpaths.mod' ],
    'extra' : [ 'deterministic_simulations/homotopy_linearization.inc' ] },
  { 'test' : [ 'deterministic_simulations/homotopy_exclude_varexo.mod' ] },
  { 'test' : [ 'deterministic_simulations/rbc_det_exo_lag_2a.mod',
               'deterministic_simulations/rbc_det_exo_lag_2b.mod',
               'deterministic_simulations/rbc_det_exo_lag_2c.mod' ] },
  { 'test' : [ 'deterministic_simulations/multiple_lead_lags/sim_base.mod',
               'deterministic_simulations/multiple_lead_lags/sim_exo_lead_lag_aux_vars.mod',
               'deterministic_simulations/multiple_lead_lags/sim_exo_lead_lag.mod',
               'deterministic_simulations/multiple_lead_lags/sim_endo_lead_lag_aux_vars.mod',
               'deterministic_simulations/multiple_lead_lags/sim_endo_lead_lag.mod',
               'deterministic_simulations/multiple_lead_lags/sim_lead_lag_aux_vars.mod',
               'deterministic_simulations/multiple_lead_lags/sim_lead_lag.mod' ] },
  { 'test' : [ 'deterministic_simulations/ramst_block_mfs1.mod' ] },
  { 'test' : [ 'deterministic_simulations/linear_approximation/sw.mod' ] },
  { 'test' : [ 'deterministic_simulations/multiple_lead_lags/AR2.mod' ] },
  { 'test' : [ 'deterministic_simulations/multiple_lead_lags/AR2_forward.mod' ] },
  { 'test' : [ 'deterministic_simulations/ramst.mod',
               'deterministic_simulations/multiple_lead_lags/ramst_augmented_histval.mod' ] },
  { 'test' : [ 'deterministic_simulations/rbc_det.mod',
               'deterministic_simulations/rbc_det_stack_solve_algo_7.mod',
               'deterministic_simulations/rbc_det_stack_solve_algo_7_exo_lag.mod',
               'deterministic_simulations/rbc_det_stack_solve_algo_7_exo_lead.mod' ] },
  { 'test' : [ 'deterministic_simulations/pfwee.mod' ],
    'extra' : [ 'deterministic_simulations/pfwee.csv' ] },
  { 'test' : [ 'deterministic_simulations/pfwee_constant_sim_length.mod' ],
    'extra' : [ 'deterministic_simulations/pfwee.csv' ] },
  { 'test' : [ 'deterministic_simulations/pfwee_learnt_in.mod' ],
    'extra' : [ 'deterministic_simulations/pfwee_learnt_in.inc' ] },
  { 'test' : [ 'deterministic_simulations/pfwee_learnt_in_bare_1st_period.mod' ],
    'extra' : [ 'deterministic_simulations/pfwee_learnt_in.inc' ] },
  { 'test' : [ 'deterministic_simulations/pfwee_learnt_in_dseries.mod' ],
    'extra' : [ 'deterministic_simulations/pfwee_learnt_in.inc' ] },
  { 'test' : [ 'deterministic_simulations/pfwee_learnt_in_dseries_bare_1st_period.mod' ],
    'extra' : [ 'deterministic_simulations/pfwee_learnt_in.inc' ] },
  { 'test' : [ 'deterministic_simulations/pfwee_multiple_shocks.mod' ] },
  { 'test' : [ 'deterministic_simulations/pfwee_homotopy_linearization.mod' ],
    'extra' : [ 'deterministic_simulations/pfwee_homotopy.inc' ] },
  { 'test' : [ 'deterministic_simulations/pfwee_homotopy_marginal_linearization.mod' ],
    'extra' : [ 'deterministic_simulations/pfwee_homotopy.inc' ] },
  { 'test' : [ 'deterministic_simulations/pfwee_homotopy_linearization_condpaths.mod' ],
    'extra' : [ 'deterministic_simulations/pfwee_homotopy.inc' ] },
  { 'test' : [ 'deterministic_simulations/pfwee_homotopy_marglin_condpaths.mod' ],
    'extra' : [ 'deterministic_simulations/pfwee_homotopy.inc' ] },
  { 'test' : [ 'deterministic_simulations/pfwee_marginal_linearization.mod' ] },
  { 'test' : [ 'deterministic_simulations/ramst_shock_paths.mod' ] },
  { 'test' : [ 'deterministic_simulations/ramst_shock_paths_pfwee.mod' ] },
  { 'test' : [ 'deterministic_simulations/ramst_shock_paths_controlled.mod' ] },
  { 'test' : [ 'deterministic_simulations/lbj/rbc.mod' ] },
  { 'test' : [ 'lmmcp/rbc.mod' ] },
  { 'test' : [ 'lmmcp/sw_lmmcp.mod',
               'lmmcp/sw_newton.mod' ],
    'extra' : [ 'lmmcp/sw-common-header.inc',
                'lmmcp/sw-common-footer.inc' ] },
  { 'test' : [ 'trend_var/fs2000_nonstationary.mod' ] },
  { 'test' : [ 'trend_var/fs2000_log_nonstationary.mod' ] },
  { 'test' : [ 'decision_rules/third_order/FV2011.mod' ],
    'extra' : [ 'decision_rules/third_order/comparison_policy_functions_dynare_mathematica.m',
                'decision_rules/third_order/FV_2011_policyfunctions.mat' ] },
  { 'test' : [ 'decision_rules/example1.mod',
               'decision_rules/example1_use_dll.mod' ],
    'extra' : [ 'decision_rules/example1_results_dyn_432.mat' ] },
  { 'test' : [ 'decision_rules/ar_qz_test.mod' ] },
  { 'test' : [ 'shock_decomposition/example1_calib_shock_decomp.mod' ],
    'extra' : [ 'shock_decomposition/example1_calib_shock_decomp_data.mat' ] },
  { 'test' : [ 'shock_decomposition/fs2000_est.mod' ],
    'extra' : [ 'shock_decomposition/fsdat_simul.m',
                'shock_decomposition/fs2000.inc' ] },
  { 'test' : [ 'shock_decomposition/fs2000_est_varlist.mod' ],
    'extra' : [ 'shock_decomposition/fsdat_simul.m',
                'shock_decomposition/fs2000.inc' ] },
  { 'test' : [ 'shock_decomposition/fs2000_cal_groups.mod' ],
    'extra' : [ 'shock_decomposition/fsdat_simul.m',
                'shock_decomposition/fs2000.inc' ] },
  { 'test' : [ 'shock_decomposition/fs2000_forecasts.mod' ],
    'extra' : [ 'shock_decomposition/fsdat_simul.m',
                'shock_decomposition/fs2000.inc' ] },
  { 'test' : [ 'shock_decomposition/ls2003_plot.mod' ],
    'extra' : [ 'ls2003/data_ca1.m' ] },
  { 'test' : [ 'shock_decomposition/shock_decomp_backward.mod' ] },
  { 'test' : [ 'stochastic_purely_forward/stochastic_purely_forward.mod' ] },
  { 'test' : [ 'stochastic_purely_forward/stochastic_purely_forward_with_static.mod' ] },
  { 'test' : [ 'forecast/forecast_test_AR1.mod' ] },
  { 'test' : [ 'forecast/forecast_test_AR1_order_2.mod' ] },
  { 'test' : [ 'forecast/forecast_test_AR1_varexo_det.mod' ] },
  { 'test' : [ 'forecast/fs2000_forecast.mod' ] },
  { 'test' : [ 'forecast/fs2000_forecast_varexo_det.mod' ] },
  { 'test' : [ 'forecast/forecast_test_AR1_varexo_det_order_2.mod' ] },
  { 'test' : [ 'forecast/fs2000_forecast.mod' ] },
  { 'test' : [ 'forecast/fs2000_forecast_varexo_det.mod' ] },
  { 'test' : [ 'forecast/example1_varexo_det.mod' ] },
  { 'test' : [ 'forecast/Hansen_exo_det_forecast.mod' ] },
  { 'test' : [ 'forecast/rbc_basic_forecast_varexo_det.mod' ] },
  { 'test' : [ 'forecast/linear_exo_det_forecast.mod' ] },
  { 'test' : [ 'forecast/ls2003_rolling_window_forecast.mod' ],
    'extra' : [ 'recursive/data_ca1.m' ] },
  { 'test' : [ 'gradient/fs2000_numgrad_13.mod' ],
    'extra' : [ 'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'gradient/fs2000_numgrad_15.mod' ],
    'extra' : [ 'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'gradient/fs2000_numgrad_2.mod' ],
    'extra' : [ 'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'gradient/fs2000_numgrad_3.mod' ],
    'extra' : [ 'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'gradient/fs2000_numgrad_5.mod' ],
    'extra' : [ 'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'filter_step_ahead/fs2000_filter_step_ahead_bayesian.mod' ],
    'extra' : [ 'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'filter_step_ahead/fs2000_filter_step_ahead_ML.mod' ],
    'extra' : [ 'fs2000/fsdat_simul.m' ] },
  { 'test' : [ 'filter_step_ahead/trend_cycle_decomposition.mod' ],
    'extra' : [ 'filter_step_ahead/trend_cycle_decomposition_data.m' ] },
  { 'test' : [ 'loglinear/example4_exp.mod',
               'loglinear/example4_loglinear.mod' ] },
  { 'test' : [ 'loglinear/example4_loglinear_lagged_exogenous.mod',
               'loglinear/example4_loglinear_lagged_exogenous_a.mod' ] },
  { 'test' : [ 'loglinear/example4_exp_histval.mod',
               'loglinear/example4_loglinear_histval.mod' ] },
  { 'test' : [ 'smoother2histval/fs2000_smooth.mod',
               'smoother2histval/fs2000_simul.mod' ],
    'extra' : [ 'smoother2histval/fsdat_simul.m' ] },
  { 'test' : [ 'smoother2histval/fs2000_smooth_stoch_simul.mod' ],
    'extra' : [ 'smoother2histval/fsdat_simul.m' ] },
  { 'test' : [ 'smoother2histval/fs2000_smooth_forecast.mod' ],
    'extra' : [ 'smoother2histval/fsdat_simul.m' ] },
  { 'test' : [ 'smoother2histval/fs2000_calib_smoother.mod' ],
    'extra' : [ 'smoother2histval/fsdat_simul.m' ] },
  { 'test' : [ 'epilogue/example1.mod' ] },
  { 'test' : [ 'differentiate_forward_vars/RBC_differentiate_forward.mod' ] },
  { 'test' : [ 'partitioning/partitions.mod' ] },
  { 'test' : [ 'prior_posterior_function/fs2000_prior_posterior_function.mod' ],
    'extra' : [ 'fs2000/fsdat_simul.m',
                'prior_posterior_function/posterior_function_demo.m' ] },
  { 'test' : [ 'observation_trends_and_prefiltering/ML/Trend_loglinear_no_prefilter.mod',
               'observation_trends_and_prefiltering/ML/Trend_loglinear_no_prefilter_first_obs.mod' ],
    'extra' : [ 'observation_trends_and_prefiltering/Trend_exp_model_no_prefilter_common.inc',
                'observation_trends_and_prefiltering/Trend_load_data_common.inc',
                'observation_trends_and_prefiltering/Trend_diagnostics_ML_common.inc',
                'observation_trends_and_prefiltering/Trend_no_prefilter_conditional_forecast.inc',
                'observation_trends_and_prefiltering/generate_trend_stationary_AR1.m' ] },
  { 'test' : [ 'observation_trends_and_prefiltering/ML/Trend_loglinear_prefilter.mod' ],
    'extra' : [ 'observation_trends_and_prefiltering/Trend_exp_model_prefilter_common.inc',
                'observation_trends_and_prefiltering/Trend_load_data_common.inc',
                'observation_trends_and_prefiltering/Trend_diagnostics_ML_common.inc',
                'observation_trends_and_prefiltering/generate_trend_stationary_AR1.m' ] },
  { 'test' : [ 'observation_trends_and_prefiltering/ML/Trend_loglinear_prefilter_first_obs.mod' ],
    'extra' : [ 'observation_trends_and_prefiltering/Trend_exp_model_prefilter_common.inc',
                'observation_trends_and_prefiltering/Trend_load_data_common.inc',
                'observation_trends_and_prefiltering/Trend_diagnostics_ML_common.inc',
                'observation_trends_and_prefiltering/generate_trend_stationary_AR1.m' ] },
  { 'test' : [ 'observation_trends_and_prefiltering/ML/Trend_no_prefilter.mod',
               'observation_trends_and_prefiltering/ML/Trend_no_prefilter_first_obs.mod' ],
    'extra' : [ 'observation_trends_and_prefiltering/Trend_model_no_prefilter_common.inc',
                'observation_trends_and_prefiltering/Trend_load_data_common.inc',
                'observation_trends_and_prefiltering/Trend_diagnostics_ML_common.inc',
                'observation_trends_and_prefiltering/Trend_no_prefilter_conditional_forecast.inc',
                'observation_trends_and_prefiltering/generate_trend_stationary_AR1.m' ] },
  { 'test' : [ 'observation_trends_and_prefiltering/ML/Trend_prefilter.mod' ],
    'extra' : [ 'observation_trends_and_prefiltering/Trend_model_prefilter_common.inc',
                'observation_trends_and_prefiltering/Trend_load_data_common.inc',
                'observation_trends_and_prefiltering/Trend_diagnostics_ML_common.inc',
                'observation_trends_and_prefiltering/generate_trend_stationary_AR1.m' ] },
  { 'test' : [ 'observation_trends_and_prefiltering/ML/Trend_prefilter_first_obs.mod' ],
    'extra' : [ 'observation_trends_and_prefiltering/Trend_model_prefilter_common.inc',
                'observation_trends_and_prefiltering/Trend_load_data_common.inc',
                'observation_trends_and_prefiltering/Trend_diagnostics_ML_common.inc',
                'observation_trends_and_prefiltering/generate_trend_stationary_AR1.m' ] },
  { 'test' : [ 'observation_trends_and_prefiltering/ML/Trend_no_prefilter_selected_var.mod' ],
    'extra' : [ 'observation_trends_and_prefiltering/generate_trend_stationary_AR1.m' ] },
  { 'test' : [ 'observation_trends_and_prefiltering/calib_smoother/Tr_no_prefilter_calib_smoother.mod',
               'observation_trends_and_prefiltering/calib_smoother/Tr_no_prefil_f_obs_loglin_cal_smoother.mod',
               'observation_trends_and_prefiltering/calib_smoother/Tr_no_prefilter_loglin_calib_smoother.mod',
               'observation_trends_and_prefiltering/calib_smoother/Tr_no_prefilt_first_obs_cal_smooth.mod',
               'observation_trends_and_prefiltering/calib_smoother/Tr_prefil_f_obs_loglin_cal_smoother.mod',
               'observation_trends_and_prefiltering/calib_smoother/Tr_prefilter_calib_smoother.mod',
               'observation_trends_and_prefiltering/calib_smoother/Tr_prefilter_loglin_calib_smoother.mod',
               'observation_trends_and_prefiltering/calib_smoother/Tr_prefilt_first_obs_cal_smooth.mod' ],
    'extra' : [ 'observation_trends_and_prefiltering/Trend_exp_model_calib_no_prefilter_common.inc',
                'observation_trends_and_prefiltering/Trend_exp_model_calib_prefilter_common.inc',
                'observation_trends_and_prefiltering/Trend_model_calib_no_prefilter_common.inc',
                'observation_trends_and_prefiltering/Trend_model_calib_prefilter_common.inc',
                'observation_trends_and_prefiltering/Trend_load_data_common.inc',
                'observation_trends_and_prefiltering/Trend_diagnostics_calib_common.inc',
                'observation_trends_and_prefiltering/generate_trend_stationary_AR1.m' ] },
  { 'test' : [ 'bgp/solow-1/solow.mod' ] },
  { 'test' : [ 'bgp/nk-1/nk.mod' ] },
  { 'test' : [ 'bgp/ramsey-1/ramsey.mod' ] },
  { 'test' : [ 'model-inversion/bk-1/solow_ces.mod' ] },
  { 'test' : [ 'model-inversion/bk-2/solow_ces.mod' ] },
  { 'test' : [ 'model-inversion/bk-3a/varmodel.mod' ] },
  { 'test' : [ 'model-inversion/bk-3a-nlsolve/varmodel.mod' ] },
  { 'test' : [ 'model-inversion/bk-3b/varmodel.mod' ] },
  { 'test' : [ 'model-inversion/bk-3b-nlsolve/varmodel.mod' ] },
  { 'test' : [ 'model-inversion/bk-4/varmodel.mod' ] },
  { 'test' : [ 'model-inversion/bk-4-nlsolve/varmodel.mod' ] },
  { 'test' : [ 'model-inversion/bk-5/varmodel.mod' ] },
  { 'test' : [ 'model-inversion/bk-5-nlsolve/varmodel.mod' ] },
  { 'test' : [ 'model-inversion/bk-6/varmodel.mod' ] },
  { 'test' : [ 'model-inversion/bk-7/varmodel.mod' ] },
  { 'test' : [ 'model-inversion/bk-8/varmodel.mod' ] },
  { 'test' : [ 'model-inversion/bk-9/msm1.mod' ] },
  { 'test' : [ 'model-inversion/bk-9/msm2.mod' ] },
  { 'test' : [ 'model-inversion/bk-diff-1/varmodel.mod' ] },
  { 'test' : [ 'model-inversion/bk-diff-2/varmodel.mod' ] },
  { 'test' : [ 'model-inversion/bk-diff-3/varmodel.mod' ] },
  { 'test' : [ 'model-inversion/bk-diff-4/varmodel.mod' ] },
  { 'test' : [ 'model-inversion/bk-diff-5/varmodel.mod' ] },
  { 'test' : [ 'model-inversion/nk-1/simulate.mod',
               'model-inversion/nk-1/invert.mod' ] },
  { 'test' : [ 'model-inversion/nk-2/simulate.mod',
               'model-inversion/nk-2/invert.mod',
               'model-inversion/nk-2/z_check_inversion.mod' ] },
  { 'test' : [ 'dynare-command-options/ramst.mod' ] },
  { 'test' : [ 'log_transform/example1.mod' ] },
  { 'test' : [ 'deterministic_simulations/ramst.mod',
               'log_transform/ramst.mod' ] },
  { 'test' : [ 'log_transform/fs2000_nonstationary.mod' ] },
  { 'test' : [ 'log_transform/nk_ramsey.mod' ] },
  { 'test' : [ 'solve_algo_12_14/simul_backward_reference.mod',
               'solve_algo_12_14/simul_backward_12.mod',
               'solve_algo_12_14/simul_backward_14.mod' ],
    'extra' : [ 'solve_algo_12_14/simul_backward_common.inc',
                'solve_algo_12_14/backward_model.inc' ] },
  { 'test' : [ 'solve_algo_12_14/purely_backward_reference.mod',
               'solve_algo_12_14/purely_backward_12.mod',
               'solve_algo_12_14/purely_backward_14.mod' ],
    'extra' : [ 'solve_algo_12_14/purely_backward_common.inc',
                'solve_algo_12_14/backward_model.inc' ] },
  { 'test' : [ 'solve_algo_12_14/purely_static_reference.mod',
               'solve_algo_12_14/purely_static_12.mod',
               'solve_algo_12_14/purely_static_14.mod' ],
    'extra' : [ 'solve_algo_12_14/purely_static_common.inc' ] },
  { 'test' : [ 'solve_algo_12_14/purely_forward_reference.mod',
               'solve_algo_12_14/purely_forward_12.mod',
               'solve_algo_12_14/purely_forward_14.mod' ],
    'extra' : [ 'solve_algo_12_14/purely_forward_common.inc' ] },
  { 'test' : [ 'optimal_policy/Ramsey/ramsey_ex_initval.mod',
               'deprecated/ramsey_ex.mod' ] },
  { 'test' : [ 'deterministic_simulations/ramst.mod',
               'deprecated/ramst.mod' ] },
  { 'test' : [ 'deterministic_simulations/ramst_pf_controlled_paths.mod',
               'deterministic_simulations/ramst_pf_controlled_paths_lbj.mod',
               'deterministic_simulations/ramst_pf_controlled_paths_ssa7.mod',
               'deterministic_simulations/ramst_pf_controlled_paths_check.mod' ],
    'extra' : [ 'deterministic_simulations/ramst_pf_controlled_paths_common.inc'] },
  { 'test' : [ 'deterministic_simulations/ramst_pfwee_controlled_paths.mod' ],
    'extra' : [ 'deterministic_simulations/ramst_pf_controlled_paths_common.inc'] },

  # ECB modfiles
  { 'test' : [ 'var-expectations/1/example1.mod' ] },
  { 'test' : [ 'var-expectations/1-with-time-shift-a/example1.mod' ] },
  { 'test' : [ 'var-expectations/2/example1.mod' ] },
  { 'test' : [ 'var-expectations/2-with-time-shift/example1.mod' ] },
  { 'test' : [ 'var-expectations/3/example1.mod',
               'var-expectations/4/example1.mod' ] },
  { 'test' : [ 'var-expectations/3-with-time-shift/example1.mod',
               'var-expectations/4-with-time-shift/example1.mod' ] },
  { 'test' : [ 'var-expectations/5/example1.mod' ] },
  { 'test' : [ 'var-expectations/5-with-time-shift/example1.mod' ] },
  { 'test' : [ 'var-expectations/6/example1.mod',
               'var-expectations/6/substitution.mod' ] },
  { 'test' : [ 'var-expectations/6-with-time-shift/example1.mod',
               'var-expectations/6-with-time-shift/substitution.mod' ] },
  { 'test' : [ 'var-expectations/7/example1.mod',
               'var-expectations/7/substitution.mod' ] },
  { 'test' : [ 'var-expectations/8/example1.mod',
               'var-expectations/8/substitution.mod' ] },
  { 'test' : [ 'var-expectations/9/example1.mod' ] },
  { 'test' : [ 'var-expectations/10/example1.mod' ] },
  { 'test' : [ 'var-expectations/11/example1.mod' ] },
  { 'test' : [ 'var-expectations/12/example1.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/vm1.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/vm2.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/vm3.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/vm4.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/vm5.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/tcm1.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/tcm2.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/tcm3.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/tcm4.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/tcm5.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/tcm6.mod',
               'trend-component-and-var-models/tcm7.mod',
               'trend-component-and-var-models/tcm8.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/tcm9.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/tcm10.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/tcm11.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/tcm12.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/legacy/vm1.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/legacy/vm2.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/legacy/vm3.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/legacy/tcm1.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/legacy/tcm2.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/legacy/tcm3.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/legacy/tcm4.mod' ] },
  { 'test' : [ 'trend-component-and-var-models/legacy/tcm5.mod' ] },
  { 'test' : [ 'pac/var-1/example1.mod' ] },
  { 'test' : [ 'pac/var-2/example1.mod' ] },
  { 'test' : [ 'pac/var-3/example1.mod' ] },
  { 'test' : [ 'pac/var-4/example1.mod' ] },
  { 'test' : [ 'pac/var-5/example1.mod',
               'pac/var-5/substitution.mod' ] },
  { 'test' : [ 'pac/var-6/example1.mod',
               'pac/var-6/substitution.mod' ] },
  { 'test' : [ 'pac/var-7/example1.mod',
               'pac/var-7/substitution.mod' ] },
  { 'test' : [ 'pac/var-8/example1.mod',
               'pac/var-8/substitution.mod' ] },
  { 'test' : [ 'pac/var-8e/example1.mod' ] },
  { 'test' : [ 'pac/var-9/example1.mod',
               'pac/var-9/substitution.mod' ] },
  { 'test' : [ 'pac/var-9e/example1.mod' ] },
  { 'test' : [ 'pac/var-10e/example1.mod',
               'pac/var-10e/example2.mod' ] },
  { 'test' : [ 'pac/var-11e/example1.mod' ] },
  { 'test' : [ 'pac/var-12/example1.mod' ] },
  { 'test' : [ 'pac/var-12/example2.mod' ] },
  { 'test' : [ 'pac/var-12/example3.mod' ] },
  { 'test' : [ 'pac/var-12/example4.mod' ] },
  { 'test' : [ 'pac/var-12/example5.mod' ] },
  { 'test' : [ 'pac/var-12/example11.mod' ] },
  { 'test' : [ 'pac/var-12/example12.mod' ] },
  { 'test' : [ 'pac/trend-component-1/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-2/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-3/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-4/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-5/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-6/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-7/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-9/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-10/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-11/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-12/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-13a/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-13b/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-14/example1.mod',
               'pac/trend-component-14/substitution.mod' ] },
  { 'test' : [ 'pac/trend-component-15/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-16/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-17/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-18/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-19/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-19/example2.mod' ] },
  { 'test' : [ 'pac/trend-component-19/example3.mod' ] },
  { 'test' : [ 'pac/trend-component-19-growth-lin-comb/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-1-mce/example_det.mod' ] },
  { 'test' : [ 'pac/trend-component-1-mce/example_sto.mod' ] },
  { 'test' : [ 'pac/trend-component-2-mce/example_det.mod' ] },
  { 'test' : [ 'pac/trend-component-2-mce/example_sto.mod' ] },
  { 'test' : [ 'pac/trend-component-20-1/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-20-2/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-20-3/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-20-4/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-21/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-22/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-23/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-24/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-25/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-26/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-27/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-28/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-28/example2.mod' ] },
  { 'test' : [ 'pac/trend-component-28/example3.mod' ] },
  { 'test' : [ 'pac/trend-component-28/example4.mod' ] },
  { 'test' : [ 'pac/trend-component-28/example5.mod' ] },
  { 'test' : [ 'pac/trend-component-28/example6.mod' ] },
  { 'test' : [ 'pac/trend-component-29/example1.mod' ] },
  { 'test' : [ 'pac/trend-component-29/example2.mod' ] },
  { 'test' : [ 'pac/trend-component-30/example1.mod' ] },
  { 'test' : [ 'write/example1.mod' ] },
  { 'test' : [ 'ecb/backward-models/irf/solow_1.mod' ] },
  { 'test' : [ 'ecb/backward-models/irf/solow_2.mod' ] },
  { 'test' : [ 'ecb/cherrypick/test1.mod' ] },
  { 'test' : [ 'ecb/cherrypick/test2.mod' ] },
  { 'test' : [ 'ecb/contribution-plots/contrib1.mod' ] },
  { 'test' : [ 'ecb/contribution-plots/contrib2.mod' ] },

  # OLS modfiles
  { 'test' : [ 'estimation/univariate/bayesian.mod' ],
    'extra' : [ 'estimation/univariate/data.csv' ] },
  { 'test' : [ 'estimation/univariate/bayesian_param_names.mod' ],
    'extra' : [ 'estimation/univariate/data.csv' ] },
  { 'test' : [ 'ecb/SURGibbs/fulton_fish.mod' ],
    'extra' : [ 'ecb/SURGibbs/fishdata.csv' ] },
  { 'test' : [ 'estimation/univariate/ols/ols_wc_1.mod' ],
    'extra' : [ 'estimation/univariate/ols/mc-ols.inc' ] },
  { 'test' : [ 'estimation/univariate/ols/ols_wc_2.mod' ],
    'extra' : [ 'estimation/univariate/ols/mc-ols.inc' ] },
  { 'test' : [ 'estimation/univariate/ols/ols_wc_3.mod' ],
    'extra' : [ 'estimation/univariate/ols/mc-ols.inc' ] },
  { 'test' : [ 'estimation/univariate/ols/ols_date_range.mod' ],
    'extra' : [ 'estimation/univariate/ols/mc-ols.inc' ] },
  { 'test' : [ 'estimation/univariate/ols/ols_param_names.mod' ],
    'extra' : [ 'estimation/univariate/ols/mc-ols.inc' ] },
  { 'test' : [ 'estimation/univariate/ols/ols_base.mod' ],
    'extra' : [ 'estimation/univariate/ols/mc-ols.inc' ] },
  { 'test' : [ 'ecb/SUR/sur_noniterative.mod' ] },
  { 'test' : [ 'ecb/SUR/sur_params_noniterative.mod' ] },
  { 'test' : [ 'ecb/SUR/panel_var_diff_NB_simulation_test.mod' ] },
  { 'test' : [ 'ecb/SUR/panel_var_diff_NB_simulation_zero_eq.mod' ] },
  { 'test' : [ 'ecb/pooled_ols/panel_var_diff_NB_simulation_test.mod' ] },
  { 'test' : [ 'ecb/pooled_ols/test_param_names.mod' ] },
  { 'test' : [ 'ecb/pooled_fgls/panel_var_diff_NB_simulation_test.mod' ] },
  { 'test' : [ 'ecb/pooled_fgls/test_param_names.mod' ] },

  # Heterogeneity
  { 'test' : [ 'heterogeneity/ks.mod' ],
    'extra': [ 'heterogeneity/ks.mat',
               'heterogeneity/expect_error.m',
               'heterogeneity/print_test_summary.m',
               'heterogeneity/run_test.m',
               'heterogeneity/test_solve_functionality.m',
               'heterogeneity/test_simulate_stochastic.m',
               'heterogeneity/test_check_steady_state.m',
               'heterogeneity/test_field_validation.m',
               'heterogeneity/test_permutation.m' ] },
  { 'test' : [ 'heterogeneity/hank_1a.mod' ],
    'extra': [ 'heterogeneity/hank_1a.mat',
               'heterogeneity/print_test_summary.m',
               'heterogeneity/run_test.m',
               'heterogeneity/test_simulate_stochastic.m' ] },
  { 'test' : [ 'heterogeneity/hank_1a.mod' ],
    'dynare_args': [ '-DDET_SHOCKS=1' ],
    'extra': [ 'heterogeneity/hank_1a.mat',
               'heterogeneity/print_test_summary.m',
               'heterogeneity/run_test.m',
               'heterogeneity/test_simulate_news_shocks.m' ] },
  { 'test' : [ 'heterogeneity/hank_2a.mod' ],
    'extra': [ 'heterogeneity/hank_2a.mat',
               'heterogeneity/print_test_summary.m',
               'heterogeneity/run_test.m',
               'heterogeneity/test_simulate_stochastic.m',
               'heterogeneity/test_permutation.m' ] },
  { 'test' : [ 'heterogeneity/ks_commands.mod' ],
    'extra' : [ 'heterogeneity/ks.mat' ] },
  { 'test' : [ 'heterogeneity/ks_commands.mod' ],
    'dynare_args': [ '-DDET_SHOCKS=1' ],
    'extra' : [ 'heterogeneity/ks.mat' ] },
  { 'test' : [ 'heterogeneity/ks_news_shocks_stochastic_options_xfail.mod' ],
    'dynare_args': [ '-DTEST_OPTION="periods"' ],
    'extra' : [ 'heterogeneity/ks.mat' ],
    'should_fail' : true },
  { 'test' : [ 'heterogeneity/ks_news_shocks_stochastic_options_xfail.mod' ],
    'dynare_args': [ '-DTEST_OPTION="irf"' ],
    'extra' : [ 'heterogeneity/ks.mat' ],
    'should_fail' : true },
  { 'test' : [ 'heterogeneity/ks_news_shocks_stochastic_options_xfail.mod' ],
    'dynare_args': [ '-DTEST_OPTION="irf_shocks"' ],
    'extra' : [ 'heterogeneity/ks.mat' ],
    'should_fail' : true },
  { 'test' : [ 'heterogeneity/ks_news_shocks_stochastic_options_xfail.mod' ],
    'dynare_args': [ '-DTEST_OPTION="relative_irf"' ],
    'extra' : [ 'heterogeneity/ks.mat' ],
    'should_fail' : true },
  { 'test' : [ 'heterogeneity/ks_news_shocks_stochastic_options_xfail.mod' ],
    'dynare_args': [ '-DTEST_OPTION="combined"' ],
    'extra' : [ 'heterogeneity/ks.mat' ],
    'should_fail' : true },
  { 'test' : [ 'heterogeneity/nonseparable_het_lead_lag_xfail.mod' ],
    'should_fail' : true },
  { 'test' : [ 'heterogeneity/het_exo_lag_xfail.mod' ],
    'should_fail' : true },
  { 'test' : [ 'heterogeneity/het_exo_lead_xfail.mod' ],
    'should_fail' : true },
  { 'test' : [ 'heterogeneity/het_endo_lag_deep_xfail.mod' ],
    'should_fail' : true },
  { 'test' : [ 'heterogeneity/het_endo_lead_deep_xfail.mod' ],
    'should_fail' : true },
  { 'test' : [ 'heterogeneity/analytical_rbc.mod' ] },
  { 'test' : [ 'heterogeneity/analytical_rbc_2capital.mod' ] },
  { 'test' : [ 'heterogeneity/ks_two_shocks.mod' ],
    'extra': [ 'heterogeneity/ks_two_shocks.mat' ] },
  { 'test' : [ 'heterogeneity/ks_steady_state.mod' ],
    'extra': [ 'heterogeneity/ks_sp.mat' ] },
  { 'test' : [ 'heterogeneity/ks_steady_state_dense.mod' ],
    'extra': [ 'heterogeneity/ks_sp_dense.mat' ] },
  { 'test' : [ 'heterogeneity/hank_1a_steady_state.mod' ],
    'extra': [ 'heterogeneity/hank_1a_sp.mat' ] },
  { 'test' : [ 'heterogeneity/hank_2a_steady_state.mod' ],
    'extra': [ 'heterogeneity/hank_2a_sp.mat' ] },

  # Heterogenous agents tests requiring large file downloads
  # They are excluded from the default test setup, see add_test_setup() call below
  { 'test' : [ 'heterogeneity/large_file_download/ks_ssj.mod' ],
    'extra': [ 'heterogeneity/large_file_download/print_test_summary.m',
               'heterogeneity/large_file_download/run_test.m',
               'heterogeneity/large_file_download/test_solve_with_custom_impulse_responses.m' ] },
  { 'test' : [ 'heterogeneity/large_file_download/hank_1a_ssj.mod' ],
    'extra': [ 'heterogeneity/large_file_download/print_test_summary.m',
               'heterogeneity/large_file_download/run_test.m',
               'heterogeneity/large_file_download/test_solve_functionality.m',
               'heterogeneity/large_file_download/test_solve_with_custom_impulse_responses.m',
               'heterogeneity/large_file_download/test_simulate_stochastic.m' ] },
  { 'test' : [ 'heterogeneity/large_file_download/hank_2a_ssj.mod' ],
    'extra': [ 'heterogeneity/large_file_download/print_test_summary.m',
               'heterogeneity/large_file_download/run_test.m',
               'heterogeneity/large_file_download/test_solve_functionality.m',
               'heterogeneity/large_file_download/test_solve_with_custom_impulse_responses.m' ] },

  # Expected failures
  { 'test' : [ 'preprocessor_checks/ramst_all_values_required_xfail.mod' ],
    'should_fail' : true },
  { 'test' : [ 'preprocessor_checks/example1_undeclared_vars_xfail.mod' ],
    'should_fail' : true },
  { 'test' : [ 'estimation/estim_param_in_shock_value_xfail.mod' ],
    'extra' : [ 'estimation/fsdat_simul.m' ],
    'should_fail' : true },
  { 'test' : [ 'optimal_policy/Ramsey/ramsey_ex_wrong_ss_file_xfail.mod' ],
    'extra' : [ 'optimal_policy/Ramsey/find_c.m' ],
    'should_fail' : true },
  { 'test' : [ 'estimation/fs2000_mixed_ML_xfail.mod' ],
    'extra' : [ 'estimation/fsdat_simul.m' ],
    'should_fail' : true },
  { 'test' : [ 'estimation/fs2000_stochastic_singularity_xfail.mod' ],
    'extra' : [ 'estimation/fsdat_simul.m' ],
    'should_fail' : true },
  { 'test' : [ 'estimation/heteroskedastic_shocks/fs2000_het_XFAIL.mod' ],
    'extra' : [ 'estimation/heteroskedastic_shocks/fs2000_het_model.inc',
                'estimation/fsdat_simul.m' ],
    'should_fail' : true },
  { 'test' : [ 'estimation/heteroskedastic_shocks/fs2000_het_XFAIL2.mod' ],
    'extra' : [ 'estimation/heteroskedastic_shocks/fs2000_het_model.inc',
                'estimation/fsdat_simul.m' ],
    'should_fail' : true },
  { 'test' : [ 'identification/ident_unit_root/ident_unit_root_xfail.mod' ],
    'should_fail' : true },
  { 'test' : [ 'identification/LindeTrabandt/LindeTrabandt2019_xfail.mod' ],
    'should_fail' : true },
  { 'test' : [ 'steady_state/Linear_steady_state_xfail.mod' ],
    'should_fail' : true },
  { 'test' : [ 'optimal_policy/Ramsey/ramsey_histval_xfail.mod' ],
    'extra' : [ 'optimal_policy/Ramsey/find_c.m' ],
    'should_fail' : true },
  { 'test' : [ 'occbin/model_irrcap_twoconstraints/dynrbc_token_xfail.mod' ],
    'should_fail' : true },
  { 'test' : [ 'particle/first_spec_xfail_0.mod' ],
    'extra' : [ 'particle/first_spec_common.inc' ],
    'should_fail' : true },
  { 'test' : [ 'particle/first_spec_xfail_1.mod' ],
    'extra' : [ 'particle/first_spec_common.inc' ],
    'should_fail' : true },
  { 'test' : [ 'kalman_initial_state/fs2000_kalman_initial_xfail.mod' ],
    'extra' : [ 'kalman_initial_state/fs2000_common.inc',
                'fs2000/fsdat_simul.m' ],
    'should_fail' : true },
  { 'test' : [ 'preprocessor_checks/example1_extra_exo_xfail.mod' ],
    'should_fail' : true },
  { 'test' : [ 'estimation/tune_mh_jscale/fs2000_1_xfail.mod' ],
    'extra' : [ 'estimation/tune_mh_jscale/fs2000.inc',
                'estimation/fsdat_simul.m' ],
    'should_fail' : true },
  { 'test' : [ 'estimation/tune_mh_jscale/fs2000_2_xfail.mod' ],
    'extra' : [ 'estimation/tune_mh_jscale/fs2000.inc',
                'estimation/fsdat_simul.m' ],
    'should_fail' : true },
  { 'test' : [ 'estimation/no_init_estimation_check_first_obs/fs2000_init_check_XFAIL.mod' ],
    'extra' : [ 'estimation/no_init_estimation_check_first_obs/fsdat_mat.m',
                'estimation/no_init_estimation_check_first_obs/fsdat_simul.m' ],
    'should_fail' : true },
  { 'test' : [ 'estimation/no_init_estimation_check_first_obs/fs2000_init_check_XFAIL2.mod' ],
    'extra' : [ 'estimation/no_init_estimation_check_first_obs/fsdat_mat_XFAIL.m',
                'estimation/no_init_estimation_check_first_obs/fsdat_simul.m' ],
    'should_fail' : true },
  { 'test' : [ 'var-expectations/1-with-time-shift-b/example1.mod' ],
    'should_fail' : true },
  { 'test' : [ 'estimation/univariate/nls/staticmodelx.mod' ],
    'should_fail' : true },

  # Misc tests with .m files
  { 'test' : [ 'run_all_unit_tests.m' ] },
  { 'test' : [ 'histval_initval_file_unit_tests.m' ],
    'extra' : [ 'histval_initval_file/my_assert.m' ] },
  { 'test' : [ 'run_kronecker_tests.m' ],
    'extra' : [ 'kronecker/test_kron.m',
                'kronecker/nash_matrices.mat' ] },
  { 'test' : [ 'nonlinearsolvers.m' ],
    'extra' : [ 'solver-test-functions/brown.m',
                'solver-test-functions/broydenbanded.m',
                'solver-test-functions/broydentridiagonal.m',
                'solver-test-functions/chebyquad.m',
                'solver-test-functions/discreteboundaryvalue.m',
                'solver-test-functions/discreteintegralequation.m',
                'solver-test-functions/helicalvalley.m',
                'solver-test-functions/powell1.m',
                'solver-test-functions/powell2.m',
                'solver-test-functions/rosenbrock.m',
                'solver-test-functions/trigonometric.m',
                'solver-test-functions/variablydimensioned.m',
                'solver-test-functions/watson.m',
                'solver-test-functions/wood.m',] },
  { 'test' : [ 'cyclereduction.m' ] },
  { 'test' : [ 'logarithmicreduction.m' ] },
  { 'test' : [ 'kalman/likelihood/test_kalman_mex.m' ] },
  { 'test' : [ 'contribs.m' ],
    'extra' : [ 'sandbox.mod',
                'simulateddata.m' ] },
  { 'test' : [ 'run_block_bytecode_tests.m' ],
    'extra' : [ 'block_bytecode/run_ls2003.m',
                'block_bytecode/ls2003.mod' ],
    'env' : { 'HAS_PARU' : (get_option('suitesparse_src_path') != '').to_string(),
              # Only enable PARDISO test if presence of the library is explicitly required,
              # to avoid cases where the library might be there but no license is configured.
              'HAS_PARDISO' : get_option('pardiso').enabled().to_string() } },
  { 'test' : [ 'run_reporting_tests.m' ],
    'extra' : [ 'reporting/AnnualTable.m',
                'reporting/CommResidTablePage.m',
                'reporting/CountryGraphPage.m',
                'reporting/CountryTablePage.m',
                'reporting/db_a.csv',
                'reporting/db_q.csv',
                'reporting/dc_a.csv',
                'reporting/dc_q.csv',
                'reporting/ResidTablePage.m',
                'reporting/runDynareReport.m' ] },

  # Examples — Stochastic simulations
  { 'test' : [ '../examples/stochastic_simulations/aguiar_gopinath_2007_trend.mod' ] },
  { 'test' : [ '../examples/stochastic_simulations/collard_2001_theoretical_moments.mod' ] },
  { 'test' : [ '../examples/stochastic_simulations/collard_2001_simulated_moments.mod' ] },
  { 'test' : [ '../examples/stochastic_simulations/collard_2001_analytical_steady_state.mod' ],
    'extra' : [ '../examples/stochastic_simulations/collard_2001_analytical_steady_state_helper.m' ] },
  { 'test' : [ '../examples/stochastic_simulations/schorfheide_2000_nonstationary.mod' ] },
  { 'test' : [ '../examples/stochastic_simulations/nk_baseline.mod' ],
    'extra' : [ '../examples/stochastic_simulations/nk_baseline_steadystate.m' ] },

  # Examples — Estimation
  { 'test' : [ '../examples/estimation/schorfheide_2000.mod' ],
    'extra' : [ '../examples/estimation/schorfheide_2000_data.m' ] },
  { 'test' : [ '../examples/estimation/gali_2015.mod' ],
    'extra' : [ '../examples/estimation/gali_2015_phillips_curve_slope.m',
                '../examples/estimation/gali_2015_prior_restrictions.m' ] },
  { 'test' : [ '../examples/estimation/rbc_irf_matching.mod' ],
    'extra' : [ '../examples/estimation/rbc_irf_matching_data.csv',
                '../examples/estimation/rbc_irf_matching_transformations.m' ] },

  # Examples — Perfect foresight
  { 'test' : [ '../examples/perfect_foresight/perfect_foresight_rbc.mod' ] },
  { 'test' : [ '../examples/perfect_foresight/perfect_foresight_expectation_errors.mod' ] },

  # Examples — Optimal policy
  { 'test' : [ '../examples/optimal_policy/nk_ramsey_osr.mod' ] },
  { 'test' : [ '../examples/optimal_policy/nk_ramsey_steady_file.mod' ],
    'extra' : [ '../examples/optimal_policy/nk_ramsey_steady_file_steadystate.m' ] },

  # Examples — Heterogeneity
  { 'test' : [ '../examples/heterogeneity/hank_one_asset.mod' ],
    'extra' : [ '../examples/heterogeneity/hank_one_asset.mat' ] },
  { 'test' : [ '../examples/heterogeneity/hank_one_asset_steady_state.mod' ],
    'extra' : [ '../examples/heterogeneity/hank_one_asset_sp.mat' ] },
  { 'test' : [ '../examples/heterogeneity/hank_two_assets.mod' ],
    'extra' : [ '../examples/heterogeneity/hank_two_assets.mat' ] },
  { 'test' : [ '../examples/heterogeneity/hank_two_assets_steady_state.mod' ],
    'extra' : [ '../examples/heterogeneity/hank_two_assets_sp.mat' ] },
  { 'test' : [ '../examples/heterogeneity/krusell_smith_1998.mod' ],
    'extra' : [ '../examples/heterogeneity/krusell_smith_1998.mat' ] },
  { 'test' : [ '../examples/heterogeneity/krusell_smith_1998_steady_state.mod' ],
    'extra' : [ '../examples/heterogeneity/krusell_smith_1998_sp.mat' ] },

  # Examples — OccBin
  { 'test' : [ '../examples/occbin/rbc_occbin.mod' ] },

  # Examples — Macroprocessor
  { 'test' : [ '../examples/macroprocessor/bkk_1992.mod' ] },

  # Examples — Reporting
  { 'test' : [ '../examples/reporting/collard_2001_reporting.mod' ] },

  # Examples — Semistructural
  { 'test' : [ '../examples/semistructural/pac_model.mod' ] }
]

base_test_driver_args = [ get_option('build_for') ]
if get_option('build_for') == 'matlab'
  base_test_driver_args += [ matlab_exe.full_path(), matlab_version, matlab_arch ]
else
  base_test_driver_args += [ octave_exe.full_path(), octave_version, '' ]
endif
base_test_driver_args += [ meson.current_source_dir(), meson.current_build_dir() ]

test_driver_exe = files('scripts/test-driver')

foreach t : mod_and_m_tests
  test_driver_args = base_test_driver_args + t['test']
  if t.has_key('extra')
    test_driver_args += [ '--' ] + t['extra']
  endif
  test_name = t['test'][-1] # Use the last test file as the test name

  test_env = {}

  # Handle optional dynare_args field for passing -D options to preprocessor
  if t.has_key('dynare_args')
    # Pass each argument as a separate environment variable
    args_list = t['dynare_args']
    test_env += { 'dynare_arg_count': args_list.length().to_string() }
    foreach i : range(args_list.length())
      test_env += { f'dynare_arg_@i@': args_list[i] }
    endforeach
  endif

  # Handle optional environment variables
  if t.has_key('env')
    test_env += t['env']
  endif

  # Construct the suite name(s) using the path to the test file
  # If test file is a/b/c/test.mod, then the suites will be:
  # [ 'a', 'a/b', 'a/b/c' ]
  test_suite = []
  foreach it : t['test']
    it_path_components = it.split('/')
    suite_name = ''
    foreach i : range(it_path_components.length()-1)
      if suite_name != ''
        suite_name += '/'
      endif
      suite_name += it_path_components[i]
      if not test_suite.contains(suite_name)
        test_suite += suite_name
      endif

      # Workaround for https://github.com/mesonbuild/meson/issues/11359
      # Do not add the subsuites of that specific suite, to avoid requiring too
      # much horizontal width of the terminal screen
      if suite_name == 'observation_trends_and_prefiltering'
        break
      endif
    endforeach
  endforeach

  test(test_name, test_driver_exe, args : test_driver_args, suite : test_suite,
       env : test_env, should_fail : t.get('should_fail', false), timeout : 0)
endforeach

# Exclude some tests from the default test setup
# They can still be run manually using the --suite option
add_test_setup('default', is_default : true, exclude_suites : [ 'heterogeneity/large_file_download' ])

### Developper stuff

## Tag file for Emacs (created in *source* directory)

git_exe = find_program('git', required : false)
etags_exe = find_program('etags', required : false)

fs = import('fs')

if fs.is_dir('.git') and git_exe.found() and etags_exe.found()
  all_files = run_command(git_exe,
                          [ '--git-dir=@0@/.git'.format(meson.project_source_root()),
                            'ls-files', '--recurse-submodules',
                            ':/*.cc', ':/*.hh', ':/*.[fF]08', ':/*.[ch]' ],
                          check : true)
  all_files = files(all_files.stdout().split())

  custom_target('tags', output : 'tags', # Dummy output argument to make Meson happy
                command : [etags_exe, '-o', '@0@/TAGS'.format(meson.project_source_root())] + all_files)
endif
