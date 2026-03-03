<a name="logo"/>
<div align="center">
<a href="https://www.dynare.org/" target="_blank">
<img src="https://www.dynare.org/assets/images/logo/dlogo.svg" alt="Dynare Logo"></img>
</a>
</div>

# Dynare

Described on the homepage: <https://www.dynare.org/>

Most users should use the precompiled package available for their OS, also
available via the Dynare homepage: <https://www.dynare.org/download/>.

# Contributions

To contribute to Dynare and participate in the Dynare community, please see: [CONTRIBUTING.md](https://git.dynare.org/Dynare/dynare/blob/master/CONTRIBUTING.md)

# License

Most of the source files are covered by the GNU General Public Licence version
3 or later (there are some exceptions to this, see [license.txt](license.txt) in
Dynare distribution for specifics).

# Building Dynare From Source

Here, we explain how to build from source:
- Dynare, including preprocessor and MEX files for MATLAB and Octave
- all the associated documentation (PDF and HTML)

This source can be retrieved in three forms:
- via git, at <https://git.dynare.org/Dynare/dynare.git>
- using the stable source archive of the latest Dynare version from <https://www.dynare.org/download/>
- using a source snapshot of the unstable version, also from <https://www.dynare.org/download/>

The first section of this page provides general instructions that apply to all platforms. Then, some specific platforms are discussed.

**Contents**

1. [**General Instructions**](#general-instructions)
1. [**Debian or Ubuntu**](#debian-or-ubuntu)
1. [**Fedora, RHEL, or derivatives**](#fedora-rhel-or-derivatives)
1. [**Arch Linux**](#arch-linux)
1. [**Windows**](#windows)
1. [**macOS**](#macos)
1. [**Docker**](#docker)

## General Instructions

### Prerequisites

A number of tools and libraries are needed in order to recompile everything. You don't necessarily need to install everything, depending on what you want to compile.

- The [GNU Compiler Collection](https://gcc.gnu.org/), version 13 or later, with
  gcc, g++ and gfortran
- [MATLAB](https://mathworks.com), version R2020a or later (if you want to compile the MEX for MATLAB)
- [GNU Octave](https://www.octave.org), version 8.4.0 or later, with
  - the development headers (if you want to compile the MEX for Octave)
  - the development libraries corresponding to the [UMFPACK](https://people.engr.tamu.edu/davis/suitesparse.html) packaged with Octave (if you want to compile the MEX for Octave)
  - the [statistics](https://gnu-octave.github.io/packages/statistics/) package and, optionally, the [control](https://gnu-octave.github.io/packages/control/), [datatypes](https://gnu-octave.github.io/packages/datatypes/), [io](https://gnu-octave.github.io/packages/io/) and [optim](https://gnu-octave.github.io/packages/optim/) packages
- [Meson](https://mesonbuild.com), version 1.3.0 or later
- [Pkgconf](http://pkgconf.org/), or another pkg-config implementation
- [Bash](https://www.gnu.org/software/bash/)
- [Boost libraries](https://www.boost.org), version 1.36 or later (only the
  headers are needed, not any of the compiled libraries)
- [Bison](https://www.gnu.org/software/bison/), version 3.2 or later
- [Flex](https://github.com/westes/flex), version 2.5.4 or later
- [SLICOT](https://www.slicot.org)
- A decent LaTeX distribution (if you want to compile PDF documentation),
  ideally with Beamer
- For building the reference manual:
  - [Sphinx](https://www.sphinx-doc.org/), with the [extension for BibTeX style
    citations](https://github.com/mcmtroffaes/sphinxcontrib-bibtex)
  - [MathJax](https://www.mathjax.org/)
- Optionally, the [X-13ARIMA-SEATS Seasonal Adjustment Program](https://www.census.gov/data/software/x13as.html)
- Optionally, the unpacked source tree of
  [SuiteSparse](https://www.suitesparse.com), for compiling the ParU MEX file
  (needed for `stack_solve_algo=8`); see the `-Dsuitesparse_src_path` option
  documented below
- Optionally, the [Panua PARDISO](https://www.panua.ch) library, for compiling
  the PARDISO MEX files (needed for `stack_solve_algo=9` and
  `stack_solve_algo=10`); see the `-Dpardiso` option documented below

### Renaming MATLAB's bundled GCC libraries on Linux

MATLAB ships its own copies of GCC runtime libraries (`libgcc_s`, `libstdc++`,
`libgfortran`, `libquadmath`) that conflict with the newer versions provided
by most Linux distributions. On Debian or Ubuntu this is handled automatically by
installing the `matlab-support` package, but on other distributions (Fedora,
RHEL and derivatives, Arch Linux, etc.) you need to rename the conflicting
libraries by appending `.bak`:
```sh
MATLAB_LIB=/usr/local/MATLAB/R2025b/sys/os/glnxa64
for f in $MATLAB_LIB/libgcc_s.so.1 \
         $MATLAB_LIB/libstdc++.so.6 \
         $MATLAB_LIB/libgfortran.so.5 \
         $MATLAB_LIB/libquadmath.so.0; do
    [ -e "$f" ] && sudo mv "$f" "$f.bak" && echo "Renamed $f to $f.bak"
done
```
The SONAMEs listed above are the same ones handled by [Debian's matlab-support package](https://packages.debian.org/bullseye/matlab-support).
Note that you need to adapt the MATLAB path to your installation.
Effectively, this forces MATLAB to use the system GCC libraries instead of its own
(potentially outdated) copies, which prevents crashes and link-time errors when
compiling MEX files. On Windows and macOS, you don't need to rename the libraries
as they are in different locations and dynamic linking works differently.

### Octave and Octave add-on packages
Octave and/or some Octave add-on packages (e.g. `datatypes`) may not be available
through the default package manager depending on the operating system and distribution.
The [Octave wiki](https://wiki.octave.org/Category:Installation) provides a list
of packaged versions and instructions to compile from source.
Alternatively, it is available as a [snap package](https://snapcraft.io/octave).
Add-on packages can be installed from *within Octave* using:
```octave
pkg install io datatypes statistics control struct optim
```

**Note:** On Octave 10 and older, the `-forge` flag is required:
```octave
pkg install -forge io datatypes statistics control struct optim
```

### Preparing the sources

If you have downloaded the sources from an official source archive or the source snapshot, just unpack them.

If you want to use Git, do the following from a terminal (note that you must
have the [Git LFS](https://git-lfs.github.com/) extension installed):
```sh
git clone --recurse-submodules https://git.dynare.org/Dynare/dynare.git
cd dynare
```
If you want a certain version (e.g., 7.0), then add `--single-branch --branch 7.0` to the git clone command.

### Configuring the build directory

If you want to compile for MATLAB, please run the following (after adapting the path to MATLAB):
```sh
meson setup -Dmatlab_path=/usr/local/MATLAB/R2025b --buildtype=debugoptimized build-matlab
```
The build directory will thus be `build-matlab`.

Or for Octave:
```sh
meson setup -Dbuild_for=octave --buildtype=debugoptimized build-octave
```
The build directory will thus be `build-octave`.

Note that if you do not choose `build-matlab` (under MATLAB) or `build-octave`
(under Octave) as the build directory, you will need to set the environment
variable `DYNARE_BUILD_DIR` to the full path of your build tree, before running
MATLAB or Octave, if you want Dynare to be able to find the preprocessor and
the MEX files.

Optionally, you can set the `DYNARE_PROFILE_SUITE` environment variable to `true` to
enable profiling during test runs. This will
generate a profiling summary of the top 20 functions by total execution time,
which can be useful for performance analysis and optimization.

The `-Dsuitesparse_src_path=…` option can be used to pass the path to the
unpacked SuiteSparse source tree, so that the ParU MEX file can be compiled.

The `-Dpardiso=enabled` option can be set to force the detection of Panua
PARDISO in the library path and compile the corresponding MEX files.
Conversely, `-Dpardiso=disabled` can be used to ignore it. The default,
`-Dpardiso=auto` will compile the PARDISO MEX files if the PARDISO library is
detected, otherwise it will ignore it.

It is possible to specify various Meson options, see the Meson documentation
for more details. Modifying the options of an existing build directory can be
done using the `meson configure` command.

### Building

For compiling the preprocessor and the MEX files:
```sh
meson compile -C <builddir>
```
where `<builddir>` is the build directory, typically either `build-matlab` or `build-octave`.

PDF and HTML documentation can be built with:
```sh
meson compile -C <builddir> doc
```

### Check

Dynare includes unit tests (in the MATLAB functions) and integration tests (in the `tests` subfolder). All the tests can be run with:
```sh
meson test -C <builddir>
```

Depending on your machine's performance, this can take several hours.

Note that running the testsuite with Octave requires the additional packages `pstoedit`, `epstool`, `xfig`, and `gnuplot`.

Often, it does not make sense to run the complete testsuite. For instance, if you modify codes only related to the perfect foresight model solver, you can decide to run only a subset of the integration tests, with:
```sh
meson test -C <builddir> --suite deterministic_simulations
```
This will run all the integration tests in `tests/deterministic_simulations`.
This syntax also works with a nested directory (e.g., `--suite deterministic_simulations/purely_forward`).

Finally, if you want to run a single integration test, e.g., `deterministic_simulations/lbj/rbc.mod`:
```sh
meson test -C <builddir> deterministic_simulations/lbj/rbc.mod
```
NB: Some individual tests cannot be run using that syntax, if they are a dependency in a chain of tests (see the `mod_and_m_tests` variable `meson.build`); in that case, you should use the name of the last `.mod` file in the chain as the test name to be passed to `meson test`.

## Debian or Ubuntu

### Prerequisites
All the prerequisites are packaged:

- `gcc`
- `g++`
- `gfortran`
- `octave-dev`
- `libboost-dev`
- `libslicot-dev` and `libslicot64-pic` (the latter is not available in Debian ⩽ 13 and Ubuntu ⩽ 25.10; `libslicot-pic` should be used instead)
- `libsuitesparse-dev`
- `flex` and `libfl-dev`
- `bison`
- `meson`
- `pkgconf`
- `texlive`
- `texlive-publishers` (for Econometrica bibliographic style)
- `texlive-latex-extra` (for fullpage.sty)
- `texlive-fonts-extra` (for ccicons)
- `texlive-science` (for amstex)
- `lmodern` (for macroprocessor PDF)
- `cm-super` (for BVAR à la Sims and GSA PDFs)
- `python3-sphinx`
- `python3-sphinxcontrib.bibtex`
- `tex-gyre`
- `latexmk`
- `libjs-mathjax`
- `x13as`

You can install them all at once with:
```sh
apt install gcc g++ gfortran octave-dev libboost-dev libslicot-dev libslicot-pic libslicot64-pic libsuitesparse-dev flex libfl-dev bison meson pkgconf texlive texlive-publishers texlive-latex-extra texlive-fonts-extra texlive-science lmodern cm-super python3-sphinx python3-sphinxcontrib.bibtex make tex-gyre latexmk libjs-mathjax x13as
```

>>> [!note] ⚠️ MATLAB GCC libraries
If you use MATLAB on Debian or Ubuntu, you need to run `apt install matlab-support` and confirm to rename the GCC libraries shipped with MATLAB to avoid conflicts with GCC libraries shipped by your distribution.
>>>

### Compile Dynare from source

- Prepare the Dynare sources:
```sh
export DYNAREDIR=$HOME/dynare
git clone --recurse-submodules https://git.dynare.org/Dynare/dynare.git $DYNAREDIR/unstable
cd $DYNAREDIR/unstable
```

- Configure Dynare from the source directory for MATLAB (adapt the path to your MATLAB installation):
```sh
meson setup -Dmatlab_path=/usr/local/MATLAB/R2025b \
      --buildtype=debugoptimized \
      build-matlab
```
To compile for Octave (additionally or instead):
```sh
meson setup -Dbuild_for=octave \
      --buildtype=debugoptimized \
      build-octave
```

- Compile:
```sh
meson compile -C build-matlab
meson compile -C build-octave
```

- Optionally, run the testsuite either for MATLAB or Octave:
```sh
meson test -C build-matlab --num-processes=$(nproc)
meson test -C build-octave --num-processes=$(nproc)
```

- The documentation can be compiled directly:
```sh
meson compile -C build-matlab doc
```

## Fedora, RHEL, or derivatives

These instructions apply to Fedora, RHEL, and RHEL-based distributions (AlmaLinux, CentOS Stream, Rocky Linux). 

### Prerequisites on Fedora

All prerequisites are available from the default Fedora repositories:

- `gcc`, `gcc-c++`, `gcc-gfortran`
- `boost-devel`
- `suitesparse-devel`
- `flex`
- `bison`
- `meson`
- `pkgconf`
- `redhat-rpm-config`
- `cmake` (for building SLICOT from source)
- `git-lfs`
- `octave`, `octave-devel`, `octave-statistics`, `octave-io`, `octave-optim`, `octave-control`, `octave-datatypes` (only needed for Octave MEX files; see also [Octave and Octave add-on packages](#octave-and-octave-add-on-packages))
- `texlive-scheme-minimal`, `texlive-collection-publishers`, `texlive-collection-latexextra`, `texlive-collection-fontsextra`, `texlive-collection-fontsrecommended`, `texlive-collection-latexrecommended`, `texlive-collection-science`, `texlive-collection-plaingeneric` (only needed for documentation)
- `python3-sphinx`, `python3-sphinxcontrib-bibtex`, `python3-sphinx-latex` (only needed for documentation)
- `latexmk` (only needed for documentation)
- `mathjax` (only needed for documentation)

You can install them all at once with:
```sh
sudo dnf install -y gcc gcc-c++ gcc-gfortran boost-devel suitesparse-devel flex bison meson pkgconf redhat-rpm-config cmake git-lfs octave octave-devel octave-statistics octave-io octave-optim octave-control octave-datatypes texlive-scheme-minimal texlive-collection-publishers texlive-collection-latexextra texlive-collection-fontsextra texlive-collection-fontsrecommended texlive-collection-latexrecommended texlive-collection-science texlive-collection-plaingeneric python3-sphinx python3-sphinxcontrib-bibtex python3-sphinx-latex latexmk mathjax
```
If you are installing `git-lfs` for the first time, you also need to run
`git lfs install` once after installing it.

>>> [!note] ⚠️ MATLAB GCC libraries
If you use MATLAB on Fedora, you need to [rename MATLAB's bundled GCC libraries](#renaming-matlabs-bundled-gcc-libraries-on-linux) manually.
>>>

### Prerequisites on RHEL, AlmaLinux, Rocky Linux, or CentOS Stream

On RHEL and derivatives, you first need to enable two additional repositories
that provide packages not included in the base distribution:

- [**EPEL**](https://docs.fedoraproject.org/en-US/epel/) (Extra Packages for
  Enterprise Linux) — a Fedora project that provides high-quality additional
  packages for RHEL-based distributions.
- [**CRB**](https://docs.fedoraproject.org/en-US/epel/getting-started/)
  (CodeReady Builder) — a repository of development packages (headers, static
  libraries) that are not enabled by default.

On **AlmaLinux, CentOS Stream, or Rocky Linux**, enable both with:
```sh
sudo dnf install -y epel-release
sudo dnf config-manager --set-enabled crb
```

On **RHEL** with a subscription, replace `RHEL_VERSION` with your major release:
```sh
RHEL_VERSION=10
sudo subscription-manager repos --enable codeready-builder-for-rhel-${RHEL_VERSION}-$(arch)-rpms
sudo dnf install -y https://dl.fedoraproject.org/pub/epel/epel-release-latest-${RHEL_VERSION}.noarch.rpm
```

Install the core build dependencies:

- `gcc`, `gcc-c++`, `gcc-gfortran`
- `boost-devel`
- `suitesparse-devel`
- `flex`
- `bison`
- `meson`
- `pkgconf`
- `redhat-rpm-config`
- `cmake` (for building SLICOT from source)
- `git-lfs`

```sh
sudo dnf install -y gcc gcc-c++ gcc-gfortran boost-devel suitesparse-devel flex bison meson pkgconf redhat-rpm-config cmake git-lfs
```
If you are installing `git-lfs` for the first time, you also need to run
`git lfs install` once after installing it.

>>> [!note] ⚠️ MATLAB GCC libraries
If you use MATLAB on RHEL-based distributions, you need to [rename MATLAB's bundled GCC libraries](#renaming-matlabs-bundled-gcc-libraries-on-linux) manually.
>>>

>>> [!note] ⚠️ Octave MEX files on RHEL-based distributions
Octave is currently not available via EPEL for RHEL-based distributions
and needs to be [compiled from source](https://wiki.octave.org/Octave_for_Red_Hat_Linux_systems).
Alternatively, it is available as a [snap package](https://snapcraft.io/octave).
Please install add-on packages using `pkg install`; see [Octave and Octave add-on packages](#octave-and-octave-add-on-packages).
>>>

>>> [!note] ⚠️ Building the documentation on RHEL-based distributions
Many of the TeX Live collection packages available on Fedora are not packaged
on RHEL and derivatives. We recommend installing TeX Live manually from
[TUG](https://www.tug.org/texlive/) using the upstream installer, which gives
you access to all packages and the `tlmgr` package manager.
For Sphinx and extensions, create a Python virtual environment:
```sh
cd $DYNAREDIR
python3 -m venv sphinx-env
sphinx-env/bin/pip install sphinx sphinxcontrib-bibtex
```
This environment only needs to be **created once** and can be **reused** for all
Dynare source trees within `$DYNAREDIR`.
>>>

### Compile dependencies

- Set a working directory for Dynare and its dependencies:
```sh
export DYNAREDIR=$HOME/dynare
```

- Compile and install SLICOT (not packaged on Fedora or RHEL):
```sh
mkdir -p $DYNAREDIR/slicot/lib
cd $DYNAREDIR/slicot
wget -O slicot-5.9.1.tar.gz https://github.com/SLICOT/SLICOT-Reference/archive/refs/tags/v5.9.1.tar.gz
tar xf slicot-5.9.1.tar.gz
cd SLICOT-Reference-5.9.1
sed -i 's/FIND_PACKAGE(\(.*\) REQUIRED)/FIND_PACKAGE(\1)/g' CMakeLists.txt
cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DSLICOT_TESTING=OFF -DSLICOT_INTEGER8=ON
cmake --build build
cp build/lib/libslicot64.a $DYNAREDIR/slicot/lib/libslicot64_pic.a
```
If you also want to compile for Octave, additionally build the non-integer8 variant:
```sh
rm -rf build
cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DSLICOT_TESTING=OFF
cmake --build build
cp build/lib/libslicot.a $DYNAREDIR/slicot/lib/libslicot_pic.a
```

- Optionally, compile and install the X-13ARIMA-SEATS Seasonal Adjustment
  Program (not packaged on Fedora or RHEL):
```sh
mkdir -p $DYNAREDIR/x13as
cd $DYNAREDIR/x13as
wget https://www2.census.gov/software/x-13arima-seats/x13as/unix-linux/program-archives/x13as_asciisrc-v1-1-b62.tar.gz
tar xf x13as_asciisrc-v1-1-b62.tar.gz
cd x13as_asciisrc-v1-1-b62
sed -i "s|-static| |" makefile.gf
make -f makefile.gf FFLAGS="-O2 -std=legacy" PROGRAM=x13as
sudo cp x13as /usr/local/bin/
```
Alternatively, if you don't have admin privileges, you can install it into
`$HOME/.local/bin` (make sure this directory is in your `PATH`):
```sh
mkdir -p $HOME/.local/bin
cp x13as $HOME/.local/bin/
```

### Compile Dynare from source

- Prepare the Dynare sources:
```sh
git clone --recurse-submodules https://git.dynare.org/Dynare/dynare.git $DYNAREDIR/unstable
cd $DYNAREDIR/unstable
```

- Configure Dynare from the source directory for MATLAB (adapt the path to your MATLAB installation):
```sh
meson setup -Dmatlab_path=/usr/local/MATLAB/R2025b \
      --buildtype=debugoptimized \
      -Dfortran_args="['-B','$DYNAREDIR/slicot/lib']" \
      build-matlab
```
To compile for Octave (additionally or instead):
```sh
meson setup -Dbuild_for=octave \
      --buildtype=debugoptimized \
      -Dfortran_args="['-B','$DYNAREDIR/slicot/lib']" \
      build-octave
```

- Compile:
```sh
meson compile -C build-matlab
meson compile -C build-octave
```

- Optionally, run the testsuite either for MATLAB or Octave:
```sh
meson test -C build-matlab --num-processes=$(nproc)
meson test -C build-octave --num-processes=$(nproc)
```

### Compile the documentation from source

The documentation packages listed above provide `sphinx-build`,
`sphinxcontrib-bibtex`, and the required LaTeX distribution. If they are
installed and `sphinx-build` was in PATH when you initially ran `meson setup`,
you can compile the documentation directly:
```sh
meson compile -C build-matlab doc
```

If `sphinx-build` was not in PATH during initial setup, reconfigure the build directory.
On RHEL-based distributions where Sphinx was installed in a virtual environment,
the reconfigure command would be:
```sh
env PATH="$DYNAREDIR/sphinx-env/bin:$PATH" meson setup --reconfigure build-matlab
```

### 🛠️ Troubleshooting

>>> [!warning] MATLAB fails to start on RHEL and RHEL-derivatives (Alma Linux, CentOS Stream, Rocky Linux)
MATLAB may fail to start due to missing libraries, e.g.:
```sh
Unexpected exception: 'N7mwboost10wrapexceptINS_16exception_detail39current_exception_std_exception_wrapperISt13runtime_errorEEEE: Error loading /usr/local/MATLAB/R2025b/bin/glnxa64/matlab_startup_plugins/matlab_graphics_ui/mwuixloader.so. libXt.so.6: cannot open shared object file: No such file or directory: Success: Success' in createMVMAndWaitForExit phase 'Creating local MVM'
```
Install them with:
```sh
sudo dnf install -y libXt libXext libXrender libXrandr libXcursor libXinerama libXi libXtst
```
>>>

>>> [!warning] Static TLS block error with MATLAB MEX files
You may encounter the following error when MATLAB tries to load Dynare's MEX files:

> Invalid MEX-file '…/sparse_hessian_times_B_kronecker_C.mexa64':
> 'libgomp.so.1: cannot allocate memory in static TLS block'

This happens because MATLAB loads MEX files at runtime via `dlopen`, and
libraries like `libgomp.so.1` (OpenMP) require static Thread-Local Storage
(TLS) that must be reserved at program startup. When the static TLS block is
exhausted, loading fails.
The fix is to preload `libgomp` so it gets its TLS slot early, and optionally
increase the available static TLS space. Start MATLAB with:

```sh
LD_PRELOAD=/usr/lib64/libgomp.so.1 matlab
```
>>>

## Arch Linux

These instructions show how to compile Dynare from source on Arch Linux.

### Prerequisites

All core prerequisites are available from the official Arch repositories:

- `gcc`, `gcc-fortran`
- `boost`
- `suitesparse`
- `flex`
- `bison`
- `meson`
- `pkgconf`
- `cmake` (for building SLICOT from source)
- `git-lfs`

You can install them all at once with:
```sh
sudo pacman -S gcc gcc-fortran boost suitesparse flex bison meson pkgconf cmake git-lfs
```
If you are installing `git-lfs` for the first time, you also need to run
`git lfs install` once after installing it.

If you want to compile the MEX files for Octave (instead of or in addition to
MATLAB), also install:
```sh
sudo pacman -S octave
```
The required Octave add-on packages can be installed from *within Octave*, see
[Octave and Octave add-on packages](#octave-and-octave-add-on-packages).

If you want to compile the documentation, install:
```sh
sudo pacman -S texlive-basic texlive-binextra texlive-latexextra texlive-publishers texlive-fontsextra texlive-fontsrecommended texlive-latexrecommended texlive-plaingeneric texlive-bibtexextra texlive-mathscience python-sphinx python-sphinxcontrib-bibtex mathjax
```

>>> [!note] ⚠️ MATLAB GCC libraries
If you use MATLAB on Arch Linux, you need to [rename MATLAB's bundled GCC libraries](#renaming-matlabs-bundled-gcc-libraries-on-linux) manually.
>>>

### Compile dependencies

- Set a working directory for Dynare and its dependencies:
```sh
export DYNAREDIR=$HOME/dynare
```

- Compile and install SLICOT (not packaged on Arch Linux):
```sh
mkdir -p $DYNAREDIR/slicot/lib
cd $DYNAREDIR/slicot
wget -O slicot-5.9.1.tar.gz https://github.com/SLICOT/SLICOT-Reference/archive/refs/tags/v5.9.1.tar.gz
tar xf slicot-5.9.1.tar.gz
cd SLICOT-Reference-5.9.1
sed -i 's/FIND_PACKAGE(\(.*\) REQUIRED)/FIND_PACKAGE(\1)/g' CMakeLists.txt
cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DSLICOT_TESTING=OFF -DSLICOT_INTEGER8=ON
cmake --build build
cp build/lib/libslicot64.a $DYNAREDIR/slicot/lib/libslicot64_pic.a
```
If you also want to compile for Octave, additionally build the non-integer8 variant:
```sh
rm -rf build
cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DSLICOT_TESTING=OFF
cmake --build build
cp build/lib/libslicot.a $DYNAREDIR/slicot/lib/libslicot_pic.a
```

- Optionally, compile and install the X-13ARIMA-SEATS Seasonal Adjustment
  Program. You can compile it from source:
```sh
mkdir -p $DYNAREDIR/x13as
cd $DYNAREDIR/x13as
wget https://www2.census.gov/software/x-13arima-seats/x13as/unix-linux/program-archives/x13as_asciisrc-v1-1-b62.tar.gz
tar xf x13as_asciisrc-v1-1-b62.tar.gz
cd x13as_asciisrc-v1-1-b62
sed -i "s|-static| |" makefile.gf
make -f makefile.gf FFLAGS="-O2 -std=legacy" PROGRAM=x13as
sudo cp x13as /usr/local/bin/
```
Alternatively, if you don't have admin privileges, you can install it into
`$HOME/.local/bin` (make sure this directory is in your `PATH`):
```sh
mkdir -p $HOME/.local/bin
cp x13as $HOME/.local/bin/
```

### Compile Dynare from source

- Prepare the Dynare sources:
```sh
git clone --recurse-submodules https://git.dynare.org/Dynare/dynare.git $DYNAREDIR/unstable
cd $DYNAREDIR/unstable
```

- Configure Dynare from the source directory for MATLAB (adapt the path to your MATLAB installation):
```sh
meson setup -Dmatlab_path=/usr/local/MATLAB/R2025b \
      --buildtype=debugoptimized \
      -Dfortran_args="['-B','$DYNAREDIR/slicot/lib']" \
      build-matlab
```
To compile for Octave (additionally or instead):
```sh
meson setup -Dbuild_for=octave \
      --buildtype=debugoptimized \
      -Dfortran_args="['-B','$DYNAREDIR/slicot/lib']" \
      build-octave
```

- Compile:
```sh
meson compile -C build-matlab
meson compile -C build-octave
```

- Optionally, run the testsuite either for MATLAB or Octave:
```sh
meson test -C build-matlab --num-processes=$(nproc)
meson test -C build-octave --num-processes=$(nproc)
```

- The documentation can be compiled directly:
```sh
meson compile -C build-matlab doc
```

## Windows

These instructions target Windows 11 `x86_64` with a [MSYS2](https://www.msys2.org) toolchain.
Compilation for `arm64` is currently not supported.

>>> [!note] ℹ️ Note
Building Octave MEX files and the documentation under MSYS2 is not supported.
The official Windows packages (which include Octave MEX files) are cross-compiled from Linux.
>>>

- Install [MSYS2](https://www.msys2.org)
- Run a MSYS2 UCRT64 shell (make sure to run the one with "UCRT" in
  the name)
- Update the system:
```sh
pacman -Syu
```
  You may be asked to close the window at the end of the
  first upgrade batch, in which case you should rerun the upgrade in a new
  window to complete the upgrade.
- Install all needed dependencies:
```sh
pacman -S git bison flex make cmake tar mingw-w64-ucrt-x86_64-meson mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-gcc-fortran mingw-w64-ucrt-x86_64-boost mingw-w64-ucrt-x86_64-pkgconf
```
- Compile and install SLICOT:
```sh
wget -O slicot-5.9.1.tar.gz https://github.com/SLICOT/SLICOT-Reference/archive/refs/tags/v5.9.1.tar.gz
tar xf slicot-5.9.1.tar.gz
cd SLICOT-Reference-5.9.1
sed -i 's/FIND_PACKAGE(\(.*\) REQUIRED)/FIND_PACKAGE(\1)/g' CMakeLists.txt
FFLAGS="-fno-underscoring" cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DSLICOT_TESTING=OFF -DSLICOT_INTEGER8=ON
cmake --build build
mkdir -p /usr/local/lib
cp build/lib/libslicot64.a /usr/local/lib/libslicot64_pic.a
cd ..
```
- Prepare the Dynare sources, either by unpacking the source tarball or with:
```sh
git clone --recurse-submodules https://git.dynare.org/Dynare/dynare.git
cd dynare
```
- Configure Dynare from the source directory for MATLAB:
```sh
meson setup -Dmatlab_path=/c/Progra~1/MATLAB/R2025b --buildtype=debugoptimized --prefer-static -Dfortran_args="['-B','C:/msys64/usr/local/lib']" build-matlab
```
where the path of MATLAB is specified. Note that you should use
the MSYS2 notation and not put spaces in the MATLAB path.
Alternatively, if your filesystem does not have short filenames (8dot3),
then you can run `mkdir -p /usr/local/MATLAB && mount c:/Program\ Files/MATLAB /usr/local/MATLAB`,
and then pass `/usr/local/MATLAB/R2025b` as MATLAB path to the `meson setup` command.
- Compile:
```sh
meson compile -C build-matlab
```
- Optionally, run the testsuite:
```sh
meson test -C build-matlab
```

## macOS

These instructions are for Apple Silicon Macs (`arm64`) using a
[Homebrew](https://brew.sh/) toolchain.
If you need to compile for Intel (`x86_64`) using Rosetta 2, see the
[note at the end of this section](#note-compiling-for-x86_64-using-rosetta-2).

>>> [!note] ℹ️ Note
The official macOS packages are compiled for both `arm64` and `x86_64` versions.
>>>

### Prerequisites

- Install the Xcode Command Line Tools:
```sh
xcode-select --install
```

- Install [Homebrew](https://brew.sh/) (the default installation on Apple Silicon
  is for `arm64` under `/opt/homebrew`):
```sh
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
Don't forget to run the displayed commands (**Next steps**) in the terminal to
add Homebrew to your PATH.

- Install the required Homebrew packages:
```sh
brew install gcc meson bison flex boost pkg-config cmake git-lfs
```
If you are installing `git-lfs` for the first time, you also need to run
`git lfs install` once after installing it.

  If you want to compile the MEX files for Octave (instead of or in addition to
  MATLAB), also install:
  ```sh
  brew install octave
  ```
  Then install the required Octave add-on packages from *within Octave*, see
  [Octave and Octave add-on packages](#octave-and-octave-add-on-packages).

  If you want to compile the documentation, install
  [MacTeX](https://www.tug.org/mactex/index.html). It is advised to symlink
  `pdflatex`, `bibtex`, and `latexmk` into `/usr/local/bin`:
  ```sh
  sudo ln -s /Library/TeX/texbin/pdflatex /usr/local/bin/pdflatex
  sudo ln -s /Library/TeX/texbin/bibtex /usr/local/bin/bibtex
  sudo ln -s /Library/TeX/texbin/latexmk /usr/local/bin/latexmk
  ```
  If you don't have admin privileges, you can symlink them into
  `$HOME/.local/bin` instead (make sure this directory is in your `PATH`).

- Install MATLAB (R2023b or later for `arm64`) and any additional toolboxes.
See [the official instructions](https://www.mathworks.com/support/requirements/apple-silicon.html)
for installing MATLAB on Apple Silicon.
Run MATLAB at least once to make sure you have a valid license.

### Compile dependencies

- Set a working directory for Dynare and its dependencies:
```sh
export DYNAREDIR=$HOME/dynare
```

- Compile and install SLICOT:
```sh
mkdir -p $DYNAREDIR/slicot/lib
cd $DYNAREDIR/slicot
curl -L -o slicot-5.9.1.tar.gz https://github.com/SLICOT/SLICOT-Reference/archive/refs/tags/v5.9.1.tar.gz
tar xf slicot-5.9.1.tar.gz
cd SLICOT-Reference-5.9.1
sed -i .bak 's/FIND_PACKAGE(\(.*\) REQUIRED)/FIND_PACKAGE(\1)/g' CMakeLists.txt
cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DSLICOT_TESTING=OFF -DSLICOT_INTEGER8=ON
cmake --build build
cp build/lib/libslicot64.a $DYNAREDIR/slicot/lib/libslicot64_pic.a
```
If you also want to compile for Octave, additionally build the non-integer8 variant:
```sh
rm -rf build
cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DSLICOT_TESTING=OFF
cmake --build build
cp build/lib/libslicot.a $DYNAREDIR/slicot/lib/libslicot_pic.a
```

- Optionally, compile and install the X-13ARIMA-SEATS Seasonal Adjustment Program:
```sh
mkdir -p $DYNAREDIR/x13as
cd $DYNAREDIR/x13as
curl -O https://www2.census.gov/software/x-13arima-seats/x13as/unix-linux/program-archives/x13as_asciisrc-v1-1-b62.tar.gz
tar xf x13as_asciisrc-v1-1-b62.tar.gz
cd x13as_asciisrc-v1-1-b62
sed -i '' 's/-static//g' makefile.gf
make -j$(sysctl -n hw.perflevel0.physicalcpu) -f makefile.gf FC=/opt/homebrew/bin/gfortran LINKER=/opt/homebrew/bin/gfortran FFLAGS="-O2 -std=legacy" LDFLAGS="-static-libgcc -static-libgfortran -static-libquadmath" PROGRAM=x13as
mkdir -p $HOME/.local/bin
cp x13as $HOME/.local/bin/x13as
```
Alternatively, if you have admin privileges: `sudo cp x13as /usr/local/bin/x13as`.
Make sure the install directory is in your `PATH`.

- Optionally, if you want to build the documentation, create a Python virtual
  environment with Sphinx and the BibTeX extension (the `sphinxcontrib-bibtex`
  package is not available via Homebrew):
```sh
cd $DYNAREDIR
python3 -m venv sphinx-env
sphinx-env/bin/pip install sphinx sphinxcontrib-bibtex
```
This environment only needs to be **created once** and can be **reused** for all
Dynare source trees within `$DYNAREDIR` (e.g., `$DYNAREDIR/unstable`,
`$DYNAREDIR/7.0`, etc.).

### Compile Dynare from source

- Prepare the Dynare sources:
```sh
git clone --recurse-submodules https://git.dynare.org/Dynare/dynare.git $DYNAREDIR/unstable
cd $DYNAREDIR/unstable
```

- Configure Dynare from the source directory for MATLAB (adapt the path to your MATLAB installation):
```sh
meson setup --native-file macOS/homebrew-native-arm64.ini \
      -Dmatlab_path=/Applications/MATLAB_R2025b.app \
      --buildtype=debugoptimized \
      -Dfortran_args="['-B','$DYNAREDIR/slicot/lib']" \
      build-matlab
```
To compile for Octave (additionally or instead):
```
meson setup --native-file macOS/homebrew-native-arm64.ini \
      -Dbuild_for=octave \
      --buildtype=debugoptimized \
      -Dfortran_args="['-B','$DYNAREDIR/slicot/lib']" \
      build-octave
```

- Compile:
```sh
meson compile -C build-matlab
meson compile -C build-octave
```

- Optionally, run the testsuite either for MATLAB or Octave:
```sh
meson test -C build-matlab --num-processes=$(sysctl -n hw.perflevel0.physicalcpu)
meson test -C build-octave --num-processes=$(sysctl -n hw.perflevel0.physicalcpu)
```
where `--num-processes` is set to the number of performance cores on your Mac.

### Compile the documentation from source

If you set up the Sphinx virtual environment as described above, make sure
`sphinx-build` is in PATH. If it was not in PATH when you initially ran
`meson setup`, reconfigure the build directory:
```sh
cd $DYNAREDIR/unstable
env PATH="$DYNAREDIR/sphinx-env/bin:$PATH" meson setup --reconfigure build-matlab
```

Then compile the documentation:
```sh
meson compile -C build-matlab doc
```

### 🛠️ Troubleshooting

>>> [!warning] Font Awesome font not found during PDF build
If the PDF build fails with an error like
`font fa7free1solid at 600 not found`, the Font Awesome 7 font map shipped with
MacTeX is not registered. Fix it by running:
```sh
sudo updmap-sys --enable Map=fontawesome7.map
```
>>>

### Optional: pass the full PATH to MATLAB for system commands

When starting MATLAB by clicking its application icon (rather than from a
terminal), it does not inherit your shell's PATH, 
which does not include `/usr/local/bin`, `$HOME/.local/bin`, or `/Library/TeX/texbin`.
You may check this by running in MATLAB:
```matlab
system('echo $PATH')
% /usr/bin:/bin:/usr/sbin:/sbin`)
```
To use system commands like `pdflatex`, `latexmk`, or `x13as` from within
MATLAB, you can either call them by their full path
(e.g., `/Library/TeX/texbin/pdflatex`) or extend the PATH by running in MATLAB:
```matlab
setenv('PATH', [getenv('PATH') ':/usr/local/bin:' getenv('HOME') '/.local/bin:/Library/TeX/texbin']);
```
You can place this in a `startup.m` file, at the top of your `.mod` file, or modify the
system default PATH via `/etc/paths`.

### Note: compiling for x86_64 using Rosetta 2

Apple has announced that Rosetta 2 will be deprecated in a future version of
macOS. The following is only needed if you must compile Dynare for Intel
(`x86_64`), e.g., for an older MATLAB version that does not support `arm64`.

Install Rosetta 2:
```sh
softwareupdate --install-rosetta --agree-to-license
```

Install the `x86_64` Homebrew (installs under `/usr/local`) alongside your
`arm64` installation:
```sh
arch -x86_64 /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Before running any build command for `x86_64`, prepend the `x86_64` Homebrew to
your PATH so that its packages take precedence:
```sh
export PATH="/usr/local/bin:$PATH"
```

Install the required packages under the `x86_64` Homebrew:
```sh
arch -x86_64 /usr/local/bin/brew install gcc meson bison flex boost pkg-config cmake git-lfs
```

All dependency compilation (SLICOT, x13as) must also be run under Rosetta 2,
replacing `/opt/homebrew` paths with `/usr/local` where applicable.

When configuring Dynare, use the `x86_64` native file and prefix all `meson`
commands with `arch -x86_64`:
```sh
arch -x86_64 meson setup --native-file macOS/homebrew-native-x86_64.ini \
      -Dmatlab_path=/Applications/MATLAB_R2025b.app \
      --buildtype=debugoptimized \
      -Dfortran_args="['-B','$DYNAREDIR/slicot/lib']" \
      build-matlab
arch -x86_64 meson compile -C build-matlab
```

Note that MATLAB for `arm64` (Apple Silicon) is available from R2023b onwards,
while `x86_64` supports R2020a and later.


## Docker
We offer a variety of pre-configured Docker containers for Dynare, pre-configured with Octave and MATLAB, including all recommended toolboxes.
These are readily available for your convenience on [Docker Hub](https://hub.docker.com/r/dynare/dynare).
The `scripts/docker` folder contains [information and instructions](scripts/docker/README.md) for interacting with, building, and customizing the containers.
