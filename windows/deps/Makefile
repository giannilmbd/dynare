# Copyright © 2017-2026 Dynare Team
#
# This file is part of Dynare.
#
# Dynare is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Dynare is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

include versions.mk

ROOT_PATH = $(realpath .)

WGET_OPTIONS := --no-verbose --no-use-server-timestamps --retry-connrefused --retry-on-host-error

.PHONY: all clean clean-all

all: sysroot-matlab sysroot-octave gcc x13as src/suitesparse

# Clean everything except downloaded tarballs
clean:
	rm -rf src octave sysroot-matlab sysroot-octave gcc x13as panua-pardiso

# Clean everything (including downloaded tarballs)
clean-all: clean
	rm -rf tarballs

#
# Slicot
#

tarballs/slicot-$(SLICOT_VERSION).tar.gz:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@ https://github.com/SLICOT/SLICOT-Reference/archive/refs/tags/v$(SLICOT_VERSION).tar.gz

src/slicot-matlab/build/lib/libslicot64.a: tarballs/slicot-$(SLICOT_VERSION).tar.gz
	rm -rf $(realpath $(dir $@)/../..)
	mkdir -p $(dir $@)/../..
	tar xf $< --directory $(dir $@)/../.. --strip-components=1
	patch -d $(dir $@)/../.. -p1 < slicot-no-blas-lapack-lookup.patch
	# NB: -fno-underscoring is passed in FFLAGS, and is prepended to CMAKE_Fortran_FLAGS_INIT from the toolchain file.
	# We don’t put it in that latter variable, so that the toolchain file can be shared with Octave.
	FFLAGS="-fno-underscoring" cmake -S $(dir $@)/../.. -B $(dir $@)/.. -DCMAKE_TOOLCHAIN_FILE=$(CURDIR)/toolchain-mingw64.cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DSLICOT_TESTING=OFF -DSLICOT_INTEGER8=ON
	cmake --build $(dir $@)/.. --verbose

src/slicot-octave/build/lib/libslicot.a: tarballs/slicot-$(SLICOT_VERSION).tar.gz
	rm -rf $(realpath $(dir $@)/../..)
	mkdir -p $(dir $@)/../..
	tar xf $< --directory $(dir $@)/../.. --strip-components=1
	patch -d $(dir $@)/../.. -p1 < slicot-no-blas-lapack-lookup.patch
	cmake -S $(dir $@)/../.. -B $(dir $@)/.. -DCMAKE_TOOLCHAIN_FILE=$(CURDIR)/toolchain-mingw64.cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DSLICOT_TESTING=OFF -DSLICOT_INTEGER8=OFF
	cmake --build $(dir $@)/.. --verbose


#
# Octave
#

tarballs/octave-${OCTAVE_VERSION}${OCTAVE_W64_BUILD}-w64.7z:
	mkdir -p $(dir $@)
	wget $(WGET_OPTIONS) -O $@ https://ftpmirror.gnu.org/gnu/octave/windows/$(notdir $@)

octave: tarballs/octave-${OCTAVE_VERSION}${OCTAVE_W64_BUILD}-w64.7z
	rm -rf $@
	7zr x $< > /dev/null
	mv octave-${OCTAVE_VERSION}-w64 $@
	touch $@

#
# MATLAB
#

tarballs/matlab%.tar.xz:
	mkdir -p $(dir $@)
	wget $(WGET_OPTIONS) -O $@.gpg https://www.dynare.org/pkg-build/$(notdir $@).gpg
	gpg --output $@ --decrypt $@.gpg
	rm $@.gpg

#
# The three root directories:
# - for building MATLAB MEX
# - for builing Octave MEX
# - the embedded GCC compiler (for use_dll)
#

MSYS2_EXCLUDES := --exclude .MTREE --exclude .BUILDINFO --exclude .PKGINFO

tarballs/mingw-w64-ucrt-x86_64-%-any.pkg.tar.zst:
	mkdir -p $(dir $@)
	wget $(WGET_OPTIONS) -O $@ https://mirror.msys2.org/mingw/ucrt64/$(notdir $@)

MATLAB_MSYS2_DEPS := \
	tarballs/mingw-w64-ucrt-x86_64-boost-$(MSYS2_BOOST_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-gsl-$(MSYS2_GSL_VERSION)-any.pkg.tar.zst

sysroot-matlab: tarballs/matlab64-$(MATLAB64_VERSION).tar.xz $(MATLAB_MSYS2_DEPS) src/slicot-matlab/build/lib/libslicot64.a panua-pardiso
	rm -rf $@
	mkdir $@
	for f in $(MATLAB_MSYS2_DEPS); do tar xf $$f --directory $@ $(MSYS2_EXCLUDES); done
	mkdir $@/MATLAB
	tar xf $< --directory $@/MATLAB/
	cp src/slicot-matlab/build/lib/libslicot64.a $@/ucrt64/lib/libslicot64_pic.a
	cp panua-pardiso/lib/libpardiso.dll $@/ucrt64/lib
	touch $@

OCTAVE_MSYS2_DEPS := \
	tarballs/mingw-w64-ucrt-x86_64-boost-$(MSYS2_BOOST_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-matio-$(MSYS2_MATIO_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-zlib-$(MSYS2_ZLIB_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-hdf5-$(MSYS2_HDF5_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-libaec-$(MSYS2_LIBAEC_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-openssl-$(MSYS2_OPENSSL_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-curl-$(MSYS2_CURL_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-zstd-$(MSYS2_ZSTD_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-brotli-$(MSYS2_BROTLI_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-libidn2-$(MSYS2_LIBIDN2_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-libpsl-$(MSYS2_LIBPSL_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-libssh2-$(MSYS2_LIBSSH2_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-nghttp2-$(MSYS2_NGHTTP2_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-nghttp3-$(MSYS2_NGHTTP3_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-ngtcp2-$(MSYS2_NGTCP2_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-libiconv-$(MSYS2_LIBICONV_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-libunistring-$(MSYS2_LIBUNISTRING_VERSION)-any.pkg.tar.zst

sysroot-octave: octave $(OCTAVE_MSYS2_DEPS) src/slicot-octave/build/lib/libslicot.a panua-pardiso
	rm -rf $@
	mkdir $@
	for f in $(OCTAVE_MSYS2_DEPS); do tar xf $$f --directory $@ $(MSYS2_EXCLUDES); done
	cp src/slicot-octave/build/lib/libslicot.a $@/ucrt64/lib/libslicot_pic.a
	cp -R octave/mingw64/include/octave-* octave/mingw64/include/suitesparse octave/mingw64/include/gsl $@/ucrt64/include
	cp -R octave/mingw64/lib/octave $@/ucrt64/lib
	cp octave/mingw64/lib/libblas.dll.a octave/mingw64/lib/liblapack.dll.a octave/mingw64/lib/libumfpack.dll.a octave/mingw64/lib/libgsl.dll.a octave/mingw64/lib/libgslcblas.dll.a octave/mingw64/lib/libopenblas.dll.a $@/ucrt64/lib
	cp octave/mingw64/lib/pkgconfig/gsl.pc octave/mingw64/lib/pkgconfig/UMFPACK.pc octave/mingw64/lib/pkgconfig/AMD.pc octave/mingw64/lib/pkgconfig/SuiteSparse_config.pc octave/mingw64/lib/pkgconfig/openblas.pc $@/ucrt64/lib/pkgconfig
	cp panua-pardiso/lib/libpardiso.dll $@/ucrt64/lib
	patch -p0 < octave-pkgconfig.patch
	touch $@

GCC_MSYS2_DEPS := \
	tarballs/mingw-w64-ucrt-x86_64-gcc-$(MSYS2_GCC_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-gcc-libs-$(MSYS2_GCC_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-gmp-$(MSYS2_GMP_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-binutils-$(MSYS2_BINUTILS_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-headers-git-$(MSYS2_HEADERS_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-crt-git-$(MSYS2_CRT_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-winpthreads-$(MSYS2_WINPTHREADS_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-libwinpthread-$(MSYS2_WINPTHREADS_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-zlib-$(MSYS2_ZLIB_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-zstd-$(MSYS2_ZSTD_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-isl-$(MSYS2_ISL_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-mpc-$(MSYS2_MPC_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-mpfr-$(MSYS2_MPFR_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-libiconv-$(MSYS2_LIBICONV_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-windows-default-manifest-$(MSYS2_WINDOWS_DEFAULT_MANIFEST_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-ucrt-x86_64-gettext-runtime-$(MSYS2_GETTEXT_RUNTIME_VERSION)-any.pkg.tar.zst

gcc: $(GCC_MSYS2_DEPS)
	rm -rf $@
	mkdir $@
	for f in $(GCC_MSYS2_DEPS); do tar xf $$f --directory $@ $(MSYS2_EXCLUDES); done
	touch $@

#
# X13AS
#

tarballs/x13as_ascii-v$(X13AS_VERSION).zip:
	mkdir -p $(dir $@)
	wget $(WGET_OPTIONS) -O $@ https://www2.census.gov/software/x-13arima-seats/x13as/windows/program-archives/$(notdir $@)

x13as: tarballs/x13as_ascii-v$(X13AS_VERSION).zip
	rm -rf $@
	unzip -o $<
	mv $@/x13as_ascii.exe $@/x13as.exe
	touch $@

#
# Panua PARDISO
#

tarballs/panua-pardiso-%-win.zip:
	mkdir -p $(dir $@)
	wget $(WGET_OPTIONS) -O $@.gpg https://www.dynare.org/pkg-build/$(notdir $@).gpg
	gpg --output $@ --decrypt $@.gpg
	rm $@.gpg

panua-pardiso: tarballs/panua-pardiso-$(PANUA_PARDISO_VERSION)-win.zip
	rm -rf $@
	unzip -o $<
	mv $(basename $(notdir $<)) $@
	touch $@

#
# SuiteSparse
#

tarballs/SuiteSparse-$(SUITESPARSE_VERSION).tar.gz:
	wget $(WGET_OPTIONS) -O $@ https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v$(SUITESPARSE_VERSION).tar.gz

src/suitesparse: tarballs/SuiteSparse-$(SUITESPARSE_VERSION).tar.gz
	rm -rf $@
	mkdir -p $@
	tar xf $< --directory $@ --strip-components=1
	touch $@
