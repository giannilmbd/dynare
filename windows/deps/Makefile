# Copyright © 2017-2025 Dynare Team
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

# See comments in windows/mingw-cross.ini
FFLAGS := -O3 -flto=auto -march=nocona -msahf -mtune=generic

.PHONY: all clean clean-all

all: sysroot-matlab sysroot-octave gcc x13as

# Clean everything except downloaded tarballs
clean:
	rm -rf src octave sysroot-matlab sysroot-octave gcc x13as

# Clean everything (including downloaded tarballs)
clean-all: clean
	rm -rf tarballs

#
# Slicot
#

tarballs/slicot-$(SLICOT_VERSION).tar.gz:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@ https://github.com/SLICOT/SLICOT-Reference/archive/refs/tags/v$(SLICOT_VERSION).tar.gz

src/slicot-matlab/libslicot64_pic.a: tarballs/slicot-$(SLICOT_VERSION).tar.gz
	rm -rf src/slicot-matlab
	mkdir -p src/slicot-matlab
	tar xf $< --directory src/slicot-matlab --strip-components=1
	make -C src/slicot-matlab -f makefile_Unix lib SLICOTLIB=../libslicot64_pic.a OPTS="$(FFLAGS) -fno-underscoring -fdefault-integer-8" FORTRAN=x86_64-w64-mingw32-gfortran LOADER=x86_64-w64-mingw32-gfortran ARCH=x86_64-w64-mingw32-gcc-ar

src/slicot-octave/libslicot_pic.a: tarballs/slicot-$(SLICOT_VERSION).tar.gz
	rm -rf src/slicot-octave
	mkdir -p src/slicot-octave
	tar xf $< --directory src/slicot-octave --strip-components=1
	make -C src/slicot-octave -f makefile_Unix lib SLICOTLIB=../libslicot_pic.a OPTS="$(FFLAGS)" FORTRAN=x86_64-w64-mingw32-gfortran LOADER=x86_64-w64-mingw32-gfortran ARCH=x86_64-w64-mingw32-gcc-ar


#
# Octave
#

tarballs/octave-${OCTAVE_VERSION}${OCTAVE_W64_BUILD}-w64.7z:
	mkdir -p tarballs
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
	mkdir -p tarballs
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

tarballs/mingw-w64-x86_64-%-any.pkg.tar.zst:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@ http://repo.msys2.org/mingw/x86_64/$(notdir $@)

MATLAB_MSYS2_DEPS := \
	tarballs/mingw-w64-x86_64-boost-$(MINGW64_BOOST_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-gsl-$(MINGW64_GSL_VERSION)-any.pkg.tar.zst

sysroot-matlab: tarballs/matlab64-$(MATLAB64_VERSION).tar.xz $(MATLAB_MSYS2_DEPS) src/slicot-matlab/libslicot64_pic.a
	rm -rf $@
	mkdir $@
	for f in $(MATLAB_MSYS2_DEPS); do tar xf $$f --directory $@ $(MSYS2_EXCLUDES); done
	mkdir $@/MATLAB
	tar xf $< --directory $@/MATLAB/
	cp src/slicot-matlab/libslicot64_pic.a $@/mingw64/lib/
	touch $@

OCTAVE_MSYS2_DEPS := \
	tarballs/mingw-w64-x86_64-boost-$(MINGW64_BOOST_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-matio-$(MINGW64_MATIO_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-zlib-$(MINGW64_ZLIB_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-hdf5-$(MINGW64_HDF5_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-libaec-$(MINGW64_LIBAEC_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-openssl-$(MINGW64_OPENSSL_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-curl-$(MINGW64_CURL_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-zstd-$(MINGW64_ZSTD_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-brotli-$(MINGW64_BROTLI_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-libidn2-$(MINGW64_LIBIDN2_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-libpsl-$(MINGW64_LIBPSL_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-libssh2-$(MINGW64_LIBSSH2_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-nghttp2-$(MINGW64_NGHTTP2_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-nghttp3-$(MINGW64_NGHTTP3_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-libiconv-$(MINGW64_LIBICONV_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-libunistring-$(MINGW64_LIBUNISTRING_VERSION)-any.pkg.tar.zst

sysroot-octave: octave $(OCTAVE_MSYS2_DEPS) src/slicot-octave/libslicot_pic.a
	rm -rf $@
	mkdir $@
	for f in $(OCTAVE_MSYS2_DEPS); do tar xf $$f --directory $@ $(MSYS2_EXCLUDES); done
	cp src/slicot-octave/libslicot_pic.a $@/mingw64/lib/
	cp -R octave/mingw64/include/octave-* octave/mingw64/include/suitesparse octave/mingw64/include/gsl $@/mingw64/include
	cp -R octave/mingw64/lib/octave $@/mingw64/lib
	cp octave/mingw64/lib/libblas.dll.a octave/mingw64/lib/liblapack.dll.a octave/mingw64/lib/libumfpack.dll.a octave/mingw64/lib/libgsl.dll.a octave/mingw64/lib/libgslcblas.dll.a $@/mingw64/lib
	cp octave/mingw64/lib/pkgconfig/gsl.pc octave/mingw64/lib/pkgconfig/UMFPACK.pc octave/mingw64/lib/pkgconfig/AMD.pc octave/mingw64/lib/pkgconfig/SuiteSparse_config.pc $@/mingw64/lib/pkgconfig
	patch -p0 < octave-pkgconfig.patch
	touch $@

GCC_MSYS2_DEPS := \
	tarballs/mingw-w64-x86_64-gcc-$(MINGW64_GCC_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-gcc-libs-$(MINGW64_GCC_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-gmp-$(MINGW64_GMP_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-binutils-$(MINGW64_BINUTILS_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-headers-git-$(MINGW64_HEADERS_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-crt-git-$(MINGW64_CRT_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-winpthreads-git-$(MINGW64_WINPTHREADS_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-libwinpthread-git-$(MINGW64_WINPTHREADS_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-zlib-$(MINGW64_ZLIB_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-zstd-$(MINGW64_ZSTD_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-isl-$(MINGW64_ISL_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-mpc-$(MINGW64_MPC_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-mpfr-$(MINGW64_MPFR_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-libiconv-$(MINGW64_LIBICONV_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-windows-default-manifest-$(MINGW64_WINDOWS_DEFAULT_MANIFEST_VERSION)-any.pkg.tar.zst \
	tarballs/mingw-w64-x86_64-gettext-runtime-$(MINGW64_GETTEXT_RUNTIME_VERSION)-any.pkg.tar.zst

gcc: $(GCC_MSYS2_DEPS)
	rm -rf $@
	mkdir $@
	for f in $(GCC_MSYS2_DEPS); do tar xf $$f --directory $@ $(MSYS2_EXCLUDES); done
	touch $@

#
# X13AS
#

tarballs/x13as_ascii-v$(X13AS_VERSION).zip:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@ https://www2.census.gov/software/x-13arima-seats/x13as/windows/program-archives/$(notdir $@)

x13as: tarballs/x13as_ascii-v$(X13AS_VERSION).zip
	rm -rf $@
	unzip -o $<
	mv $@/x13as_ascii.exe $@/x13as.exe
	touch $@
