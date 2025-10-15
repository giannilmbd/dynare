# Copyright © 2019-2025 Dynare Team
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

# settings for different architectures
DEPS_ARCH ?= x86_64 # use x86_64 by default

BREWDIR := $(if $(filter arm64,$(DEPS_ARCH)),/opt/homebrew,/usr/local)

GCC_VERSION = $(shell sed -En "/^c[[:space:]]*=/s/c[[:space:]]*=[[:space:]]*'.*gcc-([0-9]+)'/\1/p" ../homebrew-native-$(DEPS_ARCH).ini)

ROOT_PATH = $(realpath .)

WGET_OPTIONS := --no-verbose --no-use-server-timestamps --retry-connrefused --retry-on-host-error

NTHREADS=$(shell sysctl -n hw.perflevel0.physicalcpu)

.PHONY: all build build-slicot build-x13as clean clean-all

all: build

build: build-slicot build-x13as

# Clean everything except downloaded tarballs
clean:
	rm -rf $(DEPS_ARCH)/sources64 $(DEPS_ARCH)/lib64

# Clean everything (including downloaded tarballs)
clean-all: clean
	rm -rf tarballs

#
# slicot
#
tarballs/slicot-$(SLICOT_VERSION).tar.gz:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@ https://github.com/SLICOT/SLICOT-Reference/archive/refs/tags/v$(SLICOT_VERSION).tar.gz

$(DEPS_ARCH)/sources64/slicot-$(SLICOT_VERSION): tarballs/slicot-$(SLICOT_VERSION).tar.gz
	rm -rf $(DEPS_ARCH)/sources64/slicot-*
	mkdir -p $@
	tar xf $< --directory $@ --strip-components=1
	touch $@

$(DEPS_ARCH)/lib64/slicot/libslicot64_pic.a: $(DEPS_ARCH)/sources64/slicot-$(SLICOT_VERSION)
	gmake -C $< -f makefile_Unix FORTRAN=$(BREWDIR)/bin/gfortran LOADER=$(BREWDIR)/bin/gfortran SLICOTLIB=../libslicot64_pic.a OPTS="-O3 -fdefault-integer-8" lib -j$(NTHREADS)
	strip -S $</libslicot64_pic.a
	mkdir -p $(dir $@)
	cp $</libslicot64_pic.a $@

build-slicot: $(DEPS_ARCH)/lib64/slicot/libslicot64_pic.a


#
# X13AS
#
tarballs/x13as_asciisrc-v$(X13AS_VERSION).tar.gz:
	mkdir -p tarballs
	wget $(WGET_OPTIONS) -O $@ https://www2.census.gov/software/x-13arima-seats/x13as/unix-linux/program-archives/x13as_asciisrc-v$(X13AS_VERSION).tar.gz

$(DEPS_ARCH)/sources64/x13as-$(X13AS_VERSION): tarballs/x13as_asciisrc-v$(X13AS_VERSION).tar.gz
	rm -rf $(DEPS_ARCH)/sources64/x13as-*
	mkdir -p $@
	tar xf $< --directory $@ --strip-components=1

$(DEPS_ARCH)/lib64/x13as/x13as: $(DEPS_ARCH)/sources64/x13as-$(X13AS_VERSION)
	# Statically link x13as (see #1865).
	# Using -static is not possible, it does not work under Darwin.
	# Implement an ugly workaround in the absence of -static-libquadmath flag.
	# In particular, gcc must be used as the linker, because gfortran does not
	# honour static linking of libquadmath.
	# Once the -static-libquadmath flag is added (GCC 13? ; see
	# https://gcc.gnu.org/bugzilla/show_bug.cgi?id=46539), go back to using
	# gfortran as the linker with -static-libgfortran and
	# -static-libquadmath flags, and drop the GCC_VERSION variable.
	cd $< && sed -i '' 's/-static//g' makefile.gf
	gmake -C $< -f makefile.gf FC=$(BREWDIR)/bin/gfortran LINKER=$(BREWDIR)/bin/gcc-$(GCC_VERSION) FFLAGS="-O3 -std=legacy" LDFLAGS=-static-libgcc LIBS="$(BREWDIR)/lib/gcc/$(GCC_VERSION)/libgfortran.a $(BREWDIR)/lib/gcc/$(GCC_VERSION)/libquadmath.a" PROGRAM=x13as -j$(NTHREADS)

	strip $</x13as
	mkdir -p $(dir $@)
	cp $</x13as $@

build-x13as: $(DEPS_ARCH)/lib64/x13as/x13as
