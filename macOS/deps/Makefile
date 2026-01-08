# Copyright © 2019-2026 Dynare Team
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
PKG_ARCH ?= x86_64 # use x86_64 by default

BREWDIR := $(if $(filter arm64,$(PKG_ARCH)),/opt/homebrew,/usr/local)

ROOT_PATH = $(realpath .)

WGET_OPTIONS := --no-verbose --no-use-server-timestamps --retry-connrefused --retry-on-host-error

NTHREADS=$(shell sysctl -n hw.perflevel0.physicalcpu)

.PHONY: all clean clean-all

all: $(PKG_ARCH)/src/slicot-matlab/build/lib/libslicot64_pic.a $(PKG_ARCH)/src/x13as/x13as $(PKG_ARCH)/panua-pardiso $(PKG_ARCH)/src/suitesparse

# Clean everything except downloaded tarballs
clean:
	rm -rf $(PKG_ARCH)/src $(PKG_ARCH)/panua-pardiso

# Clean everything (including downloaded tarballs)
clean-all: clean
	rm -rf tarballs

#
# Slicot
#

tarballs/slicot-$(SLICOT_VERSION).tar.gz:
	mkdir -p $(dir $@)
	wget $(WGET_OPTIONS) -O $@ https://github.com/SLICOT/SLICOT-Reference/archive/refs/tags/v$(SLICOT_VERSION).tar.gz

$(PKG_ARCH)/src/slicot-matlab/build/lib/libslicot64_pic.a: tarballs/slicot-$(SLICOT_VERSION).tar.gz
	rm -rf $(realpath $(dir $@)/../..)
	mkdir -p $(dir $@)/../..
	tar xf $< --directory $(dir $@)/../.. --strip-components=1
	patch -d $(dir $@)/../.. -p1 < ../../windows/deps/slicot-no-blas-lapack-lookup.patch
	cmake -S $(dir $@)/../.. -B $(dir $@)/.. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DSLICOT_TESTING=OFF -DSLICOT_INTEGER8=ON
	cmake --build $(dir $@)/.. --verbose
	mv $(dir $@)/libslicot64.a $@


#
# X13AS
#

tarballs/x13as_asciisrc-v$(X13AS_VERSION).tar.gz:
	mkdir -p $(dir $@)
	wget $(WGET_OPTIONS) -O $@ https://www2.census.gov/software/x-13arima-seats/x13as/unix-linux/program-archives/$(notdir $@)

$(PKG_ARCH)/src/x13as/x13as: tarballs/x13as_asciisrc-v$(X13AS_VERSION).tar.gz
	rm -rf $(dir $@)
	mkdir -p $(dir $@)
	tar xf $< --directory $(dir $@) --strip-components=1
	# Statically link x13as (see #1865).
	# Using -static is not possible, it does not work under Darwin.
	cd $(dir $@) && sed -i '' 's/-static//g' makefile.gf
	gmake -C $(dir $@) -f makefile.gf FC=$(BREWDIR)/bin/gfortran LINKER=$(BREWDIR)/bin/gfortran FFLAGS="-O3 -std=legacy" LDFLAGS="-static-libgcc -static-libgfortran -static-libquadmath" PROGRAM=x13as -j$(NTHREADS)
	strip $@

#
# Panua PARDISO
#

$(PKG_ARCH)/panua-pardiso: tarballs/panua-pardiso-$(PANUA_PARDISO_VERSION)-mac_$(if $(filter x86_64,$(PKG_ARCH)),x86,arm64).zip
	rm -rf $@
	unzip -o $<
	mv $(basename $(notdir $<)) $@
	touch $@

#
# SuiteSparse
#

tarballs/SuiteSparse-$(SUITESPARSE_VERSION).tar.gz:
	wget $(WGET_OPTIONS) -O $@ https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v$(SUITESPARSE_VERSION).tar.gz

$(PKG_ARCH)/src/suitesparse: tarballs/SuiteSparse-$(SUITESPARSE_VERSION).tar.gz
	rm -rf $@
	mkdir -p $@
	tar xf $< --directory $@ --strip-components=1
	touch $@
