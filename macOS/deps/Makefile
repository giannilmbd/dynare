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
PKG_ARCH ?= x86_64 # use x86_64 by default

BREWDIR := $(if $(filter arm64,$(PKG_ARCH)),/opt/homebrew,/usr/local)

ROOT_PATH = $(realpath .)

WGET_OPTIONS := --no-verbose --no-use-server-timestamps --retry-connrefused --retry-on-host-error

NTHREADS=$(shell sysctl -n hw.perflevel0.physicalcpu)

.PHONY: all clean clean-all

all: $(PKG_ARCH)/src/slicot-matlab/libslicot64_pic.a $(PKG_ARCH)/src/x13as/x13as

# Clean everything except downloaded tarballs
clean:
	rm -rf $(PKG_ARCH)/src

# Clean everything (including downloaded tarballs)
clean-all: clean
	rm -rf tarballs

#
# Slicot
#

tarballs/slicot-$(SLICOT_VERSION).tar.gz:
	mkdir -p $(dir $@)
	wget $(WGET_OPTIONS) -O $@ https://github.com/SLICOT/SLICOT-Reference/archive/refs/tags/v$(SLICOT_VERSION).tar.gz

$(PKG_ARCH)/src/slicot-matlab/libslicot64_pic.a: tarballs/slicot-$(SLICOT_VERSION).tar.gz
	rm -rf $(dir $@)
	mkdir -p $(dir $@)
	tar xf $< --directory $(dir $@) --strip-components=1
	gmake -C $(dir $@) -f makefile_Unix lib SLICOTLIB=../$(notdir $@) OPTS="-O3 -fdefault-integer-8" FORTRAN=$(BREWDIR)/bin/gfortran LOADER=$(BREWDIR)/bin/gfortran -j$(NTHREADS)


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
