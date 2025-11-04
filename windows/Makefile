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

.PHONY: build deps clean-deps clean-all-deps clean-pkgs clean clean-all install

NTHREADS=$(shell nproc)

build: deps
	./build.sh

deps:
	make -C deps -j$(NTHREADS)

clean-deps:
	make -C deps clean

clean-all-deps:
	make -C deps clean-all

clean-pkgs:
	rm -rf exe/ 7z/

# Clean everything except downloaded tarballs of dependencies
clean: clean-pkgs clean-deps

# Clean everything (including downloaded tarballs of dependencies)
clean-all: clean-pkgs clean-all-deps

install:
	./install-packages.sh
