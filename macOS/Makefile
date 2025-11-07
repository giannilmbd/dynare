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

ARCHS = arm64 x86_64
BUILDS = $(addprefix build-, $(ARCHS))
DEPS = $(addprefix deps-, $(ARCHS))
CLEAN_DEPS = $(addprefix clean-deps-, $(ARCHS))
CLEAN_DEPS_ALL = $(addprefix clean-deps-all-, $(ARCHS))

.PHONY: build $(BUILDS) $(DEPS) $(CLEAN_DEPS) $(CLEAN_DEPS_ALL) clean-pkgs clean clean-all

build: $(BUILDS)

$(BUILDS): build-%: deps-%
	./build.sh $*

$(DEPS): deps-%:
	gmake -C deps PKG_ARCH=$*

$(CLEAN_DEPS): clean-deps-%:
	gmake -C deps clean PKG_ARCH=$*

$(CLEAN_DEPS_ALL): clean-deps-all-%:
	gmake -C deps clean-all PKG_ARCH=$*

clean-pkgs:
	rm -rf pkg

# Clean everything except downloaded tarballs of dependencies
clean: clean-pkgs $(CLEAN_DEPS)

# Clean everything (including downloaded tarballs of dependencies)
clean-all: clean-pkgs $(CLEAN_DEPS_ALL)
