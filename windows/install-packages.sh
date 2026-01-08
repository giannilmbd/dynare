#!/bin/bash

# On a Debian system, install the packages needed for Windows
# cross-compilation, and also setup the cross-compiler alternatives.

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

[[ $(id -u) == 0 ]] || { echo "You must be root" >&2; exit 1; }

PACKAGES=(make 7zip zstd wget meson pkg-config-mingw-w64-ucrt64
          gcc-mingw-w64-ucrt64 g++-mingw-w64-ucrt64
          gfortran-mingw-w64-ucrt64 flex libfl-dev bison texlive
          texlive-publishers texlive-latex-extra texlive-science
          texlive-fonts-extra lmodern cm-super python3-sphinx latexmk nsis
          cmake)

apt install "${PACKAGES[@]}"
