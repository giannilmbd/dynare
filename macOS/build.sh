#!/usr/bin/env bash

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

set -ex
#exec > >(tee build-logfile.log) 2>&1 # uncomment for debugging

# Set root directory
ROOTDIR=$(dirname "$(readlink -f "$0")")/..

# Check that build directories do not already exist
[[ -d "$ROOTDIR"/build-macOS-matlab ]] && { echo "Please remove the build-macOS-matlab directory" 2>&1; exit 1; }
[[ -z $CI ]] && [[ -d "$ROOTDIR"/build-doc ]] && { echo "Please remove the build-doc directory" 2>&1; exit 1; }

cleanup()
{
    rm -rf "$ROOTDIR"/build-macOS-matlab
    [[ -n $CI ]] || rm -rf "$ROOTDIR"/build-doc
}
trap cleanup EXIT

##
## Set settings based on architecture
##
path_remove ()  { export "$1"="$(echo -n "${!1}" | awk -v RS=: -v ORS=: '$1 != "'"$2"'"' | sed 's/:$//')"; }
path_prepend () { path_remove "$1" "$2"; export "$1"="$2:${!1}"; }
PKG_ARCH=${1:-x86_64} # default to x86_64
if [[ "$PKG_ARCH" == arm64 ]]; then
    BREWDIR=/opt/homebrew
    # Make sure /opt/homebrew/bin is set first in PATH (as it might come last)
    path_prepend PATH /opt/homebrew/bin
    MATLAB_ARCH=maca64
else
    BREWDIR=/usr/local
    # Remove /opt/homebrew/bin from PATH, so it does not interfere with the x86_64 compilations
    path_remove PATH /opt/homebrew/bin
    MATLAB_ARCH=maci64
fi
MATLAB_PATH=/Applications/"$PKG_ARCH"/MATLAB_R2025b.app

# Append texbin to PATH to access latexmk and friends
path_prepend PATH /Library/TeX/texbin

# Set dependency directory
DEPS_DIR="$ROOTDIR"/macOS/deps/"$PKG_ARCH"/

GCC_VERSION=$(sed -En "/^c[[:space:]]*=/s/c[[:space:]]*=[[:space:]]*'.*gcc-([0-9]+)'/\1/p" "$ROOTDIR"/macOS/homebrew-native-"$PKG_ARCH".ini)

##
## Compile Dynare
##
cd "$ROOTDIR"

common_meson_opts=(-Dbuild_for=matlab --buildtype=release --prefer-static \
                   -Dfortran_args="[ '-B', '$DEPS_DIR/src/slicot-matlab/', '-B', '$DEPS_DIR/panua-pardiso/lib/' ]" \
                   --native-file macOS/homebrew-native-$PKG_ARCH.ini \
                   -Dpardiso=enabled -Dsuitesparse_src_path="macOS/deps/$PKG_ARCH/src/suitesparse")

# Build for MATLAB ⩾ R2020a (x86_64) and MATLAB ⩾ R2023b (arm64)
arch -"$PKG_ARCH" meson setup "${common_meson_opts[@]}" -Dmatlab_path="$MATLAB_PATH" build-macOS-matlab
arch -"$PKG_ARCH" meson compile -v -C build-macOS-matlab

# If not in CI, build the docs
if [[ -z $CI ]]; then
    arch -"$PKG_ARCH" meson compile -v -C build-macOS-matlab doc
    ln -s build-macOS-matlab build-doc
fi

## Strip binaries
strip build-macOS-matlab/preprocessor/src/dynare-preprocessor
strip -x -S build-macOS-matlab/*.mex"$MATLAB_ARCH"

##
## Create package
##

# Determine Dynare version if not passed by an environment variable as in the CI
if [[ -z $VERSION ]]; then
    cd build-macOS-matlab
    VERSION=$(meson introspect --projectinfo | sed -En 's/^.*"version": "([^"]*)".*$/\1/p')
    cd ..
fi

# Other useful variables
DATE=$(date +%Y-%m-%d-%H%M)
DATELONG=$(date '+%d %B %Y')

# Install location must not be too long for gcc.
# Otherwise, the headers of the compiled libraries cannot be modified
# obliging recompilation on the user's system.
# If VERSION is not a official release number, then do some magic
# to get something not too long, still more or less unique.
if [[ "$VERSION" =~ ^[0-9\.]+$ ]]; then
    LOCATION=$VERSION
else
    # Get the first component, truncate it to 5 characters, and add the date
    LOCATION=$(echo "$VERSION" | cut -f1 -d"-" | cut -c 1-5)-"$DATE"
fi
# Add architecture to LOCATION and VERSION
VERSION="$VERSION"-"$PKG_ARCH"
LOCATION="$LOCATION"-"$PKG_ARCH"

NAME=dynare-"$VERSION"
PKGFILES="$ROOTDIR"/macOS/pkg/"$NAME"
mkdir -p \
      "$PKGFILES"/preprocessor \
      "$PKGFILES"/doc \
      "$PKGFILES"/scripts \
      "$PKGFILES"/contrib/ms-sbvar/TZcode
if [[ "$PKG_ARCH" == x86_64 ]]; then
    mkdir -p "$PKGFILES"/mex/matlab/"$MATLAB_ARCH"-9.8-25.2
else
    mkdir -p "$PKGFILES"/mex/matlab/"$MATLAB_ARCH"-23.2-25.2
fi      

cp -p  "$ROOTDIR"/NEWS.md                                            "$PKGFILES"
cp -p  "$ROOTDIR"/COPYING                                            "$PKGFILES"
cp -p  "$ROOTDIR"/license.txt                                        "$PKGFILES"

cp -pr "$ROOTDIR"/matlab                                             "$PKGFILES"
cp -p  "$ROOTDIR"/build-macOS-matlab/dynare_version.m                "$PKGFILES"/matlab

cp -pr "$ROOTDIR"/examples                                           "$PKGFILES"

cp -p  "$ROOTDIR"/build-macOS-matlab/preprocessor/src/dynare-preprocessor "$PKGFILES"/preprocessor

# Create backward-compatibility symlink
mkdir -p                                                             "$PKGFILES"/matlab/preprocessor64
ln -sf ../../preprocessor/dynare-preprocessor                        "$PKGFILES"/matlab/preprocessor64/dynare_m

if [[ "$PKG_ARCH" == x86_64 ]]; then
    cp -L  "$ROOTDIR"/build-macOS-matlab/*.mex"$MATLAB_ARCH"         "$PKGFILES"/mex/matlab/"$MATLAB_ARCH"-9.8-25.2
else
    cp -L  "$ROOTDIR"/build-macOS-matlab/*.mex"$MATLAB_ARCH"         "$PKGFILES"/mex/matlab/"$MATLAB_ARCH"-23.2-25.2
fi

cp -p  "$ROOTDIR"/scripts/dynare.el                                  "$PKGFILES"/scripts
cp -pr "$ROOTDIR"/contrib/ms-sbvar/TZcode/MatlabFiles                "$PKGFILES"/contrib/ms-sbvar/TZcode

cp     "$ROOTDIR"/build-doc/*.pdf                                    "$PKGFILES"/doc
cp     "$ROOTDIR"/build-doc/preprocessor/doc/*.pdf                   "$PKGFILES"/doc
cp -r  "$ROOTDIR"/build-doc/dynare-manual.html                       "$PKGFILES"/doc

mkdir -p                                                             "$PKGFILES"/matlab/dseries/externals/x13/macOS/64
cp -p  "$DEPS_DIR"/src/x13as/x13as                                   "$PKGFILES"/matlab/dseries/externals/x13/macOS/64


cd "$ROOTDIR"/macOS/pkg

# Dynare option
arch -"$PKG_ARCH" pkgbuild --root "$PKGFILES" --identifier org.dynare."$VERSION" --version "$VERSION" --install-location /Applications/Dynare/"$LOCATION" "$NAME".pkg

# Create distribution.xml by replacing variables in distribution_template.xml
sed -e "s/VERSION_NO_SPACE/$VERSION/g" \
    -e "s/LOCATION/$LOCATION/g" \
    "$ROOTDIR"/macOS/distribution_template.xml > distribution.xml

# Create welcome.html by replacing variables in welcome_template.html
sed -e "s/VERSION_NO_SPACE/$VERSION/g" \
    -e "s/DATE/$DATELONG/g" \
    "$ROOTDIR"/macOS/welcome_template.html > "$ROOTDIR"/macOS/welcome.html

# Create conclusion.html by replacing variables in conclusion_template.html
sed -e "s/GCC_VERSION/$GCC_VERSION/g" \
    "$ROOTDIR"/macOS/conclusion_template.html > "$ROOTDIR"/macOS/conclusion.html

# Create installer
arch -"$PKG_ARCH" productbuild --distribution distribution.xml --resources "$ROOTDIR"/macOS --package-path ./"$NAME".pkg "$NAME"-productbuild.pkg

# Cleanup
rm -f ./distribution.xml
rm -rf "$PKGFILES"
rm -f "$ROOTDIR"/macOS/welcome.html
rm -f "$ROOTDIR"/macOS/conclusion.html

# Final pkg
mv "$NAME"-productbuild.pkg "$NAME".pkg
