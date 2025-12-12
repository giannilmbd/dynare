#!/bin/bash
set -exo pipefail

# Creates a dynare-X.Y.mltbx in the current repository, using the settings
# below. Needs to be run from Ubuntu, with the needed packages installed.
# The required Ubuntu version can be obtained by running “!lsb_release -a” in
# MATLAB Online.

X13AS_VERSION=1-1-b62
MATLAB_VERSION=R2025b

MATLAB_PATH=/opt/MATLAB/${MATLAB_VERSION}
# TODO: change size and put white background for better rendering in MATLAB Add-Ons browser
DYNARE_PNG_LOGO=../../preprocessor/doc/logos/dlogo.png

# Prepare temporary workspace and setup cleanup function
tmpdir=$(mktemp -d)
cleanup ()
{
    rm -rf "$tmpdir"
}
trap cleanup EXIT

# Get X13 binary from the Census Bureau website if not cached, and unpack it into the temporary workspace
# The binary from Ubuntu has some shared library dependencies, so it is safer to use a static binary
X13AS_TARBALL=x13as_ascii-v${X13AS_VERSION}.tar.gz
if [[ ! -f ${X13AS_TARBALL} ]]; then
    wget --no-verbose --retry-connrefused --retry-on-host-error https://www2.census.gov/software/x-13arima-seats/x13as/unix-linux/program-archives/${X13AS_TARBALL}
fi
tar -xf ${X13AS_TARBALL} -C "${tmpdir}"

pushd ../..
meson setup -Dbuild_for=matlab -Dmatlab_path="$MATLAB_PATH" --buildtype=release -Db_lto=true --prefer-static "$tmpdir"/build-matlab-online

cd "$tmpdir"/build-matlab-online
meson compile -v
meson install --destdir "$tmpdir"
DYNARE_VERSION=$(meson introspect --projectinfo | jq -r '.version')

# Sanitize the version number so that is corresponds to MATLAB toolbox
# requirements: the version must be a scalar string or character vector of the
# form Major.Minor.Bug.Build, where Bug and Build are optional.
# Hence remove any character which is not a number or a dot, and ensure that we
# have at least a minor version number.
DYNARE_VERSION_SANITIZED=${DYNARE_VERSION//[^0-9.]/}
if [[ ${DYNARE_VERSION_SANITIZED} != *.* ]]; then
    DYNARE_VERSION_SANITIZED=${DYNARE_VERSION_SANITIZED}.0
fi

cd ..
strip usr/local/bin/dynare-preprocessor
strip usr/local/lib/dynare/mex/matlab/*.mexa64

# Populate staging area for the .mltbx archive
cp -pRL usr/local/lib/dynare dynare # -L is needed to dereference the preprocessor symlink
mkdir -p dynare/matlab/dseries/externals/x13/linux/64
cp -p x13as/x13as_ascii dynare/matlab/dseries/externals/x13/linux/64/x13as

# make toolbox
popd
"$MATLAB_PATH/bin/matlab" -batch "packageDynare('$tmpdir/dynare', '$DYNARE_VERSION', '$DYNARE_VERSION_SANITIZED', '$DYNARE_PNG_LOGO')"
