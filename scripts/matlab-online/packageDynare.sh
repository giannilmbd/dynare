#!/bin/bash
set -exo pipefail

# Creates a dynare-X.Y.mltbx in the current repository, using the settings
# below. Needs to be run from Ubuntu, with the needed packages installed.
# The required Ubuntu version can be obtained by running “!lsb_release -a” in
# MATLAB Online.

X13ASVER=1-1-b62
MATLABVER=R2024b

MATLABPATH=/opt/MATLAB/${MATLABVER}
# TODO: change size and put white background for better rendering in MATLAB Add-Ons browser
DYNARE_PNG_LOGO=../../preprocessor/doc/logos/dlogo.png

# Prepare temporary workspace and setup cleanup function
tmpdir=$(mktemp -d)
cleanup ()
{
    rm -rf "$tmpdir"
}
trap cleanup EXIT

pushd ../..
meson setup -Dbuild_for=matlab -Dmatlab_path="$MATLABPATH" -Dbuildtype=release -Dprefer_static=true "$tmpdir"/build-matlab-online

cd "$tmpdir"/build-matlab-online
meson compile -v
meson install --destdir "$tmpdir"
DYNAREVER=$(meson introspect --projectinfo | jq -r '.version')

# Sanitize the version number so that is corresponds to MATLAB toolbox
# requirements: the version must be a scalar string or character vector of the
# form Major.Minor.Bug.Build, where Bug and Build are optional.
# Hence remove any character which is not a number or a dot, and ensure that we
# have at least a minor version number.
DYNAREVER_SANITIZED=${DYNAREVER//[^0-9.]/}
if [[ ${DYNAREVER_SANITIZED} != *.* ]]; then
    DYNAREVER_SANITIZED=${DYNAREVER_SANITIZED}.0
fi

cd ..
strip usr/local/bin/dynare-preprocessor
strip usr/local/lib/dynare/mex/matlab/*.mexa64

# Get X13 binary from the Census Bureau website
# The binary from Ubuntu has some shared library dependencies, so it is safer to use a static binary
wget --no-verbose --retry-connrefused --retry-on-host-error https://www2.census.gov/software/x-13arima-seats/x13as/unix-linux/program-archives/x13as_ascii-v${X13ASVER}.tar.gz
tar xf x13as_ascii-v${X13ASVER}.tar.gz

# Populate staging area for the zip
cp -pRL usr/local/lib/dynare dynare # -L is needed to dereference the preprocessor symlink
mkdir -p dynare/matlab/dseries/externals/x13/linux/64
cp -p x13as/x13as_ascii dynare/matlab/dseries/externals/x13/linux/64/x13as

# zip dynare
cd dynare
zip -q -r "$tmpdir"/dynare.zip *

# make toolbox
popd
"$MATLABPATH/bin/matlab" -batch "packageDynare('$tmpdir/dynare.zip', '$DYNAREVER', '$DYNAREVER_SANITIZED', '$DYNARE_PNG_LOGO')"
