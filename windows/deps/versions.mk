SLICOT_VERSION = 5.9.1
X13AS_VERSION = 1-1-b62

OCTAVE_VERSION = 10.3.0
OCTAVE_W64_BUILD =

MATLAB64_VERSION = 20250618

PANUA_PARDISO_VERSION = 20240630
SUITESPARSE_VERSION = 7.11.0

### MSYS2 packages
# Determine the versions by:
# - first running: pacman -Sy
# - and then with appropriate queries using: pacman -Ss <regex>
# Dependencies can be determined using: pacman -Si <pkg>
# File lists can be determined using: pacman -Fl <pkg>
# The same information can be gathered from: https://packages.msys2.org/search

## Build dependencies

# pacman -Ss mingw-w64-ucrt-x86_64-boost
# (also used by the preprocessor CI, so it can be updated there also)
MSYS2_BOOST_VERSION = 1.90.0-2

# pacman -Ss mingw-w64-ucrt-x86_64-gsl
MSYS2_GSL_VERSION = 2.8-1

# pacman -Ss mingw-w64-ucrt-x86_64-matio
MSYS2_MATIO_VERSION = 1.5.29-1

# Dependency of matio and libssh2 (and of the MinGW compiler)
# pacman -Ss mingw-w64-ucrt-x86_64-zlib
MSYS2_ZLIB_VERSION = 1.3.1-1

# Dependency of matio
# pacman -Ss mingw-w64-ucrt-x86_64-hdf5
MSYS2_HDF5_VERSION = 1.14.6-3

# Dependency of HDF5 (provides szip library)
# pacman -Ss mingw-w64-ucrt-x86_64-libaec
MSYS2_LIBAEC_VERSION = 1.1.4-1

# Dependency of HDF5 and libssh2
# pacman -Ss mingw-w64-ucrt-x86_64-openssl
MSYS2_OPENSSL_VERSION = 3.6.0-1

# Dependency of HDF5
# pacman -Ss mingw-w64-ucrt-x86_64-curl
MSYS2_CURL_VERSION = 8.17.0-1

# Dependency of curl (and of the MinGW compiler)
# pacman -Ss mingw-w64-ucrt-x86_64-zstd
MSYS2_ZSTD_VERSION = 1.5.7-1

# Dependency of curl
# pacman -Ss mingw-w64-ucrt-x86_64-brotli
MSYS2_BROTLI_VERSION = 1.2.0-1

# Dependency of curl
# pacman -Ss mingw-w64-ucrt-x86_64-libpsl
MSYS2_LIBPSL_VERSION = 0.21.5-3

# Dependency of curl and of libpsl
# pacman -Ss mingw-w64-ucrt-x86_64-libidn2
MSYS2_LIBIDN2_VERSION = 2.3.8-4

# Dependency of curl
# pacman -Ss mingw-w64-ucrt-x86_64-libssh2
MSYS2_LIBSSH2_VERSION = 1.11.1-1

# Dependency of curl
# pacman -Ss mingw-w64-ucrt-x86_64-nghttp2
MSYS2_NGHTTP2_VERSION = 1.68.0-1

# Dependency of curl
# pacman -Ss mingw-w64-ucrt-x86_64-nghttp3
MSYS2_NGHTTP3_VERSION = 1.13.1-1

# Dependency of curl
# pacman -Ss mingw-w64-ucrt-x86_64-ngtcp2
MSYS2_NGTCP2_VERSION = 1.18.0-1

# Dependency of libpsl and libunistring (and of gettext-runtime for the MinGW compiler)
# pacman -Ss mingw-w64-ucrt-x86_64-libiconv
MSYS2_LIBICONV_VERSION = 1.18-1

# Dependency of libpsl and libidn2
# pacman -Ss mingw-w64-ucrt-x86_64-libunistring
MSYS2_LIBUNISTRING_VERSION = 1.3-1

## Packages for the embedded compiler

# pacman -Ss mingw-w64-ucrt-x86_64-gcc$
MSYS2_GCC_VERSION = 15.2.0-9

# Dependency of gcc, isl, mpc and mpfr
# pacman -Ss mingw-w64-ucrt-x86_64-gmp
MSYS2_GMP_VERSION = 6.3.0-2

# pacman -Ss mingw-w64-ucrt-x86_64-binutils
MSYS2_BINUTILS_VERSION = 2.45.1-1

# pacman -Ss mingw-w64-ucrt-x86_64-headers
MSYS2_HEADERS_VERSION = 13.0.0.r391.g848cce552-1

# pacman -Ss mingw-w64-ucrt-x86_64-crt
MSYS2_CRT_VERSION = 13.0.0.r391.g848cce552-1

# pacman -Ss mingw-w64-ucrt-x86_64-winpthreads
MSYS2_WINPTHREADS_VERSION = 13.0.0.r391.g848cce552-1

# pacman -Ss mingw-w64-ucrt-x86_64-isl
MSYS2_ISL_VERSION = 0.27-1

# pacman -Ss mingw-w64-ucrt-x86_64-mpc
MSYS2_MPC_VERSION = 1.3.1-2

# Dependency of mpc
# pacman -Ss mingw-w64-ucrt-x86_64-mpfr
MSYS2_MPFR_VERSION = 4.2.2-1

# pacman -Ss mingw-w64-ucrt-x86_64-windows-default-manifest
MSYS2_WINDOWS_DEFAULT_MANIFEST_VERSION = 6.4-4

# Dependency of binutils
# pacman -Ss mingw-w64-ucrt-x86_64-gettext-runtime
MSYS2_GETTEXT_RUNTIME_VERSION = 0.26-2
