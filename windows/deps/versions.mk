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
MSYS2_BOOST_VERSION = 1.90.0-3

## Packages for the embedded compiler

# pacman -Ss mingw-w64-ucrt-x86_64-gcc$
MSYS2_GCC_VERSION = 15.2.0-11

# Dependency of gcc, isl, mpc and mpfr
# pacman -Ss mingw-w64-ucrt-x86_64-gmp
MSYS2_GMP_VERSION = 6.3.0-2

# pacman -Ss mingw-w64-ucrt-x86_64-binutils
MSYS2_BINUTILS_VERSION = 2.46-2

# pacman -Ss mingw-w64-ucrt-x86_64-headers
MSYS2_HEADERS_VERSION = 13.0.0.r505.g7d006b2ea-1

# pacman -Ss mingw-w64-ucrt-x86_64-crt
MSYS2_CRT_VERSION = 13.0.0.r505.g7d006b2ea-1

# pacman -Ss mingw-w64-ucrt-x86_64-winpthreads
MSYS2_WINPTHREADS_VERSION = 13.0.0.r505.g7d006b2ea-1

# pacman -Ss mingw-w64-ucrt-x86_64-isl
MSYS2_ISL_VERSION = 0.27-1

# pacman -Ss mingw-w64-ucrt-x86_64-mpc
MSYS2_MPC_VERSION = 1.3.1-2

# Dependency of mpc
# pacman -Ss mingw-w64-ucrt-x86_64-mpfr
MSYS2_MPFR_VERSION = 4.2.2-1

# pacman -Ss mingw-w64-ucrt-x86_64-windows-default-manifest
MSYS2_WINDOWS_DEFAULT_MANIFEST_VERSION = 6.4-4

# pacman -Ss mingw-w64-ucrt-x86_64-zlib
MSYS2_ZLIB_VERSION = 1.3.1-1

# pacman -Ss mingw-w64-ucrt-x86_64-zstd
MSYS2_ZSTD_VERSION = 1.5.7-1

# Dependency of binutils
# pacman -Ss mingw-w64-ucrt-x86_64-gettext-runtime
MSYS2_GETTEXT_RUNTIME_VERSION = 1.0-1

# Dependency of gettext-runtime
# pacman -Ss mingw-w64-ucrt-x86_64-libiconv
MSYS2_LIBICONV_VERSION = 1.18-1

