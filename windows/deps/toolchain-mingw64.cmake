# Inspired from https://www.mingw-w64.org/build-systems/cmake/

set(CMAKE_SYSTEM_NAME Windows)
set(CMAKE_SYSTEM_PROCESSOR x86_64)

# Specify the cross compiler
set(CMAKE_C_COMPILER x86_64-w64-mingw32ucrt-gcc)
set(CMAKE_CXX_COMPILER x86_64-w64-mingw32ucrt-g++)
set(CMAKE_Fortran_COMPILER x86_64-w64-mingw32ucrt-gfortran)

# See comments in windows/mingw-cross.ini
set(CMAKE_C_FLAGS_INIT "-march=nocona -msahf -mtune=generic")
set(CMAKE_CXX_FLAGS_INIT "-march=nocona -msahf -mtune=generic")
set(CMAKE_Fortran_FLAGS_INIT "-march=nocona -msahf -mtune=generic")

# Where is the target environment
set(CMAKE_FIND_ROOT_PATH /usr/x86_64-w64-mingw32ucrt)

# Search for programs in the build host directories
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
# For libraries and headers in the target directories
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
