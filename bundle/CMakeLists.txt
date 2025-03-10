# ------------------------------------------------------------------------- #
# JEDI GDAS Bundle             #
# ------------------------------------------------------------------------- #

# Check for minimim cmake requirement
cmake_minimum_required( VERSION 3.20 FATAL_ERROR )

find_package(ecbuild 3.5 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild)

project(DA-Utils-bundle VERSION 1.0.0 LANGUAGES C CXX Fortran )

include(GNUInstallDirs)
enable_testing()

# Build type.
if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()

# Include ecbuild_bundle macro
include( ecbuild_bundle )

# Enable MPI
set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )

# Depend path for non-ecbuild packages
set(DEPEND_LIB_ROOT ${CMAKE_CURRENT_BINARY_DIR}/Depends)
list(APPEND CMAKE_PREFIX_PATH ${DEPEND_LIB_ROOT})

# Library path for non-ecbuild packages
link_directories(${CMAKE_CURRENT_BINARY_DIR}/lib)

include( GNUInstallDirs )
if(APPLE)
  list( APPEND CMAKE_INSTALL_RPATH $ENV{llvm_openmp_ROOT}/lib )
endif()
list( APPEND CMAKE_INSTALL_RPATH ${CMAKE_CURRENT_BINARY_DIR}/fv3 )

# Initialize bundle
# -----------------
ecbuild_bundle_initialize()

ecbuild_bundle( PROJECT da-utils SOURCE "../" )
# Finalize bundle
# ---------------
ecbuild_bundle_finalize()
