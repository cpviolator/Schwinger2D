cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# set the project name and version
project(Schwinger2D VERSION 1.0.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add a directory for cmake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(ExternalProject)
# Link OpempMP CMake automatically adds 
# _OPENMP to the CXX flags
option(ENABLE_OPENMP "enable OpenMP" ON)
if(ENABLE_OPENMP)
  find_package(OpenMP)
endif()

find_package(MPFR REQUIRED)

# Add Eigen for small, dense algebra
set(EIGEN_VERSION 3.3.7)
set(EIGEN_DOWNLOAD_LOCATION ${CMAKE_SOURCE_DIR}/externals/eigen/${EIGEN_VERSION}.tar.bz2)
set(EIGEN_URL https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/eigen-${EIGEN_VERSION}.tar.bz2)
set(EIGEN_SHA 685adf14bd8e9c015b78097c1dc22f2f01343756f196acdc76a678e1ae352e11)
if(NOT EXISTS ${EIGEN_DOWNLOAD_LOCATION})
  message(STATUS "Checking for Eigen tarball and downloading if necessary.")
endif()
file(
  DOWNLOAD ${EIGEN_URL} ${EIGEN_DOWNLOAD_LOCATION}
  EXPECTED_HASH SHA256=${EIGEN_SHA}
  STATUS EIGEN_DOWNLOADED
  )
list(GET EIGEN_DOWNLOADED 0 EIGEN_DOWNLOADED_CODE)
list(GET EIGEN_DOWNLOADED 1 EIGEN_DOWNLOADED_MSG)
if(${EIGEN_DOWNLOADED_CODE})
  message( SEND_ERROR "Could not download Eigen automatically (${EIGEN_DOWNLOADED_MSG}). "
    "Please download eigen from ${EIGEN_URL} and save it to ${EIGEN_DOWNLOAD_LOCATION} "
    "and try running cmake again.")
endif()

ExternalProject_Add(
  Eigen
  URL ${CMAKE_SOURCE_DIR}/externals/eigen/${EIGEN_VERSION}.tar.bz2
  URL_HASH SHA256=${EIGEN_SHA}
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/externals/eigen/
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_Property(Eigen source_dir)
set(EIGEN_INCLUDE_DIRS ${source_dir})

# Add HDF5 for experimental binary IO
option(ENABLE_HDF5 "enable HDF5" OFF)
option(DOWNLOAD_HDF5 "download HDF5" OFF)

if(ENABLE_HDF5)
  if(NOT DOWNLOAD_HDF5)
    find_package(HDF5 REQUIRED COMPONENTS CXX)
    include_directories(${HDF5_INCLUDE_DIR})
  else()
    include(GNUInstallDirs)
    ExternalProject_Add(
      HDF5
      GIT_REPOSITORY https://github.com/HDFGroup/hdf5.git
      GIT_TAG hdf5-1_12_0
      GIT_SHALLOW YES
      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/externals/hdf5
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
      -DHDF5_BUILD_EXAMPLES=OFF
      -DBUILD_TESTING=OFF
      -DONLY_SHARED_LIBS=ON
      -DHDF5_BUILD_CPP_LIB=ON
      -DHDF5_BUILD_C_LIB=ON
      #HDF5_ENABLE_SANITIZERS
      CMAKE_GENERATOR "Unix Makefiles")
    ExternalProject_Get_Property(HDF5 INSTALL_DIR)
    set(HDF5_HOME ${INSTALL_DIR})
    add_library(hdf5 STATIC IMPORTED)
    add_dependencies(hdf5 HDF5)
    set_target_properties(hdf5 PROPERTIES IMPORTED_LOCATION ${HDF5_HOME}/${CMAKE_INSTALL_LIBDIR}/libhdf5.a)
    
  endif()
endif()

option(ENABLE_ALG_REMEZ "enable AlgRemez" OFF)
if(ENABLE_ALG_REMEZ)
  ExternalProject_Add(
    ALG_REMEZ
    GIT_REPOSITORY https://github.com/maddyscientist/AlgRemez.git
    GIT_TAG ceee4425139bc507ed3ae34ced847c544a8e851f
    GIT_SHALLOW YES
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/externals/AlgRemez
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    CMAKE_GENERATOR "Unix Makefiles")
  
  add_library(algremez STATIC IMPORTED)
  add_dependencies(algremez ALG_REMEZ)
  set_target_properties(algremez PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/externals/AlgRemez/lib/libalgremez.a)
  
  ExternalProject_Get_Property(ALG_REMEZ source_dir)
  set(ALG_REMEZ_INCLUDE_DIRS ${source_dir}/include)
endif()

add_subdirectory(lib)
add_subdirectory(wilson)
