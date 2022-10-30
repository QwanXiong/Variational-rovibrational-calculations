cmake_minimum_required(VERSION 3.10)
project(co2 CXX Fortran)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED 11)
set(CMAKE_CXX_FLAGS "-g -O1")
set(CMAKE_CXX_COMPILER mpic++)

set(EIGEN_DIR "/home/include/eigen-3.3.8")
set(SPECTRA_DIR "/home/include/spectra-1.0.1/include/")
set(MPI_DIR "/home/qq_ubuntu/include")
set(CVODE_DIR "/home/include")
set(LIBDIR "/home/lib")
include_directories(${EIGEN_DIR} ${CVODE_DIR} ${SPECTRA_DIR})

add_executable(
    co2
    co2.cpp
    VarSolver.cpp
    constants.cpp
    ames-2/co2potlongrange.f90
    )

target_link_libraries(co2
          gfortran
          gsl
          gslcblas
          )

target_link_directories(co2 PUBLIC ${LIBDIR})
