cmake_minimum_required(VERSION 3.20)
project(FSIproject)

set(CMAKE_CXX_STANDARD 20)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Build)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS_DEBUG "-O3")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}) #setting the output directory to FSIproject folder instead of cmake-build-debug

add_executable(FSIproject
        main.cpp
        cpp/fluid/fvm_utilities.hpp cpp/fluid/fvm_utilities.cpp
        cpp/fluid/fvm_solver.cpp cpp/fluid/fvm_solver.hpp
        cpp/fsi/fsi_solver.cpp cpp/fsi/fsi_solver.hpp
        includes.hpp
        cpp/fluid/HLLC.cpp
        cpp/solid/solid_utilities.cpp cpp/solid/solid_utilities.hpp cpp/solid/solid_body.cpp cpp/solid/solid_body.hpp cpp/fsi/simulate.hpp cpp/fsi/setup_cases.cpp cpp/fsi/setup_cases.hpp cpp/fsi/simulate.cpp cpp/solid/dynamic_rigid.cpp cpp/solid/dynamic_rigid.hpp)

SET( EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
SET( CMAKE_MODULE_PATH "$ENV{EIGEN3_DIR}" )

find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

find_package(Boost)
if(Boost_FOUND)
    message(STATUS "Boost found")
    include_directories(${Boost_INCLUDE_DIRS})
endif()


