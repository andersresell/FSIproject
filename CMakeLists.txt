cmake_minimum_required(VERSION 3.21)
project(FSIproject)

set(CMAKE_CXX_STANDARD 14)

add_executable(FSIproject
        main.cpp
        cpp/fluid/fvm_utilities.hpp cpp/fluid/fvm_utilities.cpp
        cpp/fluid/fvm_solver.cpp cpp/fluid/fvm_solver.hpp
        cpp/fsi/fsi_solver.cpp cpp/fsi/fsi_solver.hpp includes.hpp cpp/fluid/fvm_test.cpp cpp/fluid/fvm_test.hpp)
