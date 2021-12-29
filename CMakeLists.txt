cmake_minimum_required(VERSION 3.21)
project(FSIproject)

set(CMAKE_CXX_STANDARD 14)

add_executable(FSIproject
        main.cpp
        cpp/fluid/fvm_utilities.hpp cpp/fluid/fvm_utilities.cpp
        cpp/fluid/fvm_step.cpp cpp/fluid/fvm_step.hpp
        cpp/fsi/fsi.cpp cpp/fsi/fsi.h)
