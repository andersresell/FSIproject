cmake_minimum_required(VERSION 3.21)
project(FSIproject)

set(CMAKE_CXX_STANDARD 14)

add_executable(FSIproject main.cpp cpp/fluid/structuredFVM.hpp cpp/fluid/structuredFVM.cpp)
