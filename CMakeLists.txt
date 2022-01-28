cmake_minimum_required(VERSION 3.16.3)
project(FSIproject)

set(CMAKE_CXX_STANDARD 20)


add_executable(FSIproject
        main.cpp
        cpp/fluid/fvm_utilities.hpp cpp/fluid/fvm_utilities.cpp
        cpp/fluid/fvm_solver.cpp cpp/fluid/fvm_solver.hpp
        cpp/fsi/fsi_solver.cpp cpp/fsi/fsi_solver.hpp includes.hpp cpp/fluid/fvm_test.cpp cpp/fluid/fvm_test.hpp cpp/fluid/HLLC.cpp)


find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
add_subdirectory(pybind11)
pybind11_add_module(fluid_solver_py main.cpp)
pybind11_add_module(module_name main.cpp)
target_link_libraries(FSIproject ${PYTHON_LIBRARIES})
