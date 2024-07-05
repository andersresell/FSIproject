
cmake_minimum_required(VERSION 3.20)
set(CMAKE_CXX_STANDARD 20)

#set(CMAKE_CXX_COMPILER g++)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE "Release")
endif()


# Add source files
set(SRC
main.cpp
simulate.cpp
input_parser.cpp
fluid/fvm_solver.cpp
fluid/fvm_utilities.cpp
fluid/fvm_bcs.cpp
fsi/fsi_solver.cpp
fsi/setup_cases.cpp
solid/dynamic_rigid.cpp
solid/solid_body.cpp
solid/solid_utilities.cpp
)

project(FSIproject)
# Compile sources
add_executable(${PROJECT_NAME} ${SRC})
# Set compiler flags
target_compile_options(${PROJECT_NAME} PRIVATE
    -Wall -Wno-reorder -Wno-unused-function -fmax-errors=1
)
if (CMAKE_BUILD_TYPE STREQUAL "release")
target_compile_options(${PROJECT_NAME} PRIVATE
    -DNDEBUG -O3 -march=native #-ffast-math #-flto 
)
elseif (CMAKE_BUILD_TYPE STREQUAL "debug")
target_compile_options(${PROJECT_NAME} PRIVATE
    -g -D_GLIBCXX_ASSERTIONS
)
else()
    message(FATAL_ERROR "Illegal build type")
endif()

# Find yaml-cpp package
find_package(yaml-cpp REQUIRED)

# Link yaml-cpp
target_link_libraries(${PROJECT_NAME} PRIVATE yaml-cpp)

# Link OpenMP
#if(USE_OPENMP)
find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
#endif()

# Set include directories
include_directories(${FSI_PROJECT_HOME}/include)

# Set output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

