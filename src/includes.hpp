#pragma once

#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <omp.h>
#include <set>
#include <sstream>
#include <vector>

#include <eigen3/Eigen/Dense>

using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::make_shared;
using std::make_unique;
using std::map;
using std::move;
using std::pair;
using std::runtime_error;
using std::shared_ptr;
using std::string;
using std::to_string;
using std::unique_ptr;
using std::vector;

#ifndef _OPENMP
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#define omp_set_num_threads(num_threads)
#endif
