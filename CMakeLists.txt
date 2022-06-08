# Copyright 2022 Jij Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(cimod)

include(CheckLanguage)

check_language(C)
check_language(CXX)

if(CMAKE_C_COMPILER)
  enable_language(C)
endif() 
if(CMAKE_CXX_COMPILER)
  enable_language(CXX)
endif() 

if(UNIX AND NOT APPLE)
  check_language(Fortran)
  if(CMAKE_Fortran_COMPILER)
      enable_language(Fortran)
  endif()
endif()

# detect operating system and host processor
message(STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
message(STATUS "The host processor is ${CMAKE_HOST_SYSTEM_PROCESSOR}")
message(STATUS "CMAKE_GENERATOR = ${CMAKE_GENERATOR}")

message(STATUS "Started CMake for ${PROJECT_NAME} \n")

if (APPLE)
    set(CMAKE_MACOSX_RPATH TRUE)
    add_compile_options(-Wno-unused-command-line-argument)
elseif (MSVC)
    add_compile_options(/bigobj)
    set(CMAKE_Fortran_COMPILER OFF CACHE BOOL "" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_VERBOSE_MAKEFILE TRUE)

##### Set default behavior #####
set(DEFAULT_USE_OMP Yes)

if(NOT APPLE) 
  if(NOT DEFINED USE_OMP) 
    set(USE_OMP ${DEFAULT_USE_OMP})
  endif() 
else() 
  if(FORCE_USE_OMP) 
    set(USE_OMP ${DEFAULT_USE_OMP})
    option(FORCE_USE_OMP "Enable Use OMP FORCE." OFF)
    message(STATUS "APPLE_FORCE_USE_OMP = ${FORCE_USE_OMP}")
  endif()
endif()

if(USE_OMP)
  if(WIN32)
    set(CMAKE_INSTALL_OPENMP_LIBRARIES ON)
  endif()
  find_package(OpenMP REQUIRED)
  message(STATUS "OMP_FOUND = YES")
  message(STATUS "USE_OMP = ${USE_OMP}")
endif()

message(STATUS "CIBUILDWHEEL = $ENV{CIBUILDWHEEL}")
set(CIMOD_BLA_STATIC $ENV{CIBUILDWHEEL})
if(CIMOD_BLA_STATIC)
  message(STATUS "ENABLE BLAS STATIC")
  set(BLA_STATIC ON CACHE BOOL "" FORCE)
endif()

find_package(BLAS)
find_package(LAPACK)

if(BLAS_FOUND) 
  set(EIGEN_USE_BLAS ON) 
endif()

if(LAPACK_FOUND) 
  set(EIGEN_USE_LAPACKE ON)
endif() 

list(APPEND CMAKE_MODULE_PATH external)

find_package(Eigen3 CONFIG)

if(EIGEN3_FOUND)
  message(STATUS "Skip Download eigen3")
else() 
  message(STATUS "Downlod eigen3")
  include(external/eigen.cmake)
endif() 

find_package(nlohmann_json)
if(TARGET nlohmann_json::nlohmann_json)
  message(STATUS "Skip Download nlohmann_json")
else()
  message(STATUS "Downlod nlohmann_json")
  set(JSON_BuildTests OFF CACHE INTERNAL "")
  include(external/json.cmake)
endif()

if(ENABLE_COVERAGE)
  list(APPEND CMAKE_MODULE_PATH cmake)
  message(STATUS "Use Codecov")
  find_package(codecov)
endif()

add_subdirectory(src)

# Build Python Extension 
if(SKBUILD)
  message(STATUS "Build Python Extension.")
  find_package(Python COMPONENTS Interpreter Development)
  # Scikit-Build does not add your site-packages to the search path
  # automatically, so we need to add it _or_ the pybind11 specific directory
  # here.
  execute_process(
      COMMAND "${PYTHON_EXECUTABLE}" -c
              "import pybind11; print(pybind11.get_cmake_dir())"
      OUTPUT_VARIABLE _tmp_dir
      OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
  list(APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")
  find_package(pybind11 CONFIG REQUIRED)
  
  find_package(pybind11_json) 
  if(TARGET pybind11_json)
    message(STATUS "Skip Download pybind11_json")
  else()
    message(STATUS "Downlod pybind11_json")
    include(external/pybind11-json.cmake)
  endif()
  
  add_subdirectory(cimod)
elseif(${CMAKE_BUILD_TYPE} MATCHES Debug) 
  
  message(STATUS "Use Googletest")
  
  find_package(GTest)
  if(GTest_FOUND)
    message(STATUS "Found googletest")
  else()
    message(STATUS "Downlod Googletest")
    include(external/googletest.cmake)
  endif()
  
  enable_testing()
  
  add_subdirectory(tests)
endif()

if(ENABLE_COVERAGE)
  coverage_evaluate()
endif()
