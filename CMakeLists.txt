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

cmake_minimum_required(VERSION 3.20)
project(cimod)

enable_language(C CXX)

set(CMAKE_VERBOSE_MAKEFILE TRUE)

if (APPLE)
    set(CMAKE_MACOSX_RPATH TRUE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-command-line-argument")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wno-unused-command-line-argument")
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
    message(STATUS "APPLE_FORCE_USE_OMP = ${FORCE_USE_OMP}")
  endif()
endif()

if(USE_OMP)
  set(CMAKE_INSTALL_OPENMP_LIBRARIES ON)
  find_package(OpenMP REQUIRED)
  message(STATUS "OMP_FOUND = YES")
endif()

#message(STATUS "Apple Native Build")
#set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp") 
#set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
#find_package(OpenMP REQUIRED)
#message(STATUS "OMP_FOUND = YES")
#add_definitions(-DUSE_OMP)


if (APPLE)
    find_package(BLAS)
    find_package(LAPACK)
endif()

message(STATUS "USE_OMP = ${USE_OMP}")

list(APPEND CMAKE_MODULE_PATH external)

find_package(Eigen3 3.3 NO_MODULE)

if(EIGEN3_FOUND)
  message(STATUS "Found eigen3")
else() 
  message(STATUS "Downlod eigen3")
  include(external/eigen.cmake)
endif() 

include(external/json.cmake)

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
  
  include(external/pybind11-json.cmake)
  
  add_subdirectory(cimod)
elseif(${CMAKE_BUILD_TYPE} MATCHES Debug) 
  
  message(STATUS "Use Googletest")
  message(STATUS "Downlod googletest")
  include(external/googletest.cmake)
  find_package(GTest)
  
  if(GTest_FOUND)
    message(STATUS "Found googletest")
  endif() 
  
  enable_testing()
  
  add_subdirectory(tests)
endif()
