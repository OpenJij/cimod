cmake_minimum_required(VERSION 3.11)
project(openjij)


set(CMAKE_CXX_STANDARD 14)
# Set default compile flags for GCC
if(CMAKE_COMPILER_IS_GNUCXX)
    message(STATUS "GCC detected, adding compile flags")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-long-long -pedantic")
else(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif(CMAKE_COMPILER_IS_GNUCXX)


# add fPIC option to all object files
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
 

##### Set default behavior #####
set(DEFAULT_USE_OMP No)
set(DEFAULT_USE_TEST No)
option(USE_TEST "Use test code" No)


# Use OpenMP as default behavior
if(NOT DEFINED USE_OMP)
    set(USE_OMP ${DEFAULT_USE_OMP})
endif()

if(NOT DEFINED USE_TEST)
    set(USE_TEST ${DEFAULT_USE_TEST})
endif()

message(STATUS "USE_OMP = ${USE_OMP}")
message(STATUS "USE_TEST = ${USE_TEST}")

#TODO: OMP settings

list(APPEND CMAKE_MODULE_PATH external)
if(USE_TEST)
    include(${CMAKE_SOURCE_DIR}/external/googletest.cmake)
else()
	message(STATUS "Skip downloding googletest")
endif()
include(${CMAKE_SOURCE_DIR}/external/pybind11.cmake)
include(${CMAKE_SOURCE_DIR}/external/eigen.cmake)
include(${CMAKE_SOURCE_DIR}/external/json.cmake)
include(${CMAKE_SOURCE_DIR}/external/pybind11-json.cmake)

add_subdirectory(src)
add_subdirectory(cimod)
if(USE_TEST)
    add_subdirectory(tests)
endif()
