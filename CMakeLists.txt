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

CMAKE_MINIMUM_REQUIRED (VERSION 3.20 FATAL_ERROR)
PROJECT (cimod
         DESCRIPTION "C++ library for a binary (and polynomial) quadratic model."
         HOMEPAGE_URL "https://www.openjij.org"
         )

# Standard includes
INCLUDE (CheckLanguage)

SET (CIMOD_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

IF (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    IF (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
        SET (lines "You are building in-place.")
        MESSAGE (AUTHOR_WARNING ${lines})
    ENDIF ()
    SET (CIMOD_MAIN_PROJECT ON)
    MESSAGE (STATUS "CMake ${CMAKE_VERSION}")
    IF (CMAKE_CXX_STANDARD)
        SET (CMAKE_CXX_EXTENSIONS OFF)
        SET (CMAKE_CXX_STANDARD_REQUIRED ON)
    ENDIF ()
ELSE ()
    SET (CIMOD_MAIN_PROJECT OFF)
ENDIF ()

# Options
OPTION (CIMOD_INSTALL "Install cimod header files?" ${CIMOD_MAIN_PROJECT})
OPTION (CIMOD_TEST "Build cimod test suite?" ${CIMOD_MAIN_PROJECT})
OPTION (CIMOD_DOCS "Build cimod docs?" ${CIMOD_MAIN_PROJECT})
OPTION (BUILD_DOCS "Enable Doxygen support." OFF)
OPTION (CMAKE_REQUIRE_FAILE "If CMake could not find dependencies, build will faile." OFF)

LIST (APPEND CMAKE_MODULE_PATH external)
LIST (APPEND CMAKE_MODULE_PATH cmake)

CHECK_LANGUAGE (C)
CHECK_LANGUAGE (CXX)

IF (CMAKE_C_COMPILER)
    ENABLE_LANGUAGE (C)
ENDIF ()
IF (CMAKE_CXX_COMPILER)
    ENABLE_LANGUAGE (CXX)
ENDIF ()

IF (UNIX AND NOT APPLE)
    CHECK_LANGUAGE (Fortran)
    IF (CMAKE_Fortran_COMPILER)
        ENABLE_LANGUAGE (Fortran)
    ENDIF ()
ENDIF ()

IF (CMAKE_REQUIRE_FAILE)
    SET (CMAKE_REQUIRE_FIND_PACKAGE_Eigen3 ON)
    SET (CMAKE_REQUIRE_FIND_PACKAGE_nlohmann_json ON)
    SET (CMAKE_REQUIRE_FIND_PACKAGE_pybind11 ON)
    SET (CMAKE_REQUIRE_FIND_PACKAGE_pybind11_json ON)
    SET (CMAKE_REQUIRE_FIND_PACKAGE_GTest ON)
ENDIF ()

# detect operating system and host processor
MESSAGE (STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
MESSAGE (STATUS "The host processor is ${CMAKE_HOST_SYSTEM_PROCESSOR}")
MESSAGE (STATUS "CMAKE_GENERATOR = ${CMAKE_GENERATOR}")
MESSAGE (STATUS "CIBUILDWHEEL = $ENV{CIBUILDWHEEL}")
MESSAGE (STATUS "CMAKE_REQUIRE_FIND_PACKAGE_Eigen3 = ${CMAKE_REQUIRE_FIND_PACKAGE_Eigen3}")
MESSAGE (STATUS "CMAKE_REQUIRE_FIND_PACKAGE_nlohmann_json = ${CMAKE_REQUIRE_FIND_PACKAGE_nlohmann_json}")
MESSAGE (STATUS "CMAKE_REQUIRE_FIND_PACKAGE_pybind11 = ${CMAKE_REQUIRE_FIND_PACKAGE_pybind11}")
MESSAGE (STATUS "CMAKE_REQUIRE_FIND_PACKAGE_pybind11_json = ${CMAKE_REQUIRE_FIND_PACKAGE_pybind11_json}")
MESSAGE (STATUS "CMAKE_REQUIRE_FIND_PACKAGE_GTest = ${CMAKE_REQUIRE_FIND_PACKAGE_GTest}")
MESSAGE (STATUS "Started CMake for ${PROJECT_NAME} \n")

IF (APPLE)
    SET (CMAKE_MACOSX_RPATH TRUE)
    MESSAGE (STATUS "CMAKE_OSX_DEPLOYMENT_TARGET = ${CMAKE_OSX_DEPLOYMENT_TARGET}")
    #include(cmake/PythonAutoDetectOSX.cmake)
    MESSAGE (STATUS "CMAKE_OSX_DEPLOYMENT_TARGET = ${CMAKE_OSX_DEPLOYMENT_TARGET}")
    ADD_COMPILE_OPTIONS (-Wno-unused-command-line-argument)
ELSEIF (MSVC)
    ADD_COMPILE_OPTIONS (/bigobj /EHsc)
    SET (CMAKE_Fortran_COMPILER OFF CACHE BOOL "" FORCE)
ENDIF ()

SET (CMAKE_CXX_STANDARD 17)
SET (CXX_STANDARD_REQUIRED TRUE)
SET (CMAKE_VERBOSE_MAKEFILE TRUE)

IF (NOT APPLE)
    IF ($ENV{CIBUILDWHEEL})
        MESSAGE (STATUS "ENABLE BLAS STATIC")
        SET (BLA_STATIC ON CACHE BOOL "" FORCE)
    ENDIF ()
ENDIF ()

FIND_PACKAGE (BLAS)
FIND_PACKAGE (LAPACK)

##### Set default behavior #####
SET (DEFAULT_USE_OMP ON)

IF (NOT APPLE)
    IF (NOT DEFINED USE_OMP)
        SET (USE_OMP ${DEFAULT_USE_OMP})
    ENDIF ()
ELSEIF ($ENV{CIBUILDWHEEL})
    MESSAGE (STATUS "OpenMP is disabled on macos for distribution binaries.")
ELSE ()
    SET (USE_OMP ${DEFAULT_USE_OMP})
ENDIF ()

IF (USE_OMP)
    FIND_PACKAGE (OpenMP)
    IF (OpenMP_FOUND)
        MESSAGE (STATUS "OMP_FOUND = YES")
        MESSAGE (STATUS "USE_OMP = ON")
        IF (WIN32)
            SET (CMAKE_INSTALL_OPENMP_LIBRARIES ON)
        ENDIF ()
        ADD_DEFINITIONS (-DUSE_OMP)
    ELSE ()
        MESSAGE (STATUS "OMP_FOUND = NO")
        MESSAGE (STATUS "USE_OMP = OFF")
        SET (USE_OMP OFF)
    ENDIF ()
ENDIF ()

FIND_PACKAGE (Eigen3 CONFIG)

IF (EIGEN3_FOUND)
    MESSAGE (STATUS "Skip Download eigen3")
ELSE ()
    MESSAGE (STATUS "Downlod eigen3")
    INCLUDE (external/eigen.cmake)
ENDIF ()

FIND_PACKAGE (nlohmann_json CONFIG)
IF (TARGET nlohmann_json::nlohmann_json)
    MESSAGE (STATUS "Skip Download nlohmann_json")
ELSE ()
    MESSAGE (STATUS "Downlod nlohmann_json")
    SET (JSON_BuildTests OFF CACHE INTERNAL "")
    INCLUDE (external/json.cmake)
ENDIF ()


ADD_SUBDIRECTORY (include)

# Build Python Extension 
IF (SKBUILD AND CIMOD_MAIN_PROJECT)
    MESSAGE (STATUS "Build Python Extension.")
    OPTION (WITH_THREAD "Compile in rudimentary thread support" ON)
    OPTION (BUILD_TESTING "Enable CTest support." OFF)
    FIND_PACKAGE (Python COMPONENTS Interpreter Development)
    # Scikit-Build does not add your site-packages to the search path
    # automatically, so we need to add it _or_ the pybind11 specific directory
    # here.
    EXECUTE_PROCESS (
            COMMAND "${PYTHON_EXECUTABLE}" -c
            "import pybind11; print(pybind11.get_cmake_dir())"
            OUTPUT_VARIABLE _tmp_dir
            OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
    LIST (APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")
    FIND_PACKAGE (pybind11 CONFIG REQUIRED)
    
    FIND_PACKAGE (pybind11_json CONFIG)
    IF (TARGET pybind11_json)
        MESSAGE (STATUS "Skip Download pybind11_json")
    ELSE ()
        MESSAGE (STATUS "Downlod pybind11_json")
        INCLUDE (external/pybind11-json.cmake)
    ENDIF ()
    
    ADD_SUBDIRECTORY (cimod)

ELSEIF (CIMOD_MAIN_PROJECT AND ((${CMAKE_BUILD_TYPE} MATCHES Debug) OR CIMOD_TEST))
    INCLUDE (CTest)
    IF (ENABLE_COVERAGE AND (NOT MSVC))
        MESSAGE (STATUS "Use Codecov")
        FIND_PACKAGE (codecov)
    ENDIF ()
    
    MESSAGE (STATUS "Use Googletest")
    MESSAGE (STATUS "Downlod Googletest")
    INCLUDE (external/googletest.cmake)
    
    ENABLE_TESTING ()
    
    ADD_SUBDIRECTORY (tests)
    
    IF (ENABLE_COVERAGE AND (NOT MSVC))
        COVERAGE_EVALUATE ()
    ENDIF ()
ENDIF ()

IF (CIMOD_MAIN_PROJECT AND CIMOD_DOCS AND BUILD_DOCS)
    FIND_PACKAGE (pybind11 CONFIG)
    IF (TARGET pybind11)
        MESSAGE (STATUS "Skip Download pybind11")
    ELSE ()
        MESSAGE (STATUS "Downlod pybind11_json")
        INCLUDE (external/pybind11.cmake)
    ENDIF ()
    
    FIND_PACKAGE (pybind11_json CONFIG)
    IF (TARGET pybind11_json)
        MESSAGE (STATUS "Skip Download pybind11_json")
    ELSE ()
        MESSAGE (STATUS "Downlod pybind11_json")
        INCLUDE (external/pybind11-json.cmake)
    ENDIF ()
    
    ADD_SUBDIRECTORY (cimod)
    INCLUDE (cmake/GenerateDocs.cmake)
ENDIF ()
