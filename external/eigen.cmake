include(FetchContent)

message(CHECK_START "Fetching Eigen3")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

set(FETCHCONTENT_QUIET OFF)

set(BUILD_TESTING OFF) 
set(TEST_LIB OFF)

set(EIGEN_MPL2_ONLY ON CACHE BOOL "" FORCE)
set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_Fortran_COMPILER_WORKS OFF CACHE BOOL "" FORCE)

if(BLAS_FOUND) 
  set(EIGEN_USE_BLAS ON) 
endif()

if(LAPACK_FOUND) 
  set(EIGEN_USE_LAPACKE ON)
endif()

#### Eigen ####
FetchContent_Declare(
    eigen
    GIT_REPOSITORY  https://gitlab.com/libeigen/eigen
    GIT_TAG         3.4.0
    GIT_SHALLOW     TRUE
    CMAKE_ARGS 
    -DEIGEN_MPL2_ONLY=ON 
    -DEIGEN_BUILD_PKGCONFIG=OFF
    -DEIGEN_BUILD_DOC=OFF
    -DEIGEN_BUILD_TESTING=OFF
    -DEIGEN_Fortran_COMPILER_WORKS=OFF
    )

FetchContent_MakeAvailable(eigen)

add_library(cimod-eigen_lib INTERFACE)
target_include_directories(cimod-eigen_lib INTERFACE ${eigen_SOURCE_DIR})
target_compile_definitions(cimod-eigen_lib INTERFACE EIGEN_MPL2_ONLY)

list(POP_BACK CMAKE_MESSAGE_INDENT)
message(CHECK_PASS "fetched")
