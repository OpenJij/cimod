include(FetchContent)

message(CHECK_START "Fetching Eigen3")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

set(FETCHCONTENT_QUIET OFF)

#### Eigen ####
FetchContent_Declare(
    Eigen
    GIT_REPOSITORY  https://gitlab.com/libeigen/eigen
    GIT_TAG          3.4.0
    GIT_SHALLOW TRUE
)

set(BUILD_TESTING OFF)

set(EIGEN_BUILD_PKGCONFIG OFF)
set(EIGEN_BUILD_DOC OFF)
set(EIGEN_BUILD_TESTING OFF)

set(EIGEN_MPL2_ONLY ON)

set( OFF)

if(BLAS_FOUND AND LAPACK_FOUND) 
  set(EIGEN_USE_BLAS ON)
  set(EIGEN_USE_LAPACKE ON)
endif()

FetchContent_MakeAvailable(Eigen)
#FetchContent_GetProperties(Eigen)

#find_package(Eigen3)

#if(Eigen3_FOUND)
#    message(STATUS "Found Eigen3")
#endif()

message(STATUS "eigen_SOURCE_DIR = ${eigen_SOURCE_DIR}")

list(POP_BACK CMAKE_MESSAGE_INDENT)
message(CHECK_PASS "fetched")
