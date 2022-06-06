include(FetchContent)

set(FETCHCONTENT_QUIET OFF)
set(EIGEN_MPL2_ONLY ON)
set(BUILD_TESTING OFF)
#set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
#set(EIGEN_LEAVE_TEST_IN_ALL_TARGET OFF) 
#set(EIGEN_BUILD_PKGCONFIG OFF) 
#set(EIGEN_BUILD_BTL OFF)
#set(EIGEN_BUILD_DOC OFF)
#set(EIGEN_TEST_NOQT OFF)
#set(EIGEN_SPLIT_LARGE_TESTS OFF)
#set(EIGEN_NO_ASSERTION_CHECKING OFF)
#set( OFF)

if (APPLE)
    if(BLAS_FOUND AND LAPACK_FOUND) 
      set(EIGEN_USE_BLAS ON)
      set(EIGEN_USE_LAPACKE O)
    endif()
endif()

#### Eigen ####
FetchContent_Declare(
    Eigen
    GIT_REPOSITORY  https://gitlab.com/libeigen/eigen
    GIT_TAG          3.4.0
    GIT_SHALLOW TRUE
)

FetchContent_MakeAvailable(Eigen)
FetchContent_GetProperties(Eigen)

find_package(Eigen3)

if(Eigen3_FOUND)
    message(STATUS "Found Eigen3")
endif()

message(STATUS "eigen_SOURCE_DIR = ${eigen_SOURCE_DIR}")

