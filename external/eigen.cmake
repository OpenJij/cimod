include(FetchContent)

#### Eigen ####
FetchContent_Declare(
    Eigen
    GIT_REPOSITORY  https://gitlab.com/libeigen/eigen
    GIT_TAG          3.4.0
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
)
 
set(EIGEN_MPL2_ONLY ON)
set(BUILD_TESTING OFF)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_LEAVE_TEST_IN_ALL_TARGET OFF) 
set(EIGEN_BUILD_PKGCONFIG OFF) 
set(EIGEN_BUILD_BTL OFF)
set(EIGEN_BUILD_DOC OFF)
set(EIGEN_TEST_NOQT OFF)
set(EIGEN_SPLIT_LARGE_TESTS OFF)
set(EIGEN_NO_ASSERTION_CHECKING OFF)
set( OFF)

#set_target_properties(eigen PROPERTIES
#     EIGEN_MPL2_ONLY ON
#     BUILD_TESTING OFF
#     EIGEN_BUILD_TESTING OFF
#     EIGEN_BUILD_TEST OFF
#) 

#add_library(cimod-eigen_lib INTERFACE)
#set_target_properties(cimod-eigen_lib PROPERTIES
#    EIGEN_MPL2_ONLY ON
#    BUILD_TESTING OFF
#    EIGEN_LEAVE_TEST_IN_ALL_TARGET OFF
#    EIGEN_SPLIT_LARGE_TESTS OFF
#    EIGEN_BUILD_PKGCONFIG OFF
#    EIGEN_BUILD_DOC OFF
#    EIGEN_TEST_NOQT OFF
#    EIGEN_BUILD_BTL OFF
#    EIGEN_NO_ASSERTION_CHECKING OFF
#)
#target_include_directories(cimod-eigen_lib INTERFACE ${eigen_SOURCE_DIR})
#target_compile_definitions(cimod-eigen_lib INTERFACE EIGEN_MPL2_ONLY)
if (APPLE)
    if(BLAS_FOUND AND LAPACK_FOUND) 
      set(EIGEN_USE_BLAS ON)
      set(EIGEN_USE_LAPACKE O)
    endif()
endif()

FetchContent_MakeAvailable(Eigen)
FetchContent_GetProperties(Eigen)

find_package(eigen3)
find_package(Eigen3)

message(STATUS "eigen_SOURCE_DIR = ${eigen_SOURCE_DIR}")

