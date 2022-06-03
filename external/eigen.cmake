include(FetchContent)

#### Eigen ####
FetchContent_Declare(
    eigen
    GIT_REPOSITORY  https://gitlab.com/libeigen/eigen
    GIT_TAG          3.4.0
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
)
 
set(EIGEN_MPL2_ONLY ON)
set(BUILD_TESTING OFF) 
set(EIGEN_LEAVE_TEST_IN_ALL_TARGET OFF) 
set(EIGEN_BUILD_PKGCONFIG OFF) 
set(EIGEN_BUILD_BTL OFF)
set(EIGEN_TEST_NOQT OFF)

FetchContent_MakeAvailable(eigen)

#set_target_properties(eigen PROPERTIES
#     EIGEN_MPL2_ONLY ON
#     BUILD_TESTING OFF
#     EIGEN_BUILD_TESTING OFF
#     EIGEN_BUILD_TEST OFF
#) 

add_library(cimod-eigen_lib INTERFACE)
set_target_properties(cimod-eigen_lib PROPERTIES
    EIGEN_MPL2_ONLY ON
    BUILD_TESTING OFF
    EIGEN_LEAVE_TEST_IN_ALL_TARGET OFF
    EIGEN_BUILD_PKGCONFIG OFF
    EIGEN_TEST_NOQT OFF
    EIGEN_BUILD_BTL OFF
)
target_include_directories(cimod-eigen_lib INTERFACE ${eigen_SOURCE_DIR})
#target_compile_definitions(cimod-eigen_lib INTERFACE EIGEN_MPL2_ONLY)
if (APPLE)
    if(BLAS_FOUND AND LAPACK_FOUND) 
      set_target_properties(cimod-eigen_lib PROPERTIES
        EIGEN_USE_BLAS ON
        EIGEN_USE_LAPACKE O
      )
    endif()
endif()



