include(FetchContent)

#### Eigen ####
FetchContent_Declare(
    eigen
    GIT_REPOSITORY  https://gitlab.com/libeigen/eigen
    GIT_TAG         3.3.9
    CMAKE_ARGS -DEIGEN_MPL2_ONLY
    )
    
set(EIGEN_BUILD_DOC OFF)
set(BUILD_TESTING OFF)
set(EIGEN_BUILD_PKGCONFIG OFF)
set(EIGEN_CPP_STANDARD 11)
set(EIGEN_MPL2_ONLY ON)
FetchContent_MakeAvailable(eigen)

add_library(cimod-eigen_lib INTERFACE)
set_target_properties(cimod-eigen_lib PROPERTIES
    EIGEN_BUILD_DOC OFF
    BUILD_TESTING OFF
    EIGEN_BUILD_PKGCONFIG OFF
    EIGEN_CPP_STANDARD 11
    EIGEN_MPL2_ONLY ON
)
target_include_directories(cimod-eigen_lib INTERFACE ${eigen_SOURCE_DIR})
target_compile_definitions(cimod-eigen_lib INTERFACE EIGEN_MPL2_ONLY)
if (APPLE)
    if(BLAS_FOUND AND LAPACK_FOUND) 
      target_compile_definitions(cimod-eigen_lib INTERFACE EIGEN_USE_BLAS=ON)
      target_compile_definitions(cimod-eigen_lib INTERFACE EIGEN_USE_LAPACKE=ON)
    endif()
endif()

if(OpenMP_FOUND)
  target_link_libraries(cimod-eigen_lib INTERFACE OpenMP::OpenMP_CXX)
endif()
