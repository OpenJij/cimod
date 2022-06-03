include(FetchContent)

#### Eigen ####
FetchContent_Declare(
    eigen
    GIT_REPOSITORY  https://gitlab.com/libeigen/eigen
    GIT_TAG          3.4.0
    CMAKE_ARGS -DEIGEN_MPL2_ONLY=ON -DBUILD_TESTING=OFF
    )
    
set(EIGEN_MPL2_ONLY ON)

set_target_properties(eigen PROPERTIES
     EIGEN_MPL2_ONLY ON
     BUILD_TESTING OFF
     EIGEN_BUILD_TESTING OFF
     EIGEN_BUILD_TEST OFF
) 

#add_library(cimod-eigen_lib INTERFACE)
#set_target_properties(cimod-eigen_lib PROPERTIES
#    EIGEN_MPL2_ONLY ON
#)
#target_include_directories(cimod-eigen_lib INTERFACE ${eigen_SOURCE_DIR})
#target_compile_definitions(cimod-eigen_lib INTERFACE EIGEN_MPL2_ONLY)
if (APPLE)
    if(BLAS_FOUND AND LAPACK_FOUND) 
      target_compile_definitions(eigen INTERFACE EIGEN_USE_BLAS=ON)
      target_compile_definitions(eigen INTERFACE EIGEN_USE_LAPACKE=ON)
    endif()
endif()

if(OpenMP_FOUND)
  target_link_libraries(eigen INTERFACE OpenMP::OpenMP_CXX)
endif()

FetchContent_MakeAvailable(eigen)
