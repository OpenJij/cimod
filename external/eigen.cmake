include(ExternalProject)

message(CHECK_START "Fetching Eigen3")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

set(BUILD_TESTING OFF)

set(EIGEN_Fortran_COMPILER_WORKS OFF)
set(EIGEN_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen-install/")

#### Eigen ####

ExternalProject_Add(
    eigen
    URL  https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen-src"
    BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen-build"
    INSTALL_DIR "${EIGEN_INSTALL_DIR}"
    CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    -DCMAKE_BUILD_TYPE=Release
    -DEIGEN_MPL2_ONLY=ON
    -DEIGEN_BUILD_PKGCONFIG=OFF 
    -DEIGEN_BUILD_DOC=OFF 
    -DEIGEN_BUILD_TESTING=OFF
)

file(MAKE_DIRECTORY ${EIGEN_INSTALL_DIR}/include)  # avoid race condition

add_library(eigenlib INTERFACE IMPORTED GLOBAL)
add_dependencies(eigenlib eigen)
target_compile_features(eigenlib INTERFACE cxx_std_14)

set_target_properties(eigenlib PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES ${EIGEN_INSTALL_DIR}/include/eigen3
  EIGEN_MPL2_ONLY ON
  EIGEN_BUILD_PKGCONFIG OFF 
  EIGEN_BUILD_DOC OFF 
  EIGEN_BUILD_TESTING OFF
)

if(BLAS_FOUND AND LAPACK_FOUND) 
  set_target_properties(eigenlib PROPERTIES 
    EIGEN_USE_BLAS ON
    EIGEN_USE_LAPACKE ON
  )
endif()

list(POP_BACK CMAKE_MESSAGE_INDENT)
message(CHECK_PASS "fetched")
