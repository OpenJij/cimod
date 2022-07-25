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

add_library(cxxcimod_header_only INTERFACE)

target_include_directories(cxxcimod_header_only INTERFACE 
  $<BUILD_INTERFACE:${CIMOD_INCLUDE_DIR}>
  $<INSTALL_INTERFACE:include>
)

target_link_libraries(cxxcimod_header_only INTERFACE 
    nlohmann_json::nlohmann_json
    $<IF:$<TARGET_EXISTS:Eigen3::Eigen>,Eigen3::Eigen,cimod-eigen_lib>
    $<$<TARGET_EXISTS:OpenMP::OpenMP_CXX>:OpenMP::OpenMP_CXX>
    $<$<TARGET_EXISTS:BLAS::BLAS>:BLAS::BLAS>
    $<$<TARGET_EXISTS:LAPACK::LAPACK>:LAPACK::LAPACK>
)

target_compile_definitions(cxxcimod_header_only INTERFACE 
    EIGEN_MPL2_ONLY
    BUILD_TESTING=OFF
    TEST_LIB=OFF
    EIGEN_BUILD_PKGCONFIG=OFF
    EIGEN_BUILD_DOC=OFF
    EIGEN_DOC_USE_MATHJAX=OFF 
    EIGEN_BUILD_TESTING=OFF 
    EIGEN_TEST_NOQT=OFF 
    EIGEN_LEAVE_TEST_IN_ALL_TARGET=OFF 
    $<$<TARGET_EXISTS:BLAS::BLAS>:EIGEN_USE_BLAS>
    $<$<TARGET_EXISTS:LAPACK::LAPACK>:EIGEN_USE_LAPACKE>
    $<$<CXX_COMPILER_ID:MSVC>:EIGEN_Fortran_COMPILER_WORKS=OFF>   
)
