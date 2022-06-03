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

include(GoogleTest)

add_executable(cimod_test
    test.cpp
)

target_include_directories(cimod_test PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(cimod_test PRIVATE
    cxxcimod_header_only
    $<IF:$<TARGET_EXISTS:Eigen3::Eigen>,Eigen3::Eigen,cimod-eigen_lib>
    $<$<TARGET_EXISTS:OpenMP::OpenMP_CXX>:OpenMP::OpenMP_CXX>
    nlohmann_json::nlohmann_json
    $<IF:$<TARGET_EXISTS:GTest::gtest>,GTest::gtest,cimod-gtest>
    $<IF:$<TARGET_EXISTS:GTest::gtest_main>,GTest::gtest_main,cimod-gtest_main>
    $<IF:$<TARGET_EXISTS:GTest::gmock_main>,GTest::gmock_main,cimod-gmock_main>
    $<IF:$<TARGET_EXISTS:GTest::gmock>,GTest::gmock,cimod-gmock>
)

gtest_discover_tests(cimod_test)
