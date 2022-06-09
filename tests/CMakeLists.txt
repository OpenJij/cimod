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

#include(GoogleTest)

add_executable(cimod_test
    test.cpp
)

if(ENABLE_COVERAGE AND (NOT MSVC))
  add_coverage(cimod_test)
endif()


target_include_directories(cimod_test PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_definitions(cimod_test PRIVATE
    BUILD_TESTING=OFF
)

target_link_libraries(cimod_test PRIVATE
    cxxcimod_header_only
    nlohmann_json::nlohmann_json
    $<IF:$<TARGET_EXISTS:GTest::gtest>,GTest::gtest,gtest>
    $<IF:$<TARGET_EXISTS:GTest::gtest_main>,GTest::gtest_main,gtest_main>
    $<IF:$<TARGET_EXISTS:GTest::gmock>,GTest::gmock,gmock>
    $<IF:$<TARGET_EXISTS:GTest::gmock_main>,GTest::gmock_main,gmock_main>
)

gtest_discover_tests(cimod_test)
