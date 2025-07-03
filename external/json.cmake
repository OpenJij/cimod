# Copyright 2025 Jij Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

include(FetchContent)

message(CHECK_START "Fetching JSON")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

set(BUILD_TESTING OFF)

#### nlohmann_json ####
FetchContent_Declare(
     nlohmann_json
     GIT_REPOSITORY  https://github.com/nlohmann/json
     GIT_TAG         v3.12.0
     GIT_SHALLOW     TRUE
     )
     
FetchContent_MakeAvailable(nlohmann_json)

# Since the git repository of nlohmann/json is huge, we store only a single-include file json.hpp in our project.
#set(BUILD_TESTING OFF)
#add_library(nlohmann_json INTERFACE)
#add_library(nlohmann_json::nlohmann_json ALIAS nlohmann_json)
#target_include_directories(nlohmann_json INTERFACE ${CMAKE_SOURCE_DIR}/external/nlohmann_json)

list(POP_BACK CMAKE_MESSAGE_INDENT)
message(CHECK_PASS "fetched")
