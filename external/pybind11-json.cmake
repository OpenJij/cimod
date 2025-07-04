# Copyright 2020-2025 Jij Inc.

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

message(CHECK_START "Fetching Pybind11 Json")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

set(FETCHCONTENT_QUIET OFF)
#### pybind11_json ####
FetchContent_Declare(
    pybind11_json
    GIT_REPOSITORY  https://github.com/pybind/pybind11_json
    GIT_TAG         0.2.15
    GIT_SHALLOW     TRUE
)
set(BUILD_TESTING OFF)
FetchContent_MakeAvailable(pybind11_json)

list(POP_BACK CMAKE_MESSAGE_INDENT)
message(CHECK_PASS "fetched")
