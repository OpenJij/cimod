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

set(python_module_name cxxcimod)

add_library(${python_module_name} MODULE main.cpp)

target_include_directories(${python_module_name} PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
)

pybind11_extension(${python_module_name})

if(NOT MSVC AND NOT ${CMAKE_BUILD_TYPE} MATCHES Debug|RelWithDebInfo)
    # Strip unnecessary sections of the binary on Linux/macOS
    pybind11_strip(${python_module_name})
endif()

#pybind11_add_module(${python_module_name} MODULE THIN_LTO main.cpp)

set_target_properties(${python_module_name} PROPERTIES
    INTERPROCEDURAL_OPTIMIZATION ON
    CXX_VISIBILITY_PRESET "hidden"
    VISIBILITY_INLINES_HIDDEN ON
    POSITION_INDEPENDENT_CODE ON
)

target_link_libraries(${python_module_name} PRIVATE
    pybind11::module
    pybind11::thin_lto
    pybind11::windows_extras
    pybind11_json
    nlohmann_json::nlohmann_json
    cxxcimod_header_only
)

target_compile_definitions(${python_module_name} PRIVATE 
    $<$<CONFIG:Debug>:PYBIND11_DETAILED_ERROR_MESSAGES>  
    PYTHON_VERSION_INFO=${PYTHON_VERSION_STRING}
)

install(TARGETS ${python_module_name} LIBRARY DESTINATION .)
