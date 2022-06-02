set(python_module_name cxxcimod)
#WITH_SOABI
add_library(${python_module_name} MODULE main.cpp)

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
    $<IF:$<TARGET_EXISTS:Eigen3::Eigen>,Eigen3::Eigen,cimod-eigen_lib>
    $<$<TARGET_EXISTS:OpenMP::OpenMP_CXX>:OpenMP::OpenMP_CXX>
    cxxcimod_header_only
)

add_custom_target(python
    DEPENDS ${python_module_name}
)
