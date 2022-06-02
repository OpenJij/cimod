add_library(cxxcimod MODULE main.cpp)

set(python_module_name cxxcimod)

if(NOT MSVC AND NOT ${CMAKE_BUILD_TYPE} MATCHES Debug|RelWithDebInfo)
    # Strip unnecessary sections of the binary on Linux/macOS
    pybind11_strip(${python_module_name})
endif()

set_target_properties(${python_module_name} PROPERTIES
    INTERPROCEDURAL_OPTIMIZATION ON
    CXX_VISIBILITY_PRESET ON
    VISIBILITY_INLINES_HIDDEN ON
    POSITION_INDEPENDENT_CODE ON
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
    VERSION_INFO "${CIMOD_VERSION_INFO}"
)


target_link_libraries(${python_module_name}
    pybind11::module
    pybind11::thin_lto
    pybind11_json
    nlohmann_json::nlohmann_json
    $<IF:$<TARGET_EXISTS:Eigen3::Eigen>,Eigen3::Eigen,cimod-eigen_lib>
    $<$<TARGET_EXISTS:OpenMP::OpenMP_CXX>:OpenMP::OpenMP_CXX>
    cxxcimod_header_only
)

add_custom_target(python
    DEPENDS ${python_module_name}
)
