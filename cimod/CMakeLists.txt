add_library(cxxcimod MODULE main.cpp)

set(python_module_name cxxcimod)

set_target_properties(${python_module_name} PROPERTIES
    CXX_VISIBILITY_PRESET "hidden"
    VISIBLITY_INLINES_HIDDEN ON
    POSITION_INDEPENDENT_CODE ON
    INTERPROCEDURAL_OPTIMIZATION TRUE
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
