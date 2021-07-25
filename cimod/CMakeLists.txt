add_library(cxxcimod MODULE main.cpp)


set_target_properties(cxxcimod PROPERTIES
    CXX_VISIBILITY_PRESET "hidden"
    VISIBLITY_INLINES_HIDDEN ON
    POSITION_INDEPENDENT_CODE ON
    INTERPROCEDURAL_OPTIMIZATION TRUE
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
    VERSION_INFO "${CIMOD_VERSION_INFO}"
)


target_link_libraries(cxxcimod
    pybind11::module
    pybind11::thin_lto
    pybind11_json
    nlohmann_json::nlohmann_json
    $<IF:$<TARGET_EXISTS:Eigen3::Eigen>,Eigen3::Eigen,cimod-eigen_lib>
    $<$<TARGET_EXISTS:OpenMP::OpenMP_CXX>:OpenMP::OpenMP_CXX>
    cxxcimod_header_only
)

add_custom_target(python
    DEPENDS cxxcimod
)
