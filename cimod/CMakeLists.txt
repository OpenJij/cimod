add_library(cxxcimod MODULE main.cpp)

set_target_properties(cxxcimod PROPERTIES
    CXX_VISIBILITY_PRESET "hidden"
    INTERPROCEDURAL_OPTIMIZATION TRUE
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
)


target_link_libraries(cxxcimod
    pybind11::module
    eigen_lib
    cxxcimod_header_only
    pybind11_json
    nlohmann_json::nlohmann_json
)

add_custom_target(python
    DEPENDS cxxcimod
)
