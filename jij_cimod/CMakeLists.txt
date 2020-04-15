add_library(cimod MODULE main.cpp)

set_target_properties(cimod PROPERTIES
    CXX_VISIBILITY_PRESET "hidden"
    INTERPROCEDURAL_OPTIMIZATION TRUE
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
)


target_link_libraries(cimod
    pybind11::module
    eigen_lib
    cimod_header_only
    pybind11_json
    nlohmann_json::nlohmann_json
)

add_custom_target(python
    DEPENDS cimod
)
