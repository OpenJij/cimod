if (${CIMOD_BUILD_SHARED_LIBS})
  add_library(cxxjij SHARED main.cpp)
else()
  add_library(cxxjij MODULE main.cpp)
endif()


set_target_properties(cxxcimod PROPERTIES
    CXX_VISIBILITY_PRESET "hidden"
    VISIBLITY_INLINES_HIDDEN ON
    POSITION_INDEPENDENT_CODE ON
    INTERPROCEDURAL_OPTIMIZATION TRUE
    PREFIX "${PYTHON_MODULE_PREFIX}"
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
)


target_link_libraries(cxxcimod
    pybind11::module
    pybind11::lto
    pybind11_json
    nlohmann_json::nlohmann_json
    $<IF:$<TARGET_EXISTS:Eigen3::Eigen>,Eigen3::Eigen,eigen_lib>
    $<$<TARGET_EXISTS:OpenMP::OpenMP_CXX>:OpenMP::OpenMP_CXX>
    cxxcimod_header_only
)

add_custom_target(python
    DEPENDS cxxcimod
)
