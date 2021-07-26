add_executable(cimod_test
    EXCLUDE_FROM_ALL
    test.cpp
)

set_target_properties(cimod_test PROPERTIES
    CXX_VISIBILITY_PRESET "hidden"
    INTERPROCEDURAL_OPTIMIZATION TRUE
)

target_link_libraries(cimod_test PRIVATE
    cxxcimod_header_only
    $<IF:$<TARGET_EXISTS:Eigen3::Eigen>,Eigen3::Eigen,cimod-eigen_lib>
    $<$<TARGET_EXISTS:OpenMP::OpenMP_CXX>:OpenMP::OpenMP_CXX>
    nlohmann_json::nlohmann_json
    gmock_main
    gmock
    gtest_main
    gtest
    pthread
)
