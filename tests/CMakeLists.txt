cmake_minimum_required(VERSION 3.11)

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
    eigen_lib
    nlohmann_json::nlohmann_json
    gmock_main
    gmock
    gtest_main
    gtest
    pthread
)
