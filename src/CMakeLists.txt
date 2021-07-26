add_library(cxxcimod_header_only INTERFACE)

target_include_directories(cxxcimod_header_only INTERFACE 
  ${CMAKE_CURRENT_SOURCE_DIR}
)

  
target_link_libraries(cxxcimod_header_only INTERFACE 
    $<IF:$<TARGET_EXISTS:Eigen3::Eigen>,Eigen3::Eigen,cimod-eigen_lib>
    $<$<TARGET_EXISTS:OpenMP::OpenMP_CXX>:OpenMP::OpenMP_CXX>
)

