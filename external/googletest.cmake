include(FetchContent)

message(CHECK_START "Fetching GoogleTest")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

set(CMAKE_CXX_STANDARD 17)
set(FETCHCONTENT_QUIET OFF)

#### Google test ####
FetchContent_Declare(
    googletest
    GIT_REPOSITORY  https://github.com/google/googletest.git
    GIT_TAG         release-1.11.0
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
FetchContent_GetProperties(googletest)

find_package(GTest)
if(GTest_FOUND)
    message(STATUS "Found googletest")
endif()

message(STATUS "gtest_SOURCE_DIR = ${gtest_SOURCE_DIR}")
message(STATUS "gmock_SOURCE_DIR = ${gmock_SOURCE_DIR}")


list(POP_BACK CMAKE_MESSAGE_INDENT)
message(CHECK_PASS "fetched")
