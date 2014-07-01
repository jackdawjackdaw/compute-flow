enable_testing()
find_package(Check REQUIRED)
include_directories(${CHECK_INCLUDE_DIR})

## add the lib we want to test in here
set(libs ${libs} ${CHECK_LIBRARIES} flow)
include_directories(. ../src)
## now build the test executable
add_executable(test-flow test-flow)
target_link_libraries(test-flow ${libs})

add_test(test_flow ${CMAKE_CURRENT_BINARY_DIR}/test-flow)