 enable_testing()
 find_package(Check REQUIRED)
 include_directories(${CHECK_INCLUDE_DIRS})

 ## add the lib we want to test in here
# set(LIBS ${LIBS} ${CHECK_LIBRARIES} flow-lib)
 set(LIBS ${LIBS} ${CHECK_LIBRARIES})
 include_directories(. ../src)
 ## need to actually 
 add_executable(test-flow test-flow.c)
 target_link_libraries(test-flow ${LIBS})
 add_test(test_flow ${CMAKE_CURRENT_BINARY_DIR}/test-flow)