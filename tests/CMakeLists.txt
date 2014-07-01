 enable_testing()
 find_package(Check REQUIRED)
 include_directories(${CHECK_INCLUDE_DIRS})

 ## add the lib we want to test in here
# set(LIBS ${LIBS} ${CHECK_LIBRARIES} isti_clib)
 set(LIBS ${LIBS} ${CHECK_LIBRARIES})
 include_directories(. ../src)
 add_executable(test_flow test_flow.c)
 target_link_libraries(test_flow ${LIBS})
 add_test(test_flow ${CMAKE_CURRENT_BINARY_DIR}/test_flow)