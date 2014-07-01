cmake_minimum_required(VERSION 2.8)

## the event by event 
add_executable(flow-evt-plane compute-flow)
target_link_libraries(flow-evt-plane ${libs})

## the running means version
add_executable(flow-evt-plane-run compute-v2v3)
target_link_libraries(flow-evt-plane-run ${libs})

add_executable(flow-naiive compute-v2-naiive)
target_link_libraries(flow-naiive ${libs})

install_targets( /bin flow-evt-plane flow-evt-plane-run flow-naiive)



