include_directories ("${PROJECT_SOURCE_DIR}/src/library")
include_directories ("${PROJECT_SOURCE_DIR}/src/extra/sample-moments")

add_executable(run-box-events run-box-events.cpp include/run-box-events-aux.h)
set_target_properties(run-box-events PROPERTIES CUDA_SEPARABLE_COMPILATION OFF)
target_link_libraries(run-box-events PUBLIC mdsystem)

add_executable(process-box-events-coord process-box-events-coord.cpp)
set_target_properties(process-box-events-coord PROPERTIES CUDA_SEPARABLE_COMPILATION OFF)
target_link_libraries(process-box-events-coord PUBLIC mdsystem)

add_executable(process-box-events-coord-timedep process-box-events-coord-timedep.cpp)
set_target_properties(process-box-events-coord-timedep PROPERTIES CUDA_SEPARABLE_COMPILATION OFF)
target_link_libraries(process-box-events-coord-timedep PUBLIC mdsystem)
