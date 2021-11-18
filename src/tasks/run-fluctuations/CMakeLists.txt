include_directories ("${PROJECT_SOURCE_DIR}/src/library")
include_directories ("${PROJECT_SOURCE_DIR}/src/extra/sample-moments")

add_executable(run-fluctuations run-fluctuations.cpp include/run-fluctuations-aux.h ../auxiliary/time-average-aux.h)

set_target_properties(run-fluctuations PROPERTIES CUDA_SEPARABLE_COMPILATION OFF)

target_link_libraries(run-fluctuations PUBLIC mdsystem)
