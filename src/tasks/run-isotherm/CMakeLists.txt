include_directories ("${PROJECT_SOURCE_DIR}/src/library")
include_directories ("${PROJECT_SOURCE_DIR}/src/extra/sample-moments")

add_executable(run-isotherm run-isotherm.cpp include/run-isotherm-aux.h ../auxiliary/time-average-aux.h)

set_target_properties(run-isotherm PROPERTIES CUDA_SEPARABLE_COMPILATION OFF)

target_link_libraries(run-isotherm PUBLIC mdsystem)
