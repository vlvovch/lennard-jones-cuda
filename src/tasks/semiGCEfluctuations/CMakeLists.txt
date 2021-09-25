include_directories ("${PROJECT_SOURCE_DIR}/src/library")
include_directories ("${PROJECT_SOURCE_DIR}/src/extra/sample-moments")

add_executable(semiGCEfluctuations semiGCEfluctuations.cpp)

set_target_properties(semiGCEfluctuations PROPERTIES CUDA_SEPARABLE_COMPILATION OFF)

target_link_libraries(semiGCEfluctuations PUBLIC mdsystem)
