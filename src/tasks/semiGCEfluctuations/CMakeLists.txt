include_directories ("${PROJECT_SOURCE_DIR}/src/library")

add_executable(semiGCEfluctuations semiGCEfluctuations.cpp NumberStatistics.h NumberStatistics.cpp)

set_target_properties(semiGCEfluctuations PROPERTIES CUDA_SEPARABLE_COMPILATION OFF)

target_link_libraries(semiGCEfluctuations PUBLIC mdsystem)
