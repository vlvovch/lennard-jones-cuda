# Base part of the library  
set(SRCS_MDSystem
MDSystem.h
MDSystem.cu
MDSystem.cpp
splinefunction.h
splinefunction.cpp
)

add_library(mdsystem STATIC ${SRCS_MDSystem})

target_include_directories(mdsystem PUBLIC
"${PROJECT_SOURCE_DIR}/thirdparty")

if(NOT CUDA_SAMPLES_DIR)
    target_include_directories(mdsystem PUBLIC "${PROJECT_BINARY_DIR}/thirdparty/CUDA_SAMPLES")
else()
    target_include_directories(mdsystem PUBLIC "${CUDA_SAMPLES_DIR}/common/inc")
endif()

# We need to explicitly state that we need all CUDA files in the particle
# library to be built with -dc as the member functions could be called by
# other libraries and executables
set_target_properties(mdsystem
                      PROPERTIES  CUDA_RESOLVE_DEVICE_SYMBOLS ON
                      CUDA_SEPARABLE_COMPILATION OFF
                      )
