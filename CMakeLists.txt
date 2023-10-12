cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(LennardJonesCUDA LANGUAGES CXX)


include(CheckLanguage)
check_language(CUDA)


set (CMAKE_CXX_STANDARD 11)
#
#if(NOT CUDA_SAMPLES_DIR)
#   set(CUDA_SAMPLES_DIR C:/ProgramData/NVIDIA\ Corporation/CUDA\ Samples/v10.2 CACHE STRING
#      "Set the path of CUDA SDK/samples folder."
#      FORCE)
#endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebug RelWithDebInfo MinSizeRel."
      FORCE)
   # Set the possible values of build type for cmake-gui
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
             "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()


if(CMAKE_CUDA_COMPILER)
   enable_language(CUDA)
   include_directories("${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")

   add_definitions(-DUSE_CUDA_TOOLKIT)

   #if(NOT CUDA_TARGET_ARCH)
   #   set(CUDA_TARGET_ARCH "-arch=sm_30" CACHE STRING
   #      "CUDA compute capability target. Default is -arch=sm_30. Other option: -gencode arch=compute_XX,code=sm_XX"
   #      FORCE)
   #endif(NOT CUDA_TARGET_ARCH)

   #set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS} ${CUDA_TARGET_ARCH})

   if(NOT CUDA_SAMPLES_DIR)
   file(DOWNLOAD https://raw.githubusercontent.com/NVIDIA/cuda-samples/master/Common/helper_cuda.h
           ${PROJECT_BINARY_DIR}/thirdparty/CUDA_SAMPLES/helper_cuda.h)

   file(DOWNLOAD https://raw.githubusercontent.com/NVIDIA/cuda-samples/master/Common/helper_string.h
           ${PROJECT_BINARY_DIR}/thirdparty/CUDA_SAMPLES/helper_string.h)
   endif()

endif(CMAKE_CUDA_COMPILER)

add_subdirectory(src/library)
add_subdirectory(src/gui)
add_subdirectory(src/tasks)
