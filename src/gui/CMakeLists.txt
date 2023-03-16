find_package(Qt5Widgets)
find_package(Qt5PrintSupport)
find_package(Qt5OpenGL)
find_package(Qt5OpenGLExtensions)
find_package(OpenGL)

if (Qt5Widgets_FOUND AND Qt5PrintSupport_FOUND AND Qt5OpenGL_FOUND AND OpenGL_FOUND)

#    find_package(GLEW REQUIRED)
#    if (GLEW_FOUND)
#        include_directories(${GLEW_INCLUDE_DIRS})
#        link_libraries(${GLEW_LIBRARIES})
#    endif()

# Find includes in corresponding build directories
#set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS} ${Qt5PrintSupport_EXECUTABLE_COMPILE_FLAGS} ${Qt5OpenGL_EXECUTABLE_COMPILE_FLAGS}" )
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS} ${Qt5PrintSupport_EXECUTABLE_COMPILE_FLAGS}" )

# Base part of the library  
set(SRCS_QtLennardJones
constants.h
constants.cpp
glwidget.h
glwidget.cpp
main.cpp
mainwindow.h
mainwindow.cpp
MDSystemGL.h
MDSystemGL.cpp
${PROJECT_SOURCE_DIR}/thirdparty/QCustomPlot/qcustomplot.h
${PROJECT_SOURCE_DIR}/thirdparty/QCustomPlot/qcustomplot.cpp
)

include_directories (
        "${PROJECT_SOURCE_DIR}/src/library"
        "${PROJECT_SOURCE_DIR}/thirdparty/QCustomPlot"
        )

if(CMAKE_CUDA_COMPILER)
    include_directories ("${CUDA_INCLUDE_DIRS}")
endif()

add_executable(QtLennardJones ${SRCS_QtLennardJones} qttficon.rc)

set_target_properties(QtLennardJones PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#set_target_properties(QtLennardJones PROPERTIES POSITION_INDEPENDENT_CODE ON)
#set_target_properties(QtLennardJones PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
target_link_libraries(QtLennardJones PUBLIC mdsystem
Qt5::Core Qt5::Widgets Qt5::PrintSupport Qt5::OpenGL Qt5::OpenGLExtensions
${OPENGL_gl_LIBRARY} )

if(APPLE AND CMAKE_CUDA_COMPILER)
    # We need to add the default path to the driver (libcuda.dylib) as an rpath,
    # so that the static cuda runtime can find it at runtime.
    set_property(TARGET QtLennardJones PROPERTY BUILD_RPATH ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
endif()

endif (Qt5Widgets_FOUND AND Qt5PrintSupport_FOUND AND Qt5OpenGL_FOUND AND OpenGL_FOUND)
