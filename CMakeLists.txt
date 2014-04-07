# CMakeList.txt for nonlinear code
cmake_minimum_required (VERSION 2.8.5)
project(nonlinear_tmtubes)
# Make the scripts available in the 'cmake' directory available for the
# 'include()' command, 'find_package()' command.
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/cmake )
include(UseCython)
add_definitions(-DEIGEN_NO_DEBUG -DARMA_NO_DEBUG)
#enable_testing()
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_VERBOSE_MAKEFILE on)
set (CMAKE_CXX_FLAGS "-std=c++11 -fopenmp -O2 -fPIC -Wall -Wno-unused-function  -Wno-unused-variable -Wno-unused-but-set-variable -Wno-return-local-addr -Wno-cpp -Wno-address")
#configure_file (logger.h.in  "${CMAKE_CURRENT_SOURCE_DIR}/logger.h")

AUX_SOURCE_DIRECTORY(tube tubesources)
AUX_SOURCE_DIRECTORY(var varsources)

set (nonlinear_SRC globalconf.cpp ${tubesources} ${varsources})

#  tube/continuityeq.cpp tube/drag.cpp
#    tube/energyeq.cpp tube/geom.cpp tube/momentumeq.cpp
#    tube/solidenergyeq.cpp tube/stateeq.cpp tube/tube.cpp
#    tube/tubeequation.cpp tube/vertex.cpp 

include_directories(
  /usr/include/python2.7
  /usr/include/eigen3
  /usr/local/include
  )

# math_common should be linked to armadillo
set(PYTHON_LIBRARIES libpython2.7.so)
set(PYTHON_LIBRARY libpython2.7.so)
add_executable(test test.cpp ${nonlinear_SRC})
# If the pyx file is a C++ file, we should specify that here.
set_source_files_properties(
#  ${CYTHON_CMAKE_EXAMPLE_SOURCE_DIR}/src/rect.pyx
#  ${CYTHON_CMAKE_EXAMPLE_SOURCE_DIR}/src/pi_calculator.pyx
  tube.pyx
  PROPERTIES CYTHON_IS_CXX TRUE )
set(CYTHON_CXX_EXTENSION cpp)
set(CYTHON_FLAGS -I/usr/local/include --timestamps)
cython_add_module(tube tube.pyx ${nonlinear_SRC})
target_link_libraries(tube math_common armadillo)
#add_library(pytube SHARED pytube.cpp)
#set_target_properties(pytube PROPERTIES PREFIX "")
target_link_libraries(test math_common armadillo)
#target_link_libraries(pytube python2.7)
# The headers which have to be installed
#set (math_common_HEADERS math_common.h logger.h vtypes.h material.h bessel.h)
#enable_testing()
#add_test(test_cpp test COMMAND test)
#add_custom_target(testtest COMMAND test)
add_custom_command(TARGET test POST_BUILD COMMAND test)
#Install library to library directory
#install(TARGETS math_common DESTINATION lib)
#add_library(math_common SHARED ${math_common_SRC})
# Install headers to include directory
#install(FILES ${math_common_HEADERS}  DESTINATION include)


#set_target_properties(math_common PROPERTIES VERSION ${math_common_VERSION_STR#ING} SOVERSION ${math_common_VERSION_MAJOR})









