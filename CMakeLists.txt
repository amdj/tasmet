# CMakeList.txt for nonlinear code
cmake_minimum_required (VERSION 2.6)
#enable_testing()
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_VERBOSE_MAKEFILE on)
set (CMAKE_CXX_FLAGS "-std=c++11 -fPIC -Wall -Wno-unused-variable -Wno-unused-but-set-variable -Wno-return-local-addr")
#configure_file (logger.h.in  "${CMAKE_CURRENT_SOURCE_DIR}/logger.h")

AUX_SOURCE_DIRECTORY(tube tubesources)
AUX_SOURCE_DIRECTORY(var varsources)

set (nonlinear_SRC globalconf.cpp system.cpp ${tubesources} ${varsources})

#  tube/continuityeq.cpp tube/drag.cpp
#    tube/energyeq.cpp tube/geom.cpp tube/momentumeq.cpp
#    tube/solidenergyeq.cpp tube/stateeq.cpp tube/tube.cpp
#    tube/tubeequation.cpp tube/vertex.cpp 

include_directories(
  /usr/include/python2.7
  )

# math_common should be linked to armadillo
add_executable(test test.cpp ${nonlinear_SRC})
add_library(pytube SHARED pytube.cpp)
set_target_properties(pytube PROPERTIES PREFIX "")
target_link_libraries(test math_common armadillo)
target_link_libraries(pytube python2.7)
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









