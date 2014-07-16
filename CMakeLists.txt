# CMakeList.txt for nonlinear code

#CMakeLists.txt for branch master, producing optimized code (slower)
# math_common should be linked to armadillo
set(PYVER 3.4m)
set(PYTHON_LIBRARY libpython${PYVER}.so)

cmake_minimum_required (VERSION 2.8.5)
project(nonlinear_tmtubes)
# Make the scripts available in the 'cmake' directory available for the
# 'include()' command, 'find_package()' command.
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/cmake )
include(UseCython)
#add_definitions(-DEIGEN_NO_DEBUG -DARMA_NO_DEBUG)
add_definitions(-DLOGLEVEL=nonlinearlogger)
#enable_testing()

set(CMAKE_C_COMPILER glang)
set(CMAKE_CXX_COMPILER clang++)

# set(CMAKE_CXX_COMPILER g++)
######################### VERBOSE MAKEFILE
#set(CMAKE_VERBOSE_MAKEFILE on)
##########################################
 # -pipe
 #           Use pipes rather than temporary files for communication between the various stages of compilation.  This fails to work on some systems where the
 #           assembler is unable to read from a pipe; but the GNU assembler has no trouble.


 
 #Clang++ flags
set (CMAKE_CXX_FLAGS "-std=c++11 -pipe -march=native -fPIC -Wall -Wno-unused-function  -Wno-unused-variable -Wno-unused-but-set-variable -Wno-return-local-addr -Wno-cpp  -Wno-address")
#Optimized clang++ flags
set(CMAKE_CXX_FLAGS $(CMAKE_CXX_FLAGS) "-O2 -march=x86_64 ")


# set (CMAKE_CXX_FLAGS "-std=c++11 -march=native -fPIC -Wall -Wno-unused-function -Wno-unknown-warning-option  -Wno-unused-variable -Wno-unused-but-set-variable -Wno-return-local-addr -Wno-cpp  -Wno-address")

#set (CMAKE_CXX_FLAGS "-std=c++11 -fopenmp -O3 -march=native -fPIC -Wall -Wno-unused-function  -Wno-unused-variable -Wno-unused-but-set-variable -Wno-return-local-addr -Wno-cpp  -Wno-address  -ffunction-sections -fdata-sections")
set(CMAKE_CXX_FLAGS "-std=c++11 -fPIC")
#configure_file (logger.h.in  "${CMAKE_CURRENT_SOURCE_DIR}/logger.h")

AUX_SOURCE_DIRECTORY(tube tubesources)
AUX_SOURCE_DIRECTORY(var varsources)
AUX_SOURCE_DIRECTORY(seg segsources)
AUX_SOURCE_DIRECTORy(sys syssources)
AUX_SOURCE_DIRECTORy(sol solsources)
AUX_SOURCE_DIRECTORy(common commonsources)
set (nonlinear_SRC ${tubesources} ${varsources} ${segsources}
${syssources} ${solsources}
${varsources} ${commonsources})

#  tube/continuityeq.cpp tube/drag.cpp
#    tube/energyeq.cpp tube/geom.cpp tube/momentumeq.cpp
#    tube/solidenergyeq.cpp tube/stateeq.cpp tube/tube.cpp
#    tube/tubeequation.cpp tube/vertex.cpp 

include_directories(
  /usr/include/eigen3
  /usr/local/include
  /usr/local/include/eigen3
  sol
  sys
  seg
  var
  tube
  common
  models
  )


add_library(nonlin ${nonlinear_SRC})

# add_executable(multitube multitube.cpp)
# target_link_libraries(multitube nonlin math_common armadillo)

# If the pyx file is a C++ file, we should specify that here.
set_source_files_properties(
#  ${CYTHON_CMAKE_EXAMPLE_SOURCE_DIR}/src/rect.pyx
#  ${CYTHON_CMAKE_EXAMPLE_SOURCE_DIR}/src/pi_calculator.pyx
  pytube.pyx doubletube.pyx
  PROPERTIES CYTHON_IS_CXX TRUE )

set(CYTHON_CXX_EXTENSION cpp)
set(CYTHON_FLAGS -I/usr/local/include --timestamps)
add_executable(test123 test123.cpp)
cython_add_module(pytube pytube.pyx)
# cython_add_module(doubletube doubletube.pyx)
set_target_properties(pytube PROPERTIES COMPILE_FLAGS "-Wno-cpp -Wno-strict-aliasing")

# target_link_libraries(doubletube nonlin math_common armadillo)
target_link_libraries(test123 nonlin math_common armadillo)
target_link_libraries(pytube nonlin math_common armadillo)


# To execute the test directly
# add_custom_command(TARGET test123 POST_BUILD COMMAND test123 )

# add_library(pytube SHARED pytube.cpp)
# set_target_properties(pytube PROPERTIES PREFIX "")

#target_link_libraries(pytube python2.7)
# The headers which have to be installed
#set (math_common_HEADERS math_common.h logger.h vtypes.h material.h bessel.h)
#enable_testing()
#add_test(test_cpp test COMMAND test)
#add_custom_target(testtest COMMAND test)

#Install library to library directory
#install(TARGETS math_common DESTINATION lib)
#add_library(math_common SHARED ${math_common_SRC})
# Install headers to include directory
#install(FILES ${math_common_HEADERS}  DESTINATION include)


#set_target_properties(math_common PROPERTIES VERSION ${math_common_VERSION_STR#ING} SOVERSION ${math_common_VERSION_MAJOR})









