# CMakeList.txt for nonlinear code
set(CMAKE_CXX_FLAGS "-march=native -mtune=native ")
############################## Optimized code flags
# set(CMAKE_CXX_FLAGS " -O3 ")
# add_definitions(-DLOGGER=0 -DTRACER=0 -DNDEBUG)
# add_definitions(-DEIGEN_NO_DEBUG -DARMA_NO_DEBUG)
############################## End optimized code flags

#PYthon version
set(PYVER 3.4m)

######################### VERBOSE MAKEFILE
#set(CMAKE_VERBOSE_MAKEFILE on)
##########################################





############################## Not so often changed
set(PYTHON_LIBRARY libpython${PYVER}.so)
add_definitions(-DLOGLEVEL=nonlinearlogger)
cmake_minimum_required (VERSION 2.8.5)
project(nonlinear_tmtubes)
# Make the scripts available in the 'cmake' directory available for
# the 'include()' command, 'find_package()' command.

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
${CMAKE_CURRENT_LIST_DIR}/cmake )
set(CLANG_FLAGS "-Wno-unknown-warning-option")

set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CLANG_FLAGS} ")
# set(CMAKE_CXX_COMPILER g++)

# -pipe Use pipes rather than temporary files for communication
#           between the various stages of compilation.  This fails to
#           work on some systems where the assembler is unable to
#           read from a pipe; but the GNU assembler has no trouble.
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pipe -fPIC -Wall \
 -Wno-unused-function -Woverloaded-virtual  \
 -Wno-unused-variable -Wno-unused-but-set-variable \
 -Wno-return-local-addr -Wno-cpp  -Wno-address")

include(UseCython)
AUX_SOURCE_DIRECTORY(tube tubesources)
AUX_SOURCE_DIRECTORY(var varsources)
AUX_SOURCE_DIRECTORY(seg segsources)
AUX_SOURCE_DIRECTORy(sys syssources)
AUX_SOURCE_DIRECTORy(sol solsources)
AUX_SOURCE_DIRECTORy(common commonsources)
set (nonlinear_SRC ${tubesources} ${varsources} ${segsources}
  ${syssources} ${solsources}
  ${varsources} ${commonsources})


include_directories(
/usr/include/python${PYVER}
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
add_executable(threetubes threetubes.cpp)
cython_add_module(pytube pytube.pyx)
# cython_add_module(doubletube doubletube.pyx)
set_target_properties(pytube PROPERTIES COMPILE_FLAGS "-Wno-cpp -Wno-strict-aliasing")

# target_link_libraries(doubletube nonlin math_common armadillo)
target_link_libraries(test123 nonlin math_common armadillo)
target_link_libraries(threetubes nonlin math_common armadillo)
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
