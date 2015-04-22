# CMakeList.txt for linear code
cmake_minimum_required (VERSION 3.1)
project(TATwente)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/common/cmake_tools)

# CMakeList.txt for nonlinear code

# Python
find_package(PythonLibs REQUIRED)

# Initialize swig
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
SET(CMAKE_SWIG_FLAGS -py3 -Wall -DSWIG_PYTHON)

add_subdirectory(${PROJECT_SOURCE_DIR}/common)

# It is of importance that this definition is done AFTER subdirectory common is built
add_definitions(-DTRACERNAME=nltracer)


############################## Optimized code flags
# set(CMAKE_CXX_FLAGS " -O3 ")
set(CMAKE_CXX_FLAGS " -O2 ")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native -fopenmp")
# set(CMAKE_CLANG "${CMAKE_GCC} -march=native -mtune=native -fopenmp")

#add_definitions(-DNDEBUG)	#This is for BOOST functions, and to disable assertions!
# add_definitions(-DTRACER=0)

#Debug mode
# set(CMAKE_GCC "${CMAKE_GCC} -g -ggdb")
############################## End optimized code flags

############################## Compiler settings
# set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_COMPILER "clang++")
# set(CMAKE_CXX_COMPILER_ARG1="clang++")

set(CMAKE_CLANG "${CMAKE_CLANG} -pipe -Woverloaded-virtual")
# set(CMAKE_GCC "${CMAKE_GCC} -Wno-unused-function  -ffunction-sections\
#  -fdata-sections -Wno-unused-local-typedefs -Wno-empty-body\
# -Wno-unused-but-set-variable\
# -Wno-return-local-addr -Wno-cpp \
# ")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pipe -fPIC -Wall \
 -Wextra -Wno-unused-parameter \
 -Wno-unused-variable  \
 -Wno-address")

# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_GCC}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CLANG}")
add_definitions(-DNPY_NO_DEPRECATED_API=NPY_1_4_API_VERSION)

SET_PROPERTY(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
SET_PROPERTY(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache) 


############################## End compiler settings

AUX_SOURCE_DIRECTORY(${PROJECT_SOURCE_DIR}/src/ src)
AUX_SOURCE_DIRECTORY(${PROJECT_SOURCE_DIR}/src/seg seg)
AUX_SOURCE_DIRECTORy(${PROJECT_SOURCE_DIR}/src/sol sol)
AUX_SOURCE_DIRECTORy(${PROJECT_SOURCE_DIR}/src/sys sys)
AUX_SOURCE_DIRECTORY(${PROJECT_SOURCE_DIR}/src/var var)
AUX_SOURCE_DIRECTORy(${PROJECT_SOURCE_DIR}/src/tube/geom/grid grid)
AUX_SOURCE_DIRECTORy(${PROJECT_SOURCE_DIR}/src/tube/eq eq)
AUX_SOURCE_DIRECTORY(${PROJECT_SOURCE_DIR}/src/tube tube)
AUX_SOURCE_DIRECTORy(${PROJECT_SOURCE_DIR}/src/tube/geom geom)
AUX_SOURCE_DIRECTORy(${PROJECT_SOURCE_DIR}/src/tube/drag drag)
AUX_SOURCE_DIRECTORy(${PROJECT_SOURCE_DIR}/src/tube/cell cell)
AUX_SOURCE_DIRECTORy(${PROJECT_SOURCE_DIR}/src/tube/heat heat)
AUX_SOURCE_DIRECTORy(${PROJECT_SOURCE_DIR}/src/tube/con con)

set (nonlinear_SRC
  ${src}
  ${seg}
  ${sol}
  ${sys}
  ${var}
  ${tube}
  ${con}
  ${drag}
  ${eq}
  ${geom}
  ${grid}
  ${heat}
  ${cell}
  )

include_directories(
  ${PYTHON_INCLUDE_PATH}
  ${ARMADILLO_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/common/src
  ${PROJECT_SOURCE_DIR}/common/src/swig
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_SOURCE_DIR}/src/sys
  ${PROJECT_SOURCE_DIR}/src/sol
  ${PROJECT_SOURCE_DIR}/src/seg
  ${PROJECT_SOURCE_DIR}/src/var
  ${PROJECT_SOURCE_DIR}/src/tube
  ${PROJECT_SOURCE_DIR}/src/tube/geom
  ${PROJECT_SOURCE_DIR}/src/tube/geom/grid
  ${PROJECT_SOURCE_DIR}/src/tube/drag
  ${PROJECT_SOURCE_DIR}/src/tube/con
  ${PROJECT_SOURCE_DIR}/src/tube/cell
  ${PROJECT_SOURCE_DIR}/src/tube/eq
  ${PROJECT_SOURCE_DIR}/src/tube/heat
  )
link_directories(${PROJECT_SOURCE_DIR}/common)
# add_library(nonlin_lib SHARED ${nonlinear_SRC})
add_library(nonlin_lib  ${nonlinear_SRC})

target_link_libraries(nonlin_lib common ${ARMADILLO_LIBRARY})

set_source_files_properties(${PROJECT_SOURCE_DIR}/src/swig/nonlin.i
  PROPERTIES CPLUSPLUS ON)

swig_add_module(TATwente python src/swig/nonlin.i)

set(SWIG_MODULE_nonlin_FILE_nonlin_EXTRA_DEPS
  ${PROJECT_SOURCE_DIR}/src/common/swig_files/arma_numpy.i
)

swig_link_libraries(TATwente nonlin_lib common ${PYTHON_LIBRARIES}
  ${ARMADILLO_LIBRARIES})
set_source_files_properties( ${swig_generated_file_fullname}
      PROPERTIES COMPILE_FLAGS "-Wno-missing-field-initializers -Wno-self-assign -Wno-#warnings")

# To execute the test directly
# add_custom_command(TARGET test123 POST_BUILD COMMAND test123 )

#enable_testing()
#add_test(test_cpp test COMMAND test)
#add_custom_target(testtest COMMAND test)

#Install library to library directory
#install(TARGETS common DESTINATION lib)
# Install headers to include directory
#install(FILES ${common_HEADERS}  DESTINATION include)

# add_executable(test123 test123.cpp)
# target_link_libraries(test123 nonlin_lib common armadillo)

# add_definitions(-DNODRAG)

# Does this work
