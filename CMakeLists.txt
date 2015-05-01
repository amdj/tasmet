# CMakeList.txt for nonlinear code
#=============
# General info
#=============

cmake_minimum_required (VERSION 3.2)
project(TATwente)

# Tracer name (name of the variable)
add_definitions(-DTRACERNAME=TATwenteTracer)
add_subdirectory(${PROJECT_SOURCE_DIR}/common)
#==================
# Compiler settings
#==================
SET_PROPERTY(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
SET_PROPERTY(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache) 

set(CMAKE_CLANG "${CMAKE_CLANG} -Woverloaded-virtual ")
# set(CMAKE_GCC "${CMAKE_GCC} -Wno-unused-function  -ffunction-sections\
#  -fdata-sections -Wno-unused-local-typedefs -Wno-empty-body\
# -Wno-unused-but-set-variable\
# -Wno-return-local-addr -Wno-cpp \
# ")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pipe -fPIC -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -Wno-address ")


# SET THE COMPILER
# set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_COMPILER "clang++")

# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_GCC}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CLANG}")
add_definitions(-DNPY_NO_DEPRECATED_API=NPY_1_4_API_VERSION)

#=====================
# Optimized code flags
#=====================

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native -fopenmp")
# set(CMAKE_CLANG "${CMAKE_GCC} -march=native -mtune=native -fopenmp")

# Pre-allocation size for matrices. Very important setting to tune the
# code in case we know that we are going to run with certain sizes of matrices and vectors. For Nf=6, set this to (2*Nf+1)^2=169
# For Nf=1: 9
# For Nf=2: 25
# For Nf=3: 49
# For Nf=4: 81
# For Nf=5: 121
# For Nf=6: 169
# For Nf=7: 225
# For Nf=8: 289
# For Nf=9: 361
# For Nf=10: 441
# For Nf=11: 529

# Watch out! Setting prealloc too high can give too much overhead for smaller Nf's
add_definitions(-DARMA_MAT_PREALLOC=169)


# To disable bound checking on std::vector, and to disable assertions
# add_definitions(-DNDEBUG)

# To increase speed on Armadillo
# add_definitions(-DARMA_NO_DEBUG)

# Disable traces
# add_definitions(-DTRACER=0)

#Debug mode
# set(CMAKE_GCC "${CMAKE_GCC} -g -ggdb")


#==================
# Required packages
#==================

# For importing find directives for Cmake
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/common/cmake_tools)

# Python
find_package(PythonLibs REQUIRED)

# Initialize swig
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
SET(CMAKE_SWIG_FLAGS -py3 -Wall -DSWIG_PYTHON)

# Armadillo
find_package(Armadillo REQUIRED)
add_definitions(-DARMA_USE_SUPERLU -DARMA_USE_CXX11 -DARMA_USE_BLAS -DARMA_USE_LAPACK)


#============
# Make common
#============

# add_subdirectory(${PROJECT_SOURCE_DIR}/common)


#=========
# The code
#=========

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
# For a shared lib:
# add_library(nonlin_lib SHARED ${nonlinear_SRC})
add_library(nonlin_lib  ${nonlinear_SRC})

target_link_libraries(nonlin_lib common_src)

set_source_files_properties(${PROJECT_SOURCE_DIR}/src/swig/nonlin.i
  PROPERTIES CPLUSPLUS ON)
swig_add_module(TATwente python src/swig/nonlin.i)

set(SWIG_MODULE_nonlin_FILE_nonlin_EXTRA_DEPS
  ${PROJECT_SOURCE_DIR}/src/common/swig_files/arma_numpy.i
)
swig_link_libraries(TATwente nonlin_lib common_src ${PYTHON_LIBRARIES}
  ${ARMADILLO_LIBRARIES})
set_source_files_properties( ${swig_generated_file_fullname}
      PROPERTIES COMPILE_FLAGS "-Wno-missing-field-initializers -Wno-self-assign -Wno-#warnings")
