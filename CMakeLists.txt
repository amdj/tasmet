# CMakeList.txt for TaSMET
cmake_minimum_required (VERSION 3.2)
project(TaSMET)

# Tracer name (name of the variable)
add_definitions(-DTRACERNAME=TaSMETTracer)

# This is the common code (gas and solid libs, etc)
add_subdirectory(${PROJECT_SOURCE_DIR}/common)

#==================
# Compiler settings
#==================

# If we would like to cache the compilation process, please activate these settings
SET_PROPERTY(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
SET_PROPERTY(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache) 

# SET THE COMPILER
# set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_COMPILER "clang++")

# Choose some compilation settings
set(CMAKE_CLANG "${CMAKE_CLANG} -Woverloaded-virtual ")
# set(CMAKE_GCC "${CMAKE_GCC} -Wno-unused-function  -ffunction-sections\
#  -fdata-sections -Wno-unused-local-typedefs -Wno-empty-body\
# -Wno-unused-but-set-variable\
# -Wno-return-local-addr -Wno-cpp \
# ")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pipe -fPIC -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -Wno-address ")

# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_GCC}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CLANG}")
add_definitions(-DNPY_NO_DEPRECATED_API=NPY_1_4_API_VERSION)

#=====================
# Optimized code flags
#=====================

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -march=native -mtune=native -fopenmp")
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
# For Nf=12: 625
# Watch out! Setting prealloc too high can give too much overhead for smaller Nf's
# add_definitions(-DARMA_MAT_PREALLOC=625)

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


# Add the code subdirectory
add_subdirectory(${PROJECT_SOURCE_DIR}/src)
