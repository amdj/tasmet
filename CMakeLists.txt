# CMakeList.txt for TaSMET
cmake_minimum_required (VERSION 2.8)
project(TaSMET)
set(PACKAGE_VERSION 0.1)
message("Running Cmake for TaSMET version ${PACKAGE_VERSION}")

# Set the Python version to link against. If this is not set, Cmake tries to find it automatically.
# set(TaSMET_PY_VERSION "2.7")
# set(TaSMET_PY_VERSION "3.5m")

# Tracer name (name of the variable)
add_definitions(-DTRACERNAME=TaSMETTracer)
# add_definitions(-DTRACER_IN_COMMON)

#====================================================
# Compiler settings ************************************************
#====================================================

# SET THE COMPILER
# set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_COMPILER "clang++")

# Always required make flags in case of both compilers
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pipe -fPIC -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable  ")


if(${CMAKE_CXX_COMPILER} STREQUAL "g++")
  message("Compiler chosen: " ${CMAKE_CXX_COMPILER})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -ffunction-sections -fdata-sections -Wno-unused-local-typedefs")

elseif(${CMAKE_CXX_COMPILER} STREQUAL "clang++")
  message("Compiler chosen: " ${CMAKE_CXX_COMPILER})
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Woverloaded-virtual" )  

endif(${CMAKE_CXX_COMPILER} STREQUAL "g++")

# Stop letting Numpy complain about its API
add_definitions(-DNPY_NO_DEPRECATED_API=NPY_1_4_API_VERSION)

#==================================================
# Optimized code flags *******************************************
#==================================================

#Debug mode
# set(CMAKE_CXX_FLAGS  "${CMAKE_CXX)FLAGS} -g -ggdb")

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -march=native -mtune=native")
# set(CMAKE_CLANG "${CMAKE_GCC} -march=native -mtune=native -fopenmp")

# To disable bound checking on std::vector, and to disable assertions
# add_definitions(-DNDEBUG)

# To increase speed on Armadillo
# add_definitions(-DARMA_NO_DEBUG)

# Disable traces
# add_definitions(-DTRACER=0)

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


# ##########################
# Finding the presence of the prerequisites
# ##########################


# For importing find directives for Cmake
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/common/cmake_tools)

# ##########################
# Python #####################
# ##########################

if(TaSMET_PY_VERSION)
  # Find major version from version string
  set(PYTHON_LIBRARY "/usr/lib/libpython${TaSMET_PY_VERSION}.so")
  set(PYTHON_INCLUDE_DIR "/usr/include/python${TaSMET_PY_VERSION}")  
  set(PYTHON_INCLUDE_DIRS "/usr/include/python${TaSMET_PY_VERSION}")  

endif(TaSMET_PY_VERSION)
message("Python include dirs: ${PYTHON_INCLUDE_DIRS}")

find_package(PythonLibs REQUIRED)
string(REGEX MATCH "^." TaSMET_PY_MAJOR_VERSION ${PYTHONLIBS_VERSION_STRING})
MESSAGE("Python major version: ${TaSMET_PY_MAJOR_VERSION}")  

# Find the site_packages directory of python
execute_process(COMMAND python${TaSMET_PY_MAJOR_VERSION} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)

# ################################
# Initialize swig
# ################################

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

SET(CMAKE_SWIG_FLAGS  -Wall -DSWIG_PYTHON)
if(${TaSMET_PY_MAJOR_VERSION}=="3")
  SET(CMAKE_SWIG_FLAGS  ${CMAKE_SWIG_FLAGS} -py${TaSMET_PY_MAJOR_VERSION})  
endif(${TaSMET_PY_MAJOR_VERSION}=="3")

# Armadillo
find_package(Armadillo REQUIRED)
add_definitions(-DARMA_USE_SUPERLU -DARMA_USE_CXX11)

# ==================== Compile the code in common and src

# This is the common code (gas and solid libs, etc)
add_subdirectory(${PROJECT_SOURCE_DIR}/common)

# Add the code subdirectory
add_subdirectory(${PROJECT_SOURCE_DIR}/src)

# ================================== Installation

# Install common files
install(FILES ${PROJECT_SOURCE_DIR}/common/_common.so
  DESTINATION ${PYTHON_SITE_PACKAGES}/${PROJECT_NAME}/common)
install(FILES ${PROJECT_SOURCE_DIR}/common/__init__.py
  DESTINATION ${PYTHON_SITE_PACKAGES}/${PROJECT_NAME}/common)
install(FILES ${PROJECT_SOURCE_DIR}/common/common.py
  DESTINATION ${PYTHON_SITE_PACKAGES}/${PROJECT_NAME}/common)

# Rest of the files is installed from src/CMakeLists.txt
