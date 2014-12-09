# CMakeList.txt for nonlinear code
add_definitions(-DTRACERNAME=nltracer -DARMA_USE_BLAS -DARMA_USE_LAPACK)

############################## Not so often changed



# -pipe Use pipes rather than temporary files for communication
#           between the various stages of compilation.  This fails to
#           work on some systems where the assembler is unable to
#           read from a pipe; but the GNU assembler has no troubl.


AUX_SOURCE_DIRECTORY(src/tube tubesources)
AUX_SOURCE_DIRECTORY(src/var varsources)
AUX_SOURCE_DIRECTORY(src/seg segsources)
AUX_SOURCE_DIRECTORy(src/sys syssources)
AUX_SOURCE_DIRECTORy(src/tube/geom/grid gridsources)
AUX_SOURCE_DIRECTORy(src/sol solsources)
AUX_SOURCE_DIRECTORy(src/gui guisources)

AUX_SOURCE_DIRECTORy(src/tube/eq tubeeq)
AUX_SOURCE_DIRECTORy(src/tube/geom geomsources)
AUX_SOURCE_DIRECTORy(src/tube/drag tubedrag)
AUX_SOURCE_DIRECTORy(src/tube/heat tubeheat)
AUX_SOURCE_DIRECTORy(src/tube/vertex tubeheat)
AUX_SOURCE_DIRECTORy(src/tube/con tubecon)
# set (nonlinear_SRC ${tubesources} ${varsources} ${segsources}
#   ${syssources}  ${geomsources} ${solsources}
#   ${varsources} ${commonsources}
#   ${tubeeq}  ${tubebc} ${tubedrag} ${gridsources}
#   )
set (nonlinear_SRC ${varsources} ${segsources} ${tubeeq}
  ${geomsources} ${syssources} ${tubeheat} ${tubedrag}
  ${varsources} ${commonsources} ${gridsources}
  ${gridsources} ${tubesources} ${tubecon}
  )

# ${tubedrag}

include_directories(
  include
  src/sys
  src/sol
  src/seg
  src/var
  src/tube
  src/gui
  src/tube/geom
  src/tube/geom/grid
  src/tube/drag
  src/tube/con
  src/tube/vertex
  src/tube/eq
  src/tube/heat
  ../common
  ../common/gas
  ../common/interpolate
  ../common/special
  ../common/fsolve  
  ../common/solid
  ../common/rottfuncs
  /usr/include/umfpack
  )
add_library(gui ${guisources})
add_library(nonlin SHARED ${nonlinear_SRC})

target_link_libraries(nonlin math_common boost_thread boost_system blas lapack boost_iostreams boost_serialization)
# add_executable(multitube multitube.cpp)
# target_link_libraries(multitube nonlin math_common armadillo)

# If the pyx file is a C++ file, we should specify that here.
set_source_files_properties(
  #  ${CYTHON_CMAKE_EXAMPLE_SOURCE_DIR}/src/rect.pyx
  #  ${CYTHON_CMAKE_EXAMPLE_SOURCE_DIR}/src/pi_calculator.pyx
  pytube.pyx doubletube.pyx
  PROPERTIES CYTHON_IS_CXX TRUE )
find_package(wxWidgets COMPONENTS core base REQUIRED)
include("${wxWidgets_USE_FILE}")

set(CYTHON_CXX_EXTENSION cpp)
set(CYTHON_FLAGS -I/usr/local/include --timestamps)
add_executable(test123 test123.cpp)
# add_executable(gui_test gui_test.cpp)
# add_executable(enginetest enginetest.cpp)
# add_executable(threetubes threetubes.cpp)
# cython_add_module(pytube pytube.pyx)
# cython_add_module(doubletube doubletube.pyx)
# set_target_properties(pytube PROPERTIES COMPILE_FLAGS "-Wno-cpp  -Wno-strict-aliasing")


# target_link_libraries(doubletube nonlin math_common armadillo)
# target_link_libraries(gui_test gui ${wxWidgets_LIBRARIES} nonlin math_common blas lapack)
target_link_libraries(test123 nonlin math_common blas lapack boost_iostreams boost_serialization)
# target_link_libraries(enginetest nonlin  math_common blas lapack)
# target_link_libraries(threetubes nonlin math_common blas lapack)
# target_link_libraries(pytube models nonlin math_common blas lapack)


# To execute the test directly
# add_custom_command(TARGET test123 POST_BUILD COMMAND test123 )

# add_library(pytube SHARED pytube.cpp)
# set_target_properties(pytube PROPERTIES PREFIX "")

#enable_testing()
#add_test(test_cpp test COMMAND test)
#add_custom_target(testtest COMMAND test)

#Install library to library directory
#install(TARGETS math_common DESTINATION lib)
#add_library(math_common SHARED ${math_common_SRC})
# Install headers to include directory
#install(FILES ${math_common_HEADERS}  DESTINATION include)


#set_target_properties(math_common PROPERTIES VERSION ${math_common_VERSION_STR#ING} SOVERSION ${math_common_VERSION_MAJOR})














