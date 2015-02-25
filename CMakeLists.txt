# CMakeList.txt for nonlinear code
add_definitions(-DTRACERNAME=nltracer -DARMA_USE_BLAS -DARMA_USE_LAPACK)

############################## Not so often changed



# -pipe Use pipes rather than temporary files for communication
#           between the various stages of compilation.  This fails to
#           work on some systems where the assembler is unable to
#           read from a pipe; but the GNU assembler has no troubl.


AUX_SOURCE_DIRECTORY(src/seg seg)
AUX_SOURCE_DIRECTORy(src/sol sol)
AUX_SOURCE_DIRECTORy(src/sys sys)
AUX_SOURCE_DIRECTORY(src/var var)
AUX_SOURCE_DIRECTORY(src/tube tube)

AUX_SOURCE_DIRECTORy(src/tube/geom/grid grid)
AUX_SOURCE_DIRECTORy(src/tube/eq eq)
AUX_SOURCE_DIRECTORy(src/tube/geom geom)
AUX_SOURCE_DIRECTORy(src/tube/drag drag)
AUX_SOURCE_DIRECTORy(src/tube/heat heat)
AUX_SOURCE_DIRECTORy(src/tube/vertex vertex)
AUX_SOURCE_DIRECTORy(src/tube/con con)
# set (nonlinear_SRC ${tubesources} ${varsources} ${segsources}
#   ${syssources}  ${geomsources} ${solsources}
#   ${varsources} ${commonsources}
#   ${tubeeq}  ${tubebc} ${tubedrag} ${gridsources}
#   )
set (nonlinear_SRC ${seg} ${sol} ${sys} ${var} ${tube}
  #bc
  ${con} ${drag} ${eq} ${geom} ${grid} ${heat} ${vertex}
  )

include_directories(
  include
  src/sys
  src/sol
  src/seg
  src/var
  src/tube
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

add_library(nonlin SHARED ${nonlinear_SRC})
target_link_libraries(nonlin math_common boost_thread boost_system boost_iostreams armadillo)

# If the pyx file is a C++ file, we should specify that here.
# find_package(wxWidgets COMPONENTS core base REQUIRED)
# include("${wxWidgets_USE_FILE}")

add_executable(test123 test123.cpp)
# set_target_properties(pytube PROPERTIES COMPILE_FLAGS "-Wno-cpp  -Wno-strict-aliasing")

# target_link_libraries(gui_test gui ${wxWidgets_LIBRARIES} nonlin math_common blas lapack)
target_link_libraries(test123 nonlin math_common armadillo)

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

