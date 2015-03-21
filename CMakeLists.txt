# CMakeList.txt for nonlinear code
add_definitions(-DTRACERNAME=nltracer)

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
  swig
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
  ../common/swig_files
  ../common/gas
  ../common/interpolate
  ../common/special
  ../common/fsolve  
  ../common/solid
  ../common/rottfuncs
  )

# add_library(nonlin_lib SHARED ${nonlinear_SRC})
add_library(nonlin_lib  ${nonlinear_SRC})
target_link_libraries(nonlin_lib math_common ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${ARMADILLO_LIBRARY})

add_executable(test123 test123.cpp)
target_link_libraries(test123 nonlin_lib math_common armadillo)


set_source_files_properties( ${swig_generated_file_fullname}
      PROPERTIES COMPILE_FLAGS "-Wno-missing-field-initializers")

# Does this work
set_source_files_properties(nonlin.i PROPERTIES CPLUSPLUS ON)

swig_add_module(nonlin python nonlin.i)
set(SWIG_MODULE_nonlin_FILE_nonlin_EXTRA_DEPS
  ${PROJECT_SOURCE_DIR}/src/common/swig_files/arma_numpy.i
)
swig_link_libraries(nonlin nonlin_lib math_common ${PYTHON_LIBRARIES} armadillo)

# To execute the test directly
# add_custom_command(TARGET test123 POST_BUILD COMMAND test123 )

#enable_testing()
#add_test(test_cpp test COMMAND test)
#add_custom_target(testtest COMMAND test)

#Install library to library directory
#install(TARGETS math_common DESTINATION lib)
# Install headers to include directory
#install(FILES ${math_common_HEADERS}  DESTINATION include)

