# CMakeList.txt for nonlinear code
add_definitions(-DTRACERNAME=nltracer)

############################## Not so often changed

# -pipe Use pipes rather than temporary files for communication
#           between the various stages of compilation.  This fails to
#           work on some systems where the assembler is unable to
#           read from a pipe; but the GNU assembler has no troubl.

AUX_SOURCE_DIRECTORY(src/ src)
AUX_SOURCE_DIRECTORY(src/seg seg)
AUX_SOURCE_DIRECTORy(src/sol sol)
AUX_SOURCE_DIRECTORy(src/sys sys)
AUX_SOURCE_DIRECTORY(src/var var)

AUX_SOURCE_DIRECTORy(src/tube/geom/grid grid)
# AUX_SOURCE_DIRECTORy(src/tube/eq eq)
set(eqdir src/tube/eq)
set(eq ${eqdir}/mu.cpp ${eqdir}/tubeequation.cpp ${eqdir}/continuity.cpp  ${eqdir}/state.cpp
  ${eqdir}/isentropic.cpp ${eqdir}/weightfactors.cpp ${eqdir}/solidenergy.cpp ${eqdir}/momentum.cpp)

# AUX_SOURCE_DIRECTORY(src/tube tube)
set(tubedir src/tube)
set(tube ${tubedir}/tube.cpp ${tubedir}/isentropictube.cpp)

AUX_SOURCE_DIRECTORy(src/tube/geom geom)
# AUX_SOURCE_DIRECTORy(src/tube/drag drag)
set(dragdir src/tube/drag)
set(drag ${dragdir}/drag.cpp)

# AUX_SOURCE_DIRECTORy(src/tube/heat heat)
set(heatdir src/tube/heat)
set(heat ${heatdir}/heat.cpp)
AUX_SOURCE_DIRECTORy(src/tube/cell cell)

# AUX_SOURCE_DIRECTORy(src/tube/con con)
set(condir src/tube/con)
set(con ${condir}/tubebc.cpp ${condir}/adiabaticwall.cpp ${condir}/prescribeqty.cpp
  ${condir}/pressurebc.cpp)
# set (nonlinear_SRC ${tubesources} ${varsources} ${segsources}
#   ${syssources}  ${geomsources} ${solsources}
#   ${varsources} ${commonsources}
#   ${tubeeq}  ${tubebc} ${tubedrag} ${gridsources}
#   )
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
  swig
  src
  src/sys
  src/sol
  src/seg
  src/var
  src/tube
  src/tube/geom
  src/tube/geom/grid
  src/tube/drag
  src/tube/con
  src/tube/cell
  src/tube/eq
  src/tube/heat
  ../common
  ../common/exception
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

# add_executable(test123 test123.cpp)
# target_link_libraries(test123 nonlin_lib math_common armadillo)

add_definitions(-DNODRAG)

# Does this work
set_source_files_properties(nonlin.i PROPERTIES CPLUSPLUS ON)

swig_add_module(nonlin python nonlin.i)
set(SWIG_MODULE_nonlin_FILE_nonlin_EXTRA_DEPS
  ${PROJECT_SOURCE_DIR}/src/common/swig_files/arma_numpy.i
)
swig_link_libraries(nonlin nonlin_lib math_common ${PYTHON_LIBRARIES} armadillo)
set_source_files_properties( ${swig_generated_file_fullname}
      PROPERTIES COMPILE_FLAGS "-Wno-missing-field-initializers -Wno-self-assign -Wno-#warnings")


# To execute the test directly
# add_custom_command(TARGET test123 POST_BUILD COMMAND test123 )

#enable_testing()
#add_test(test_cpp test COMMAND test)
#add_custom_target(testtest COMMAND test)

#Install library to library directory
#install(TARGETS math_common DESTINATION lib)
# Install headers to include directory
#install(FILES ${math_common_HEADERS}  DESTINATION include)

