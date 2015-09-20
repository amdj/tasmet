
#=========
# The TaSMET code
#=========

# AUX_SOURCE_DIRECTORY(${PROJECT_SOURCE_DIR}/src/ src)
AUX_SOURCE_DIRECTORY(${PROJECT_SOURCE_DIR}/src/volume volume)
AUX_SOURCE_DIRECTORY(${PROJECT_SOURCE_DIR}/src/mech mech)
AUX_SOURCE_DIRECTORY(${PROJECT_SOURCE_DIR}/src/seg seg)
AUX_SOURCE_DIRECTORy(${PROJECT_SOURCE_DIR}/src/sol sol)
AUX_SOURCE_DIRECTORy(${PROJECT_SOURCE_DIR}/src/sys sys)
AUX_SOURCE_DIRECTORY(${PROJECT_SOURCE_DIR}/src/var var)
AUX_SOURCE_DIRECTORy(${PROJECT_SOURCE_DIR}/src/duct/geom/grid grid)
AUX_SOURCE_DIRECTORy(${PROJECT_SOURCE_DIR}/src/duct/eq eq)
AUX_SOURCE_DIRECTORy(${PROJECT_SOURCE_DIR}/src/duct/eq/solid solideq)
AUX_SOURCE_DIRECTORY(${PROJECT_SOURCE_DIR}/src/duct duct)
AUX_SOURCE_DIRECTORy(${PROJECT_SOURCE_DIR}/src/duct/geom geom)
AUX_SOURCE_DIRECTORy(${PROJECT_SOURCE_DIR}/src/duct/drag drag)
AUX_SOURCE_DIRECTORy(${PROJECT_SOURCE_DIR}/src/duct/cell cell)
AUX_SOURCE_DIRECTORy(${PROJECT_SOURCE_DIR}/src/duct/heat heat)
AUX_SOURCE_DIRECTORy(${PROJECT_SOURCE_DIR}/src/duct/connectors connectors)

set (nonlinear_SRC
  ${src}
  ${seg}
  ${volume}
  ${mech}
  ${sol}
  ${sys}
  ${var}
  ${duct}
  ${connectors}
  ${drag}
  ${eq}
  ${solideq}
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
  ${PROJECT_SOURCE_DIR}/src/mech
  ${PROJECT_SOURCE_DIR}/src/sol
  ${PROJECT_SOURCE_DIR}/src/seg
  ${PROJECT_SOURCE_DIR}/src/var
  ${PROJECT_SOURCE_DIR}/src/volume
  ${PROJECT_SOURCE_DIR}/src/duct
  ${PROJECT_SOURCE_DIR}/src/duct/geom
  ${PROJECT_SOURCE_DIR}/src/duct/geom/grid
  ${PROJECT_SOURCE_DIR}/src/duct/drag
  ${PROJECT_SOURCE_DIR}/src/duct/connectors
  ${PROJECT_SOURCE_DIR}/src/duct/cell
  ${PROJECT_SOURCE_DIR}/src/duct/eq
  ${PROJECT_SOURCE_DIR}/src/duct/eq/solid
  ${PROJECT_SOURCE_DIR}/src/duct/heat
  )
link_directories(${PROJECT_SOURCE_DIR}/common)
# For a shared lib:
# add_library(nonlin_lib SHARED ${nonlinear_SRC})
add_library(nonlin_lib  ${nonlinear_SRC})

target_link_libraries(nonlin_lib common_src)

set_source_files_properties(${PROJECT_SOURCE_DIR}/src/swig/nonlin.i
  PROPERTIES CPLUSPLUS ON)
swig_add_module(tasmet python ${PROJECT_SOURCE_DIR}/src/swig/nonlin.i)

set(SWIG_MODULE_nonlin_FILE_nonlin_EXTRA_DEPS
  ${PROJECT_SOURCE_DIR}/src/common/swig_files/arma_numpy.i
)
swig_link_libraries(tasmet nonlin_lib common_src ${PYTHON_LIBRARIES}
  ${ARMADILLO_LIBRARIES})
set_source_files_properties( ${swig_generated_file_fullname}
      PROPERTIES COMPILE_FLAGS "-Wno-missing-field-initializers -Wno-self-assign -Wno-#warnings")
