
#=========
# The TaSMET code
#=========

AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/ src)
AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/volume volume)
AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/mech mech)
AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/seg seg)
AUX_SOURCE_DIRECTORy(${CMAKE_CURRENT_SOURCE_DIR}/sol sol)
AUX_SOURCE_DIRECTORy(${CMAKE_CURRENT_SOURCE_DIR}/sys sys)
AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/var var)
AUX_SOURCE_DIRECTORy(${CMAKE_CURRENT_SOURCE_DIR}/duct/geom/grid grid)
AUX_SOURCE_DIRECTORy(${CMAKE_CURRENT_SOURCE_DIR}/duct/eq eq)
AUX_SOURCE_DIRECTORy(${CMAKE_CURRENT_SOURCE_DIR}/duct/eq/solid solideq)
AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/duct duct)
AUX_SOURCE_DIRECTORy(${CMAKE_CURRENT_SOURCE_DIR}/duct/geom geom)
AUX_SOURCE_DIRECTORy(${CMAKE_CURRENT_SOURCE_DIR}/duct/drag drag)
AUX_SOURCE_DIRECTORy(${CMAKE_CURRENT_SOURCE_DIR}/duct/cell cell)
AUX_SOURCE_DIRECTORy(${CMAKE_CURRENT_SOURCE_DIR}/duct/heat heat)
AUX_SOURCE_DIRECTORy(${CMAKE_CURRENT_SOURCE_DIR}/duct/connectors connectors)

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
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/sys
  ${CMAKE_CURRENT_SOURCE_DIR}/mech
  ${CMAKE_CURRENT_SOURCE_DIR}/sol
  ${CMAKE_CURRENT_SOURCE_DIR}/seg
  ${CMAKE_CURRENT_SOURCE_DIR}/var
  ${CMAKE_CURRENT_SOURCE_DIR}/volume
  ${CMAKE_CURRENT_SOURCE_DIR}/duct
  ${CMAKE_CURRENT_SOURCE_DIR}/duct/geom
  ${CMAKE_CURRENT_SOURCE_DIR}/duct/geom/grid
  ${CMAKE_CURRENT_SOURCE_DIR}/duct/drag
  ${CMAKE_CURRENT_SOURCE_DIR}/duct/connectors
  ${CMAKE_CURRENT_SOURCE_DIR}/duct/cell
  ${CMAKE_CURRENT_SOURCE_DIR}/duct/eq
  ${CMAKE_CURRENT_SOURCE_DIR}/duct/eq/solid
  ${CMAKE_CURRENT_SOURCE_DIR}/duct/heat
  )
link_directories(${PROJECT_SOURCE_DIR}/common)
# For a shared lib:
# add_library(nonlin_lib SHARED ${nonlinear_SRC})
add_library(nonlin_lib  ${nonlinear_SRC})

target_link_libraries(nonlin_lib common_src)

set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/swig/nonlin.i
  PROPERTIES CPLUSPLUS ON)
swig_add_module(TaSMET python ${CMAKE_CURRENT_SOURCE_DIR}/swig/nonlin.i)

set(SWIG_MODULE_TaSMET_FILE_nonlin_EXTRA_DEPS
  ${CMAKE_CURRENT_SOURCE_DIR}/common/swig_files/arma_numpy.i
)
swig_link_libraries(TaSMET nonlin_lib common_src ${PYTHON_LIBRARIES}
  ${ARMADILLO_LIBRARIES})

# To stop letting the compiler complain about the generated files by Swig

if(${CMAKE_CXX_COMPILER} STREQUAL "g++")

  message("Using gcc compiler")

  set(SWIG_COMMON_COMPILE_FLAGS "-Wno-missing-field-initializers")
  # Swig wrapper file compile flags

elseif(${CMAKE_CXX_COMPILER} STREQUAL "clang++")

  message("Using clang++ compiler")

  # Swig wrapper file compile flags
  set(SWIG_COMMON_COMPILE_FLAGS "-Wno-missing-field-initializers -Wno-self-assign -Wno-#warnings")

endif(${CMAKE_CXX_COMPILER} STREQUAL "g++")

set_source_files_properties( ${swig_generated_file_fullname}
  PROPERTIES COMPILE_FLAGS "${SWIG_COMMON_COMPILE_FLAGS} ")
  
# Install swig files to the right place
install(TARGETS _TaSMET
  DESTINATION ${PYTHON_SITE_PACKAGES}/${PROJECT_NAME})
install(FILES ${CMAKE_BINARY_DIR}/src/TaSMET.py
  DESTINATION ${PYTHON_SITE_PACKAGES}/${PROJECT_NAME})
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/python/__init__.py
  DESTINATION ${PYTHON_SITE_PACKAGES}/${PROJECT_NAME})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python/post
  DESTINATION ${PYTHON_SITE_PACKAGES}/${PROJECT_NAME})
# install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python/gui
#   DESTINATION ${PYTHON_SITE_PACKAGES}/${PROJECT_NAME})

