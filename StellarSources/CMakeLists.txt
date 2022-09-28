project(singlestar C Fortran)

message(STATUS "Making Stellar Sources...")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# List source files
set ( VH1_SRC_MOD_FILES
    "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/global.f90"
    "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/zonemod.f90"
)
set ( VH1_SRC_FILES
    "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/dtcon.f90"
    "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/dump.f90"
    "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/init.f90"
    "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/prin.f90"
    "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/sweepx.f90"
    "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/sweepy.f90"
    "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/sweepz.f90"
    "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/images.f90"
    "${CMAKE_SOURCE_DIR}/VH1/src/Patch/PPMLR/ppmlr.f90"
    "${CMAKE_SOURCE_DIR}/VH1/src/Patch/PPMLR/forces.f90"
    "${CMAKE_SOURCE_DIR}/VH1/src/Patch/PPMLR/flatten.f90"
    "${CMAKE_SOURCE_DIR}/VH1/src/Patch/PPMLR/evolve.f90"
    "${CMAKE_SOURCE_DIR}/VH1/src/Patch/PPMLR/remap.f90"
    "${CMAKE_SOURCE_DIR}/VH1/src/Patch/PPMLR/states.f90"
    "${CMAKE_SOURCE_DIR}/VH1/src/Patch/PPMLR/boundary.f90"
    "${CMAKE_SOURCE_DIR}/VH1/src/Patch/PPMLR/volume.f90"
    "${CMAKE_SOURCE_DIR}/VH1/src/Patch/PPMLR/riemann.f90"
    "${CMAKE_SOURCE_DIR}/VH1/src/Patch/PPMLR/parabola.f90"
)
set ( VH1_OBJ_MOD_FILES
    global
    zonemod
    )
set ( VH1_OBJ_FILES
    dtcon
    dump
    init
    prin
    sweepx
    sweepy
    sweepz
    images
    ppmlr
    forces
    flatten
    evolve
    remap
    states
    boundary
    volume
    riemann
    parabola
)


foreach(modsrcfile modobjfile IN ZIP_LISTS VH1_SRC_MOD_FILES VH1_OBJ_MOD_FILES)
    add_library(${modobjfile} OBJECT ${modsrcfile})
endforeach()
foreach(srcfile objfile IN ZIP_LISTS VH1_SRC_FILES VH1_OBJ_FILES)
    add_library(${objfile} OBJECT ${srcfile})
    foreach(modsrcfile modobjfile IN ZIP_LISTS VH1_SRC_MOD_FILES VH1_OBJ_MOD_FILES)
        add_dependencies(${objfile} ${modobjfile})
    endforeach()
endforeach()

set(f2py_module_name "vhone")
set(fortran_src_file "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/vhone.f90")
set(generated_module_file ${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})

add_custom_target(${f2py_module_name} ALL
  DEPENDS "${fortran_src_file}"
  )

add_custom_command(
  OUTPUT "${f2py_module_name}module.c" "${f2py_module_name}-f2pywrappers2.f90"
  COMMAND ${F2PY_EXECUTABLE}
    -m ${f2py_module_name}
    ${fortran_src_file}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${fortran_src_file}
  )

add_library(${generated_module_file} MODULE
            "${f2py_module_name}module.c"
            "${F2PY_INCLUDE_DIR}/fortranobject.c"
            "${f2py_module_name}-f2pywrappers2.f90"
             ${fortran_src_file})

target_link_libraries(${generated_module_file} ${VH1_OBJ_FILES} ${VH1_OBJ_MOD_FILES})
target_include_directories(${generated_module_file} PUBLIC
                           ${F2PY_INCLUDE_DIRS}
                           ${PYTHON_INCLUDE_DIRS})
set_target_properties(${generated_module_file} PROPERTIES SUFFIX "")
set_target_properties(${generated_module_file} PROPERTIES PREFIX "")


if (UNIX)
  if (APPLE)
    set_target_properties(${generated_module_file} PROPERTIES
    LINK_FLAGS  '-Wl,-dylib,-undefined,dynamic_lookup')
  else()
    set_target_properties(${generated_module_file} PROPERTIES
  LINK_FLAGS  '-Wl,--allow-shlib-undefined')
  endif()
endif()

if (SKBUILD)
  install(TARGETS ${generated_module_file} DESTINATION vhone)
else()
  install(TARGETS ${generated_module_file} DESTINATION ${CMAKE_SOURCE_DIR}/vhone)
endif()