project(VH1 C Fortran)

message(STATUS "Making VH1...")


# build linked files for VH1



function(add_fortran_library LIB)
    add_library(${LIB} ${ARGN})

    # set module path to LIB_DIR/mod
    get_target_property(LIB_DIR ${LIB} BINARY_DIR)
    set_target_properties(${LIB} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_DIR}/)

    # making LIB_DIR/mod available for libraries linking LIB 
    target_include_directories(${LIB} INTERFACE ${LIB_DIR}/)
endfunction(add_fortran_library)


add_fortran_library(global "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/global.f90")
add_fortran_library(zonemod "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/zonemod.f90")

add_fortran_library(dtcon OBJECT "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/dtcon.f90")
add_fortran_library(dump OBJECT "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/dump.f90")
add_fortran_library(init OBJECT "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/init.f90")
add_fortran_library(prin OBJECT "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/prin.f90")
add_fortran_library(sweepx OBJECT "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/sweepx.f90")
add_fortran_library(sweepy OBJECT "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/sweepy.f90")
add_fortran_library(sweepz OBJECT "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/sweepz.f90")
add_fortran_library(images OBJECT "${CMAKE_SOURCE_DIR}/VH1/src/Patch/f2py/images.f90")

add_fortran_library(ppmlr OBJECT "${CMAKE_SOURCE_DIR}/VH1/src/Patch/PPMLR/ppmlr.f90")
add_fortran_library(forces OBJECT "${CMAKE_SOURCE_DIR}/VH1/src/Patch/PPMLR/forces.f90")
add_fortran_library(flatten OBJECT "${CMAKE_SOURCE_DIR}/VH1/src/Patch/PPMLR/flatten.f90")
add_fortran_library(evolve OBJECT "${CMAKE_SOURCE_DIR}/VH1/src/Patch/PPMLR/evolve.f90")
add_fortran_library(remap OBJECT "${CMAKE_SOURCE_DIR}/VH1/src/Patch/PPMLR/remap.f90")
add_fortran_library(states OBJECT "${CMAKE_SOURCE_DIR}/VH1/src/Patch/PPMLR/states.f90")
add_fortran_library(boundary OBJECT "${CMAKE_SOURCE_DIR}/VH1/src/Patch/PPMLR/boundary.f90")
add_fortran_library(volume OBJECT "${CMAKE_SOURCE_DIR}/VH1/src/Patch/PPMLR/volume.f90")
add_fortran_library(riemann OBJECT "${CMAKE_SOURCE_DIR}/VH1/src/Patch/PPMLR/riemann.f90")
add_fortran_library(parabola OBJECT "${CMAKE_SOURCE_DIR}/VH1/src/Patch/PPMLR/parabola.f90")


# Build VH1
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

#target_link_libraries(${generated_module_file} global zonemod dtcon dump init prin sweepx sweepy sweepz images)
#target_link_libraries(${generated_module_file} ppmlr forces flatten evolve remap states boundary volume riemann parabola)
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