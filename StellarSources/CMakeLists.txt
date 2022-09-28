project(weltgeist C Fortran)

message(STATUS "Making Stellar Sources...")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_library(table_1d_module OBJECT "${CMAKE_SOURCE_DIR}/StellarSources/Fortran/table_1d_module.f90")
add_library(singlestar_module OBJECT "${CMAKE_SOURCE_DIR}/StellarSources/Fortran/singlestar_module.f90")
add_dependencies(singlestar_module table_1d_module)

set(f2py_module_name "singlestar")
set(fortran_src_file "${CMAKE_SOURCE_DIR}/StellarSources/Fortran/singlestar_f2py.f90")
set(generated_module_file ${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})

add_custom_target(${f2py_module_name} ALL
  DEPENDS "${fortran_src_file}"
  )

add_custom_command(
  OUTPUT "${f2py_module_name}module.c" "${f2py_module_name}-f2pywrappers.f"
  COMMAND ${F2PY_EXECUTABLE}
    -m ${f2py_module_name}
    ${fortran_src_file}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${fortran_src_file}
  )

add_library(${generated_module_file} MODULE
            "${f2py_module_name}module.c"
            "${F2PY_INCLUDE_DIR}/fortranobject.c"
            "${f2py_module_name}-f2pywrappers.f"
             ${fortran_src_file})

target_link_libraries(${generated_module_file} singlestar_module table_1d_module)
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
  install(TARGETS ${generated_module_file} DESTINATION weltgeist)
else()
  install(TARGETS ${generated_module_file} DESTINATION ${CMAKE_SOURCE_DIR}/weltgeist)
endif()