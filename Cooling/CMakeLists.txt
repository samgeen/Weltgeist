project(weltgeist C Fortran)

message(STATUS "Making Cooling...")

# Build Cooling
set(f2py_module_name "cooling_module")
set(fortran_src_file "${CMAKE_SOURCE_DIR}/Cooling/cooling_module_frig.f90")
set(generated_module_file ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})

add_custom_target(${f2py_module_name} ALL
  DEPENDS ${generated_module_file}
  )

add_custom_command(
OUTPUT ${generated_module_file}
COMMAND ${F2PY_EXECUTABLE}
  -m ${f2py_module_name}
  -c
  ${fortran_src_file}
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

install(FILES ${generated_module_file} DESTINATION weltgeist)
