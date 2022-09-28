project(weltgeist C Fortran)

message(STATUS "Making Radiation...")

# Build raytracing
set(f2py_module_name "raytracing")
set(fortran_src_file "${CMAKE_SOURCE_DIR}/Radiation/raytracing.f90")
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
