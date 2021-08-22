cmake_minimum_required(VERSION 3.14.0)
project(statislib
        LANGUAGES Fortran
        VERSION 0
        DESCRIPTION "A Statistical Distributions Library for Fortran"
)

function(fypp_f90 srcfiles tgtfiles)
  set(_tgtfiles)
  foreach(srcfile IN LISTS srcfiles)
    string(REGEX REPLACE "\.fypp$" ".f90" tgtfile ${srcfile})
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${tgtfile}
      COMMAND ${FYPP} ${CMAKE_CURRENT_SOURCE_DIR}/${srcfile} ${CMAKE_CURRENT_BINARY_DIR}/${tgtfile}
      MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/${srcfile}
      )
    list(APPEND _tgtfiles ${CMAKE_CURRENT_BINARY_DIR}/${tgtfile})
  endforeach()
  set(${tgtfiles} ${_tgtfiles} PARENT_SCOPE)
endfunction()


# --- compiler options
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 9.0)
    message(FATAL_ERROR "GCC Version 9 or newer required")
  endif()
  add_compile_options(-fimplicit-none)
  add_compile_options(-ffree-line-length-132)
  add_compile_options(-Wall)
  add_compile_options(-Wextra)
  add_compile_options(-Wimplicit-procedure)
  add_compile_options(-O3)
endif()


# --- find preprocessor
find_program(FYPP fypp)
if(NOT FYPP)
  message(FATAL_ERROR "Preprocessor fypp not found!")
endif()

add_subdirectory(src)
