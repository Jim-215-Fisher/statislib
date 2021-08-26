### Pre-process: .fypp -> .f90 via Fypp

# Create a list of the files to be preprocessed
set(fyppFiles
    statislib.fypp
    common.fypp
    special_function.fypp
    uniform_distribution.fypp
    normal_distribution.fypp
    exponential_distribution.fypp
    gamma_distribution.fypp
    beta_distribution.fypp
    binomial_distribution.fypp
    geometric_distribution.fypp
    hypergeometric_distribution.fypp
    logarithmic_distribution.fypp
    negbinomial_distribution.fypp
    poisson_distribution.fypp
)



fypp_f90("${fyppFiles}" outFiles)

set(SRC
    ${outFiles}
)

find_package(fortran_stdlib REQUIRED)

add_library(${PROJECT_NAME} ${SRC})

set(LIB_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files/)
set_target_properties(${PROJECT_NAME} PROPERTIES
    Fortran_MODULE_DIRECTORY ${LIB_MOD_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${LIB_MOD_DIR}>
)

target_link_libraries(
  ${PROJECT_NAME}
  PRIVATE
  fortran_stdlib::fortran_stdlib
)
