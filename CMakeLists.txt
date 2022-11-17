
cmake_minimum_required(VERSION 3.14)

 get_filename_component(srcdir "${CMAKE_CURRENT_SOURCE_DIR}" REALPATH)
 get_filename_component(bindir "${CMAKE_CURRENT_BINARY_DIR}" REALPATH)

  if("${srcdir}" STREQUAL "${bindir}")
    message(FATAL_ERROR
      "It is not allowed to configure and build this project from its source folder. Please, create a \
separate build directory and invoke CMake from that directory.")
  endif()

project(
  "libwannier90"
  LANGUAGES "Fortran"
  VERSION "3.1.0"
  DESCRIPTION "Wannier90 Library for use in Siesta"
)

# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# General configuration information
add_subdirectory("config")

option(WITH_MPI "Whether libwannier90 should support MPI-parallelism" FALSE)

if(WITH_MPI)
  if (NOT TARGET MPI::MPI_Fortran)
    find_package(MPI REQUIRED)
  endif()
endif()

# Collect source of the project
add_library( "${PROJECT_NAME}-lib"

 src/w90_mpi.F90
 src/constants.F90
 src/io.F90
 src/utility.F90
 src/comms.F90
 src/parameters.F90
 src/ws_distance.F90
 src/hamiltonian.F90
 src/sitesym.F90
 src/overlap.F90
 src/kmesh.F90
 src/disentangle.F90
 src/wannierise.F90
 src/plot.F90
 src/transport.F90
 src/wannier_lib.F90

)

set_target_properties(
  "${PROJECT_NAME}-lib"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME "wannier90"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)

if(WITH_MPI)
  target_compile_definitions(
  "${PROJECT_NAME}-lib"
  PRIVATE MPI
  )

  target_link_libraries("${PROJECT_NAME}-lib" PUBLIC MPI::MPI_Fortran)

endif()

target_include_directories(
  "${PROJECT_NAME}-lib"
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)

# Add example application
#add_subdirectory("app")

# Export targets for other projects
add_library("${PROJECT_NAME}" INTERFACE)
target_link_libraries("${PROJECT_NAME}" INTERFACE "${PROJECT_NAME}-lib")
install(
  TARGETS
  "${PROJECT_NAME}"
  "${PROJECT_NAME}-lib"
  EXPORT
  "${PROJECT_NAME}-targets"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
install(
  EXPORT
  "${PROJECT_NAME}-targets"
  NAMESPACE
  "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
install(
  DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}"
)
# Package license files
install(
  FILES
  "LICENSE"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
)

#add_subdirectory("docs")

# add the testsuite

#enable_testing()
#add_subdirectory("Testers")
