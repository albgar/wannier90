@PACKAGE_INIT@

if(NOT TARGET "@PROJECT_NAME@::@PROJECT_NAME@")
  include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")

  include(CMakeFindDependencyMacro)

endif()
