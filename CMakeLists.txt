cmake_minimum_required(VERSION 2.8.9)
cmake_policy(VERSION 2.8.9)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

set(LIBRARY_SOVERSION_INFO 2)
set(LIBRARY_VERSION_INFO 2.2.0)

set(LOCAL_PROJECT_NAME ANTS)
INCLUDE(InstallRequiredSystemLibraries)
## NOTE THERE SHOULD BE NO PROJECT STATEMENT HERE!
## This file acts as a simple switch to initiate
## two completely independant CMake build environments.

#-----------------------------------------------------------------------------
# Superbuild Option - Enabled by default
#                   Phase I:  ${LOCAL_PROJECT_NAME}_SUPERBUILD is set to ON, and the
#                             supporting packages defined in "SuperBuild.cmake"
#                             are built.  The last package in "SuperBuild.cmake"
#                             to be built is a recursive call to this
#                             file with ${LOCAL_PROJECT_NAME}_SUPERBUILD explicitly
#                             set to "OFF" to initiate Phase II
#
#                   Phase II: Build the ${LOCAL_PROJECT_NAME}, referencing the support
#                             packages built in Phase I.
#-----------------------------------------------------------------------------
option(${LOCAL_PROJECT_NAME}_SUPERBUILD "Build ${LOCAL_PROJECT_NAME} and the projects it depends on via SuperBuild.cmake." ON)
mark_as_advanced(${LOCAL_PROJECT_NAME}_SUPERBUILD)

option(${LOCAL_PROJECT_NAME}_USE_QT "Find and use Qt with VTK to build GUI Tools" OFF)
mark_as_advanced(${LOCAL_PROJECT_NAME}_USE_QT)

option(${LOCAL_PROJECT_NAME}_INSTALL_DEVELOPMENT "Install development support include and libraries for external packages." OFF)
mark_as_advanced(${LOCAL_PROJECT_NAME}_INSTALL_DEVELOPMENT)

option(COPY_SCRIPT_FILES_TO_BIN_DIR "Copy the script files to the ANTS bin directory." OFF)
mark_as_advanced(COPY_SCRIPT_FILES_TO_BIN_DIR)

set(CPACK_PACKAGE_NAME "ANTs")
set(CPACK_PACKAGE_VENDOR "CMake.org")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ANTs - Advanced Normalization Tools")
set(CPACK_PACKAGE_VERSION 2.2.0)
set(CPACK_PACKAGE_VERSION_MAJOR 2)
set(CPACK_PACKAGE_VERSION_MINOR 2)
set(CPACK_PACKAGE_VERSION_PATCH 0)
set(CPACK_PACKAGE_INSTALL_DIRECTORY "ANTS_Install")
set(CPACK_BINARY_GENERATORS "DragNDrop TGZ TZ")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ANTs - robust image registration, segmentation and more")
set(CPACK_RESOURCE_FILE_LICENSE  "${CMAKE_CURRENT_SOURCE_DIR}/ANTSCopyright.txt")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
#-----------------------------------------------------------------------------
# Superbuild script
#-----------------------------------------------------------------------------
if(${LOCAL_PROJECT_NAME}_SUPERBUILD)
  project(SuperBuild_${LOCAL_PROJECT_NAME})
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  include(CPack)   # This must always be last!
  return()
else()
  project(${LOCAL_PROJECT_NAME})
  include("${CMAKE_CURRENT_SOURCE_DIR}/${LOCAL_PROJECT_NAME}.cmake")
  include(CPack)   # This must always be last!
  return()
endif()

message(FATAL_ERROR "You should never reach this point !")
