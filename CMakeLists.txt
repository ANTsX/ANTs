cmake_minimum_required(VERSION 2.8.5)
cmake_policy(VERSION 2.8.5)

set(LOCAL_PROJECT_NAME ANTS)

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

#-----------------------------------------------------------------------------
# Superbuild script
#-----------------------------------------------------------------------------
if(${LOCAL_PROJECT_NAME}_SUPERBUILD)
  project(SuperBuild_${LOCAL_PROJECT_NAME})
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
else()
  project(${LOCAL_PROJECT_NAME})
  include("${CMAKE_CURRENT_SOURCE_DIR}/${LOCAL_PROJECT_NAME}.cmake")
  return()
endif()

message(FATAL_ERROR "You should never reach this point !")
