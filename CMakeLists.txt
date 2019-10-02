cmake_minimum_required(VERSION 3.10.2...3.14.2)
cmake_policy(VERSION 3.10.2...3.14.2)

if(CMAKE_CXX_STANDARD EQUAL "98" )
   message(FATAL_ERROR "CMAKE_CXX_STANDARD:STRING=98 is not supported in ITK version 5 and greater.")
endif()

#####
##  Set the default target properties for ITK
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 11) # Supported values are ``11``, ``14``, and ``17``.
endif()
if(NOT CMAKE_CXX_STANDARD_REQUIRED)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()
if(NOT CMAKE_CXX_EXTENSIONS)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

foreach(p
    ## Only policies introduced after the cmake_minimum_required
    ## version need to explicitly be set to NEW.
    CMP0070 #3.10.0 Define ``file(GENERATE)`` behavior for relative paths.
    CMP0071 #3.10.0 Let ``AUTOMOC`` and ``AUTOUIC`` process ``GENERATED`` files.
    CMP0074 #3.12.0 `find_package()`` uses ``<PackageName>_ROOT`` variables.
    )
  if(POLICY ${p})
    cmake_policy(SET ${p} NEW)
  endif()
endforeach()

set(LOCAL_PROJECT_NAME ANTS)
include(InstallRequiredSystemLibraries)
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
SET(CMAKE_INSTALL_PREFIX /opt/ANTs CACHE PATH "Default ANTs install path")
SET(SUPERBUILD_STAGING_PREFIX ${CMAKE_BINARY_DIR}/external)

if(${LOCAL_PROJECT_NAME}_SUPERBUILD)
  project(SuperBuild_${LOCAL_PROJECT_NAME})
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  include(CPack)   # This must always be last!
  return()
else()
  project(${LOCAL_PROJECT_NAME}
    DESCRIPTION "Advanced Normalization Tools (ANTs)"
    LANGUAGES CXX C)
  include("${CMAKE_CURRENT_SOURCE_DIR}/${LOCAL_PROJECT_NAME}.cmake")
  include(CPack)   # This must always be last!
  return()
endif()

message(FATAL_ERROR "You should never reach this point !")
