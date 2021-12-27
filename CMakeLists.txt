set(ANTs_OLDEST_VALIDATED_POLICIES_VERSION "3.16.3")
set(ANTs_NEWEST_VALIDATED_POLICIES_VERSION "3.19.7")
cmake_minimum_required(VERSION ${ANTs_OLDEST_VALIDATED_POLICIES_VERSION}..${ANTs_NEWEST_VALIDATED_POLICIES_VERSION} FATAL_ERROR)

#####
##  Set the default target properties for ITK
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14) # Supported values are ``14``, ``17``, and ``20``.
endif()
# ==== Define language standard configurations requiring at least c++14 standard
if(CMAKE_CXX_STANDARD EQUAL "98" OR CMAKE_CXX_STANDARD LESS "14")
   message(FATAL_ERROR "C++98 to C++11 are no longer supported in ITK version 5.3 and greater.")
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
  project(SuperBuild_${LOCAL_PROJECT_NAME}
    DESCRIPTION "Optional superbuild step to build external support libraries for ANTs"
  )
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  include(CPack)   # This must always be last!
  return()
else()
  project(${LOCAL_PROJECT_NAME}
    DESCRIPTION "Advanced Normalization Tools (ANTs) computes high-dimensional mappings to capture the statistics of brain structure and function."
    HOMEPAGE_URL "https://github.com/ANTsX/ANTs"
    LANGUAGES CXX C)
  include("${CMAKE_CURRENT_SOURCE_DIR}/${LOCAL_PROJECT_NAME}.cmake")
  include(CPack)   # This must always be last!
  return()
endif()

message(FATAL_ERROR "You should never reach this point !")
