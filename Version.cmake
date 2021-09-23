# Version info
include(ProjectSourceVersion)

set(_GIT_VERSION "${_GIT_VERSION_MAJOR}.${_GIT_VERSION_MINOR}")
if(DEFINED _GIT_VERSION_PATCH)
  set(_GIT_VERSION "${_GIT_VERSION}.${_GIT_VERSION_PATCH}")
  if(DEFINED _GIT_VERSION_TWEAK)
    set(_GIT_VERSION "${_GIT_VERSION}.${_GIT_VERSION_TWEAK}")
  endif()
endif()

# pre-release codes are defined based on suffix of most recent tags.

# a[N]+<b[N]+<c[N]+=rc[N]+
if(NOT ${_GIT_VERSION_RC} STREQUAL "" )
  set(${PROJECT_NAME}_VERSION_RC "${_GIT_VERSION_RC}")
endif()

set(${PROJECT_NAME}_VERSION "${_GIT_VERSION}")
set(${PROJECT_NAME}_VERSION_MAJOR "${_GIT_VERSION_MAJOR}")
set(${PROJECT_NAME}_VERSION_MINOR "${_GIT_VERSION_MINOR}")
set(${PROJECT_NAME}_VERSION_PATCH "${_GIT_VERSION_PATCH}")
set(${PROJECT_NAME}_VERSION_TWEAK "${_GIT_VERSION_TWEAK}")
set(${PROJECT_NAME}_VERSION_HASH "${_GIT_VERSION_HASH}")

if(DEFINED _GIT_VERSION_POST)
  set(${PROJECT_NAME}_VERSION_POST "${_GIT_VERSION_POST}")
endif()

if(DEFINED ${PROJECT_NAME}_VERSION_RC)
  set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION}${${PROJECT_NAME}_VERSION_RC}")
endif()

if(DEFINED ${PROJECT_NAME}_VERSION_POST)
  set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION}.post${${PROJECT_NAME}_VERSION_POST}")
endif()

if( (NOT ANTS_BUILD_DISTRIBUTE) AND (NOT ${PROJECT_NAME}_VERSION_HASH STREQUAL "GITDIR-NOTFOUND") )
  set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION}-g${${PROJECT_NAME}_VERSION_HASH}")
endif()

# If no version information exists and the user has not passed version info to cmake, set defaults
# This should only happen for development snapshots, releases should have a cmake file describing
# the version
if("${${PROJECT_NAME}_VERSION_MAJOR}" STREQUAL "")
  set(${PROJECT_NAME}_VERSION_MAJOR 0)
  set(${PROJECT_NAME}_VERSION_MINOR 0)
  set(${PROJECT_NAME}_VERSION_PATCH 0)
  set(${PROJECT_NAME}_VERSION_TWEAK 0)
  set(${PROJECT_NAME}_VERSION "0.0.0.0")
  message(status "No version information found - using placeholder ${PROJECT_NAME}_VERSION")
endif()