#
# ONLY MODIFY TO CHANGE VERSION
#
# The number of commits since last this file has changed it used to
# define "dev" and "post", modification of this file will reset that
# version.
#

# Version info from git
include(ProjectSourceVersion)

set(${PROJECT_NAME}_VERSION_MAJOR "2")
set(${PROJECT_NAME}_VERSION_MINOR "5")
set(${PROJECT_NAME}_VERSION_PATCH "4")
# set(${PROJECT_NAME}_VERSION_TWEAK "")

# The hash is the current git sha1 hash tag of the HEAD.
set(${PROJECT_NAME}_VERSION_HASH "${_GIT_VERSION_HASH}")

if ( NOT ${PROJECT_NAME}_VERSION_HASH STREQUAL "GITDIR-NOTFOUND" )
  # Use version information from git if the source directory is a git repository
  # Compare git version (major.minor[.patch[.tweak]]) with version number hard-coded above
  set(_${PROJECT_NAME}_VERSION_NUMBER "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}")

  if(DEFINED ${PROJECT_NAME}_VERSION_PATCH)
    set(_${PROJECT_NAME}_VERSION_NUMBER "${_${PROJECT_NAME}_VERSION_NUMBER}.${${PROJECT_NAME}_VERSION_PATCH}")
    if(DEFINED ${PROJECT_NAME}_VERSION_TWEAK)
      set(_${PROJECT_NAME}_VERSION_NUMBER "${_${PROJECT_NAME}_VERSION_NUMBER}.${${PROJECT_NAME}_VERSION_TWEAK}")
    endif()
  endif()

  if(_GIT_VERSION VERSION_EQUAL _${PROJECT_NAME}_VERSION_NUMBER)
    if(_GIT_TAG_COUNT) #ignore if 0
      set(_GIT_VERSION_POST "${_GIT_TAG_COUNT}")
    endif()
  else()
    # The first commit after a tag should increase the project version
    # number in Version.cmake and be "dev1"
    math(EXPR _GIT_VERSION_COUNT "${_GIT_VERSION_COUNT}+1")
    set(_GIT_VERSION_DEV "${_GIT_VERSION_COUNT}")
  endif()

  # DEV or POST is set to the number of commits since this file has been
  # changed. If the MAJOR.MINOR.[PATCH[.TWEAK]] matches "closest"
  # version tag then its considered in the release branch and POST is set
  # while DEV is undefined, otherwise we are considered under
  # development and DEV is set and POST is undefined.
  #
  if(DEFINED _GIT_VERSION_POST)
    set(${PROJECT_NAME}_VERSION_POST "${_GIT_VERSION_POST}")
  elseif(DEFINED _GIT_VERSION_DEV)
    set(${PROJECT_NAME}_VERSION_DEV "${_GIT_VERSION_DEV}")
  endif()

  # pre-release codes are defined based on suffix of most recent tags.
  # a[N]+<b[N]+<c[N]+=rc[N]+
  if(NOT ${_GIT_VERSION_RC} STREQUAL "" )
    set(${PROJECT_NAME}_VERSION_RC "${_GIT_VERSION_RC}")
  endif()

  if(DEFINED ${PROJECT_NAME}_VERSION_RC)
    set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION}${${PROJECT_NAME}_VERSION_RC}")
  endif()

  set(${PROJECT_NAME}_VERSION "${_${PROJECT_NAME}_VERSION_NUMBER}")

  if(DEFINED ${PROJECT_NAME}_VERSION_POST)
    set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION}.post${${PROJECT_NAME}_VERSION_POST}")
  elseif(DEFINED ${PROJECT_NAME}_VERSION_DEV)
    set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION}.dev${${PROJECT_NAME}_VERSION_DEV}")
  endif()

  if( NOT ${PROJECT_NAME}_BUILD_DISTRIBUTE AND NOT ${PROJECT_NAME}_VERSION_HASH STREQUAL "GITDIR-NOTFOUND" )
    set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION}-g${${PROJECT_NAME}_VERSION_HASH}")
  endif()
else()
  # No git version information
  # Set version information from numbers only
  set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}")

  # Set to 1 for release commits only
  set(${PROJECT_NAME}_RELEASE_VERSION 0)

  if( NOT ${PROJECT_NAME}_RELEASE_VERSION )
    set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION}.dev")
  endif()

endif()

