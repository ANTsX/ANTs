#
# ONLY MODIFY TO CHANGE VERSION
#
# The number of commits since last this file has changed it used to
# define "dev" and "post", modification of this file will reset that
# version.
#

# Version info
include(ProjectSourceVersion)

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

# The hash is the current git sha1 hash tag of the HEAD.
set(${PROJECT_NAME}_VERSION_HASH "${_GIT_VERSION_HASH}")


# DEV or POST is set to the number of commits since this file has been
# changed. If the MAJOR.MINOR.[PATCH[.TWEAK]] matches "closest"
# version tag then its consider in the release branch and POST is set
# while DEV is undefined, otherwise we are considered under
# development and DEV is set and POST is undefined.
if(DEFINED _GIT_VERSION_POST)
  set(${PROJECT_NAME}_VERSION_POST "${_GIT_VERSION_POST}")
elseif(DEFINED _GIT_VERSION_DEV)
  set(${PROJECT_NAME}_VERSION_DEV "${_GIT_VERSION_DEV}")
endif()
