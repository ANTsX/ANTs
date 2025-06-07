#
# This CMake code extracts the information from the git repository,
# and automatically causes a reconfigure if the git HEAD changes. The
# following variable may be defined after execution:
#
# _GIT_VERSION_HASH - the SHA1 hash of the current HEAD
#
# Based on the most recent tag starting with the letter "v" for
# version, which is expected to be of the form
# vN.N[.N[.N][(a|b|c|rc[N])] the following is extracted or undefined:
#
# _GIT_VERSION_MAJOR
# _GIT_VERSION_MINOR
# _GIT_VERSION_PATCH
# _GIT_VERSION_TWEAK
# _GIT_VERSION_RC
#
# If the current project's version ( defiend by
# ${CMAKE_PROJECT_NAME}_VERSION_MAJOR and MINOR and PATCH and TWEAK
# match that of the tag, then it'll be considered that the project is
# in post release mode otherwise it's considered underdevelopment.
#
# Only one of the following variables will be defined.
# _GIT_VERSION_DEV is defined as number of commits
# since the projects Version.cmake file has been modified. While
# _GIT_VERSION_POST is defined as the number of commits since the tag.
#


include(GetGitRevisionDescription)

get_git_head_revision(GIT_REFVAR _GIT_VERSION_HASH)

# if there is not git directory we should be in a distributed package
# we will use version provided in Version.cmake
if(_GIT_VERSION_HASH STREQUAL "GITDIR-NOTFOUND")
  return()
endif()

if(_GIT_VERSION_HASH MATCHES "^[a-fA-F0-9]+$")
  # Get first seven chars of hash (git default for short hash)
  # https://git-scm.com/book/en/v2/Git-Tools-Revision-Selection
  string(SUBSTRING "${_GIT_VERSION_HASH}" 0 7 _GIT_VERSION_HASH)
endif()

# find the closest anotated tag with the v prefix for version
git_describe(_GIT_TAG "--match=v*")

git_commits_since("${PROJECT_SOURCE_DIR}/Version.cmake" _GIT_VERSION_COUNT)

set(VERSION_REGEX "^v([0-9]+)\\.([0-9]+)+(\\.([0-9]+))?(\\.([0-9]+))?((a|b|c|rc)[0-9]*)?(-[0-9]+)?")

string(REGEX MATCH "${VERSION_REGEX}" _out "${_GIT_TAG}")

if("${_out}" STREQUAL "")
  message(NOTICE "Could not parse most recent version from git tag: \"${_GIT_TAG}\"")
  set(_GIT_VERSION "GITVERSION-NOTFOUND")
  return()
endif()

set(_GIT_VERSION_MAJOR "${CMAKE_MATCH_1}")
set(_GIT_VERSION_MINOR "${CMAKE_MATCH_2}")
if(NOT "${CMAKE_MATCH_4}" STREQUAL "")
  set(_GIT_VERSION_PATCH "${CMAKE_MATCH_4}")
endif()
if(NOT "${CMAKE_MATCH_6}" STREQUAL "")
  set(_GIT_VERSION_TWEAK "${CMAKE_MATCH_6}")
endif()
if(NOT "${CMAKE_MATCH_7}" STREQUAL "")
  set(_GIT_VERSION_RC "${CMAKE_MATCH_7}" ) # a,b,rc01 etc
endif()

if(NOT "${CMAKE_MATCH_9}" STREQUAL "")
  #trim leading '-'
  string(SUBSTRING "${CMAKE_MATCH_9}" 1 -1 CMAKE_MATCH_9)

  set(_GIT_TAG_COUNT "${CMAKE_MATCH_9}")
endif()


set(_GIT_VERSION "${_GIT_VERSION_MAJOR}.${_GIT_VERSION_MINOR}")
if(DEFINED _GIT_VERSION_PATCH)
  set(_GIT_VERSION "${_GIT_VERSION}.${_GIT_VERSION_PATCH}")
  if(DEFINED _GIT_VERSION_TWEAK)
    set(_GIT_VERSION "${_GIT_VERSION}.${_GIT_VERSION_TWEAK}")
  elseif(DEFINED ${CMAKE_PROJECT_NAME}_VERSION_TWEAK)
    set(_GIT_VERSION "${_GIT_VERSION}.0")
  endif()
elseif(DEFINED ${CMAKE_PROJECT_NAME}_VERSION_PATCH)
  set(_GIT_VERSION "${_GIT_VERSION}.0")
  if(DEFINED ${CMAKE_PROJECT_NAME}_VERSION_TWEAK)
    set(_GIT_VERSION "${_GIT_VERSION}.0")
  endif()
endif()

