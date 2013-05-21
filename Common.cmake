
include(CMakeDependentOption)

#-----------------------------------------------------------------------------
# Build option(s)
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
option(BUILD_SHARED_LIBS "Build ITK with shared libraries." OFF)
set(ANTS_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})

option(BUILD_EXTERNAL_APPLICATIONS "Build applications that are bundled with ANTs for convenience." ON)
set(BUILD_EXTERNAL_APPLICATIONS ${BUILD_EXTERNAL_APPLICATIONS})

set(ITK_VERSION_MAJOR 4 CACHE STRING "Choose the expected ITK major version to build ANTS (only ITKv4 supported).")
# Set the possible values of ITK major version for cmake-gui
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "4")
if(NOT ${ITK_VERSION_MAJOR} STREQUAL "4")
  message(FATAL_ERROR "ITK_VERSION_MAJOR should be 4")
endif()

set(USE_ITKv4 ON)

# With CMake 2.8.9 or later, the UPDATE_COMMAND is required for updates to occur.
# For earlier versions, we nullify the update state to prevent updates and
# undesirable rebuild.
if(CMAKE_VERSION VERSION_LESS 2.8.9)
  set(cmakeversion_external_update UPDATE_COMMAND "")
else()
  set(cmakeversion_external_update LOG_UPDATE 1)
endif()

#-----------------------------------------------------------------------------
# Update CMake module path
#------------------------------------------------------------------------------
set(CMAKE_MODULE_PATH
  ${${PROJECT_NAME}_SOURCE_DIR}/CMake
  ${${PROJECT_NAME}_BINARY_DIR}/CMake
  ${CMAKE_MODULE_PATH}
  )

#-----------------------------------------------------------------------------
# Sanity checks
#------------------------------------------------------------------------------
include(PreventInSourceBuilds)
include(PreventInBuildInstalls)

#-----------------------------------------------------------------------------
# CMake Function(s) and Macro(s)
#-----------------------------------------------------------------------------
if(CMAKE_PATCH_VERSION LESS 3)
  include(Pre283CMakeParseArguments)
else()
  include(CMakeParseArguments)
endif()

#-----------------------------------------------------------------------------
# Platform check
#-----------------------------------------------------------------------------
set(PLATFORM_CHECK true)

if(PLATFORM_CHECK)
  # See CMake/Modules/Platform/Darwin.cmake)
  #   6.x == Mac OSX 10.2 (Jaguar)
  #   7.x == Mac OSX 10.3 (Panther)
  #   8.x == Mac OSX 10.4 (Tiger)
  #   9.x == Mac OSX 10.5 (Leopard)
  #  10.x == Mac OSX 10.6 (Snow Leopard)
  if (DARWIN_MAJOR_VERSION LESS "9")
    message(FATAL_ERROR "Only Mac OSX >= 10.5 are supported !")
  endif()
endif()

#-----------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#-----------------------------------------------------------------------------
if(NOT COMMAND SETIFEMPTY)
  macro(SETIFEMPTY)
    set(KEY ${ARGV0})
    set(VALUE ${ARGV1})
    if(NOT ${KEY})
      set(${ARGV})
    endif()
  endmacro()
endif()

#-----------------------------------------------------------------------------
SETIFEMPTY(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

#-----------------------------------------------------------------------------
SETIFEMPTY(CMAKE_INSTALL_LIBRARY_DESTINATION lib)
SETIFEMPTY(CMAKE_INSTALL_ARCHIVE_DESTINATION lib)
SETIFEMPTY(CMAKE_INSTALL_RUNTIME_DESTINATION bin)

#-------------------------------------------------------------------------
SETIFEMPTY(BRAINSTools_CLI_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
SETIFEMPTY(BRAINSTools_CLI_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
SETIFEMPTY(BRAINSTools_CLI_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

#-------------------------------------------------------------------------
SETIFEMPTY(BRAINSTools_CLI_INSTALL_LIBRARY_DESTINATION ${CMAKE_INSTALL_LIBRARY_DESTINATION})
SETIFEMPTY(BRAINSTools_CLI_INSTALL_ARCHIVE_DESTINATION ${CMAKE_INSTALL_ARCHIVE_DESTINATION})
SETIFEMPTY(BRAINSTools_CLI_INSTALL_RUNTIME_DESTINATION ${CMAKE_INSTALL_RUNTIME_DESTINATION})

#-------------------------------------------------------------------------
# Augment compiler flags
#-------------------------------------------------------------------------
include(ITKSetStandardCompilerFlags)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_DEBUG_DESIRED_FLAGS} " )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_DEBUG_DESIRED_FLAGS} " )
else() # Release, or anything else
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_RELEASE_DESIRED_FLAGS} " )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_RELEASE_DESIRED_FLAGS} " )
endif()

#-----------------------------------------------------------------------------
# Add needed flag for gnu on linux like enviroments to build static common libs
# suitable for linking with shared object libs.
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
  if(NOT "${CMAKE_CXX_FLAGS}" MATCHES "-fPIC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  endif()
  if(NOT "${CMAKE_C_FLAGS}" MATCHES "-fPIC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  endif()
endif()
