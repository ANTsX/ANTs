
# Make sure this file is included only once
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

# Sanity checks
if(DEFINED KWStyle_EXE AND NOT EXISTS ${KWStyle_EXE})
  message(FATAL_ERROR "KWStyle_EXE variable is defined but corresponds to non-existing file")
endif()

# Set dependency list
set(KWStyle_DEPENDENCIES "")

# Include dependent projects if any
SlicerMacroCheckExternalProjectDependency(KWStyle)
set(proj KWStyle)

if(NOT DEFINED KWStyle_EXE AND NOT ${USE_SYSTEM_KWStyle})

  # Set CMake OSX variable to pass down the external project
  set(CMAKE_OSX_EXTERNAL_PROJECT_ARGS)
  if(APPLE)
    list(APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
      -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
      -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
      -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET})
  endif()

  ExternalProject_add(${proj}
    CVS_REPOSITORY :pserver:anoncvs@public.kitware.com:/cvsroot/KWStyle
    CVS_MODULE KWStyle
    SOURCE_DIR ${proj}
    BINARY_DIR ${proj}-build
    LOG_CONFIGURE 0  # Wrap configure in script to ignore log output from dashboards
    LOG_BUILD     0  # Wrap build in script to to ignore log output from dashboards
    LOG_TEST      0  # Wrap test in script to to ignore log output from dashboards
    LOG_INSTALL   0  # Wrap install in script to to ignore log output from dashboards
    ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
    CMAKE_ARGS
      -Wno-dev
      --no-warn-unused-cli
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/Utils
      -DCMAKE_BUILD_TYPE:STRING=Release
    DEPENDS
      ${KWStyle_DEPENDENCIES}
    )

  set(KWStyle_EXE ${CMAKE_BINARY_DIR}/Utils/bin/KWStyle)
else()
  if(${USE_SYSTEM_KWStyle})
    find_program(KWStyle_EXE KWStyle DOC "Path of KWStyle program")
    if(NOT KWStyle_EXE)
      message(FATAL_ERROR "To use the system KWStyle, set KWStyle_EXE")
    endif()
  endif()
  # The project is provided using KWStyle_EXE, nevertheless since other
  # project may depend on KWStyle, let's add an 'empty' one
  SlicerMacroEmptyExternalProject(${proj} "${KWStyle_DEPENDENCIES}")
endif()

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS KWStyle_EXE:FILEPATH)

