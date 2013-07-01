
# Make sure this file is included only once
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

# Sanity checks
if(DEFINED Uncrustify_EXE AND NOT EXISTS ${Uncrustify_EXE})
  message(FATAL_ERROR "Uncrustify_EXE variable is defined but corresponds to non-existing file")
endif()

# Set dependency list
set(Uncrustify_DEPENDENCIES "")

# Include dependent projects if any
SlicerMacroCheckExternalProjectDependency(Uncrustify)
set(proj Uncrustify)

if(NOT DEFINED Uncrustify_EXE AND NOT ${USE_SYSTEM_Uncrustify})

  set(Uncrustify_REPOSITORY git://uncrustify.git.sourceforge.net/gitroot/uncrustify/uncrustify)
  set(Uncrustify_GIT_TAG 60f3681da60462eda539b78e0c6c3eea823481e5)

  ExternalProject_add(${proj}
    GIT_REPOSITORY ${Uncrustify_REPOSITORY}
    GIT_TAG ${Uncrustify_GIT_TAG}
    LOG_CONFIGURE 0  # Wrap configure in script to ignore log output from dashboards
    LOG_BUILD     0  # Wrap build in script to to ignore log output from dashboards
    LOG_TEST      0  # Wrap test in script to to ignore log output from dashboards
    LOG_INSTALL   0  # Wrap install in script to to ignore log output from dashboards
    ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
    SOURCE_DIR ${proj}
    BINARY_DIR ${proj}-build
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${CMAKE_BINARY_DIR}/Utils
    DEPENDS
      ${Uncrustify_DEPENDENCIES}
    )

  set(Uncrustify_EXE ${CMAKE_BINARY_DIR}/Utils/bin/uncrustify)

else()
  if(${USE_SYSTEM_Uncrustify})
    find_program(Uncrustify_EXE uncrustify DOC "Path of uncrustify program")
    if(NOT Uncrustify_EXE)
      message(FATAL_ERROR "To use the system uncrustify, set Uncrustify_EXE")
    endif()
    message("USING the system ${extProjName}, set ${extProjName}_DIR=${${extProjName}_DIR}")
  endif()

  # The project is provided using Uncrustify_EXE, nevertheless since other
  # project may depend on Uncrustify, let's add an 'empty' one
  SlicerMacroEmptyExternalProject(${proj} "${Uncrustify_DEPENDENCIES}")

endif()

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS Uncrustify_EXE:FILEPATH)

