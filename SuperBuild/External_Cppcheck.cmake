
# Make sure this file is included only once
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

# Sanity checks
if(DEFINED Cppcheck_EXE AND NOT EXISTS ${Cppcheck_EXE})
  message(FATAL_ERROR "Cppcheck_EXE variable is defined but corresponds to non-existing file")
endif()

# Set dependency list
set(Cppcheck_DEPENDENCIES "")

# Include dependent projects if any
SlicerMacroCheckExternalProjectDependency(Cppcheck)
set(proj Cppcheck)

if(NOT DEFINED Cppcheck_EXE AND NOT ${USE_SYSTEM_Cppcheck})

  set(Cppcheck_REPOSITORY git://github.com/danmar/cppcheck.git)
  set(Cppcheck_GIT_TAG origin/master)

  ExternalProject_add(${proj}
    GIT_REPOSITORY ${Cppcheck_REPOSITORY}
    GIT_TAG ${Cppcheck_GIT_TAG}
    SOURCE_DIR ${proj}
    BUILD_IN_SOURCE 1
    ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND HAVE_RULES=no ${CMAKE_MAKE_PROGRAM}
    INSTALL_COMMAND HAVE_RULES=no DESTDIR=${CMAKE_BINARY_DIR}/ PREFIX=Utils ${CMAKE_MAKE_PROGRAM} install
    DEPENDS
      ${Cppcheck_DEPENDENCIES}
    )

  set(Cppcheck_EXE ${CMAKE_BINARY_DIR}/Utils/bin/cppcheck)

else()
  if(${USE_SYSTEM_Cppcheck})
    find_program(Cppcheck_EXE cppcheck DOC "Path of Cppcheck program")
    if(NOT Cppcheck_EXE)
      message(FATAL_ERROR "To use the system Cppcheck, set Cppcheck_EXE")
    endif()
  endif()
  # The project is provided using Cppcheck_EXE, nevertheless since other
  # project may depend on Cppcheck, let's add an 'empty' one
  SlicerMacroEmptyExternalProject(${proj} "${Cppcheck_DEPENDENCIES}")
endif()

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS Cppcheck_EXE:FILEPATH)
