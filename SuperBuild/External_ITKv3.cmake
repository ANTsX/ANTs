
# Make sure this file is included only once
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

#if(${USE_SYSTEM_ITK})
#  unset(ITK_DIR CACHE)
#endif()

# Sanity checks
if(DEFINED ITK_DIR AND NOT EXISTS ${ITK_DIR})
  message(FATAL_ERROR "ITK_DIR variable is defined but corresponds to non-existing directory")
endif()

# Set dependency list
set(ITKv3_DEPENDENCIES "")

# Include dependent projects if any
SlicerMacroCheckExternalProjectDependency(ITKv3)
set(proj ITKv3)

if(NOT DEFINED ITK_DIR AND NOT ${USE_SYSTEM_ITK})
  #message(STATUS "${__indent}Adding project ${proj}")

  # Set CMake OSX variable to pass down the external project
  set(CMAKE_OSX_EXTERNAL_PROJECT_ARGS)
  if(APPLE)
    list(APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
      -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
      -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
      -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET})
  endif()

  set(ITKv3_REPOSITORY git://itk.org/ITK.git)
  set(ITKv3_GIT_TAG v3.20.1)

  ExternalProject_Add(${proj}
    GIT_REPOSITORY ${ITKv3_REPOSITORY}
    GIT_TAG ${ITKv3_GIT_TAG}
    UPDATE_COMMAND ""
    SOURCE_DIR ${proj}
    BINARY_DIR ${proj}-build
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
      ${COMMON_EXTERNAL_PROJECT_ARGS}
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      -DITK_LEGACY_REMOVE:BOOL=ON
      -DKWSYS_USE_MD5:BOOL=ON # Required by SlicerExecutionModel
      -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
      -DBUILD_SHARED_LIBS:BOOL=ON
      -DITK_INSTALL_LIB_DIR:PATH=${Slicer_INSTALL_LIB_DIR}
      -DITK_USE_REVIEW:BOOL=ON
      -DITK_USE_REVIEW_STATISTICS:BOOL=ON
      -DITK_USE_OPTIMIZED_REGISTRATION_METHODS:BOOL=ON
      -DITK_USE_CENTERED_PIXEL_COORDINATES_CONSISTENTLY:BOOL=ON
      -DITK_USE_TRANSFORM_IO_FACTORIES:BOOL=ON
      -DITK_USE_PORTABLE_ROUND:BOOL=ON # Unused
      ##HACK: THE NEXT TWO LINES ARE TO AVOID CMAKE ISSUES DURING BUILD
      -DLIBRARY_OUTPUT_PATH:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
      -DEXECUTABLE_OUTPUT_PATH:PATH=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    INSTALL_COMMAND ""
    DEPENDS
      ${ITKv3_DEPENDENCIES}
    )
  set(ITK_DIR ${CMAKE_BINARY_DIR}/${proj}-build)
else()
  if(${USE_SYSTEM_ITK})
    find_package(ITK ${ITK_VERSION_MAJOR} REQUIRED)
    if(NOT ITK_DIR)
      message(FATAL_ERROR "To use the system ITK, set ITK_DIR")
    endif()
  endif()
  # The project is provided using ITK_DIR, nevertheless since other
  # project may depend on ITKv3, let's add an 'empty' one
  SlicerMacroEmptyExternalProject(${proj} "${ITKv3_DEPENDENCIES}")
endif()

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS ITK_DIR:PATH)
