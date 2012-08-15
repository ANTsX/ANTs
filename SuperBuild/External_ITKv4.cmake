
# Make sure this file is included only once
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

# Include dependent projects if any
set(extProjName ITK) #The find_package known name
set(proj ITKv4)      #This local name

#if(${USE_SYSTEM_${extProjName}})
#  unset(${extProjName}_DIR CACHE)
#endif()

# Sanity checks
if(DEFINED ${extProjName}_DIR AND NOT EXISTS ${${extProjName}_DIR})
  message(FATAL_ERROR "${extProjName}_DIR variable is defined but corresponds to non-existing directory")
endif()

# Set dependency list
set(${proj}_DEPENDENCIES "")

SlicerMacroCheckExternalProjectDependency(${proj})

if(NOT DEFINED ${extProjName}_DIR AND NOT ${USE_SYSTEM_${extProjName}})

  # Set CMake OSX variable to pass down the external project
  set(CMAKE_OSX_EXTERNAL_PROJECT_ARGS)
  if(APPLE)
    list(APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
      -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
      -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
      -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET})
  endif()

  ### --- Project specific additions here
  set(ITKv4_WRAP_ARGS)
  #if(foo)
    #set(ITKv4_WRAP_ARGS
    #  -DINSTALL_WRAP_ITK_COMPATIBILITY:BOOL=OFF
    #  -DWRAP_float:BOOL=ON
    #  -DWRAP_unsigned_char:BOOL=ON
    #  -DWRAP_signed_short:BOOL=ON
    #  -DWRAP_unsigned_short:BOOL=ON
    #  -DWRAP_complex_float:BOOL=ON
    #  -DWRAP_vector_float:BOOL=ON
    #  -DWRAP_covariant_vector_float:BOOL=ON
    #  -DWRAP_rgb_signed_short:BOOL=ON
    #  -DWRAP_rgb_unsigned_char:BOOL=ON
    #  -DWRAP_rgb_unsigned_short:BOOL=ON
    #  -DWRAP_ITK_TCL:BOOL=OFF
    #  -DWRAP_ITK_JAVA:BOOL=OFF
    #  -DWRAP_ITK_PYTHON:BOOL=ON
    #  -DPYTHON_EXECUTABLE:PATH=${${CMAKE_PROJECT_NAME}_PYTHON_EXECUTABLE}
    #  -DPYTHON_INCLUDE_DIR:PATH=${${CMAKE_PROJECT_NAME}_PYTHON_INCLUDE}
    #  -DPYTHON_LIBRARY:FILEPATH=${${CMAKE_PROJECT_NAME}_PYTHON_LIBRARY}
    #  )
  #endif()
  # HACK This code fixes a loony problem with HDF5 -- it doesn't
  #      link properly if -fopenmp is used.
  string(REPLACE "-fopenmp" "" ITK_CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  string(REPLACE "-fopenmp" "" ITK_CMAKE_CXX_FLAGS "${CMAKE_CX_FLAGS}")

  set(${proj}_CMAKE_OPTIONS
      -DITK_LEGACY_REMOVE:BOOL=ON
      -DITK_BUILD_ALL_MODULES:BOOL=ON
      -DITK_USE_REVIEW:BOOL=ON
      -DITKV3_COMPATIBILITY:BOOL=OFF
      -DKWSYS_USE_MD5:BOOL=ON # Required by SlicerExecutionModel
      -DUSE_WRAP_ITK:BOOL=OFF ## HACK:  QUICK CHANGE
    )
  ### --- End Project specific additions
  set(${proj}_REPOSITORY ${git_protocol}://itk.org/ITK.git)
  set(${proj}_GIT_TAG 827b8f9c50d43b97f1f52b1c9598de68aae3f7dd) #2012-08-15
  ExternalProject_Add(${proj}
    GIT_REPOSITORY ${${proj}_REPOSITORY}
    GIT_TAG ${${proj}_GIT_TAG}
    SOURCE_DIR ${proj}
    BINARY_DIR ${proj}-build
    UPDATE_COMMAND ""
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
      ${COMMON_EXTERNAL_PROJECT_ARGS}
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      ${${proj}_CMAKE_OPTIONS}
    INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES}
    )
  set(${extProjName}_DIR ${CMAKE_BINARY_DIR}/${proj}-build)
else()
  if(${USE_SYSTEM_${extProjName}})
    find_package(${extProjName} ${ITK_VERSION_MAJOR} REQUIRED)
    if(NOT ${extProjName}_DIR)
      message(FATAL_ERROR "To use the system ${extProjName}, set ${extProjName}_DIR")
    endif()
  endif()
  # The project is provided using ${extProjName}_DIR, nevertheless since other
  # project may depend on ${extProjName}v4, let's add an 'empty' one
  SlicerMacroEmptyExternalProject(${proj} "${${proj}_DEPENDENCIES}")
endif()

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS ${extProjName}_DIR:PATH)
