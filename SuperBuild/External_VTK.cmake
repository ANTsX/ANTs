
# Make sure this file is included only once
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

# Include dependent projects if any
set(extProjName VTK)     #The find_package known name
set(proj ${extProjName}) #This local name

#if(${USE_SYSTEM_${extProjName}})
#  unset(${extProjName}_DIR CACHE)
#endif()

# Sanity checks
if(DEFINED ${extProjName}_DIR AND NOT EXISTS ${${extProjName}_DIR})
  message(FATAL_ERROR "${extProjName}_DIR variable is defined but corresponds to non-existing directory (${${extProjName}_DIR})")
endif()

# Set dependency list
set(${proj}_DEPENDENCIES "")
if (Slicer_USE_PYTHONQT)
  list(APPEND ${proj}_DEPENDENCIES python)
endif()

SlicerMacroCheckExternalProjectDependency(${proj})

if(NOT DEFINED ${extProjName}_DIR AND NOT ${USE_SYSTEM_${extProjName}})

  # Set CMake OSX variable to pass down the external project
  set(CMAKE_OSX_EXTERNAL_PROJECT_ARGS)
  if(APPLE)
    list(APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
      -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
      -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
      -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
      -DVTK_REQUIRED_OBJCXX_FLAGS:STRING="")
  endif()

  ### --- Project specific additions here
  set(VTK_WRAP_TCL OFF)
  set(VTK_WRAP_PYTHON OFF)

  if (Slicer_USE_PYTHONQT)
    set(VTK_WRAP_PYTHON ON)
  endif()

  set(VTK_PYTHON_ARGS)
  if(Slicer_USE_PYTHONQT)
    set(VTK_PYTHON_ARGS
      -DVTK_INSTALL_PYTHON_USING_CMAKE:BOOL=ON
      -DPYTHON_EXECUTABLE:PATH=${slicer_PYTHON_EXECUTABLE}
      -DPYTHON_INCLUDE_DIR:PATH=${slicer_PYTHON_INCLUDE}
      -DPYTHON_LIBRARY:FILEPATH=${slicer_PYTHON_LIBRARY}
      )
  endif()

  set(VTK_QT_ARGS)
  if(${LOCAL_PROJECT_NAME}_USE_QT)
    if(NOT APPLE)
      set(VTK_QT_ARGS
        #-DDESIRED_QT_VERSION:STRING=4 # Unused
        -DVTK_USE_GUISUPPORT:BOOL=ON
        -DVTK_USE_QVTK_QTOPENGL:BOOL=ON
        -DVTK_USE_QT:BOOL=ON
        -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
        )
    else()
      set(VTK_QT_ARGS
        -DVTK_USE_CARBON:BOOL=OFF
        # Default to Cocoa, VTK/CMakeLists.txt will enable Carbon and disable cocoa if needed
        -DVTK_USE_COCOA:BOOL=ON
        -DVTK_USE_X:BOOL=OFF
        #-DVTK_USE_RPATH:BOOL=ON # Unused
        #-DDESIRED_QT_VERSION:STRING=4 # Unused
        -DVTK_USE_GUISUPPORT:BOOL=ON
        -DVTK_USE_QVTK_QTOPENGL:BOOL=ON
        -DVTK_USE_QT:BOOL=ON
        -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
        )
    endif()
    find_package(Qt4 REQUIRED)
  else()
    set(VTK_QT_ARGS
        -DVTK_USE_GUISUPPORT:BOOL=OFF
        -DVTK_USE_QT:BOOL=OFF
        )
  endif()

  # Disable Tk when Python wrapping is enabled
  if (Slicer_USE_PYTHONQT)
    list(APPEND VTK_QT_ARGS -DVTK_USE_TK:BOOL=OFF)
  endif()

  set(slicer_TCL_LIB)
  set(slicer_TK_LIB)
  set(slicer_TCLSH)
  set(VTK_TCL_ARGS)
  if(VTK_WRAP_TCL)
    if(WIN32)
      set(slicer_TCL_LIB ${CMAKE_BINARY_DIR}/tcl-build/lib/tcl84.lib)
      set(slicer_TK_LIB ${CMAKE_BINARY_DIR}/tcl-build/lib/tk84.lib)
      set(slicer_TCLSH ${CMAKE_BINARY_DIR}/tcl-build/bin/tclsh.exe)
    elseif(APPLE)
      set(slicer_TCL_LIB ${CMAKE_BINARY_DIR}/tcl-build/lib/libtcl8.4.dylib)
      set(slicer_TK_LIB ${CMAKE_BINARY_DIR}/tcl-build/lib/libtk8.4.dylib)
      set(slicer_TCLSH ${CMAKE_BINARY_DIR}/tcl-build/bin/tclsh84)
    else()
      set(slicer_TCL_LIB ${CMAKE_BINARY_DIR}/tcl-build/lib/libtcl8.4.so)
      set(slicer_TK_LIB ${CMAKE_BINARY_DIR}/tcl-build/lib/libtk8.4.so)
      set(slicer_TCLSH ${CMAKE_BINARY_DIR}/tcl-build/bin/tclsh84)
    endif()
    set(VTK_TCL_ARGS
      -DTCL_INCLUDE_PATH:PATH=${CMAKE_BINARY_DIR}/tcl-build/include
      -DTK_INCLUDE_PATH:PATH=${CMAKE_BINARY_DIR}/tcl-build/include
      -DTCL_LIBRARY:FILEPATH=${slicer_TCL_LIB}
      -DTK_LIBRARY:FILEPATH=${slicer_TK_LIB}
      -DTCL_TCLSH:FILEPATH=${slicer_TCLSH}
      )
  endif()

  set(VTK_BUILD_STEP "")
  if(UNIX)
    configure_file(SuperBuild/VTK_build_step.cmake.in
      ${CMAKE_CURRENT_BINARY_DIR}/VTK_build_step.cmake
      @ONLY)
    set(VTK_BUILD_STEP ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/VTK_build_step.cmake)
  endif()

  set(${proj}_CMAKE_OPTIONS
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/${proj}-install
      -DVTK_USE_PARALLEL:BOOL=ON
      -DVTK_DEBUG_LEAKS:BOOL=${Slicer_USE_VTK_DEBUG_LEAKS}
      -DVTK_LEGACY_REMOVE:BOOL=OFF
      -DVTK_WRAP_TCL:BOOL=${VTK_WRAP_TCL}
      #-DVTK_USE_RPATH:BOOL=ON # Unused
      ${VTK_TCL_ARGS}
      -DVTK_WRAP_PYTHON:BOOL=${VTK_WRAP_PYTHON}
      -DVTK_INSTALL_LIB_DIR:PATH=${Slicer_INSTALL_LIB_DIR}
      ${VTK_PYTHON_ARGS}
      ${VTK_QT_ARGS}
      ${VTK_MAC_ARGS}
      )
  ### --- End Project specific additions
  set(${proj}_REPOSITORY ${git_protocol}://vtk.org/VTK.git)
  set(${proj}_GIT_TAG "v5.10.0")
  ExternalProject_Add(${proj}
    GIT_REPOSITORY ${${proj}_REPOSITORY}
    GIT_TAG ${${proj}_GIT_TAG}
    SOURCE_DIR ${proj}
    BINARY_DIR ${proj}-build
    LOG_CONFIGURE 0  # Wrap configure in script to ignore log output from dashboards
    LOG_BUILD     0  # Wrap build in script to to ignore log output from dashboards
    LOG_TEST      0  # Wrap test in script to to ignore log output from dashboards
    LOG_INSTALL   0  # Wrap install in script to to ignore log output from dashboards
    ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      -Wno-dev
      --no-warn-unused-cli
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
      ${COMMON_EXTERNAL_PROJECT_ARGS}
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      ${${proj}_CMAKE_OPTIONS}
## We really do want to install in order to limit # of include paths INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES}
    BUILD_COMMAND ${VTK_BUILD_STEP}
    )

  set(VTKPatchScript ${CMAKE_CURRENT_LIST_DIR}/VTKPatch.cmake)
  ExternalProject_Add_Step(${proj} VTKPatch
    COMMENT "get rid of obsolete C/CXX flags"
    DEPENDEES download
    DEPENDERS configure
    COMMAND ${CMAKE_COMMAND}
    -DVTKSource=<SOURCE_DIR>
    -P ${VTKPatchScript}
    )

  set(${extProjName}_DIR ${CMAKE_BINARY_DIR}/${proj}-install/lib/vtk-5.10)
else()
  if(${USE_SYSTEM_${extProjName}})
    set(VTK_VERSION_MAJOR 5)
    find_package(${extProjName} REQUIRED)
    if(NOT ${extProjName}_DIR)
      message(FATAL_ERROR "To use the system ${extProjName}, set ${extProjName}_DIR")
    endif()
    message("USING the system ${extProjName}, set ${extProjName}_DIR=${${extProjName}_DIR}")
  endif()
  # The project is provided using ${extProjName}_DIR, nevertheless since other
  # project may depend on ${extProjName}v4, let's add an 'empty' one
  SlicerMacroEmptyExternalProject(${proj} "${${proj}_DEPENDENCIES}")
endif()

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS ${extProjName}_DIR:PATH)

