
include(${CMAKE_CURRENT_LIST_DIR}/Common.cmake)

configure_file(${CMAKE_CURRENT_LIST_DIR}/CTestCustom.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/CTestCustom.cmake COPYONLY)

set(BUILDSCRIPTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/BuildScripts)
set(CMAKE_MODULE_PATH
    ${BUILDSCRIPTS_DIR}
    ${${PROJECT_NAME}_SOURCE_DIR}/CMake
    ${${PROJECT_NAME}_BINARY_DIR}/CMake
    ${CMAKE_MODULE_PATH}
    )

set (CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

# Set up ITK
find_package(ITK 4 REQUIRED)
include(${ITK_USE_FILE})

# Set up VTK
option(USE_VTK "Use VTK Libraries" OFF)
if(USE_VTK)
  find_package(VTK)
  if(VTK_FOUND)
    include(${VTK_USE_FILE})
  else(VTK_FOUND)
     message("Cannot build some programs without VTK.  Please set VTK_DIR if you need these programs.")
  endif(VTK_FOUND)
endif(USE_VTK)

# With MS compilers on Win64, we need the /bigobj switch, else generated
# code results in objects with number of sections exceeding object file
# format.
# see http://msdn.microsoft.com/en-us/library/ms173499.aspx
if(CMAKE_CL_64 OR MSVC)
  add_definitions(/bigobj)
endif()

option(USE_FFTWD "Use double precision fftw if found" OFF)
option(USE_FFTWF "Use single precision fftw if found" OFF)
option(USE_SYSTEM_FFTW "Use an installed version of fftw" OFF)
if (USE_FFTWD OR USE_FFTWF)
  if(USE_SYSTEM_FFTW)
      find_package( FFTW )
      link_directories(${FFTW_LIBDIR})
  else(USE_SYSTEM_FFTW)
      link_directories(${ITK_DIR}/fftw/lib)
      include_directories(${ITK_DIR}/fftw/include)
  endif(USE_SYSTEM_FFTW)
endif(USE_FFTWD OR USE_FFTWF)

# These are configure time options that specify which
# subset of tests should be run
option(RUN_SHORT_TESTS    "Run the quick unit tests."                                   ON  )
option(RUN_LONG_TESTS     "Run the time consuming tests. i.e. real world registrations" OFF )
option(OLD_BASELINE_TESTS "Use reported metrics from old tests"                         OFF )
#-----------------------------------------------------------------------------
include(CTest)
enable_testing()

configure_file(${CMAKE_CURRENT_LIST_DIR}/CTestCustom.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/CTestCustom.cmake COPYONLY)

add_subdirectory(Examples)
