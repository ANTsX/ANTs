
include(${CMAKE_CURRENT_LIST_DIR}/Common.cmake)

configure_file(${CMAKE_CURRENT_LIST_DIR}/CTestCustom.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/CTestCustom.cmake COPYONLY)

set(CMAKE_MODULE_PATH
    ${${PROJECT_NAME}_SOURCE_DIR}/CMake
    ${${PROJECT_NAME}_BINARY_DIR}/CMake
    ${CMAKE_MODULE_PATH}
    )

set (CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

# Set up ITK
find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

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

option(ITK_USE_FFTWD "Use double precision fftw if found" OFF)
option(ITK_USE_FFTWF "Use single precision fftw if found" OFF)
option(ITK_USE_SYSTEM_FFTW "Use an installed version of fftw" OFF)
if (ITK_USE_FFTWD OR ITK_USE_FFTWF)
  if(ITK_USE_SYSTEM_FFTW)
      find_package( FFTW )
      link_directories(${FFTW_LIBDIR})
  else(ITK_USE_SYSTEM_FFTW)
      link_directories(${ITK_DIR}/fftw/lib)
      include_directories(${ITK_DIR}/fftw/include)
  endif(ITK_USE_SYSTEM_FFTW)
endif(ITK_USE_FFTWD OR ITK_USE_FFTWF)

# These are configure time options that specify which
# subset of tests should be run
option(RUN_SHORT_TESTS    "Run the quick unit tests."                                   ON  )
option(RUN_LONG_TESTS     "Run the time consuming tests. i.e. real world registrations" OFF )
option(OLD_BASELINE_TESTS "Use reported metrics from old tests"                         OFF )
#-----------------------------------------------------------------------------
include(CTest)
enable_testing()
#Set the global max TIMEOUT for CTest jobs.  This is very large for the moment
#and should be revisted to reduce based on "LONG/SHORT" test times, set to 1 hr for now
set(CTEST_TEST_TIMEOUT 1800 CACHE STRING "Maximum seconds allowed before CTest will kill the test." FORCE)
set(DART_TESTING_TIMEOUT ${CTEST_TEST_TIMEOUT} CACHE STRING "Maximum seconds allowed before CTest will kill the test." FORCE)

configure_file(${CMAKE_CURRENT_LIST_DIR}/CTestCustom.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/CTestCustom.cmake COPYONLY)

add_subdirectory(Examples)
if(BUILD_EXTERNAL_APPLICATIONS)
  add_subdirectory(ExternalApplications)
endif()
