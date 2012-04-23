project(ANTS)
cmake_minimum_required(VERSION 2.8.2)

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
#The following lines are required to use Dart
include(CTest)
enable_testing()
#include(${CMAKE_ROOT}/Modules/FindITK.cmake)
#if(USE_ITK_FILE)
#  include(${USE_ITK_FILE})
#endif(USE_ITK_FILE)
set(DART_TESTING_TIMEOUT 1500)
set(PICSL_INCLUDE_DIRS
../Utilities
../ImageRegistration
../ImageSegmentation
../GraphTheory
../Tensor
../Temporary
.
)
include_directories(${PICSL_INCLUDE_DIRS})
link_directories( ${ITK_LIBRARY_PATH}  )
# message("${ITK_LIBRARIES}")
# non-templated class -- this should be stored in a library and linked in...
set(UI_SOURCES "../Utilities/antsCommandLineParser" "../Utilities/antsCommandLineOption" )

## Note that the antsUtilities can always be built static.  It will then be linked
## Directly into the other libraries.
add_library(antsUtilities STATIC antsUtilities.cxx )
target_link_libraries(antsUtilities ${ITK_LIBRARIES} )


add_library(l_ANTS ANTS.cxx ${UI_SOURCES})
target_link_libraries(l_ANTS antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME ANTS )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ANTSJacobian ANTSJacobian.cxx ${UI_SOURCES})
target_link_libraries(l_ANTSJacobian antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME ANTSJacobian )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_PrintHeader PrintHeader.cxx ${UI_SOURCES})
target_link_libraries(l_PrintHeader antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME PrintHeader )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ResetDirection ResetDirection.cxx ${UI_SOURCES})
target_link_libraries(l_ResetDirection antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME ResetDirection )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

# add_library(l_ResetDirection2 ResetDirection2.cxx ${UI_SOURCES})
# target_link_libraries(l_ResetDirection2 antsUtilities ${ITK_LIBRARIES} )
# set( ANTS_FUNCTION_NAME ResetDirection2 )
# configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
# add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
# target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

# add_library(l_BoykovGraphCutFilter BoykovGraphCutFilter.cxx ${UI_SOURCES})
# target_link_libraries(l_BoykovGraphCutFilter antsUtilities ${ITK_LIBRARIES} )
# set( ANTS_FUNCTION_NAME BoykovGraphCutFilter )
# configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
# add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
# target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

# add_library(l_GeneralLinearModel GeneralLinearModel.cxx ${UI_SOURCES})
# target_link_libraries(l_GeneralLinearModel antsUtilities ${ITK_LIBRARIES} )
# set( ANTS_FUNCTION_NAME GeneralLinearModel )
# configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
# add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
# target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ANTSUseLandmarkImagesToGetAffineTransform ANTSUseLandmarkImagesToGetAffineTransform.cxx ${UI_SOURCES})
target_link_libraries(l_ANTSUseLandmarkImagesToGetAffineTransform antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME ANTSUseLandmarkImagesToGetAffineTransform )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ANTSUseLandmarkImagesToGetBSplineDisplacementField ANTSUseLandmarkImagesToGetBSplineDisplacementField.cxx ${UI_SOURCES})
target_link_libraries(l_ANTSUseLandmarkImagesToGetBSplineDisplacementField antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME ANTSUseLandmarkImagesToGetBSplineDisplacementField )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ANTSUseDeformationFieldToGetAffineTransform ANTSUseDeformationFieldToGetAffineTransform ${UI_SOURCES})
target_link_libraries(l_ANTSUseDeformationFieldToGetAffineTransform antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME ANTSUseDeformationFieldToGetAffineTransform )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

# add_library(l_FDR FDR.cxx ${UI_SOURCES})
# target_link_libraries(l_FDR antsUtilities ${ITK_LIBRARIES} )
# set( ANTS_FUNCTION_NAME FDR )
# configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
# add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
# target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ImageMath ImageMath.cxx ${UI_SOURCES})
target_link_libraries(l_ImageMath antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME ImageMath )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_SurfaceCurvature SurfaceCurvature.cxx ${UI_SOURCES})
target_link_libraries(l_SurfaceCurvature antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME SurfaceCurvature )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_SurfaceBasedSmoothing SurfaceBasedSmoothing.cxx ${UI_SOURCES})
target_link_libraries(l_SurfaceBasedSmoothing antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME SurfaceBasedSmoothing )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )


# add_library(l_test test.cxx)
# target_link_libraries(l_test)
# set( ANTS_FUNCTION_NAME test )
# configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
# add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
# target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

# add_library(l_WarpImage WarpImage.cxx ${UI_SOURCES})
# target_link_libraries(l_WarpImage antsUtilities ${ITK_LIBRARIES} )
# set( ANTS_FUNCTION_NAME WarpImage )
# configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
# add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
# target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_LaplacianThickness LaplacianThickness.cxx ${UI_SOURCES})
target_link_libraries(l_LaplacianThickness antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME LaplacianThickness )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_SetOrigin SetOrigin.cxx ${UI_SOURCES})
target_link_libraries(l_SetOrigin antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME SetOrigin )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_SetSpacing SetSpacing.cxx ${UI_SOURCES})
target_link_libraries(l_SetSpacing antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME SetSpacing )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_SetDirectionByMatrix SetDirectionByMatrix.cxx ${UI_SOURCES})
target_link_libraries(l_SetDirectionByMatrix antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME SetDirectionByMatrix )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

# add_library(l_WarpFunctionalImage WarpFunctionalImage.cxx ${UI_SOURCES})
# target_link_libraries(l_WarpFunctionalImage antsUtilities ${ITK_LIBRARIES} )
# set( ANTS_FUNCTION_NAME WarpFunctionalImage )
# configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
# add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
# target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ConvertScalarImageToRGB ConvertScalarImageToRGB.cxx)
target_link_libraries(l_ConvertScalarImageToRGB antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME ConvertScalarImageToRGB )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_CreateWarpedGridImage CreateWarpedGridImage.cxx)
target_link_libraries(l_CreateWarpedGridImage antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME CreateWarpedGridImage )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

# add_library(l_CreateJacobianDeterminantImage CreateJacobianDeterminantImage.cxx)
# target_link_libraries(l_CreateJacobianDeterminantImage antsUtilities ${ITK_LIBRARIES})
# set( ANTS_FUNCTION_NAME CreateJacobianDeterminantImage )
# configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
# add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
# target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_MeasureImageSimilarity MeasureImageSimilarity.cxx)
target_link_libraries(l_MeasureImageSimilarity antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME MeasureImageSimilarity )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ConvertToJpg ConvertToJpg.cxx)
target_link_libraries(l_ConvertToJpg antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME ConvertToJpg )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ConvertImagePixelType ConvertImagePixelType.cxx ${UI_SOURCES})
target_link_libraries(l_ConvertImagePixelType antsUtilities ${ITK_LIBRARIES}    )
set( ANTS_FUNCTION_NAME ConvertImagePixelType )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_AverageImages AverageImages.cxx ${UI_SOURCES})
target_link_libraries(l_AverageImages antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME AverageImages )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_AverageTensorImages AverageTensorImages.cxx ${UI_SOURCES})
target_link_libraries(l_AverageTensorImages antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME AverageTensorImages )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ImageSetStatistics ImageSetStatistics.cxx ${UI_SOURCES})
target_link_libraries(l_ImageSetStatistics antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME ImageSetStatistics )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ThresholdImage ThresholdImage.cxx ${UI_SOURCES})
target_link_libraries(l_ThresholdImage antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME ThresholdImage )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_MultiplyImages MultiplyImages.cxx ${UI_SOURCES})
target_link_libraries(l_MultiplyImages antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME MultiplyImages )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_SmoothImage SmoothImage.cxx ${UI_SOURCES})
target_link_libraries(l_SmoothImage antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME SmoothImage )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ClusterImageStatistics ClusterImageStatistics.cxx ${UI_SOURCES})
target_link_libraries(l_ClusterImageStatistics antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME ClusterImageStatistics )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_LabelClustersUniquely LabelClustersUniquely.cxx ${UI_SOURCES})
target_link_libraries(l_LabelClustersUniquely antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME LabelClustersUniquely )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

# add_library(l_StudentsTestOnImages StudentsTestOnImages.cxx ${UI_SOURCES})
# target_link_libraries(l_StudentsTestOnImages antsUtilities ${ITK_LIBRARIES} )
# set( ANTS_FUNCTION_NAME StudentsTestOnImages )
# configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
# add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
# target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_LabelOverlapMeasures LabelOverlapMeasures.cxx )
target_link_libraries(l_LabelOverlapMeasures antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME LabelOverlapMeasures )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_MeasureMinMaxMean MeasureMinMaxMean.cxx ${UI_SOURCES})
target_link_libraries(l_MeasureMinMaxMean antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME MeasureMinMaxMean )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_WarpImageMultiTransform WarpImageMultiTransform ${UI_SOURCES})
target_link_libraries(l_WarpImageMultiTransform antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME WarpImageMultiTransform )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ComposeMultiTransform ComposeMultiTransform ${UI_SOURCES})
target_link_libraries(l_ComposeMultiTransform antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME ComposeMultiTransform )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_StackSlices StackSlices.cxx ${UI_SOURCES})
target_link_libraries(l_StackSlices antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME StackSlices )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_MemoryTest MemoryTest.cxx ${UI_SOURCES})
target_link_libraries(l_MemoryTest antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME MemoryTest )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

# add_library(l_ANTSOrientImage ANTSOrientImage.cxx ${UI_SOURCES})
# target_link_libraries(l_ANTSOrientImage antsUtilities ${ITK_LIBRARIES} )
# set( ANTS_FUNCTION_NAME ANTSOrientImage )
# configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
# add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
# target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_PermuteFlipImageOrientationAxes PermuteFlipImageOrientationAxes.cxx ${UI_SOURCES})
target_link_libraries(l_PermuteFlipImageOrientationAxes antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME PermuteFlipImageOrientationAxes )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ImageCompare ImageCompare.cxx)
target_link_libraries(l_ImageCompare antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME ImageCompare )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ResampleImageBySpacing ResampleImageBySpacing.cxx)
target_link_libraries(l_ResampleImageBySpacing antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME ResampleImageBySpacing )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

# add_library(l_ResampleImageByOtherImage ResampleImageByOtherImage.cxx)
# target_link_libraries(l_ResampleImageByOtherImage antsUtilities ${ITK_LIBRARIES})
# set( ANTS_FUNCTION_NAME ResampleImageByOtherImage )
# configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
# add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
# target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_CopyImageHeaderInformation CopyImageHeaderInformation.cxx)
target_link_libraries(l_CopyImageHeaderInformation antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME CopyImageHeaderInformation )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_WarpTimeSeriesImageMultiTransform WarpTimeSeriesImageMultiTransform.cxx)
target_link_libraries(l_WarpTimeSeriesImageMultiTransform antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME WarpTimeSeriesImageMultiTransform )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ResampleImage ResampleImage.cxx )
target_link_libraries(l_ResampleImage antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME ResampleImage )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ExtractSliceFromImage ExtractSliceFromImage.cxx )
target_link_libraries(l_ExtractSliceFromImage antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME ExtractSliceFromImage )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ExtractRegionFromImage ExtractRegionFromImage.cxx )
target_link_libraries(l_ExtractRegionFromImage antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME ExtractRegionFromImage )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ExtractRegionFromImageByMask ExtractRegionFromImageByMask.cxx )
target_link_libraries(l_ExtractRegionFromImageByMask antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME ExtractRegionFromImageByMask )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_LabelGeometryMeasures LabelGeometryMeasures.cxx )
target_link_libraries(l_LabelGeometryMeasures antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME LabelGeometryMeasures )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_PasteImageIntoImage PasteImageIntoImage.cxx )
target_link_libraries(l_PasteImageIntoImage antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME PasteImageIntoImage )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_TileImages TileImages.cxx )
target_link_libraries(l_TileImages antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME TileImages )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_CreateImage CreateImage.cxx )
target_link_libraries(l_CreateImage antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME CreateImage )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_WarpTensorImageMultiTransform WarpTensorImageMultiTransform.cxx)
target_link_libraries(l_WarpTensorImageMultiTransform antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME WarpTensorImageMultiTransform )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ReorientTensorImage ReorientTensorImage.cxx)
target_link_libraries(l_ReorientTensorImage antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME ReorientTensorImage )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_RebaseTensorImage RebaseTensorImage.cxx)
target_link_libraries(l_RebaseTensorImage antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME RebaseTensorImage )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_N3BiasFieldCorrection N3BiasFieldCorrection.cxx)
target_link_libraries(l_N3BiasFieldCorrection antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME N3BiasFieldCorrection )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_N4BiasFieldCorrection N4BiasFieldCorrection.cxx  ${UI_SOURCES})
target_link_libraries(l_N4BiasFieldCorrection antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME N4BiasFieldCorrection )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_KellySlater KellySlater.cxx)
target_link_libraries(l_KellySlater antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME KellySlater )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_KellyKapowski KellyKapowski.cxx  ${UI_SOURCES})
target_link_libraries(l_KellyKapowski antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME KellyKapowski )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_CreateDTICohort CreateDTICohort.cxx  ${UI_SOURCES})
target_link_libraries(l_CreateDTICohort antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME CreateDTICohort )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_AverageAffineTransform AverageAffineTransform.cxx)
target_link_libraries(l_AverageAffineTransform antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME AverageAffineTransform )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_antsRegistration antsRegistration.cxx ${UI_SOURCES})
target_link_libraries(l_antsRegistration antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME antsRegistration )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_antsApplyTransforms antsApplyTransforms.cxx ${UI_SOURCES})
target_link_libraries(l_antsApplyTransforms antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME antsApplyTransforms )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_CompositeTransformUtil CompositeTransformUtil.cxx)
target_link_libraries(l_CompositeTransformUtil antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME CompositeTransformUtil )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_antsMotionCorr antsMotionCorr.cxx ${UI_SOURCES})
target_link_libraries(l_antsMotionCorr antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME antsMotionCorr )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_sccan sccan.cxx ${UI_SOURCES})
target_link_libraries(l_sccan antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME sccan )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_CreateDisplacementField CreateDisplacementField.cxx ${UI_SOURCES})
target_link_libraries(l_CreateDisplacementField antsUtilities ${ITK_LIBRARIES} )
set( ANTS_FUNCTION_NAME CreateDisplacementField )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

if(USE_VTK)
include(${CMAKE_ROOT}/Modules/FindVTK.cmake)
if(USE_VTK_FILE)
  include(${USE_VTK_FILE})
  add_library(l_ConvertVectorFieldToVTK ConvertVectorFieldToVTK.cxx )
  target_link_libraries(l_ConvertVectorFieldToVTK antsUtilities ${ITK_LIBRARIES} vtkIO)
  set( ANTS_FUNCTION_NAME ConvertVectorFieldToVTK )
  configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
  add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
  target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

  add_library(l_GetMeshAndTopology GetMeshAndTopology.cxx)
  target_link_libraries(l_GetMeshAndTopology antsUtilities ${ITK_LIBRARIES} ITKMetaIO  vtkRendering vtksys vtkIO)
  set( ANTS_FUNCTION_NAME GetMeshAndTopology )
  configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
  add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
  target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

  add_library(l_CheckTopology CheckTopology.cxx)
  target_link_libraries(l_CheckTopology antsUtilities ${ITK_LIBRARIES} ITKMetaIO  vtkRendering vtksys vtkIO)
  set( ANTS_FUNCTION_NAME CheckTopology )
  configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
  add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
  target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

  add_library(l_WarpVTKPolyDataMultiTransform WarpVTKPolyDataMultiTransform)
  target_link_libraries(l_WarpVTKPolyDataMultiTransform antsUtilities ${ITK_LIBRARIES} vtksys vtkIO )
  set( ANTS_FUNCTION_NAME WarpVTKPolyDataMultiTransform )
  configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
  add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
  target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

#   set(FLAT_SRC
#     ../Temporary/itkFEMElement3DC0LinearTriangular.cxx
#     ../Temporary/itkFEMElement3DC0LinearTriangularLaplaceBeltrami.cxx
#     ../Temporary/itkFEMElement3DC0LinearTriangularMembrane.cxx
#   )
#   # Build the library
#   add_library(l_FLATFEM ${FLAT_SRC})
#   target_link_libraries(l_FLATFEM antsUtilities ${ITK_LIBRARIES})
#   set( ANTS_FUNCTION_NAME FLATFEM )
#   configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
#   add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
#   target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

#
#
#   add_library(l_ConformalMapping ConformalMapping.cxx)
#   target_link_libraries(l_ConformalMapping antsUtilities ${ITK_LIBRARIES} FLATFEM antsUtilities ${ITK_LIBRARIES} ITKMetaIO  vtkRendering vtkIO)
#   set( ANTS_FUNCTION_NAME ConformalMapping )
#   configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
#   add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
#   target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

#
#   add_library(l_ANTSConformalMapping ANTSConformalMapping.cxx ${UI_SOURCES})
#   target_link_libraries(l_ANTSConformalMapping antsUtilities ${ITK_LIBRARIES} FLATFEM antsUtilities ${ITK_LIBRARIES} ITKMetaIO  vtkRendering vtkIO)
#   set( ANTS_FUNCTION_NAME ANTSConformalMapping )
#   configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
#   add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
#   target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

endif(USE_VTK_FILE)
endif(USE_VTK)

add_library(l_ANTSIntegrateVectorField ANTSIntegrateVectorField.cxx)
target_link_libraries(l_ANTSIntegrateVectorField antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME ANTSIntegrateVectorField )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_ANTSIntegrateVelocityField ANTSIntegrateVelocityField.cxx)
target_link_libraries(l_ANTSIntegrateVelocityField antsUtilities ${ITK_LIBRARIES})
set( ANTS_FUNCTION_NAME ANTSIntegrateVelocityField )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_AtroposMin AtroposMin.cxx ${UI_SOURCES})
target_link_libraries(l_AtroposMin antsUtilities ${ITK_LIBRARIES}   )
set( ANTS_FUNCTION_NAME AtroposMin )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

add_library(l_Atropos Atropos.cxx ${UI_SOURCES})
target_link_libraries(l_Atropos antsUtilities ${ITK_LIBRARIES}   )
set( ANTS_FUNCTION_NAME Atropos )
configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

install(TARGETS ANTS
  AverageAffineTransform
  l_AverageAffineTransform
  ANTSJacobian
  l_ANTSJacobian
  PasteImageIntoImage
  l_PasteImageIntoImage
  PrintHeader
  l_PrintHeader
  ResetDirection
  l_ResetDirection
  ANTSUseLandmarkImagesToGetAffineTransform
  l_ANTSUseLandmarkImagesToGetAffineTransform
  ANTSUseLandmarkImagesToGetBSplineDisplacementField
  l_ANTSUseLandmarkImagesToGetBSplineDisplacementField
  ANTSUseDeformationFieldToGetAffineTransform
  l_ANTSUseDeformationFieldToGetAffineTransform
  AverageAffineTransform
  l_AverageAffineTransform
  ImageMath
  l_ImageMath
  SurfaceCurvature
  l_SurfaceCurvature
  SurfaceBasedSmoothing
  l_SurfaceBasedSmoothing
  LaplacianThickness
  l_LaplacianThickness
  SetOrigin
  l_SetOrigin
  SetSpacing
  l_SetSpacing
  SetDirectionByMatrix
  l_SetDirectionByMatrix
  ConvertScalarImageToRGB
  l_ConvertScalarImageToRGB
  CreateWarpedGridImage
  l_CreateWarpedGridImage
  MeasureImageSimilarity
  l_MeasureImageSimilarity
  ConvertToJpg
  l_ConvertToJpg
  ConvertImagePixelType
  l_ConvertImagePixelType
  AverageImages
  l_AverageImages
  AverageTensorImages
  l_AverageTensorImages
  ImageSetStatistics
  l_ImageSetStatistics
  ThresholdImage
  l_ThresholdImage
  MultiplyImages
  l_MultiplyImages
  SmoothImage
  l_SmoothImage
  ClusterImageStatistics
  l_ClusterImageStatistics
  LabelClustersUniquely
  l_LabelClustersUniquely
  LabelOverlapMeasures
  l_LabelOverlapMeasures
  MeasureMinMaxMean
  l_MeasureMinMaxMean
  WarpImageMultiTransform
  l_WarpImageMultiTransform
  ComposeMultiTransform
  l_ComposeMultiTransform
  StackSlices
  l_StackSlices
  PermuteFlipImageOrientationAxes
  l_PermuteFlipImageOrientationAxes
  ImageCompare
  l_ImageCompare
  ResampleImageBySpacing
  l_ResampleImageBySpacing
  CopyImageHeaderInformation
  l_CopyImageHeaderInformation
  WarpTimeSeriesImageMultiTransform
  l_WarpTimeSeriesImageMultiTransform
  ResampleImage
  l_ResampleImage
  ExtractSliceFromImage
  l_ExtractSliceFromImage
  ExtractRegionFromImage
  l_ExtractRegionFromImage
  ExtractRegionFromImageByMask
  l_ExtractRegionFromImageByMask
  LabelGeometryMeasures
  l_LabelGeometryMeasures
  TileImages
  l_TileImages
  CreateImage
  l_CreateImage
  WarpTensorImageMultiTransform
  l_WarpTensorImageMultiTransform
  ReorientTensorImage
  l_ReorientTensorImage
  RebaseTensorImage
  l_RebaseTensorImage
  N3BiasFieldCorrection
  l_N3BiasFieldCorrection
  N4BiasFieldCorrection
  l_N4BiasFieldCorrection
  KellySlater
  l_KellySlater
  KellyKapowski
  l_KellyKapowski
  antsMotionCorr
  l_antsMotionCorr
  antsRegistration
  l_antsRegistration
  antsApplyTransforms
  l_antsApplyTransforms
  sccan
  l_sccan
  #  ConvertVectorFieldToVTK
  #  l_ConvertVectorFieldToVTK
  #  GetMeshAndTopology
  #  l_GetMeshAndTopology
  #  CheckTopology
  #  l_CheckTopology
  #  WarpVTKPolyDataMultiTransform
  #  l_WarpVTKPolyDataMultiTransform
  #  ANTSConformalMapping
  #  l_ANTSConformalMapping
  ANTSIntegrateVectorField
  l_ANTSIntegrateVectorField
  ANTSIntegrateVelocityField
  l_ANTSIntegrateVelocityField
  Atropos
  l_Atropos
  CreateDisplacementField
  l_CreateDisplacementField
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(PROGRAMS ../Scripts/ants.sh  ../Scripts/basic_ants_example.sh
                ../Scripts/antsIntroduction.sh
                ../Scripts/antsaffine.sh
        ../Scripts/waitForXGridJobs.pl
                ../Scripts/ANTSAverage2DAffine.sh ../Scripts/ANTSAverage3DAffine.sh
                ../Scripts/buildtemplateparallel.sh ../Scripts/ANTSpexec.sh
                ../Scripts/waitForSGEQJobs.pl DESTINATION bin
                PERMISSIONS  OWNER_WRITE OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
                CONFIGURATIONS  Release )
#Set any libraries that your project depends on.
#examples: ${ITK_LIBRARIES}, VTKRendering, etc
set(Libraries
antsUtilities ${ITK_LIBRARIES}
)
set(TEST_BINARY_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set(DATA_DIR ${CMAKE_CURRENT_LIST_DIR}/Data)
set(R16_IMAGE ${DATA_DIR}/r16slice.nii.gz)
set(R64_IMAGE ${DATA_DIR}/r64slice.nii.gz)

if (0) ## This test is not yet completed.
add_test(antsRegistrationTest1 ${TEST_BINARY_DIR}/antsRegistration -d 2
  -m mi[  ${R16_IMAGE}, ${R64_IMAGE} , 1 , 20, Regular , 0.05 ]
  -t affine[ 1.0 ]
  -i 2100x1200x1200x0
  -s 3x2x1x0
  -f 4x3x2x1
  -m cc[  ${R16_IMAGE}, ${R64_IMAGE} , 1, 4 ]
  -t syn[ 5.0, 3, 0.0 ]
  -i 100x100x50
  -s 2x1x0
  -f 3x2x1
  -u 1
  -o [xtest_CC,xtest_CC.nii.gz,xtest_CC_inv.nii.gz] )
endif()


#set(CurrentExe "GlamorousGlue")
#add_library(l_${CurrentExe} GlamorousGlue.cxx)
#target_link_libraries(l_${CurrentExe} ${Libraries})
##
# Copy all the shell scripts in the Scripts/ directory to the bin directory
##
option( COPY_SCRIPT_FILES_TO_BIN_DIR "Copy the script files to the ANTS bin directory." ON )
if( COPY_SCRIPT_FILES_TO_BIN_DIR )
  file( GLOB SCRIPT_FILES "${CMAKE_SOURCE_DIR}/../Scripts/*.sh" )
  file( GLOB PSCRIPT_FILES "${CMAKE_SOURCE_DIR}/../Scripts/*.pl" )
  foreach( SCRIPT_FILE ${SCRIPT_FILES} )
    add_custom_command( TARGET ANTS POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy ${SCRIPT_FILE} ${CMAKE_BINARY_DIR} )
  endforeach( SCRIPT_FILE )
  foreach( SCRIPT_FILE ${PSCRIPT_FILES} )
    add_custom_command( TARGET ANTS POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy ${SCRIPT_FILE} ${CMAKE_BINARY_DIR} )
  endforeach( SCRIPT_FILE )
endif( COPY_SCRIPT_FILES_TO_BIN_DIR )
set(CPACK_PACKAGE_NAME "ANTs")
set(CPACK_PACKAGE_VENDOR "CMake.org")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ANTs - Advanced Normalization Tools")
set(CPACK_PACKAGE_VERSION "1.9.y")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "9")
set(CPACK_PACKAGE_VERSION_PATCH "y")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "ANTS")
# This must always be last!
include(CPack)
