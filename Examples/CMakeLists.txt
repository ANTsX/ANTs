project(ANTS)
cmake_minimum_required(VERSION 2.8.2)
# set(CMAKE_BUILD_TYPE "Release")
#Change PROJECT_NAME to the name of your project
option(USE_ITK "Use ITK Libraries" ON)
option(USE_VTK "Use VTK Libraries" OFF)
set (CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)
# Set up ITK
if(USE_ITK)
  find_package(ITK)
  if(ITK_FOUND)
    include(${ITK_USE_FILE})
  else(ITK_FOUND)
    message(FATAL_ERROR
      "Cannot build without ITK.  Please set ITK_DIR.")
  endif(ITK_FOUND)
endif(USE_ITK)
#
# now check if we have the minimum required itk version number
#
# refer to this variable from ITK .... set( ITK_VERSION_MAJOR "4")
set( number ${ITK_VERSION_MAJOR}  )
set( operation LESS )
set( limit 4 )
if( number ${operation} ${limit} )
  message(FATAL_ERROR "Oops, you have itk version ${number} ... we require at  least version ${limit}." )
endif( number ${operation} ${limit} )
# Set up VTK
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
)
include_directories(${PICSL_INCLUDE_DIRS})
link_directories( ${ITK_LIBRARY_PATH}  )
# message("${ITK_LIBRARIES}")
# non-templated class -- this should be stored in a library and linked in...
set(UI_SOURCES "../Utilities/antsCommandLineParser" "../Utilities/antsCommandLineOption" )
add_executable(ANTS ANTS.cxx ${UI_SOURCES})
target_link_libraries(ANTS ${ITK_LIBRARIES} )
add_executable(ANTSJacobian ANTSJacobian.cxx ${UI_SOURCES})
target_link_libraries(ANTSJacobian ${ITK_LIBRARIES} )
add_executable(PrintHeader PrintHeader.cxx ${UI_SOURCES})
target_link_libraries(PrintHeader ${ITK_LIBRARIES} )
add_executable(ResetDirection ResetDirection.cxx ${UI_SOURCES})
target_link_libraries(ResetDirection ${ITK_LIBRARIES} )
#add_executable(ResetDirection2 ResetDirection2.cxx ${UI_SOURCES})
#target_link_libraries(ResetDirection2 ${ITK_LIBRARIES} )
#add_executable(BoykovGraphCutFilter BoykovGraphCutFilter.cxx ${UI_SOURCES})
#target_link_libraries(BoykovGraphCutFilter ${ITK_LIBRARIES} )
#add_executable(GeneralLinearModel GeneralLinearModel.cxx ${UI_SOURCES})
#target_link_libraries(GeneralLinearModel ${ITK_LIBRARIES} )
add_executable(ANTSUseLandmarkImagesToGetAffineTransform ANTSUseLandmarkImagesToGetAffineTransform.cxx ${UI_SOURCES})
target_link_libraries(ANTSUseLandmarkImagesToGetAffineTransform ${ITK_LIBRARIES} )
add_executable(ANTSUseDeformationFieldToGetAffineTransform ANTSUseDeformationFieldToGetAffineTransform ${UI_SOURCES})
target_link_libraries(ANTSUseDeformationFieldToGetAffineTransform ${ITK_LIBRARIES} )
#add_executable(FDR FDR.cxx ${UI_SOURCES})
#target_link_libraries(FDR ${ITK_LIBRARIES} )
add_executable(ImageMath ImageMath.cxx ${UI_SOURCES})
target_link_libraries(ImageMath ${ITK_LIBRARIES} )
add_executable(SurfaceCurvature SurfaceCurvature.cxx ${UI_SOURCES})
target_link_libraries(SurfaceCurvature ${ITK_LIBRARIES} )
add_executable(SurfaceBasedSmoothing SurfaceBasedSmoothing.cxx ${UI_SOURCES})
target_link_libraries(SurfaceBasedSmoothing ${ITK_LIBRARIES} )
# add_executable(test test.cxx)
#target_link_libraries(test)
#add_executable(WarpImage WarpImage.cxx ${UI_SOURCES})
#target_link_libraries(WarpImage ${ITK_LIBRARIES} )
add_executable(LaplacianThickness LaplacianThickness.cxx ${UI_SOURCES})
target_link_libraries(LaplacianThickness ${ITK_LIBRARIES} )
add_executable(SetOrigin SetOrigin.cxx ${UI_SOURCES})
target_link_libraries(SetOrigin ${ITK_LIBRARIES} )
add_executable(SetSpacing SetSpacing.cxx ${UI_SOURCES})
target_link_libraries(SetSpacing ${ITK_LIBRARIES} )
add_executable(SetDirectionByMatrix SetDirectionByMatrix.cxx ${UI_SOURCES})
target_link_libraries(SetDirectionByMatrix ${ITK_LIBRARIES} )
#add_executable(WarpFunctionalImage WarpFunctionalImage.cxx ${UI_SOURCES})
#target_link_libraries(WarpFunctionalImage ${ITK_LIBRARIES} )
add_executable(ConvertScalarImageToRGB ConvertScalarImageToRGB.cxx)
target_link_libraries(ConvertScalarImageToRGB ${ITK_LIBRARIES})
add_executable(CreateWarpedGridImage CreateWarpedGridImage.cxx)
target_link_libraries(CreateWarpedGridImage ${ITK_LIBRARIES})
# add_executable(CreateJacobianDeterminantImage CreateJacobianDeterminantImage.cxx)
# target_link_libraries(CreateJacobianDeterminantImage ${ITK_LIBRARIES})
add_executable(MeasureImageSimilarity MeasureImageSimilarity.cxx)
target_link_libraries(MeasureImageSimilarity ${ITK_LIBRARIES})
add_executable(ConvertToJpg ConvertToJpg.cxx)
target_link_libraries(ConvertToJpg ${ITK_LIBRARIES} )
add_executable(ConvertImagePixelType ConvertImagePixelType.cxx ${UI_SOURCES})
target_link_libraries(ConvertImagePixelType ${ITK_LIBRARIES}    )
add_executable(AverageImages AverageImages.cxx ${UI_SOURCES})
target_link_libraries(AverageImages ${ITK_LIBRARIES} )
add_executable(AverageTensorImages AverageTensorImages.cxx ${UI_SOURCES})
target_link_libraries(AverageTensorImages ${ITK_LIBRARIES} )
add_executable(ImageSetStatistics ImageSetStatistics.cxx ${UI_SOURCES})
target_link_libraries(ImageSetStatistics ${ITK_LIBRARIES} )
add_executable(ThresholdImage ThresholdImage.cxx ${UI_SOURCES})
target_link_libraries(ThresholdImage ${ITK_LIBRARIES} )
add_executable(MultiplyImages MultiplyImages.cxx ${UI_SOURCES})
target_link_libraries(MultiplyImages ${ITK_LIBRARIES} )
add_executable(SmoothImage SmoothImage.cxx ${UI_SOURCES})
target_link_libraries(SmoothImage ${ITK_LIBRARIES} )
add_executable(ClusterImageStatistics ClusterImageStatistics.cxx ${UI_SOURCES})
target_link_libraries(ClusterImageStatistics ${ITK_LIBRARIES} )
add_executable(LabelClustersUniquely LabelClustersUniquely.cxx ${UI_SOURCES})
target_link_libraries(LabelClustersUniquely ${ITK_LIBRARIES} )
#add_executable(StudentsTestOnImages StudentsTestOnImages.cxx ${UI_SOURCES})
#target_link_libraries(StudentsTestOnImages ${ITK_LIBRARIES} )
add_executable(LabelOverlapMeasures LabelOverlapMeasures.cxx )
target_link_libraries(LabelOverlapMeasures ${ITK_LIBRARIES})
add_executable(MeasureMinMaxMean MeasureMinMaxMean.cxx ${UI_SOURCES})
target_link_libraries(MeasureMinMaxMean ${ITK_LIBRARIES} )
add_executable(WarpImageMultiTransform WarpImageMultiTransform ${UI_SOURCES})
target_link_libraries(WarpImageMultiTransform ${ITK_LIBRARIES} )
add_executable(ComposeMultiTransform ComposeMultiTransform ${UI_SOURCES})
target_link_libraries(ComposeMultiTransform ${ITK_LIBRARIES} )
add_executable(StackSlices StackSlices.cxx ${UI_SOURCES})
target_link_libraries(StackSlices ${ITK_LIBRARIES} )
add_executable(MemoryTest MemoryTest.cxx ${UI_SOURCES})
target_link_libraries(MemoryTest ${ITK_LIBRARIES} )
#add_executable(ANTSOrientImage ANTSOrientImage.cxx ${UI_SOURCES})
#target_link_libraries(ANTSOrientImage ${ITK_LIBRARIES} )
add_executable(PermuteFlipImageOrientationAxes PermuteFlipImageOrientationAxes.cxx ${UI_SOURCES})
target_link_libraries(PermuteFlipImageOrientationAxes ${ITK_LIBRARIES} )
add_executable(ImageCompare ImageCompare.cxx)
target_link_libraries(ImageCompare ${ITK_LIBRARIES})
add_executable(ResampleImageBySpacing ResampleImageBySpacing.cxx)
target_link_libraries(ResampleImageBySpacing ${ITK_LIBRARIES})
#add_executable(ResampleImageByOtherImage ResampleImageByOtherImage.cxx)
#target_link_libraries(ResampleImageByOtherImage ${ITK_LIBRARIES})
add_executable(CopyImageHeaderInformation CopyImageHeaderInformation.cxx)
target_link_libraries(CopyImageHeaderInformation ${ITK_LIBRARIES})
add_executable(WarpTimeSeriesImageMultiTransform WarpTimeSeriesImageMultiTransform.cxx)
target_link_libraries(WarpTimeSeriesImageMultiTransform ${ITK_LIBRARIES})
add_executable(ResampleImage ResampleImage.cxx )
target_link_libraries(ResampleImage ${ITK_LIBRARIES})
add_executable(ExtractSliceFromImage ExtractSliceFromImage.cxx )
target_link_libraries(ExtractSliceFromImage ${ITK_LIBRARIES})
add_executable(ExtractRegionFromImage ExtractRegionFromImage.cxx )
target_link_libraries(ExtractRegionFromImage ${ITK_LIBRARIES})
add_executable(ExtractRegionFromImageByMask ExtractRegionFromImageByMask.cxx )
target_link_libraries(ExtractRegionFromImageByMask ${ITK_LIBRARIES})
add_executable(LabelGeometryMeasures LabelGeometryMeasures.cxx )
target_link_libraries(LabelGeometryMeasures ${ITK_LIBRARIES})
add_executable(PasteImageIntoImage PasteImageIntoImage.cxx )
target_link_libraries(PasteImageIntoImage ${ITK_LIBRARIES})
add_executable(TileImages TileImages.cxx )
target_link_libraries(TileImages ${ITK_LIBRARIES})
add_executable(CreateImage CreateImage.cxx )
target_link_libraries(CreateImage ${ITK_LIBRARIES})
add_executable(WarpTensorImageMultiTransform WarpTensorImageMultiTransform.cxx)
target_link_libraries(WarpTensorImageMultiTransform ${ITK_LIBRARIES})
add_executable(ReorientTensorImage ReorientTensorImage.cxx)
target_link_libraries(ReorientTensorImage ${ITK_LIBRARIES})
add_executable(RebaseTensorImage RebaseTensorImage.cxx)
target_link_libraries(RebaseTensorImage ${ITK_LIBRARIES})
add_executable(N3BiasFieldCorrection N3BiasFieldCorrection.cxx)
target_link_libraries(N3BiasFieldCorrection ${ITK_LIBRARIES})
add_executable(N4BiasFieldCorrection N4BiasFieldCorrection.cxx  ${UI_SOURCES})
target_link_libraries(N4BiasFieldCorrection ${ITK_LIBRARIES})
add_executable(KellySlater KellySlater.cxx)
target_link_libraries(KellySlater ${ITK_LIBRARIES})
add_executable(KellyKapowski KellyKapowski.cxx  ${UI_SOURCES})
target_link_libraries(KellyKapowski ${ITK_LIBRARIES})
add_executable(CreateDTICohort CreateDTICohort.cxx  ${UI_SOURCES})
target_link_libraries(CreateDTICohort ${ITK_LIBRARIES})
add_executable(AverageAffineTransform AverageAffineTransform.cxx)
 target_link_libraries(AverageAffineTransform ${ITK_LIBRARIES})
add_executable(antsRegistration antsRegistration.cxx ${UI_SOURCES})
 target_link_libraries(antsRegistration ${ITK_LIBRARIES})
add_executable(antsApplyTransforms antsApplyTransforms.cxx ${UI_SOURCES})
 target_link_libraries(antsApplyTransforms ${ITK_LIBRARIES})
add_executable(CompositeTransformUtil CompositeTransformUtil.cxx)
 target_link_libraries(CompositeTransformUtil ${ITK_LIBRARIES})

add_executable(antsMotionCorr antsMotionCorr.cxx ${UI_SOURCES})
target_link_libraries(antsMotionCorr ${ITK_LIBRARIES})
add_executable(sccan sccan.cxx ${UI_SOURCES})
target_link_libraries(sccan ${ITK_LIBRARIES} )
add_executable(CreateDisplacementField CreateDisplacementField.cxx ${UI_SOURCES})
target_link_libraries(CreateDisplacementField ${ITK_LIBRARIES} )
if(USE_VTK)
include(${CMAKE_ROOT}/Modules/FindVTK.cmake)
if(USE_VTK_FILE)
  include(${USE_VTK_FILE})
  add_executable(ConvertVectorFieldToVTK ConvertVectorFieldToVTK.cxx )
  target_link_libraries(ConvertVectorFieldToVTK ${ITK_LIBRARIES} vtkIO)
  add_executable(GetMeshAndTopology GetMeshAndTopology.cxx)
  target_link_libraries(GetMeshAndTopology ${ITK_LIBRARIES} ITKMetaIO  vtkRendering vtksys vtkIO)
  add_executable(CheckTopology CheckTopology.cxx)
  target_link_libraries(CheckTopology ${ITK_LIBRARIES} ITKMetaIO  vtkRendering vtksys vtkIO)
add_executable(WarpVTKPolyDataMultiTransform WarpVTKPolyDataMultiTransform)
target_link_libraries(WarpVTKPolyDataMultiTransform ${ITK_LIBRARIES} vtksys vtkIO )
#   set(FLAT_SRC
#     ../Temporary/itkFEMElement3DC0LinearTriangular.cxx
#     ../Temporary/itkFEMElement3DC0LinearTriangularLaplaceBeltrami.cxx
#     ../Temporary/itkFEMElement3DC0LinearTriangularMembrane.cxx
#   )
#   # Build the library
#   add_library(FLATFEM ${FLAT_SRC})
#   target_link_libraries(FLATFEM ${ITK_LIBRARIES})
#
#
#   add_executable(ConformalMapping ConformalMapping.cxx)
#   target_link_libraries(ConformalMapping ${ITK_LIBRARIES} FLATFEM ${ITK_LIBRARIES} ITKMetaIO  vtkRendering vtkIO)
#
#  add_executable(ANTSConformalMapping ANTSConformalMapping.cxx ${UI_SOURCES})
#   target_link_libraries(ANTSConformalMapping ${ITK_LIBRARIES} FLATFEM ${ITK_LIBRARIES} ITKMetaIO  vtkRendering vtkIO)
endif(USE_VTK_FILE)
endif(USE_VTK)
add_executable(ANTSIntegrateVectorField ANTSIntegrateVectorField.cxx)
target_link_libraries(ANTSIntegrateVectorField ${ITK_LIBRARIES})
add_executable(ANTSIntegrateVelocityField ANTSIntegrateVelocityField.cxx)
target_link_libraries(ANTSIntegrateVelocityField ${ITK_LIBRARIES})
  add_executable(AtroposMin AtroposMin.cxx ${UI_SOURCES})
  target_link_libraries(AtroposMin ${ITK_LIBRARIES}   )
  add_executable(Atropos Atropos.cxx ${UI_SOURCES})
  target_link_libraries(Atropos ${ITK_LIBRARIES}   )
install(TARGETS ANTS
  AverageAffineTransform
  ANTSJacobian
  PasteImageIntoImage
  PrintHeader
  ResetDirection
  ANTSUseLandmarkImagesToGetAffineTransform
  ANTSUseDeformationFieldToGetAffineTransform
  AverageAffineTransform
  ImageMath
  SurfaceCurvature
  SurfaceBasedSmoothing
  LaplacianThickness
  SetOrigin
  SetSpacing
  SetDirectionByMatrix
  ConvertScalarImageToRGB
  CreateWarpedGridImage
  MeasureImageSimilarity
  ConvertToJpg
  ConvertImagePixelType
  AverageImages
  AverageTensorImages
  ImageSetStatistics
  ThresholdImage
  MultiplyImages
  SmoothImage
  ClusterImageStatistics
  LabelClustersUniquely
  LabelOverlapMeasures
  MeasureMinMaxMean
  WarpImageMultiTransform
  ComposeMultiTransform
  StackSlices
  PermuteFlipImageOrientationAxes
  ImageCompare
  ResampleImageBySpacing
  CopyImageHeaderInformation
  WarpTimeSeriesImageMultiTransform
  ResampleImage
  ExtractSliceFromImage
  ExtractRegionFromImage
  ExtractRegionFromImageByMask
  LabelGeometryMeasures
  TileImages
  CreateImage
  WarpTensorImageMultiTransform
  ReorientTensorImage
  RebaseTensorImage
  N3BiasFieldCorrection
  N4BiasFieldCorrection
  KellySlater
  KellyKapowski
  antsMotionCorr
  antsRegistration
  antsApplyTransforms
  sccan
#  ConvertVectorFieldToVTK
#  GetMeshAndTopology
#  CheckTopology
#  WarpVTKPolyDataMultiTransform
#  ANTSConformalMapping
  ANTSIntegrateVectorField
  ANTSIntegrateVelocityField
  Atropos
  CreateDisplacementField
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
${ITK_LIBRARIES}
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
#add_executable(${CurrentExe} GlamorousGlue.cxx)
#target_link_libraries(${CurrentExe} ${Libraries})
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
