project(ANTS)

cmake_minimum_required(VERSION 2.4)

if(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

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

# Set up VTK
if(USE_VTK)
  find_package(VTK)
  if(VTK_FOUND)
    include(${VTK_USE_FILE})
  else(VTK_FOUND)
     message("Cannot build some programs without VTK.  Please set VTK_DIR if you need these programs.")
  endif(VTK_FOUND)
endif(USE_VTK)

#The following lines are required to use Dart
include(CTest)
enable_testing()


#INCLUDE (${CMAKE_ROOT}/Modules/FindITK.cmake)
#IF (USE_ITK_FILE)
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

LINK_DIRECTORIES( ${ITK_LIBRARY_PATH} )

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

add_executable(TileImages TileImages.cxx )
target_link_libraries(TileImages ${ITK_LIBRARIES})


add_executable(WarpTensorImageMultiTransform WarpTensorImageMultiTransform.cxx)
target_link_libraries(WarpTensorImageMultiTransform ${ITK_LIBRARIES})

add_executable(ReorientTensorImage ReorientTensorImage.cxx)
target_link_libraries(ReorientTensorImage ${ITK_LIBRARIES})

add_executable(N3BiasFieldCorrection N3BiasFieldCorrection.cxx)
target_link_libraries(N3BiasFieldCorrection ${ITK_LIBRARIES})

add_executable(N4BiasFieldCorrection N4BiasFieldCorrection.cxx  ${UI_SOURCES})
target_link_libraries(N4BiasFieldCorrection ${ITK_LIBRARIES})

add_executable(KellySlater KellySlater.cxx)
target_link_libraries(KellySlater ${ITK_LIBRARIES})

IF (USE_VTK)
INCLUDE (${CMAKE_ROOT}/Modules/FindVTK.cmake)
IF (USE_VTK_FILE)
  include(${USE_VTK_FILE})
  add_executable(ConvertVectorFieldToVTK ConvertVectorFieldToVTK.cxx )
  target_link_libraries(ConvertVectorFieldToVTK ${ITK_LIBRARIES} vtkIO)

  add_executable(GetMeshAndTopology GetMeshAndTopology.cxx)
  target_link_libraries(GetMeshAndTopology ${ITK_LIBRARIES} ITKMetaIO  vtkRendering vtksys vtkIO)

  add_executable(CheckTopology CheckTopology.cxx)
  target_link_libraries(CheckTopology ${ITK_LIBRARIES} ITKMetaIO  vtkRendering vtksys vtkIO)

add_executable(WarpVTKPolyDataMultiTransform WarpVTKPolyDataMultiTransform)
target_link_libraries(WarpVTKPolyDataMultiTransform ${ITK_LIBRARIES} vtksys vtkIO )


  set(FLAT_SRC
    ../Temporary/itkFEMElement3DC0LinearTriangular.cxx
    ../Temporary/itkFEMElement3DC0LinearTriangularLaplaceBeltrami.cxx
    ../Temporary/itkFEMElement3DC0LinearTriangularMembrane.cxx
  )
  # Build the library
  ADD_LIBRARY(FLATFEM ${FLAT_SRC})
  target_link_libraries(FLATFEM ${ITK_LIBRARIES})


  add_executable(ConformalMapping ConformalMapping.cxx)
  target_link_libraries(ConformalMapping ${ITK_LIBRARIES} FLATFEM ${ITK_LIBRARIES} ITKMetaIO  vtkRendering vtkIO)

 add_executable(ANTSConformalMapping ANTSConformalMapping.cxx ${UI_SOURCES})
  target_link_libraries(ANTSConformalMapping ${ITK_LIBRARIES} FLATFEM ${ITK_LIBRARIES} ITKMetaIO  vtkRendering vtkIO)

ENDIF (USE_VTK_FILE)
endif(USE_VTK)



add_executable(ANTSIntegrateVectorField ANTSIntegrateVectorField.cxx)
target_link_libraries(ANTSIntegrateVectorField ${ITK_LIBRARIES})

if( ITK_USE_REVIEW )
  add_executable(Atropos Atropos.cxx ${UI_SOURCES})
  target_link_libraries(Atropos ${ITK_LIBRARIES}   )
endif( ITK_USE_REVIEW )


#Set any libraries that your project depends on.
#examples: ${ITK_LIBRARIES}, VTKRendering, etc
set(Libraries
${ITK_LIBRARIES}
)

#set(CurrentExe "GlamorousGlue")
#add_executable(${CurrentExe} GlamorousGlue.cxx)
#target_link_libraries(${CurrentExe} ${Libraries})

##
# Copy all the shell scripts in the Scripts/ directory to the bin directory
##
option( COPY_SCRIPT_FILES_TO_BIN_DIR "Copy the script files to the ANTS bin directory." ON )

if( COPY_SCRIPT_FILES_TO_BIN_DIR )
  FILE( GLOB SCRIPT_FILES "${CMAKE_SOURCE_DIR}/../Scripts/*.sh" )
  FILE( GLOB PSCRIPT_FILES "${CMAKE_SOURCE_DIR}/../Scripts/*.pl" )

  FOREACH( SCRIPT_FILE ${SCRIPT_FILES} )
    ADD_CUSTOM_COMMAND( TARGET ANTS POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy ${SCRIPT_FILE} ${CMAKE_BINARY_DIR} )
  ENDFOREACH( SCRIPT_FILE )
  FOREACH( SCRIPT_FILE ${PSCRIPT_FILES} )
    ADD_CUSTOM_COMMAND( TARGET ANTS POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy ${SCRIPT_FILE} ${CMAKE_BINARY_DIR} )
  ENDFOREACH( SCRIPT_FILE )
endif( COPY_SCRIPT_FILES_TO_BIN_DIR )

###
#  Perform testing
###

set(DATA_DIR ${CMAKE_SOURCE_DIR}/Data)
set(R16_IMAGE ${DATA_DIR}/r16slice.nii)
set(R64_IMAGE ${DATA_DIR}/r64slice.nii)
set(OUTPUT_PREFIX ${CMAKE_BINARY_DIR}/TEST)
set(WARP ${OUTPUT_PREFIX}Warp.nii.gz ${OUTPUT_PREFIX}Affine.txt )
set(INVERSEWARP -i ${OUTPUT_PREFIX}Affine.txt ${OUTPUT_PREFIX}InverseWarp.nii.gz )
set(WARP_IMAGE ${CMAKE_BINARY_DIR}/warped.nii.gz)
set(INVERSEWARP_IMAGE ${CMAKE_BINARY_DIR}/inversewarped.nii.gz)

set(DEVIL_IMAGE ${DATA_DIR}/Frown.nii)
set(ANGEL_IMAGE ${DATA_DIR}/Smile.nii)

set(SEG_IMAGE ${DATA_DIR}/nslice.nii.gz)



###
#  ANTS metric testing
###
add_test(ANTS_CC_1   ANTS 2 -m  PR[ ${R16_IMAGE}, ${R64_IMAGE},  1 ,2 ] -r Gauss[ 3 , 0 ] -t SyN[ 0.5   ] -i 50x50x30 -o ${OUTPUT_PREFIX}.nii.gz)
add_test(ANTS_CC_1_WARP WarpImageMultiTransform 2 ${R64_IMAGE} ${WARP_IMAGE} ${WARP} -R ${R16_IMAGE} )
add_test(ANTS_CC_1_JPG  ConvertToJpg ${WARP_IMAGE} ANTSCC1.jpg)
add_test(ANTS_CC_1_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 11.9992 0.05)
add_test(ANTS_CC_1_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.61 0.05)
add_test(ANTS_CC_1_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.00038593 0.05)
add_test(ANTS_CC_1_INVERSEWARP WarpImageMultiTransform 2 ${R16_IMAGE} ${INVERSEWARP_IMAGE} ${INVERSEWARP}  -R ${R16_IMAGE} )
add_test(ANTS_CC_1_INVERSEWARP_METRIC_0 MeasureImageSimilarity 2 0 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12.1606 0.05)
add_test(ANTS_CC_1_INVERSEWARP_METRIC_1 MeasureImageSimilarity 2 1 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.61 0.05)
add_test(ANTS_CC_1_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000380545 0.05)

add_test(ANTS_CC_2   ANTS 2 -m  PR[${R16_IMAGE},${R64_IMAGE},1,4] -r Gauss[3,0] -t SyN[0.5] -i 50x50x30 -o ${OUTPUT_PREFIX}.nii.gz --go-faster true )
add_test(ANTS_CC_2_WARP WarpImageMultiTransform 2 ${R64_IMAGE} ${WARP_IMAGE} ${WARP}  -R ${R16_IMAGE} )
add_test(ANTS_CC_2_JPG  ConvertToJpg ${WARP_IMAGE} ANTSCC2.jpg)
add_test(ANTS_CC_2_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 11.7083 0.05)
add_test(ANTS_CC_2_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.62 0.05)
add_test(ANTS_CC_2_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000461792 0.05)
add_test(ANTS_CC_2_INVERSEWARP WarpImageMultiTransform 2 ${R16_IMAGE} ${INVERSEWARP_IMAGE} ${INVERSEWARP}   -R ${R16_IMAGE}  )
add_test(ANTS_CC_2_INVERSEWARP_METRIC_0 MeasureImageSimilarity 2 0 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 11.912 0.05)
add_test(ANTS_CC_2_INVERSEWARP_METRIC_1 MeasureImageSimilarity 2 1 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.62 0.05)
add_test(ANTS_CC_2_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000442839 0.05)
#add_test(ANTS_CC_2_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -5 0.05)

add_test(ANTS_CC_3   ANTS 2 -m  CC[${R16_IMAGE},${R64_IMAGE},1,4] -r Gauss[3,0] -t SyN[0.5] -i 50x50x30 -o ${OUTPUT_PREFIX}.nii.gz  --go-faster true )
add_test(ANTS_CC_3_WARP WarpImageMultiTransform 2 ${R64_IMAGE} ${WARP_IMAGE} ${WARP}  -R ${R16_IMAGE} )
add_test(ANTS_CC_3_JPG  ConvertToJpg ${WARP_IMAGE} ANTSCC2.jpg)
add_test(ANTS_CC_3_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 11.7083 0.05)
add_test(ANTS_CC_3_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.62 0.05)
add_test(ANTS_CC_3_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000461792 0.05)
add_test(ANTS_CC_3_INVERSEWARP WarpImageMultiTransform 2 ${R16_IMAGE} ${INVERSEWARP_IMAGE} ${INVERSEWARP}   -R ${R16_IMAGE} )
add_test(ANTS_CC_3_INVERSEWARP_METRIC_0 MeasureImageSimilarity 2 0 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 11.912 0.05)
add_test(ANTS_CC_3_INVERSEWARP_METRIC_1 MeasureImageSimilarity 2 1 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.62 0.05)
add_test(ANTS_CC_3_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000442839 0.05)
#add_test(ANTS_CC_3_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -5 0.05)

add_test(ANTS_MSQ    ANTS 2 -m MSQ[${R16_IMAGE},${R64_IMAGE},1,0] -r Gauss[3,0] -t SyN[0.25] -i 50x50x30 -o ${OUTPUT_PREFIX}.nii.gz)
add_test(ANTS_MSQ_WARP WarpImageMultiTransform 2 ${R64_IMAGE} ${WARP_IMAGE} ${WARP}  -R ${R16_IMAGE} )
add_test(ANTS_MSQ_JPG  ConvertToJpg ${WARP_IMAGE} ANTSMSQ.jpg)
add_test(ANTS_MSQ_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 11.7416 0.05)
add_test(ANTS_MSQ_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.6 0.05)
add_test(ANTS_MSQ_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -1.2 0.05)
add_test(ANTS_MSQ_INVERSEWARP WarpImageMultiTransform 2 ${R16_IMAGE} ${INVERSEWARP_IMAGE} ${INVERSEWARP}   -R ${R16_IMAGE} )
add_test(ANTS_MSQ_INVERSEWARP_METRIC_0 MeasureImageSimilarity 2 0 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 11.7845 0.05)
add_test(ANTS_MSQ_INVERSEWARP_METRIC_1 MeasureImageSimilarity 2 1 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.6 0.05)
add_test(ANTS_MSQ_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -1.2 0.05)

add_test(ANTS_MI_1   ANTS 2 -m  MI[${R16_IMAGE},${R64_IMAGE},1,32] -r Gauss[3,0] -t SyN[0.25] -i 50x50x30 -o ${OUTPUT_PREFIX}.nii.gz)
add_test(ANTS_MI_1_WARP WarpImageMultiTransform 2 ${R64_IMAGE} ${WARP_IMAGE} ${WARP}  -R ${R16_IMAGE} )
add_test(ANTS_MI_1_JPG  ConvertToJpg ${WARP_IMAGE} ANTSMI1.jpg)
add_test(ANTS_MI_1_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12.4 0.05)
add_test(ANTS_MI_1_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.45 0.05)
add_test(ANTS_MI_1_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000370686 0.05)
add_test(ANTS_MI_1_INVERSEWARP WarpImageMultiTransform 2 ${R16_IMAGE} ${INVERSEWARP_IMAGE} ${INVERSEWARP}  -R ${R16_IMAGE} )
add_test(ANTS_MI_1_INVERSEWARP_METRIC_0 MeasureImageSimilarity 2 0 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12.6 0.05)
add_test(ANTS_MI_1_INVERSEWARP_METRIC_1 MeasureImageSimilarity 2 1 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.45 0.05)
add_test(ANTS_MI_1_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000366879 0.05)

add_test(ANTS_MI_2   ANTS 2 -m  MI[${R16_IMAGE},${R64_IMAGE},1,24] -r Gauss[3,0] -t SyN[0.25] -i 50x50x20 -o ${OUTPUT_PREFIX}.nii.gz)
add_test(ANTS_MI_2_WARP WarpImageMultiTransform 2 ${R64_IMAGE} ${WARP_IMAGE} ${WARP}  -R ${R16_IMAGE} )
add_test(ANTS_MI_2_JPG  ConvertToJpg ${WARP_IMAGE} ANTSMI2.jpg)
add_test(ANTS_MI_2_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12.4 0.05)
add_test(ANTS_MI_2_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.45 0.05)
add_test(ANTS_MI_2_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000370686 0.05)
add_test(ANTS_MI_2_INVERSEWARP WarpImageMultiTransform 2 ${R16_IMAGE} ${INVERSEWARP_IMAGE} ${INVERSEWARP}  -R ${R16_IMAGE}  )
add_test(ANTS_MI_2_INVERSEWARP_METRIC_0 MeasureImageSimilarity 2 0 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12.6 0.05)
add_test(ANTS_MI_2_INVERSEWARP_METRIC_1 MeasureImageSimilarity 2 1 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.45 0.05)
add_test(ANTS_MI_2_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000366879 0.05)



###
#  ANTS transform testing
###
add_test(ANTS_ELASTIC ANTS 2 -m PR[${R16_IMAGE},${R64_IMAGE},1,2] -t Elast[1]       -i 50x50x50 -r Gauss[0,1] -o ${OUTPUT_PREFIX}.nii.gz)
add_test(ANTS_ELASTIC_WARP WarpImageMultiTransform 2 ${R64_IMAGE} ${WARP_IMAGE} ${WARP}   -R ${R16_IMAGE}  )
add_test(ANTS_ELASTIC_JPG  ConvertToJpg ${WARP_IMAGE} ANTSELASTIC.jpg)
add_test(ANTS_ELASTIC_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12.5 0.05)
add_test(ANTS_ELASTIC_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.4 0.05)
add_test(ANTS_ELASTIC_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000459869 0.05)

add_test(ANTS_GSYN    ANTS 2 -m CC[${R16_IMAGE},${R64_IMAGE},1,3] -t SyN[0.25]       -i 50x50x50 -r Gauss[3,0.] -o ${OUTPUT_PREFIX}.nii.gz)
add_test(ANTS_GSYN_WARP WarpImageMultiTransform 2 ${R64_IMAGE} ${WARP_IMAGE} ${WARP}  -R ${R16_IMAGE} )
add_test(ANTS_GSYN_JPG  ConvertToJpg ${WARP_IMAGE} ANTSGSYN.jpg)
add_test(ANTS_GSYN_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 11.7734 0.05)
add_test(ANTS_GSYN_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.6 0.05)
add_test(ANTS_GSYN_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000478672 0.05)
add_test(ANTS_GSYN_INVERSEWARP WarpImageMultiTransform 2 ${R16_IMAGE} ${INVERSEWARP_IMAGE} ${INVERSEWARP}  -R ${R16_IMAGE}  )
add_test(ANTS_GSYN_JPGINV  ConvertToJpg ${INVERSEWARP_IMAGE} ANTSGSYNINV.jpg)
add_test(ANTS_GSYN_INVERSEWARP_METRIC_0 MeasureImageSimilarity 2 0 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12.0541 0.05)
add_test(ANTS_GSYN_INVERSEWARP_METRIC_1 MeasureImageSimilarity 2 1 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.6 0.05)
add_test(ANTS_GSYN_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000475175 0.05)

add_test(ANTS_EXP     ANTS 2 -m PR[${R16_IMAGE},${R64_IMAGE},1,4] -t Exp[0.5,2,0.5]       -i 50x50x50 -r Gauss[0.5,0.25] -o ${OUTPUT_PREFIX}.nii.gz)
add_test(ANTS_EXP_WARP WarpImageMultiTransform 2 ${R64_IMAGE} ${WARP_IMAGE} ${WARP}  -R ${R16_IMAGE} )
add_test(ANTS_EXP_JPG  ConvertToJpg ${WARP_IMAGE} ANTSEXP.jpg)
add_test(ANTS_EXP_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12 0.05)
add_test(ANTS_EXP_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.5 0.05)
add_test(ANTS_EXP_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000423693 0.05)
add_test(ANTS_EXP_INVERSEWARP WarpImageMultiTransform 2 ${R16_IMAGE} ${INVERSEWARP_IMAGE} ${INVERSEWARP}   -R ${R16_IMAGE} )
add_test(ANTS_EXP_INVERSEWARP_METRIC_0 MeasureImageSimilarity 2 0 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12.349 0.05)
add_test(ANTS_EXP_INVERSEWARP_METRIC_1 MeasureImageSimilarity 2 1 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.5 0.05)
add_test(ANTS_EXP_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000339606 0.05)

#add_test(ANTS_GSYN    ANTS 2 -m PR[${R16_IMAGE},${R64_IMAGE},1,2] -t SyN[0.75]            -i 50x50x50 -r Gauss[3,0.0,32] -o ${OUTPUT_PREFIX}.nii.gz)
 add_test(ANTS_SYN     ANTS 2 -m PR[${R16_IMAGE},${R64_IMAGE},1,2] -t SyN[0.5,2,0.05] -i 50x50x50 -r Gauss[3,0.0,32] -o ${OUTPUT_PREFIX}.nii.gz)
add_test(ANTS_SYN_WARP WarpImageMultiTransform 2 ${R64_IMAGE} ${WARP_IMAGE} ${WARP}  -R ${R16_IMAGE}  )
add_test(ANTS_SYN_JPG  ConvertToJpg ${WARP_IMAGE} ANTSSYN.jpg)
add_test(ANTS_SYN_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12.0239 0.05)
add_test(ANTS_SYN_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.6 0.05)
add_test(ANTS_SYN_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000461922 0.05)
add_test(ANTS_SYN_INVERSEWARP WarpImageMultiTransform 2 ${R16_IMAGE} ${INVERSEWARP_IMAGE} ${INVERSEWARP}  -R ${R16_IMAGE}  )
add_test(ANTS_SYN_JPGINV  ConvertToJpg ${INVERSEWARP_IMAGE} ANTSSYNINV.jpg)
add_test(ANTS_SYN_INVERSEWARP_METRIC_0 MeasureImageSimilarity 2 0 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12.5104 0.05)
add_test(ANTS_SYN_INVERSEWARP_METRIC_1 MeasureImageSimilarity 2 1 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.6 0.05)
add_test(ANTS_SYN_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000444279 0.05)

###
# PSE sub-tests:  Check to see if .txt files and .vtk files also run correctly
###

set(ANGEL_IMAGE_TXT ${DATA_DIR}/Smile.txt)
set(DEVIL_IMAGE_TXT ${DATA_DIR}/Frown.txt)
set(ANGEL_IMAGE_VTK ${DATA_DIR}/Smile.vtk)
set(DEVIL_IMAGE_VTK ${DATA_DIR}/Frown.vtk)

add_test(ANTS_PSE_MSQ_TXT ANTS 2 -i 91x70x55x40x30  -r Gauss[3,0.,32] -t SyN[0.25]
                                 -m MSQ[${DEVIL_IMAGE},${ANGEL_IMAGE},1,0]
                                 -m PSE[${DEVIL_IMAGE},${ANGEL_IMAGE},${DEVIL_IMAGE_TXT},${ANGEL_IMAGE_TXT},1,0.33,11,1,25]
                                 --continue-affine 0 --number-of-affine-iterations 0
                                 -o ${OUTPUT_PREFIX}.nii.gz)
add_test(ANTS_PSE_MSQ_VTK ANTS 2 -i 91x70x55x40x30  -r Gauss[3,0.,32] -t SyN[0.25]
                                 -m MSQ[${DEVIL_IMAGE},${ANGEL_IMAGE},1,0]
                                 -m PSE[${DEVIL_IMAGE},${ANGEL_IMAGE},${DEVIL_IMAGE_VTK},${ANGEL_IMAGE_VTK},1,0.33,11,1,25]
                                 --continue-affine 0 --number-of-affine-iterations 0
                                 -o ${OUTPUT_PREFIX}.nii.gz)

###
#  ANTS labeled data testing
###
option(RUN_LONG_TESTS "Run the time consuming tests." OFF )
if(RUN_LONG_TESTS)

add_test(ANTS_PSE_MSQ_IMG ANTS 2 -i 191x170x90x90x10  -r Gauss[6,0.25] -t SyN[1,2,0.1]
                             -m PSE[${DEVIL_IMAGE},${ANGEL_IMAGE},${DEVIL_IMAGE},${ANGEL_IMAGE},0.25,0.1,100,0,10]
                             -m MSQ[${DEVIL_IMAGE},${ANGEL_IMAGE},1,0.1]
                              --continue-affine 0 --number-of-affine-iterations 0 --geodesic 2
                             -o ${OUTPUT_PREFIX}.nii.gz)
add_test(ANTS_PSE_MSQ_IMG_WARP WarpImageMultiTransform 2 ${ANGEL_IMAGE} ${WARP_IMAGE} ${WARP} -R  ${DEVIL_IMAGE} )
add_test(ANTS_PSE_JPG  ConvertToJpg ${WARP_IMAGE} ANTSPSE.jpg)
add_test(ANTS_PSE_MSQ_IMG_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${DEVIL_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 0.0116083 0.05)
add_test(ANTS_PSE_MSQ_IMG_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${DEVIL_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.884041 0.05)
add_test(ANTS_PSE_MSQ_IMG_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${DEVIL_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000710551 0.05)
add_test(ANTS_PSE_MSQ_IMG_INVERSEWARP WarpImageMultiTransform 2 ${DEVIL_IMAGE} ${INVERSEWARP_IMAGE} ${INVERSEWARP} -R  ${ANGEL_IMAGE}  )
add_test(ANTS_PSE_JPGINV  ConvertToJpg ${INVERSEWARP_IMAGE} ANTSPSEINV.jpg)
add_test(ANTS_PSE_MSQ_IMG_INVERSEWARP_METRIC_0 MeasureImageSimilarity 2 0 ${ANGEL_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 0.0109054 0.05)
add_test(ANTS_PSE_MSQ_IMG_INVERSEWARP_METRIC_1 MeasureImageSimilarity 2 1 ${ANGEL_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.901632 0.05)
add_test(ANTS_PSE_MSQ_IMG_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${ANGEL_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000704717 0.05)


###
#  ANTS images with non-trival rotation header test
###
#set(ROT_REF_IMAGE ${DATA_DIR}/ref2.nii.gz)
#set(ROT_MOV_IMAGE ${DATA_DIR}/mov2.nii.gz)
set(ROT_REF_IMAGE ${DATA_DIR}/r16roth.nii.gz)
set(ROT_MOV_IMAGE ${DATA_DIR}/r64roth.nii.gz)

set(ROT_OUTPUT_PREFIX ${CMAKE_BINARY_DIR}/RotTEST)
set(ROT_WARP ${ROT_OUTPUT_PREFIX}Warp.nii.gz ${ROT_OUTPUT_PREFIX}Affine.txt )
set(ROT_INVERSEWARP -i ${ROT_OUTPUT_PREFIX}Affine.txt ${ROT_OUTPUT_PREFIX}InverseWarp.nii.gz )
set(ROT_WARP_IMAGE ${CMAKE_BINARY_DIR}/rotwarped.nii.gz)
set(ROT_INVERSEWARP_IMAGE ${CMAKE_BINARY_DIR}/rotinversewarped.nii.gz)
set(ROT_WARP_FILES ${ROT_OUTPUT_PREFIX}Warp*vec.nii.gz ${ROT_OUTPUT_PREFIX}InverseWarp*vec.nii.gz ${ROT_OUTPUT_PREFIX}Affine.txt )

#add_test(ANTS_ROT_SYN ANTS 3 -m MSQ[${ROT_REF_IMAGE},${ROT_MOV_IMAGE},1,0] -t SyN[0.5,2,0.05] -i 50x5x0 -r Gauss[3,0.] -o ${ROT_OUTPUT_PREFIX}.nii.gz)
#add_test(ANTS_ROT_SYN_WARP WarpImageMultiTransform 3 ${ROT_MOV_IMAGE} ${ROT_WARP_IMAGE} -R ${ROT_REF_IMAGE} ${ROT_WARP} )
#add_test(ANTS_ROT_SYN_WARP_METRIC_0 MeasureImageSimilarity 3 0 ${ROT_REF_IMAGE} ${ROT_WARP_IMAGE} ${ROT_OUTPUT_PREFIX}log.txt ${ROT_OUTPUT_PREFIX}metric.nii.gz 35.4473 0.5)
#add_test(ANTS_ROT_SYN_INVERSEWARP WarpImageMultiTransform 3 ${ROT_REF_IMAGE} ${ROT_INVERSEWARP_IMAGE} -R ${ROT_MOV_IMAGE} ${ROT_INVERSEWARP})
#add_test(ANTS_ROT_SYN_INVERSEWARP_METRIC_0 MeasureImageSimilarity 3 0 ${ROT_MOV_IMAGE} ${ROT_INVERSEWARP_IMAGE} ${ROT_OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 34.8184 0.5)
#add_test(ANTS_ROT_SYN_CLEAN rm ${ROT_WARP_FILES})

add_test(ANTS_ROT_GSYN ANTS 3 -m MSQ[${ROT_REF_IMAGE},${ROT_MOV_IMAGE},1,0] -t SyN[0.25] -i 50x5x0 -r Gauss[3,0.0,32] -o ${ROT_OUTPUT_PREFIX}.nii.gz)
add_test(ANTS_ROT_GSYN_WARP WarpImageMultiTransform 3 ${ROT_MOV_IMAGE} ${ROT_WARP_IMAGE} -R ${ROT_REF_IMAGE} ${ROT_WARP} )
add_test(ANTS_ROT_GSYN_WARP_METRIC_0 MeasureImageSimilarity 3 0 ${ROT_REF_IMAGE} ${ROT_WARP_IMAGE} ${ROT_OUTPUT_PREFIX}log.txt ${ROT_OUTPUT_PREFIX}metric.nii.gz 35.4473 0.5)
add_test(ANTS_ROT_GSYN_INVERSEWARP WarpImageMultiTransform 3 ${ROT_REF_IMAGE} ${ROT_INVERSEWARP_IMAGE} -R ${ROT_MOV_IMAGE} ${ROT_INVERSEWARP})
add_test(ANTS_ROT_GSYN_INVERSEWARP_METRIC_0 MeasureImageSimilarity 3 0 ${ROT_MOV_IMAGE} ${ROT_INVERSEWARP_IMAGE} ${ROT_OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz  34.8184 0.5)
#add_test(ANTS_ROT_GSYN_CLEAN rm ${ROT_WARP_FILES})

add_test(ANTS_ROT_EXP ANTS 3 -m MSQ[${ROT_REF_IMAGE},${ROT_MOV_IMAGE},1,0] -t Exp[3,2] -i 50x5x1 -r Gauss[3,0.0,32] -o ${ROT_OUTPUT_PREFIX}.nii.gz)
add_test(ANTS_ROT_EXP_WARP WarpImageMultiTransform 3 ${ROT_MOV_IMAGE} ${ROT_WARP_IMAGE} -R ${ROT_REF_IMAGE} ${ROT_WARP} )
add_test(ANTS_ROT_EXP_WARP_METRIC_0 MeasureImageSimilarity 3 0 ${ROT_REF_IMAGE} ${ROT_WARP_IMAGE} ${ROT_OUTPUT_PREFIX}log.txt ${ROT_OUTPUT_PREFIX}metric.nii.gz 35.4473 0.5)
add_test(ANTS_ROT_EXP_WARP2 WarpImageMultiTransform 3 ${ROT_MOV_IMAGE} ${ROT_WARP_IMAGE} -R ${ROT_REF_IMAGE} --ANTS-prefix ${ROT_OUTPUT_PREFIX} )
add_test(ANTS_ROT_EXP_WARP2_METRIC_0_2 MeasureImageSimilarity 3 0 ${ROT_REF_IMAGE} ${ROT_WARP_IMAGE} ${ROT_OUTPUT_PREFIX}log.txt ${ROT_OUTPUT_PREFIX}metric.nii.gz 35.4473 0.5)
add_test(ANTS_ROT_EXP_INVERSEWARP WarpImageMultiTransform 3 ${ROT_REF_IMAGE} ${ROT_INVERSEWARP_IMAGE} -R ${ROT_MOV_IMAGE} ${ROT_INVERSEWARP})
add_test(ANTS_ROT_EXP_INVERSEWARP_METRIC_0 MeasureImageSimilarity 3 0 ${ROT_MOV_IMAGE} ${ROT_INVERSEWARP_IMAGE} ${ROT_OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz  34.8184 0.5)
add_test(ANTS_ROT_EXP_INVERSEWARP2 WarpImageMultiTransform 3 ${ROT_REF_IMAGE} ${ROT_INVERSEWARP_IMAGE} -R ${ROT_MOV_IMAGE} --ANTS-prefix-invert ${ROT_OUTPUT_PREFIX})
add_test(ANTS_ROT_EXP_INVERSEWARP2_METRIC_0_2 MeasureImageSimilarity 3 0 ${ROT_MOV_IMAGE} ${ROT_INVERSEWARP_IMAGE} ${ROT_OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz  34.8184 0.5)
#add_test(ANTS_ROT_EXP_CLEAN rm ${ROT_WARP_FILES})

###
#  Test SyN with time
###
set(CHALF_IMAGE ${DATA_DIR}/chalf.nii.gz)
set(C_IMAGE ${DATA_DIR}/c.nii.gz)

add_test(ANTS_SYN_WITH_TIME ANTS 2 -m MSQ[${CHALF_IMAGE},${C_IMAGE},1,0.] -t SyN[1,10,0.05] -i 150x100x2x2 -r Gauss[0.5,0.1] -o ${OUTPUT_PREFIX} --geodesic 1 --number-of-affine-iterations 0)
add_test(ANTS_SYN_WITH_TIME_WARP WarpImageMultiTransform 2 ${C_IMAGE} ${OUTPUT_PREFIX}.nii.gz ${OUTPUT_PREFIX}Warp.nii.gz  -R  ${CHALF_IMAGE} )
add_test(ANTS_SYN_WITH_TIME_METRIC MeasureImageSimilarity 2 0 ${CHALF_IMAGE} ${OUTPUT_PREFIX}.nii.gz ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz  0.0943736 0.1)


##
# Apocrita tests
##
set(R16_MASK ${DATA_DIR}/r16mask.nii.gz)
set(R16_PRIORS ${DATA_DIR}/r16priors.nii.gz)

#add_test(APOC_OTSU_INIT Apocrita 2 -i otsu[${R16_IMAGE},3] -x ${R16_MASK} -n 10 -m [0.3,1x1,0.2,0.1] -o ${OUTPUT_PREFIX}.nii.gz )
#add_test(APOC_OTSU_INIT_RADIUS_2x2 Apocrita 2 -i otsu[${R16_IMAGE},3] -x ${R16_MASK} -n 10 -m [0.3,2,0.2,0.1] -o ${OUTPUT_PREFIX}.nii.gz )
#add_test(APOC_KMEANS_INIT Apocrita 2 -i kmeans[${R16_IMAGE},3] -x ${R16_MASK} -n 10 -m [0.3,1x1,0.2,0.1] -o [${OUTPUT_PREFIX}.nii.gz,${OUTPUT_PREFIX}_posteriors%d.nii.gz])
#add_test(APOC_PRIORLABELIMAGE_INIT Apocrita 2 -i priorlabelimage[${R16_IMAGE},5,${R16_PRIORS},0.5] -x ${R16_MASK} -n 10 -m [0.3,1x1,0.2,0.1] -o [${OUTPUT_PREFIX}.nii.gz,${OUTPUT_PREFIX}_posteriors%d.nii.gz] -l 1[1,0.75] -l 2[1,1.0] -l 3[0.5,0.5] -l 4[1,1])

endif(RUN_LONG_TESTS)

