cmake_minimum_required(VERSION 2.4)

if(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#Change PROJECT_NAME to the name of your project
PROJECT(ANTS)


OPTION(USE_ITK "Use ITK Libraries" ON)
OPTION(USE_VTK "Use VTK Libraries" OFF)

set (CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

FIND_PATH(VTK_DIR UseVTK.cmake $ENV{VTK_DIR} )

# Set up ITK
IF(USE_ITK)
  FIND_PACKAGE(ITK)
  IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
      "Cannot build without ITK.  Please set ITK_DIR.")
  ENDIF(ITK_FOUND)
ENDIF(USE_ITK)

# Set up VTK
IF(USE_VTK)
  FIND_PACKAGE(VTK)
  IF(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
  ELSE(VTK_FOUND)
     MESSAGE("Cannot build some programs without VTK.  Please set VTK_DIR if you need these programs.")
  ENDIF(VTK_FOUND)
ENDIF(USE_VTK)

#The following lines are required to use Dart
INCLUDE(CTest)
ENABLE_TESTING()


INCLUDE (${CMAKE_ROOT}/Modules/FindITK.cmake)
IF (USE_ITK_FILE)
  INCLUDE(${USE_ITK_FILE})
ENDIF(USE_ITK_FILE)

SET(DART_TESTING_TIMEOUT 1500)

SET(PICSL_INCLUDE_DIRS
../Utilities
../ImageRegistration
../ImageSegmentation
../GraphTheory
../Tensor
../Temporary
)

INCLUDE_DIRECTORIES(${PICSL_INCLUDE_DIRS})

LINK_DIRECTORIES( ${ITK_LIBRARY_PATH} )

# non-templated class -- this should be stored in a library and linked in...
SET(UI_SOURCES "../Utilities/antsCommandLineParser" "../Utilities/antsCommandLineOption" )

ADD_EXECUTABLE(ANTS ANTS.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(ANTS ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(ANTSJacobian ANTSJacobian.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(ANTSJacobian ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(PrintHeader PrintHeader.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(PrintHeader ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(ResetDirection ResetDirection.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(ResetDirection ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

#ADD_EXECUTABLE(ResetDirection2 ResetDirection2.cxx ${UI_SOURCES})
#TARGET_LINK_LIBRARIES(ResetDirection2 ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

#ADD_EXECUTABLE(BoykovGraphCutFilter BoykovGraphCutFilter.cxx ${UI_SOURCES})
#TARGET_LINK_LIBRARIES(BoykovGraphCutFilter ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

#ADD_EXECUTABLE(GeneralLinearModel GeneralLinearModel.cxx ${UI_SOURCES})
#TARGET_LINK_LIBRARIES(GeneralLinearModel ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )


ADD_EXECUTABLE(ANTSUseLandmarkImagesToGetAffineTransform ANTSUseLandmarkImagesToGetAffineTransform.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(ANTSUseLandmarkImagesToGetAffineTransform ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )


ADD_EXECUTABLE(ANTSUseDeformationFieldToGetAffineTransform ANTSUseDeformationFieldToGetAffineTransform ${UI_SOURCES})
TARGET_LINK_LIBRARIES(ANTSUseDeformationFieldToGetAffineTransform ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )


#ADD_EXECUTABLE(FDR FDR.cxx ${UI_SOURCES})
#TARGET_LINK_LIBRARIES(FDR ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(ImageMath ImageMath.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(ImageMath ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )


ADD_EXECUTABLE(SurfaceCurvature SurfaceCurvature.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(SurfaceCurvature ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(SurfaceBasedSmoothing SurfaceBasedSmoothing.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(SurfaceBasedSmoothing ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )



# ADD_EXECUTABLE(test test.cxx)
#TARGET_LINK_LIBRARIES(test)

#ADD_EXECUTABLE(WarpImage WarpImage.cxx ${UI_SOURCES})
#TARGET_LINK_LIBRARIES(WarpImage ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(LaplacianThickness LaplacianThickness.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(LaplacianThickness ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(SetOrigin SetOrigin.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(SetOrigin ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(SetSpacing SetSpacing.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(SetSpacing ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(SetDirectionByMatrix SetDirectionByMatrix.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(SetDirectionByMatrix ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )


#ADD_EXECUTABLE(WarpFunctionalImage WarpFunctionalImage.cxx ${UI_SOURCES})
#TARGET_LINK_LIBRARIES(WarpFunctionalImage ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(ConvertScalarImageToRGB ConvertScalarImageToRGB.cxx)
TARGET_LINK_LIBRARIES(ConvertScalarImageToRGB ITKIO)

ADD_EXECUTABLE(CreateWarpedGridImage CreateWarpedGridImage.cxx)
TARGET_LINK_LIBRARIES(CreateWarpedGridImage ITKCommon  ITKBasicFilters ITKIO)

# ADD_EXECUTABLE(CreateJacobianDeterminantImage CreateJacobianDeterminantImage.cxx)
# TARGET_LINK_LIBRARIES(CreateJacobianDeterminantImage ITKCommon  ITKBasicFilters ITKIO)

ADD_EXECUTABLE(MeasureImageSimilarity MeasureImageSimilarity.cxx)
TARGET_LINK_LIBRARIES(MeasureImageSimilarity ITKCommon  ITKBasicFilters ITKIO)

ADD_EXECUTABLE(ConvertToJpg ConvertToJpg.cxx)
TARGET_LINK_LIBRARIES(ConvertToJpg ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(ConvertImagePixelType ConvertImagePixelType.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(ConvertImagePixelType ITKCommon  ITKBasicFilters ITKIO    )


ADD_EXECUTABLE(AverageImages AverageImages.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(AverageImages ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(AverageTensorImages AverageTensorImages.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(AverageTensorImages ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(ImageSetStatistics ImageSetStatistics.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(ImageSetStatistics ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(ThresholdImage ThresholdImage.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(ThresholdImage ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(MultiplyImages MultiplyImages.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(MultiplyImages ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(SmoothImage SmoothImage.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(SmoothImage ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(ClusterImageStatistics ClusterImageStatistics.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(ClusterImageStatistics ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(LabelClustersUniquely LabelClustersUniquely.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(LabelClustersUniquely ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

#ADD_EXECUTABLE(StudentsTestOnImages StudentsTestOnImages.cxx ${UI_SOURCES})
#TARGET_LINK_LIBRARIES(StudentsTestOnImages ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(LabelOverlapMeasures LabelOverlapMeasures.cxx )
TARGET_LINK_LIBRARIES(LabelOverlapMeasures ITKIO)

ADD_EXECUTABLE(MeasureMinMaxMean MeasureMinMaxMean.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(MeasureMinMaxMean ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(WarpImageMultiTransform WarpImageMultiTransform ${UI_SOURCES})
TARGET_LINK_LIBRARIES(WarpImageMultiTransform ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(ComposeMultiTransform ComposeMultiTransform ${UI_SOURCES})
TARGET_LINK_LIBRARIES(ComposeMultiTransform ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )


ADD_EXECUTABLE(StackSlices StackSlices.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(StackSlices ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(MemoryTest MemoryTest.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(MemoryTest ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

#ADD_EXECUTABLE(ANTSOrientImage ANTSOrientImage.cxx ${UI_SOURCES})
#TARGET_LINK_LIBRARIES(ANTSOrientImage ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(PermuteFlipImageOrientationAxes PermuteFlipImageOrientationAxes.cxx ${UI_SOURCES})
TARGET_LINK_LIBRARIES(PermuteFlipImageOrientationAxes ITKCommon  ITKBasicFilters ITKIO   ITKNumerics  ITKStatistics  ITKFEM )

ADD_EXECUTABLE(ImageCompare ImageCompare.cxx)
TARGET_LINK_LIBRARIES(ImageCompare ITKIO)

ADD_EXECUTABLE(ResampleImageBySpacing ResampleImageBySpacing.cxx)
TARGET_LINK_LIBRARIES(ResampleImageBySpacing ITKIO)

#ADD_EXECUTABLE(ResampleImageByOtherImage ResampleImageByOtherImage.cxx)
#TARGET_LINK_LIBRARIES(ResampleImageByOtherImage ITKIO)

ADD_EXECUTABLE(CopyImageHeaderInformation CopyImageHeaderInformation.cxx)
TARGET_LINK_LIBRARIES(CopyImageHeaderInformation ITKIO)

ADD_EXECUTABLE(WarpTimeSeriesImageMultiTransform WarpTimeSeriesImageMultiTransform.cxx)
TARGET_LINK_LIBRARIES(WarpTimeSeriesImageMultiTransform ITKCommon ITKBasicFilters ITKIO ITKNumerics  ITKStatistics)

ADD_EXECUTABLE(ResampleImage ResampleImage.cxx )
TARGET_LINK_LIBRARIES(ResampleImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(ExtractSliceFromImage ExtractSliceFromImage.cxx )
TARGET_LINK_LIBRARIES(ExtractSliceFromImage ITKNumerics ITKIO ITKCommon)

ADD_EXECUTABLE(TileImages TileImages.cxx )
TARGET_LINK_LIBRARIES(TileImages ITKNumerics ITKIO ITKCommon)


ADD_EXECUTABLE(WarpTensorImageMultiTransform WarpTensorImageMultiTransform.cxx)
TARGET_LINK_LIBRARIES(WarpTensorImageMultiTransform ITKCommon ITKBasicFilters ITKIO ITKNumerics  ITKStatistics)

ADD_EXECUTABLE(ReorientTensorImage ReorientTensorImage.cxx)
TARGET_LINK_LIBRARIES(ReorientTensorImage ITKCommon ITKBasicFilters ITKIO)

ADD_EXECUTABLE(N3BiasFieldCorrection N3BiasFieldCorrection.cxx)
TARGET_LINK_LIBRARIES(N3BiasFieldCorrection ITKIO ITKStatistics ITKNumerics)

ADD_EXECUTABLE(N4BiasFieldCorrection N4BiasFieldCorrection.cxx  ${UI_SOURCES})
TARGET_LINK_LIBRARIES(N4BiasFieldCorrection ITKIO ITKStatistics)

ADD_EXECUTABLE(KellySlater KellySlater.cxx)
TARGET_LINK_LIBRARIES(KellySlater ITKIO ITKStatistics)

IF (USE_VTK)
INCLUDE (${CMAKE_ROOT}/Modules/FindVTK.cmake)
IF (USE_VTK_FILE)
  INCLUDE(${USE_VTK_FILE})
  ADD_EXECUTABLE(ConvertVectorFieldToVTK ConvertVectorFieldToVTK.cxx )
  TARGET_LINK_LIBRARIES(ConvertVectorFieldToVTK ITKNumerics ITKIO ITKCommon vtkIO)

  ADD_EXECUTABLE(GetMeshAndTopology GetMeshAndTopology.cxx)
  TARGET_LINK_LIBRARIES(GetMeshAndTopology ITKCommon ITKNumerics ITKFEM  ITKIO ITKMetaIO  vtkRendering vtksys vtkIO)

  ADD_EXECUTABLE(CheckTopology CheckTopology.cxx)
  TARGET_LINK_LIBRARIES(CheckTopology ITKCommon ITKNumerics ITKFEM  ITKIO ITKMetaIO  vtkRendering vtksys vtkIO)

ADD_EXECUTABLE(WarpVTKPolyDataMultiTransform WarpVTKPolyDataMultiTransform)
TARGET_LINK_LIBRARIES(WarpVTKPolyDataMultiTransform ITKIO vtksys vtkIO )


  SET(FLAT_SRC
    ../Temporary/itkFEMElement3DC0LinearTriangular.cxx
    ../Temporary/itkFEMElement3DC0LinearTriangularLaplaceBeltrami.cxx
    ../Temporary/itkFEMElement3DC0LinearTriangularMembrane.cxx
  )
  # Build the library
  ADD_LIBRARY(FLATFEM ${FLAT_SRC})
  TARGET_LINK_LIBRARIES(FLATFEM ITKBasicFilters ITKIO)


  ADD_EXECUTABLE(ConformalMapping ConformalMapping.cxx)
  TARGET_LINK_LIBRARIES(ConformalMapping ITKCommon ITKNumerics ITKFEM FLATFEM ITKIO ITKMetaIO  vtkRendering vtkIO)

 ADD_EXECUTABLE(ANTSConformalMapping ANTSConformalMapping.cxx ${UI_SOURCES})
  TARGET_LINK_LIBRARIES(ANTSConformalMapping ITKCommon ITKNumerics ITKFEM FLATFEM ITKIO ITKMetaIO  vtkRendering vtkIO)

ENDIF (USE_VTK_FILE)
ENDIF(USE_VTK)



ADD_EXECUTABLE(ANTSIntegrateVectorField ANTSIntegrateVectorField.cxx)
TARGET_LINK_LIBRARIES(ANTSIntegrateVectorField ITKIO ITKStatistics)

#IF( ITK_USE_REVIEW_STATISTICS )
  ADD_EXECUTABLE(Atropos Atropos.cxx ${UI_SOURCES})
  TARGET_LINK_LIBRARIES(Atropos ITKCommon  ITKBasicFilters ITKIO ITKStatistics   )
#ENDIF( ITK_USE_REVIEW_STATISTICS )


#Set any libraries that your project depends on.
#examples: ITKCommon, VTKRendering, etc
SET(Libraries
ITKIO
)

#SET(CurrentExe "GlamorousGlue")
#ADD_EXECUTABLE(${CurrentExe} GlamorousGlue.cxx)
#TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

##
# Copy all the shell scripts in the Scripts/ directory to the bin directory
##
OPTION( COPY_SCRIPT_FILES_TO_BIN_DIR "Copy the script files to the ANTS bin directory." ON )

IF( COPY_SCRIPT_FILES_TO_BIN_DIR )
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
ENDIF( COPY_SCRIPT_FILES_TO_BIN_DIR )

###
#  Perform testing
###

SET(DATA_DIR ${CMAKE_SOURCE_DIR}/Data)
SET(R16_IMAGE ${DATA_DIR}/r16slice.nii)
SET(R64_IMAGE ${DATA_DIR}/r64slice.nii)
SET(OUTPUT_PREFIX ${CMAKE_BINARY_DIR}/TEST)
SET(WARP ${OUTPUT_PREFIX}Warp.nii.gz ${OUTPUT_PREFIX}Affine.txt )
SET(INVERSEWARP -i ${OUTPUT_PREFIX}Affine.txt ${OUTPUT_PREFIX}InverseWarp.nii.gz )
SET(WARP_IMAGE ${CMAKE_BINARY_DIR}/warped.nii.gz)
SET(INVERSEWARP_IMAGE ${CMAKE_BINARY_DIR}/inversewarped.nii.gz)

SET(DEVIL_IMAGE ${DATA_DIR}/Frown.nii)
SET(ANGEL_IMAGE ${DATA_DIR}/Smile.nii)

SET(SEG_IMAGE ${DATA_DIR}/nslice.nii.gz)



###
#  ANTS metric testing
###
ADD_TEST(ANTS_CC_1   ANTS 2 -m  PR[ ${R16_IMAGE}, ${R64_IMAGE},  1 ,2 ] -r Gauss[ 3 , 0 ] -t SyN[ 0.5   ] -i 50x50x30 -o ${OUTPUT_PREFIX}.nii.gz)
ADD_TEST(ANTS_CC_1_WARP WarpImageMultiTransform 2 ${R64_IMAGE} ${WARP_IMAGE} ${WARP} -R ${R16_IMAGE} )
ADD_TEST(ANTS_CC_1_JPG  ConvertToJpg ${WARP_IMAGE} ANTSCC1.jpg)
ADD_TEST(ANTS_CC_1_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 11.9992 0.05)
ADD_TEST(ANTS_CC_1_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.61 0.05)
ADD_TEST(ANTS_CC_1_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.00038593 0.05)
ADD_TEST(ANTS_CC_1_INVERSEWARP WarpImageMultiTransform 2 ${R16_IMAGE} ${INVERSEWARP_IMAGE} ${INVERSEWARP}  -R ${R16_IMAGE} )
ADD_TEST(ANTS_CC_1_INVERSEWARP_METRIC_0 MeasureImageSimilarity 2 0 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12.1606 0.05)
ADD_TEST(ANTS_CC_1_INVERSEWARP_METRIC_1 MeasureImageSimilarity 2 1 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.61 0.05)
ADD_TEST(ANTS_CC_1_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000380545 0.05)

ADD_TEST(ANTS_CC_2   ANTS 2 -m  PR[${R16_IMAGE},${R64_IMAGE},1,4] -r Gauss[3,0] -t SyN[0.5] -i 50x50x30 -o ${OUTPUT_PREFIX}.nii.gz --go-faster true )
ADD_TEST(ANTS_CC_2_WARP WarpImageMultiTransform 2 ${R64_IMAGE} ${WARP_IMAGE} ${WARP}  -R ${R16_IMAGE} )
ADD_TEST(ANTS_CC_2_JPG  ConvertToJpg ${WARP_IMAGE} ANTSCC2.jpg)
ADD_TEST(ANTS_CC_2_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 11.7083 0.05)
ADD_TEST(ANTS_CC_2_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.62 0.05)
ADD_TEST(ANTS_CC_2_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000461792 0.05)
ADD_TEST(ANTS_CC_2_INVERSEWARP WarpImageMultiTransform 2 ${R16_IMAGE} ${INVERSEWARP_IMAGE} ${INVERSEWARP}   -R ${R16_IMAGE}  )
ADD_TEST(ANTS_CC_2_INVERSEWARP_METRIC_0 MeasureImageSimilarity 2 0 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 11.8012 0.05)
ADD_TEST(ANTS_CC_2_INVERSEWARP_METRIC_1 MeasureImageSimilarity 2 1 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.62 0.05)
ADD_TEST(ANTS_CC_2_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000442839 0.05)
#ADD_TEST(ANTS_CC_2_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -5 0.05)

ADD_TEST(ANTS_CC_3   ANTS 2 -m  CC[${R16_IMAGE},${R64_IMAGE},1,4] -r Gauss[3,0] -t SyN[0.5] -i 50x50x30 -o ${OUTPUT_PREFIX}.nii.gz  --go-faster true )
ADD_TEST(ANTS_CC_3_WARP WarpImageMultiTransform 2 ${R64_IMAGE} ${WARP_IMAGE} ${WARP}  -R ${R16_IMAGE} )
ADD_TEST(ANTS_CC_3_JPG  ConvertToJpg ${WARP_IMAGE} ANTSCC2.jpg)
ADD_TEST(ANTS_CC_3_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 11.7083 0.05)
ADD_TEST(ANTS_CC_3_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.62 0.05)
ADD_TEST(ANTS_CC_3_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000461792 0.05)
ADD_TEST(ANTS_CC_3_INVERSEWARP WarpImageMultiTransform 2 ${R16_IMAGE} ${INVERSEWARP_IMAGE} ${INVERSEWARP}   -R ${R16_IMAGE} )
ADD_TEST(ANTS_CC_3_INVERSEWARP_METRIC_0 MeasureImageSimilarity 2 0 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 11.8012 0.05)
ADD_TEST(ANTS_CC_3_INVERSEWARP_METRIC_1 MeasureImageSimilarity 2 1 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.62 0.05)
ADD_TEST(ANTS_CC_3_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000442839 0.05)
#ADD_TEST(ANTS_CC_3_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -5 0.05)

ADD_TEST(ANTS_MSQ    ANTS 2 -m MSQ[${R16_IMAGE},${R64_IMAGE},1,0] -r Gauss[2,0] -t SyN[0.5] -i 50x50x30 -o ${OUTPUT_PREFIX}.nii.gz)
ADD_TEST(ANTS_MSQ_WARP WarpImageMultiTransform 2 ${R64_IMAGE} ${WARP_IMAGE} ${WARP}  -R ${R16_IMAGE} )
ADD_TEST(ANTS_MSQ_JPG  ConvertToJpg ${WARP_IMAGE} ANTSMSQ.jpg)
ADD_TEST(ANTS_MSQ_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 11.7416 0.05)
ADD_TEST(ANTS_MSQ_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.6 0.05)
ADD_TEST(ANTS_MSQ_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000543902 0.05)
ADD_TEST(ANTS_MSQ_INVERSEWARP WarpImageMultiTransform 2 ${R16_IMAGE} ${INVERSEWARP_IMAGE} ${INVERSEWARP}   -R ${R16_IMAGE} )
ADD_TEST(ANTS_MSQ_INVERSEWARP_METRIC_0 MeasureImageSimilarity 2 0 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 11.7845 0.05)
ADD_TEST(ANTS_MSQ_INVERSEWARP_METRIC_1 MeasureImageSimilarity 2 1 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.6 0.05)
ADD_TEST(ANTS_MSQ_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000501458 0.05)

ADD_TEST(ANTS_MI_1   ANTS 2 -m  MI[${R16_IMAGE},${R64_IMAGE},1,32] -r Gauss[3,0] -t SyN[0.25] -i 50x50x30 -o ${OUTPUT_PREFIX}.nii.gz)
ADD_TEST(ANTS_MI_1_WARP WarpImageMultiTransform 2 ${R64_IMAGE} ${WARP_IMAGE} ${WARP}  -R ${R16_IMAGE} )
ADD_TEST(ANTS_MI_1_JPG  ConvertToJpg ${WARP_IMAGE} ANTSMI1.jpg)
ADD_TEST(ANTS_MI_1_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12.4 0.05)
ADD_TEST(ANTS_MI_1_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.45 0.05)
ADD_TEST(ANTS_MI_1_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000370686 0.05)
ADD_TEST(ANTS_MI_1_INVERSEWARP WarpImageMultiTransform 2 ${R16_IMAGE} ${INVERSEWARP_IMAGE} ${INVERSEWARP}  -R ${R16_IMAGE} )
ADD_TEST(ANTS_MI_1_INVERSEWARP_METRIC_0 MeasureImageSimilarity 2 0 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12.6 0.05)
ADD_TEST(ANTS_MI_1_INVERSEWARP_METRIC_1 MeasureImageSimilarity 2 1 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.45 0.05)
ADD_TEST(ANTS_MI_1_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000366879 0.05)

ADD_TEST(ANTS_MI_2   ANTS 2 -m  SMI[${R16_IMAGE},${R64_IMAGE},1,32] -r Gauss[3,0] -t SyN[0.25] -i 50x0x0 -o ${OUTPUT_PREFIX}.nii.gz)
ADD_TEST(ANTS_MI_2_WARP WarpImageMultiTransform 2 ${R64_IMAGE} ${WARP_IMAGE} ${WARP}  -R ${R16_IMAGE} )
ADD_TEST(ANTS_MI_2_JPG  ConvertToJpg ${WARP_IMAGE} ANTSMI2.jpg)
ADD_TEST(ANTS_MI_2_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12.4 0.05)
ADD_TEST(ANTS_MI_2_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.45 0.05)
ADD_TEST(ANTS_MI_2_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000370686 0.05)
ADD_TEST(ANTS_MI_2_INVERSEWARP WarpImageMultiTransform 2 ${R16_IMAGE} ${INVERSEWARP_IMAGE} ${INVERSEWARP}  -R ${R16_IMAGE}  )
ADD_TEST(ANTS_MI_2_INVERSEWARP_METRIC_0 MeasureImageSimilarity 2 0 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12.6 0.05)
ADD_TEST(ANTS_MI_2_INVERSEWARP_METRIC_1 MeasureImageSimilarity 2 1 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.45 0.05)
ADD_TEST(ANTS_MI_2_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000366879 0.05)



###
#  ANTS transform testing
###
ADD_TEST(ANTS_ELASTIC ANTS 2 -m PR[${R16_IMAGE},${R64_IMAGE},1,2] -t Elast[1]       -i 50x50x50 -r Gauss[0,1] -o ${OUTPUT_PREFIX}.nii.gz)
ADD_TEST(ANTS_ELASTIC_WARP WarpImageMultiTransform 2 ${R64_IMAGE} ${WARP_IMAGE} ${WARP}   -R ${R16_IMAGE}  )
ADD_TEST(ANTS_ELASTIC_JPG  ConvertToJpg ${WARP_IMAGE} ANTSELASTIC.jpg)
ADD_TEST(ANTS_ELASTIC_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12.5 0.05)
ADD_TEST(ANTS_ELASTIC_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.4 0.05)
ADD_TEST(ANTS_ELASTIC_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000459869 0.05)

ADD_TEST(ANTS_GSYN    ANTS 2 -m PR[${R16_IMAGE},${R64_IMAGE},1,2] -t SyN[0.25]       -i 50x50x50 -r Gauss[3,0.0,32] -o ${OUTPUT_PREFIX}.nii.gz)
ADD_TEST(ANTS_GSYN_WARP WarpImageMultiTransform 2 ${R64_IMAGE} ${WARP_IMAGE} ${WARP}  -R ${R16_IMAGE} )
ADD_TEST(ANTS_GSYN_JPG  ConvertToJpg ${WARP_IMAGE} ANTSGSYN.jpg)
ADD_TEST(ANTS_GSYN_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 11.7734 0.05)
ADD_TEST(ANTS_GSYN_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.6 0.05)
ADD_TEST(ANTS_GSYN_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000478672 0.05)
ADD_TEST(ANTS_GSYN_INVERSEWARP WarpImageMultiTransform 2 ${R16_IMAGE} ${INVERSEWARP_IMAGE} ${INVERSEWARP}  -R ${R16_IMAGE}  )
ADD_TEST(ANTS_GSYN_JPGINV  ConvertToJpg ${INVERSEWARP_IMAGE} ANTSGSYNINV.jpg)
ADD_TEST(ANTS_GSYN_INVERSEWARP_METRIC_0 MeasureImageSimilarity 2 0 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12.0541 0.05)
ADD_TEST(ANTS_GSYN_INVERSEWARP_METRIC_1 MeasureImageSimilarity 2 1 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.6 0.05)
ADD_TEST(ANTS_GSYN_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000475175 0.05)

ADD_TEST(ANTS_EXP     ANTS 2 -m PR[${R16_IMAGE},${R64_IMAGE},1,4] -t Exp[0.5,2,0.5]       -i 50x50x50 -r Gauss[0.5,0.25] -o ${OUTPUT_PREFIX}.nii.gz)
ADD_TEST(ANTS_EXP_WARP WarpImageMultiTransform 2 ${R64_IMAGE} ${WARP_IMAGE} ${WARP}  -R ${R16_IMAGE} )
ADD_TEST(ANTS_EXP_JPG  ConvertToJpg ${WARP_IMAGE} ANTSEXP.jpg)
ADD_TEST(ANTS_EXP_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12 0.05)
ADD_TEST(ANTS_EXP_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.5 0.05)
ADD_TEST(ANTS_EXP_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000423693 0.05)
ADD_TEST(ANTS_EXP_INVERSEWARP WarpImageMultiTransform 2 ${R16_IMAGE} ${INVERSEWARP_IMAGE} ${INVERSEWARP}   -R ${R16_IMAGE} )
ADD_TEST(ANTS_EXP_INVERSEWARP_METRIC_0 MeasureImageSimilarity 2 0 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12.349 0.05)
ADD_TEST(ANTS_EXP_INVERSEWARP_METRIC_1 MeasureImageSimilarity 2 1 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.5 0.05)
ADD_TEST(ANTS_EXP_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000339606 0.05)

#ADD_TEST(ANTS_GSYN    ANTS 2 -m PR[${R16_IMAGE},${R64_IMAGE},1,2] -t SyN[0.75]            -i 50x50x50 -r Gauss[3,0.0,32] -o ${OUTPUT_PREFIX}.nii.gz)
 ADD_TEST(ANTS_SYN     ANTS 2 -m PR[${R16_IMAGE},${R64_IMAGE},1,2] -t SyN[0.5,2,0.05] -i 50x50x50 -r Gauss[3,0.0,32] -o ${OUTPUT_PREFIX}.nii.gz)
ADD_TEST(ANTS_SYN_WARP WarpImageMultiTransform 2 ${R64_IMAGE} ${WARP_IMAGE} ${WARP}  -R ${R16_IMAGE}  )
ADD_TEST(ANTS_SYN_JPG  ConvertToJpg ${WARP_IMAGE} ANTSSYN.jpg)
ADD_TEST(ANTS_SYN_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12.0239 0.05)
ADD_TEST(ANTS_SYN_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.6 0.05)
ADD_TEST(ANTS_SYN_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${R16_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000461922 0.05)
ADD_TEST(ANTS_SYN_INVERSEWARP WarpImageMultiTransform 2 ${R16_IMAGE} ${INVERSEWARP_IMAGE} ${INVERSEWARP}  -R ${R16_IMAGE}  )
ADD_TEST(ANTS_SYN_JPGINV  ConvertToJpg ${INVERSEWARP_IMAGE} ANTSSYNINV.jpg)
ADD_TEST(ANTS_SYN_INVERSEWARP_METRIC_0 MeasureImageSimilarity 2 0 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 12.5104 0.05)
ADD_TEST(ANTS_SYN_INVERSEWARP_METRIC_1 MeasureImageSimilarity 2 1 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.6 0.05)
ADD_TEST(ANTS_SYN_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${R64_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000444279 0.05)

###
# PSE sub-tests:  Check to see if .txt files and .vtk files also run correctly
###

SET(ANGEL_IMAGE_TXT ${DATA_DIR}/Smile.txt)
SET(DEVIL_IMAGE_TXT ${DATA_DIR}/Frown.txt)
SET(ANGEL_IMAGE_VTK ${DATA_DIR}/Smile.vtk)
SET(DEVIL_IMAGE_VTK ${DATA_DIR}/Frown.vtk)

ADD_TEST(ANTS_PSE_MSQ_TXT ANTS 2 -i 91x70x55x40x30  -r Gauss[3,0.,32] -t SyN[0.25]
                                 -m MSQ[${DEVIL_IMAGE},${ANGEL_IMAGE},1,0]
                                 -m PSE[${DEVIL_IMAGE},${ANGEL_IMAGE},${DEVIL_IMAGE_TXT},${ANGEL_IMAGE_TXT},1,0.33,11,1,25]
                                 --continue-affine 0 --number-of-affine-iterations 0
                                 -o ${OUTPUT_PREFIX}.nii.gz)
ADD_TEST(ANTS_PSE_MSQ_VTK ANTS 2 -i 91x70x55x40x30  -r Gauss[3,0.,32] -t SyN[0.25]
                                 -m MSQ[${DEVIL_IMAGE},${ANGEL_IMAGE},1,0]
                                 -m PSE[${DEVIL_IMAGE},${ANGEL_IMAGE},${DEVIL_IMAGE_VTK},${ANGEL_IMAGE_VTK},1,0.33,11,1,25]
                                 --continue-affine 0 --number-of-affine-iterations 0
                                 -o ${OUTPUT_PREFIX}.nii.gz)

###
#  ANTS labeled data testing
###
OPTION(RUN_LONG_TESTS "Run the time consuming tests." OFF )
IF(RUN_LONG_TESTS)

ADD_TEST(ANTS_PSE_MSQ_IMG ANTS 2 -i 191x170x90x90x10  -r Gauss[6,0.25] -t SyN[1,2,0.1]
                             -m PSE[${DEVIL_IMAGE},${ANGEL_IMAGE},${DEVIL_IMAGE},${ANGEL_IMAGE},0.25,0.1,100,0,10]
                             -m MSQ[${DEVIL_IMAGE},${ANGEL_IMAGE},1,0.1]
                              --continue-affine 0 --number-of-affine-iterations 0 --geodesic 2
                             -o ${OUTPUT_PREFIX}.nii.gz)
ADD_TEST(ANTS_PSE_MSQ_IMG_WARP WarpImageMultiTransform 2 ${ANGEL_IMAGE} ${WARP_IMAGE} ${WARP} -R  ${DEVIL_IMAGE} )
ADD_TEST(ANTS_PSE_JPG  ConvertToJpg ${WARP_IMAGE} ANTSPSE.jpg)
ADD_TEST(ANTS_PSE_MSQ_IMG_WARP_METRIC_0 MeasureImageSimilarity 2 0 ${DEVIL_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 0.0116083 0.05)
ADD_TEST(ANTS_PSE_MSQ_IMG_WARP_METRIC_1 MeasureImageSimilarity 2 1 ${DEVIL_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.884041 0.05)
ADD_TEST(ANTS_PSE_MSQ_IMG_WARP_METRIC_2 MeasureImageSimilarity 2 2 ${DEVIL_IMAGE} ${WARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000710551 0.05)
ADD_TEST(ANTS_PSE_MSQ_IMG_INVERSEWARP WarpImageMultiTransform 2 ${DEVIL_IMAGE} ${INVERSEWARP_IMAGE} ${INVERSEWARP} -R  ${ANGEL_IMAGE}  )
ADD_TEST(ANTS_PSE_JPGINV  ConvertToJpg ${INVERSEWARP_IMAGE} ANTSPSEINV.jpg)
ADD_TEST(ANTS_PSE_MSQ_IMG_INVERSEWARP_METRIC_0 MeasureImageSimilarity 2 0 ${ANGEL_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 0.0109054 0.05)
ADD_TEST(ANTS_PSE_MSQ_IMG_INVERSEWARP_METRIC_1 MeasureImageSimilarity 2 1 ${ANGEL_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.901632 0.05)
ADD_TEST(ANTS_PSE_MSQ_IMG_INVERSEWARP_METRIC_2 MeasureImageSimilarity 2 2 ${ANGEL_IMAGE} ${INVERSEWARP_IMAGE} ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz -0.000704717 0.05)


###
#  ANTS images with non-trival rotation header test
###
#SET(ROT_REF_IMAGE ${DATA_DIR}/ref2.nii.gz)
#SET(ROT_MOV_IMAGE ${DATA_DIR}/mov2.nii.gz)
SET(ROT_REF_IMAGE ${DATA_DIR}/r16roth.nii.gz)
SET(ROT_MOV_IMAGE ${DATA_DIR}/r64roth.nii.gz)

SET(ROT_OUTPUT_PREFIX ${CMAKE_BINARY_DIR}/RotTEST)
SET(ROT_WARP ${ROT_OUTPUT_PREFIX}Warp.nii.gz ${ROT_OUTPUT_PREFIX}Affine.txt )
SET(ROT_INVERSEWARP -i ${ROT_OUTPUT_PREFIX}Affine.txt ${ROT_OUTPUT_PREFIX}InverseWarp.nii.gz )
SET(ROT_WARP_IMAGE ${CMAKE_BINARY_DIR}/rotwarped.nii.gz)
SET(ROT_INVERSEWARP_IMAGE ${CMAKE_BINARY_DIR}/rotinversewarped.nii.gz)
SET(ROT_WARP_FILES ${ROT_OUTPUT_PREFIX}Warp*vec.nii.gz ${ROT_OUTPUT_PREFIX}InverseWarp*vec.nii.gz ${ROT_OUTPUT_PREFIX}Affine.txt )

#ADD_TEST(ANTS_ROT_SYN ANTS 3 -m MSQ[${ROT_REF_IMAGE},${ROT_MOV_IMAGE},1,0] -t SyN[0.5,2,0.05] -i 50x5x0 -r Gauss[3,0.] -o ${ROT_OUTPUT_PREFIX}.nii.gz)
#ADD_TEST(ANTS_ROT_SYN_WARP WarpImageMultiTransform 3 ${ROT_MOV_IMAGE} ${ROT_WARP_IMAGE} -R ${ROT_REF_IMAGE} ${ROT_WARP} )
#ADD_TEST(ANTS_ROT_SYN_WARP_METRIC_0 MeasureImageSimilarity 3 0 ${ROT_REF_IMAGE} ${ROT_WARP_IMAGE} ${ROT_OUTPUT_PREFIX}log.txt ${ROT_OUTPUT_PREFIX}metric.nii.gz 35.4473 0.5)
#ADD_TEST(ANTS_ROT_SYN_INVERSEWARP WarpImageMultiTransform 3 ${ROT_REF_IMAGE} ${ROT_INVERSEWARP_IMAGE} -R ${ROT_MOV_IMAGE} ${ROT_INVERSEWARP})
#ADD_TEST(ANTS_ROT_SYN_INVERSEWARP_METRIC_0 MeasureImageSimilarity 3 0 ${ROT_MOV_IMAGE} ${ROT_INVERSEWARP_IMAGE} ${ROT_OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz 34.8184 0.5)
#ADD_TEST(ANTS_ROT_SYN_CLEAN rm ${ROT_WARP_FILES})

ADD_TEST(ANTS_ROT_GSYN ANTS 3 -m MSQ[${ROT_REF_IMAGE},${ROT_MOV_IMAGE},1,0] -t SyN[0.25] -i 50x5x0 -r Gauss[3,0.0,32] -o ${ROT_OUTPUT_PREFIX}.nii.gz)
ADD_TEST(ANTS_ROT_GSYN_WARP WarpImageMultiTransform 3 ${ROT_MOV_IMAGE} ${ROT_WARP_IMAGE} -R ${ROT_REF_IMAGE} ${ROT_WARP} )
ADD_TEST(ANTS_ROT_GSYN_WARP_METRIC_0 MeasureImageSimilarity 3 0 ${ROT_REF_IMAGE} ${ROT_WARP_IMAGE} ${ROT_OUTPUT_PREFIX}log.txt ${ROT_OUTPUT_PREFIX}metric.nii.gz 35.4473 0.5)
ADD_TEST(ANTS_ROT_GSYN_INVERSEWARP WarpImageMultiTransform 3 ${ROT_REF_IMAGE} ${ROT_INVERSEWARP_IMAGE} -R ${ROT_MOV_IMAGE} ${ROT_INVERSEWARP})
ADD_TEST(ANTS_ROT_GSYN_INVERSEWARP_METRIC_0 MeasureImageSimilarity 3 0 ${ROT_MOV_IMAGE} ${ROT_INVERSEWARP_IMAGE} ${ROT_OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz  34.8184 0.5)
#ADD_TEST(ANTS_ROT_GSYN_CLEAN rm ${ROT_WARP_FILES})

ADD_TEST(ANTS_ROT_EXP ANTS 3 -m MSQ[${ROT_REF_IMAGE},${ROT_MOV_IMAGE},1,0] -t Exp[3,2] -i 50x5x1 -r Gauss[3,0.0,32] -o ${ROT_OUTPUT_PREFIX}.nii.gz)
ADD_TEST(ANTS_ROT_EXP_WARP WarpImageMultiTransform 3 ${ROT_MOV_IMAGE} ${ROT_WARP_IMAGE} -R ${ROT_REF_IMAGE} ${ROT_WARP} )
ADD_TEST(ANTS_ROT_EXP_WARP_METRIC_0 MeasureImageSimilarity 3 0 ${ROT_REF_IMAGE} ${ROT_WARP_IMAGE} ${ROT_OUTPUT_PREFIX}log.txt ${ROT_OUTPUT_PREFIX}metric.nii.gz 35.4473 0.5)
ADD_TEST(ANTS_ROT_EXP_WARP2 WarpImageMultiTransform 3 ${ROT_MOV_IMAGE} ${ROT_WARP_IMAGE} -R ${ROT_REF_IMAGE} --ANTS-prefix ${ROT_OUTPUT_PREFIX} )
ADD_TEST(ANTS_ROT_EXP_WARP2_METRIC_0_2 MeasureImageSimilarity 3 0 ${ROT_REF_IMAGE} ${ROT_WARP_IMAGE} ${ROT_OUTPUT_PREFIX}log.txt ${ROT_OUTPUT_PREFIX}metric.nii.gz 35.4473 0.5)
ADD_TEST(ANTS_ROT_EXP_INVERSEWARP WarpImageMultiTransform 3 ${ROT_REF_IMAGE} ${ROT_INVERSEWARP_IMAGE} -R ${ROT_MOV_IMAGE} ${ROT_INVERSEWARP})
ADD_TEST(ANTS_ROT_EXP_INVERSEWARP_METRIC_0 MeasureImageSimilarity 3 0 ${ROT_MOV_IMAGE} ${ROT_INVERSEWARP_IMAGE} ${ROT_OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz  34.8184 0.5)
ADD_TEST(ANTS_ROT_EXP_INVERSEWARP2 WarpImageMultiTransform 3 ${ROT_REF_IMAGE} ${ROT_INVERSEWARP_IMAGE} -R ${ROT_MOV_IMAGE} --ANTS-prefix-invert ${ROT_OUTPUT_PREFIX})
ADD_TEST(ANTS_ROT_EXP_INVERSEWARP2_METRIC_0_2 MeasureImageSimilarity 3 0 ${ROT_MOV_IMAGE} ${ROT_INVERSEWARP_IMAGE} ${ROT_OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz  34.8184 0.5)
#ADD_TEST(ANTS_ROT_EXP_CLEAN rm ${ROT_WARP_FILES})

###
#  Test SyN with time
###
SET(CHALF_IMAGE ${DATA_DIR}/chalf.nii.gz)
SET(C_IMAGE ${DATA_DIR}/c.nii.gz)

ADD_TEST(ANTS_SYN_WITH_TIME ANTS 2 -m MSQ[${CHALF_IMAGE},${C_IMAGE},1,0.] -t SyN[1,10,0.05] -i 150x100x2x2 -r Gauss[0.5,0.1] -o ${OUTPUT_PREFIX} --geodesic 1 --number-of-affine-iterations 0)
ADD_TEST(ANTS_SYN_WITH_TIME_WARP WarpImageMultiTransform 2 ${C_IMAGE} ${OUTPUT_PREFIX}.nii.gz ${OUTPUT_PREFIX}Warp.nii.gz  -R  ${CHALF_IMAGE} )
ADD_TEST(ANTS_SYN_WITH_TIME_METRIC MeasureImageSimilarity 2 0 ${CHALF_IMAGE} ${OUTPUT_PREFIX}.nii.gz ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz  0.0943736 0.1)


##
# Apocrita tests
##
SET(R16_MASK ${DATA_DIR}/r16mask.nii.gz)
SET(R16_PRIORS ${DATA_DIR}/r16priors.nii.gz)

#ADD_TEST(APOC_OTSU_INIT Apocrita 2 -i otsu[${R16_IMAGE},3] -x ${R16_MASK} -n 10 -m [0.3,1x1,0.2,0.1] -o ${OUTPUT_PREFIX}.nii.gz )
#ADD_TEST(APOC_OTSU_INIT_RADIUS_2x2 Apocrita 2 -i otsu[${R16_IMAGE},3] -x ${R16_MASK} -n 10 -m [0.3,2,0.2,0.1] -o ${OUTPUT_PREFIX}.nii.gz )
#ADD_TEST(APOC_KMEANS_INIT Apocrita 2 -i kmeans[${R16_IMAGE},3] -x ${R16_MASK} -n 10 -m [0.3,1x1,0.2,0.1] -o [${OUTPUT_PREFIX}.nii.gz,${OUTPUT_PREFIX}_posteriors%d.nii.gz])
#ADD_TEST(APOC_PRIORLABELIMAGE_INIT Apocrita 2 -i priorlabelimage[${R16_IMAGE},5,${R16_PRIORS},0.5] -x ${R16_MASK} -n 10 -m [0.3,1x1,0.2,0.1] -o [${OUTPUT_PREFIX}.nii.gz,${OUTPUT_PREFIX}_posteriors%d.nii.gz] -l 1[1,0.75] -l 2[1,1.0] -l 3[0.5,0.5] -l 4[1,1])

ENDIF(RUN_LONG_TESTS)
