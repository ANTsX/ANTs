## Note that the antsUtilities can always be built static.  It will then be linked
## Directly into the other libraries.
add_library(antsUtilities STATIC antsUtilities.cxx
            antsRegistrationTemplateHeader.cxx
            antsRegistration2DDouble.cxx antsRegistration2DFloat.cxx
            antsRegistration3DDouble.cxx antsRegistration3DFloat.cxx
            antsRegistration4DDouble.cxx antsRegistration4DFloat.cxx
            ../Utilities/ReadWriteData.cxx
            ../Utilities/antsCommandLineOption.cxx
            ../Utilities/antsCommandLineParser.cxx
            ANTsVersion.cxx
            )
target_link_libraries(antsUtilities ${ITK_LIBRARIES} )

macro(STANDARD_ANTS_BUILD ANTS_FUNCTION_NAME EXTRA_LIBS)
  set( ANTS_FUNCTION_NAME ${ANTS_FUNCTION_NAME} )

  add_library(l_${ANTS_FUNCTION_NAME} ${ANTS_FUNCTION_NAME}.cxx)
  SET_TARGET_PROPERTIES(l_${ANTS_FUNCTION_NAME} PROPERTIES
    SOVERSION ${LIBRARY_SOVERSION_INFO}  VERSION ${LIBRARY_VERSION_INFO})
  message(STATUS "${ANTS_FUNCTION_NAME} ${EXTRA_LIBS}")
  target_link_libraries(l_${ANTS_FUNCTION_NAME} antsUtilities ${EXTRA_LIBS} )
  configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
  add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
  target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

  install(TARGETS l_${ANTS_FUNCTION_NAME} ${ANTS_FUNCTION_NAME}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
  )

endmacro()

#set(PROG_WITH_NO_FLAG simpleSynRegistration)
#add_executable(${PROG_WITH_NO_FLAG} ${PROG_WITH_NO_FLAG}.cxx)
#target_link_libraries(${PROG_WITH_NO_FLAG} antsUtilities ${ITK_LIBRARIES})


set(CORE_ANTS_APPS
  antsAffineInitializer
  SurfaceCurvature
  SurfaceBasedSmoothing
  CreateJacobianDeterminantImage
  ThresholdImage
  ResampleImage
  sccan
  N4BiasFieldCorrection
  N3BiasFieldCorrection
  KellyKapowski
  antsRegistration
  antsMotionCorrStats
  antsMotionCorr
  antsApplyTransforms
  LabelGeometryMeasures
  LabelClustersUniquely
  Atropos
  ImageMath
  AverageImages
  StackSlices
  AverageAffineTransform
  AverageAffineTransformNoRigid
  antsApplyTransformsToPoints
)

set(MORE_ANTS_APPS
    ANTS
    ANTSJacobian
    PrintHeader
    ResetDirection
    ANTSUseLandmarkImagesToGetAffineTransform
    ANTSUseLandmarkImagesToGetBSplineDisplacementField
    ANTSUseDeformationFieldToGetAffineTransform
    LaplacianThickness
    SetOrigin
    SetSpacing
    SetDirectionByMatrix
    ConvertScalarImageToRGB
    CreateWarpedGridImage
    MeasureImageSimilarity
    ConvertToJpg
    ConvertImagePixelType
    ConvertInputImagePixelTypeToFloat
    FitBSplineToPoints
    AverageTensorImages
    ImageSetStatistics
    MultiplyImages
    SmoothImage
    ClusterImageStatistics
    LabelOverlapMeasures
    LesionFilling
    MeasureMinMaxMean
    WarpImageMultiTransform
    ComposeMultiTransform
    MemoryTest
    PermuteFlipImageOrientationAxes
    ImageCompare
    ResampleImageBySpacing
    CopyImageHeaderInformation
    WarpTimeSeriesImageMultiTransform
    ExtractSliceFromImage
    ExtractRegionFromImage
    ExtractRegionFromImageByMask
    PasteImageIntoImage
    TileImages
    CreateTiledMosaic
    CreateImage
    WarpTensorImageMultiTransform
    ReorientTensorImage
    RebaseTensorImage
    KellySlater
    CreateDTICohort
    antsAlignOrigin
    simpleSynRegistration
    CompositeTransformUtil
    antsMotionCorrDiffusionDirection
    antsSliceRegularizedRegistration
    CreateDisplacementField
    ANTSIntegrateVectorField
    ANTSIntegrateVelocityField
    ConvertTransformFile
    compareTwoTransforms
    TimeSCCAN
    TextureCooccurrenceFeatures
    TextureRunLengthFeatures
    ImageIntensityStatistics
    GetConnectedComponentsFeatureImages
    DeNrrd
    antsTransformInfo
    antsUtilitiesTesting
  )

foreach(ANTS_APP ${CORE_ANTS_APPS})
  STANDARD_ANTS_BUILD(${ANTS_APP} "")
endforeach()

if (BUILD_ALL_ANTS_APPS)
foreach(ANTS_APP ${MORE_ANTS_APPS})
  STANDARD_ANTS_BUILD(${ANTS_APP} "")
endforeach()
endif(BUILD_ALL_ANTS_APPS)


if(USE_VTK)
find_package(VTK 6.2 REQUIRED NO_MODULE)
find_package(VTK COMPONENTS
   vtkCommonCore
   vtkCommonDataModel
   vtkIOGeometry
   vtkIOXML
   vtkIOLegacy
   vtkIOPLY
   vtkRenderingOpenGL
   vtkFiltersModeling
   vtkImagingStencil
   vtkImagingGeneral
   vtkRenderingAnnotation
   )
include(${VTK_USE_FILE})
include_directories(${VTK_INCLUDE_DIRS})

set(VTK_ANTS_APPS
#    ConvertVectorFieldToVTK
    antsSurf
    GetMeshAndTopology
    CheckTopology
    WarpVTKPolyDataMultiTransform
    # FLATFEM
    # ConformalMapping
    # ANTSConformalMapping
  )

foreach(ANTS_APP ${VTK_ANTS_APPS})
  STANDARD_ANTS_BUILD(${ANTS_APP} "${VTK_LIBRARIES}")
endforeach()
endif(USE_VTK)

install(PROGRAMS ../Scripts/ANTSpexec.sh
     ../Scripts/antsASLProcessing.sh
     ../Scripts/antsAtroposN4.sh
     ../Scripts/antsBOLDNetworkAnalysis.R
     ../Scripts/antsBrainExtraction.sh
     ../Scripts/antsCorticalThickness.sh
     ../Scripts/antsIntermodalityIntrasubject.sh
     ../Scripts/antsIntroduction.sh
     ../Scripts/antsLaplacianBoundaryCondition.R
     ../Scripts/antsLongitudinalCorticalThickness.sh
     ../Scripts/antsMalfLabeling.sh
     ../Scripts/antsMultivariateTemplateConstruction.sh
     ../Scripts/antsMultivariateTemplateConstruction2.sh
     ../Scripts/antsNetworkAnalysis.R
     ../Scripts/antsNeuroimagingBattery
     ../Scripts/antsRegistrationSyN.sh
     ../Scripts/antsRegistrationSyNQuick.sh
     ../Scripts/waitForPBSQJobs.pl
     ../Scripts/waitForSGEQJobs.pl
     ../Scripts/waitForXGridJobs.pl
                DESTINATION bin
                PERMISSIONS  OWNER_WRITE OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
                CONFIGURATIONS  Release
        )


## Build test option
if(BUILD_TESTING)
    add_subdirectory(TestSuite)
endif(BUILD_TESTING)

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
set(CPACK_PACKAGE_VERSION 2.1.0)
set(CPACK_PACKAGE_VERSION_MAJOR 2)
set(CPACK_PACKAGE_VERSION_MINOR 1)
set(CPACK_PACKAGE_VERSION_PATCH 0)
set(CPACK_PACKAGE_INSTALL_DIRECTORY "ANTS")
set(CPACK_BINARY_GENERATORS "DragNDrop TGZ TZ")
# This must always be last!
include(CPack)
