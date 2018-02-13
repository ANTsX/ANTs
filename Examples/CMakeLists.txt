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
            ImageMathHelper.cxx ImageMathHelper2D.cxx  ImageMathHelper3D.cxx  ImageMathHelper4D.cxx
            )
target_link_libraries(antsUtilities ${ITK_LIBRARIES} )
install(TARGETS antsUtilities
    RUNTIME DESTINATION ${BIN_INSTALL_DIR}
    COMPONENT RUNTIME_antsUtilities
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}
    COMPONENT RUNTIME_antsUtilities
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
    COMPONENT DEVELOPMENT_antsUtilities
  )

macro(STANDARD_ANTS_BUILD ANTS_FUNCTION_NAME EXTRA_LIBS)
  set( ANTS_FUNCTION_NAME ${ANTS_FUNCTION_NAME} )

  add_library(l_${ANTS_FUNCTION_NAME} ${ANTS_FUNCTION_NAME}.cxx)
  set_target_properties(l_${ANTS_FUNCTION_NAME} PROPERTIES
    SOVERSION ${LIBRARY_SOVERSION_INFO}  VERSION ${LIBRARY_VERSION_INFO})
  message(STATUS "${ANTS_FUNCTION_NAME} ${EXTRA_LIBS}")
  target_link_libraries(l_${ANTS_FUNCTION_NAME} antsUtilities ${EXTRA_LIBS} )
  configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
  add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
  target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

  install(TARGETS l_${ANTS_FUNCTION_NAME} ${ANTS_FUNCTION_NAME}
    RUNTIME DESTINATION ${BIN_INSTALL_DIR}
    COMPONENT RUNTIME_${ANTS_FUNCTION_NAME}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}
    COMPONENT RUNTIME_${ANTS_FUNCTION_NAME}
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
    COMPONENT DEVELOPMENT_${ANTS_FUNCTION_NAME}
  )

endmacro()

#set(PROG_WITH_NO_FLAG simpleSynRegistration)
#add_executable(${PROG_WITH_NO_FLAG} ${PROG_WITH_NO_FLAG}.cxx)
#target_link_libraries(${PROG_WITH_NO_FLAG} antsUtilities ${ITK_LIBRARIES})


set(CORE_ANTS_APPS
  antsAffineInitializer
  antsJointFusion
  DenoiseImage
  SurfaceBasedSmoothing
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
  antsApplyTransformsToPoints
)

set(MORE_ANTS_APPS
    antsAI
    antsJointTensorFusion
    ImageMath
    iMath
    ANTS
    ANTSJacobian
    CreateJacobianDeterminantImage
    PrintHeader
    ResetDirection
    ANTSUseLandmarkImagesToGetAffineTransform
    ANTSUseLandmarkImagesToGetBSplineDisplacementField
    ANTSUseDeformationFieldToGetAffineTransform
    antsLandmarkBasedTransformInitializer
    LaplacianThickness
    SetOrigin
    SetSpacing
    SetDirectionByMatrix
    SurfaceCurvature
    ConvertScalarImageToRGB
    CreateWarpedGridImage
    MeasureImageSimilarity
    ConvertToJpg
    ConvertImage
    ConvertImagePixelType
    ConvertInputImagePixelTypeToFloat
    FitBSplineToPoints
    AverageTensorImages
    ImageSetStatistics
    MultiplyImages
    SmoothDisplacementField
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
    NonLocalSuperResolution
    WarpTensorImageMultiTransform
    ReorientTensorImage
    RebaseTensorImage
    KellySlater
    CreateDTICohort
    antsAlignOrigin
    antsMotionCorrDiffusionDirection
    antsSliceRegularizedRegistration
    ANTSIntegrateVectorField
    ANTSIntegrateVelocityField
    antsTransformInfo
    antsUtilitiesTesting
    AverageAffineTransform
    AverageAffineTransformNoRigid
    AverageImages
    simpleSynRegistration
    CompositeTransformUtil
    CreateDisplacementField
    ConvertTransformFile
    compareTwoTransforms
    SuperResolution
    TimeSCCAN
    TextureCooccurrenceFeatures
    TextureRunLengthFeatures
    ImageIntensityStatistics
    GetConnectedComponentsFeatureImages
    DeNrrd
    StackSlices
  )

foreach(ANTS_APP ${CORE_ANTS_APPS})
  STANDARD_ANTS_BUILD(${ANTS_APP} "")
endforeach()

if (BUILD_ALL_ANTS_APPS)
  foreach(ANTS_APP ${MORE_ANTS_APPS})
    STANDARD_ANTS_BUILD(${ANTS_APP} "")
  endforeach()
else()
  foreach(ANTS_APP ${MORE_ANTS_APPS})
    if(ANTS_BUILD_${ANTS_APP})
      STANDARD_ANTS_BUILD(${ANTS_APP} "")
    endif()
  endforeach()
endif()


if(USE_VTK)
set(VTK_ANTS_APPS
    # ConvertVectorFieldToVTK
    antsSurf
    antsVol
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
endif()

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
     ../Scripts/antsJointLabelFusion.sh
     ../Scripts/antsMultivariateTemplateConstruction.sh
     ../Scripts/antsMultivariateTemplateConstruction2.sh
     ../Scripts/antsNetworkAnalysis.R
     ../Scripts/antsNeuroimagingBattery
     ../Scripts/antsRegistrationSyN.sh
     ../Scripts/antsRegistrationSyNQuick.sh
     ../Scripts/waitForPBSQJobs.pl
     ../Scripts/waitForSGEQJobs.pl
     ../Scripts/waitForXGridJobs.pl
     ../Scripts/waitForSlurmJobs.pl
                DESTINATION ${SCRIPTS_INSTALL_DIR}
                PERMISSIONS  OWNER_WRITE OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
                CONFIGURATIONS  Release
                COMPONENT SCRIPTS
        )

if( COPY_SCRIPT_FILES_TO_BIN_DIR )
  message( "-- Copy scripts." )
  file( GLOB SCRIPT_FILES "../Scripts/*.sh" )
  file( GLOB PSCRIPT_FILES "../Scripts/*.pl" )
  foreach( SCRIPT_FILE ${SCRIPT_FILES} )
    message( "  Copying ${SCRIPT_FILE} to ${CMAKE_CURRENT_BINARY_DIR}/../../bin/" )
    configure_file( ${SCRIPT_FILE} ${CMAKE_CURRENT_BINARY_DIR}/../../bin/ COPYONLY )
  endforeach()
  foreach( SCRIPT_FILE ${PSCRIPT_FILES} )
    message( "  Copying ${SCRIPT_FILE} to ${CMAKE_CURRENT_BINARY_DIR}/../../bin/" )
    configure_file( ${SCRIPT_FILE} ${CMAKE_CURRENT_BINARY_DIR}/../../bin/ COPYONLY )
  endforeach()
endif()

## Build test option
if(BUILD_TESTING)
    add_subdirectory(TestSuite)
endif()

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
