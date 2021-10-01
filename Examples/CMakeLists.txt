SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
## Note that the antsUtilities can always be built static.  It will then be linked
## Directly into the other libraries.
add_library(antsUtilities antsUtilities.cxx
            antsRegistrationTemplateHeader.cxx
            antsRegistration2DDouble.cxx antsRegistration2DFloat.cxx
            antsRegistration3DDouble.cxx antsRegistration3DFloat.cxx
            antsRegistration4DDouble.cxx antsRegistration4DFloat.cxx
            ../Utilities/ReadWriteData.cxx
            ../Utilities/antsCommandLineOption.cxx
            ../Utilities/antsCommandLineParser.cxx
            ANTsVersion.cxx
            ImageMathHelper2D.cxx ImageMathHelper3D.cxx ImageMathHelper4D.cxx
            )
target_link_libraries(antsUtilities ${ITK_LIBRARIES} )
install(TARGETS antsUtilities
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT RUNTIME_antsUtilities
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT RUNTIME_antsUtilities
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT DEVELOPMENT_antsUtilities
  )

macro(STATIC_ANTS_BUILD ANTS_FUNCTION_NAME EXTRA_LIBS)
  set( ANTS_FUNCTION_NAME ${ANTS_FUNCTION_NAME} )

  add_library(l_${ANTS_FUNCTION_NAME} ${ANTS_FUNCTION_NAME}.cxx)
  set_target_properties(l_${ANTS_FUNCTION_NAME} PROPERTIES
    SOVERSION ${ANTS_VERSION_MAJOR}  VERSION ${${PROJECT_NAME}_VERSION})
  message(STATUS "${ANTS_FUNCTION_NAME} ${EXTRA_LIBS}")
  target_link_libraries(l_${ANTS_FUNCTION_NAME} antsUtilities ${EXTRA_LIBS} )
  configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )

  if (ANTS_INSTALL_LIBS_ONLY)
    install(TARGETS l_${ANTS_FUNCTION_NAME}
 #    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
 #    COMPONENT RUNTIME_${ANTS_FUNCTION_NAME}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      COMPONENT RUNTIME_${ANTS_FUNCTION_NAME}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      COMPONENT DEVELOPMENT_${ANTS_FUNCTION_NAME}
    )
  else()
    add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
    target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )
    install(TARGETS l_${ANTS_FUNCTION_NAME} ${ANTS_FUNCTION_NAME}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      COMPONENT RUNTIME_${ANTS_FUNCTION_NAME}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      COMPONENT RUNTIME_${ANTS_FUNCTION_NAME}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      COMPONENT DEVELOPMENT_${ANTS_FUNCTION_NAME}
    )
  endif()

endmacro()

macro(DYNAMIC_ANTS_BUILD ANTS_FUNCTION_NAME EXTRA_LIBS)
  set( ANTS_FUNCTION_NAME ${ANTS_FUNCTION_NAME} )

  add_library(l_${ANTS_FUNCTION_NAME} ${ANTS_FUNCTION_NAME}.cxx)
  set_target_properties(l_${ANTS_FUNCTION_NAME} PROPERTIES
    SOVERSION ${ANTS_VERSION_MAJOR}  VERSION ${${PROJECT_NAME}_VERSION})
  message(STATUS "${ANTS_FUNCTION_NAME} ${EXTRA_LIBS}")
  target_link_libraries(l_${ANTS_FUNCTION_NAME} antsUtilities ${EXTRA_LIBS} )
  configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )

  if (ANTS_INSTALL_LIBS_ONLY)
    install(TARGETS l_${ANTS_FUNCTION_NAME}
  #   RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  #   COMPONENT RUNTIME_${ANTS_FUNCTION_NAME}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      COMPONENT RUNTIME_${ANTS_FUNCTION_NAME}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      COMPONENT DEVELOPMENT_${ANTS_FUNCTION_NAME}
  )
  else()
    add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
    target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME} )

    install(TARGETS l_${ANTS_FUNCTION_NAME} ${ANTS_FUNCTION_NAME}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      COMPONENT RUNTIME_${ANTS_FUNCTION_NAME}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      COMPONENT RUNTIME_${ANTS_FUNCTION_NAME}
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      COMPONENT DEVELOPMENT_${ANTS_FUNCTION_NAME}
    )
  endif()

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
    ANTSIntegrateVelocityField
    antsJointTensorFusion
    ImageMath
    iMath
    AddNoiseToImage
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
    SimulateDisplacementField
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
  if (BUILD_SHARED_LIBS)
    DYNAMIC_ANTS_BUILD(${ANTS_APP} "")
  else()
    STATIC_ANTS_BUILD(${ANTS_APP} "")
  endif()
endforeach()

if (BUILD_ALL_ANTS_APPS)
  foreach(ANTS_APP ${MORE_ANTS_APPS})
    if (BUILD_SHARED_LIBS)
     DYNAMIC_ANTS_BUILD(${ANTS_APP} "")
    else()
      STATIC_ANTS_BUILD(${ANTS_APP} "")
    endif()
  endforeach()
else()
  foreach(ANTS_APP ${MORE_ANTS_APPS})
    if(ANTS_BUILD_${ANTS_APP})
      if (BUILD_SHARED_LIBS)
        DYNAMIC_ANTS_BUILD(${ANTS_APP} "")
      else()
        STATIC_ANTS_BUILD(${ANTS_APP} "")
      endif()
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

  if (BUILD_SHARED_LIBS)
    DYNAMIC_ANTS_BUILD(${ANTS_APP} "${VTK_LIBRARIES}")
  else()
    STATIC_ANTS_BUILD(${ANTS_APP} "${VTK_LIBRARIES}")
  endif()

endforeach()
endif()

## Build test option
if(BUILD_TESTING)
    add_subdirectory(TestSuite)
endif()
