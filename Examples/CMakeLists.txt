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

macro(STANDARD_ANTS_BUILD ANTS_FUNCTION_NAME EXTRA_LIBS)
  set( ANTS_FUNCTION_NAME ${ANTS_FUNCTION_NAME} )

  add_library(l_${ANTS_FUNCTION_NAME} ${ANTS_FUNCTION_NAME}.cxx ${UI_SOURCES})
  message(STATUS "${ANTS_FUNCTION_NAME} ${EXTRA_LIBS}")
  target_link_libraries(l_${ANTS_FUNCTION_NAME} antsUtilities ${ITK_LIBRARIES} ${EXTRA_LIBS} )
  configure_file( template_for_executables.cxx.in cli_${ANTS_FUNCTION_NAME}.cxx )
  add_executable( ${ANTS_FUNCTION_NAME} cli_${ANTS_FUNCTION_NAME}.cxx )
  target_link_libraries( ${ANTS_FUNCTION_NAME} l_${ANTS_FUNCTION_NAME})

  install(TARGETS l_${ANTS_FUNCTION_NAME} ${ANTS_FUNCTION_NAME}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
  )

endmacro()

set(BASE_ANTS_APPS
  ANTS
  ANTSJacobian
  PrintHeader
  ResetDirection
  # ResetDirection2
  # BoykovGraphCutFilter
  # GeneralLinearModel
  ANTSUseLandmarkImagesToGetAffineTransform
  ANTSUseLandmarkImagesToGetBSplineDisplacementField
  ANTSUseDeformationFieldToGetAffineTransform
  # FDR
  ImageMath
  SurfaceCurvature
  SurfaceBasedSmoothing
  # test
  # WarpImage
  LaplacianThickness
  SetOrigin
  SetSpacing
  SetDirectionByMatrix
  # WarpFunctionalImage
  ConvertScalarImageToRGB
  CreateWarpedGridImage
  # CreateJacobianDeterminantImage
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
  # StudentsTestOnImages
  LabelOverlapMeasures
  MeasureMinMaxMean
  WarpImageMultiTransform
  ComposeMultiTransform
  StackSlices
  MemoryTest
  # ANTSOrientImage
  PermuteFlipImageOrientationAxes
  ImageCompare
  ResampleImageBySpacing
  # ResampleImageByOtherImage
  CopyImageHeaderInformation
  WarpTimeSeriesImageMultiTransform
  ResampleImage
  ExtractSliceFromImage
  ExtractRegionFromImage
  ExtractRegionFromImageByMask
  LabelGeometryMeasures
  PasteImageIntoImage
  TileImages
  CreateImage
  WarpTensorImageMultiTransform
  ReorientTensorImage
  RebaseTensorImage
  N3BiasFieldCorrection
  N4BiasFieldCorrection
  KellySlater
  KellyKapowski
  CreateDTICohort
  AverageAffineTransform
  antsRegistration
  antsApplyTransforms
  antsApplyTransformsToPoints
  CompositeTransformUtil
  antsMotionCorr
  sccan
  CreateDisplacementField
  ANTSIntegrateVectorField
  ANTSIntegrateVelocityField
  AtroposMin
  Atropos
  ConvertTransformFile
)
foreach(ANTS_APP ${BASE_ANTS_APPS})
  STANDARD_ANTS_BUILD(${ANTS_APP} "")
endforeach()

if(USE_VTK)
include(${CMAKE_ROOT}/Modules/FindVTK.cmake)
if(USE_VTK_FILE)
  include(${USE_VTK_FILE})
set(VTK_ANTS_APPS
    ConvertVectorFieldToVTK
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
endif(USE_VTK_FILE)
endif(USE_VTK)

install(PROGRAMS ../Scripts/ants.sh
                 ../Scripts/basic_ants_example.sh
                 ../Scripts/antsIntroduction.sh
                 ../Scripts/antsaffine.sh
                 ../Scripts/waitForXGridJobs.pl
                 ../Scripts/ANTSAverage2DAffine.sh
                 ../Scripts/ANTSAverage3DAffine.sh
                 ../Scripts/buildtemplateparallel.sh
                 ../Scripts/ANTSpexec.sh
                 ../Scripts/waitForSGEQJobs.pl
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
set(CPACK_PACKAGE_VERSION "1.9.y")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "9")
set(CPACK_PACKAGE_VERSION_PATCH "y")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "ANTS")
# This must always be last!
include(CPack)
