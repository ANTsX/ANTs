###############################################################################
###############################################################################
## MakeTestDriverFromANTSbinary
## For tools made with the slicer execution model,
## This macro will build a test driver that adds the
## --compare
## --compareIntensityTolerance
## --compareRadiusTolerance
## --compareNumberOfPixelsTolerance
## to the SEM tools.
macro(MakeTestDriverFromANTSbinary ANTSbinaryName ANTSbinaryTestSourceName ANTS_FUNCTION_NAME)

  set(CMAKE_TESTDRIVER_BEFORE_TESTMAIN "#include \"itkTestDriverBeforeTest.inc\"")
  set(CMAKE_TESTDRIVER_AFTER_TESTMAIN "#include \"itkTestDriverAfterTest.inc\"")

  set(ANTS_FUNCTION_NAME ${ANTS_FUNCTION_NAME})
  configure_file( template_for_executableTestWrapper.cxx.in ${ANTSbinaryTestSourceName} )

  include_directories(${ITK_INCLUDE_DIRS})
  create_test_sourcelist(${ANTSbinaryName}   ${ANTSbinaryName}Driver.cxx ${ANTSbinaryTestSourceName}
     EXTRA_INCLUDE itkTestDriverIncludeRequiredIOFactories.h
     FUNCTION  ProcessArgumentsAndRegisterRequiredFactories
     )

  add_executable(${ANTSbinaryName}Driver ${ANTSbinaryName}Driver.cxx ${ANTSbinaryTestSourceName})
  target_link_libraries(${ANTSbinaryName}Driver l_${ANTS_FUNCTION_NAME} )
  set_target_properties(${ANTSbinaryName}Driver PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )
endmacro()

#-----------------------------------------------------------------------
# Setup locations to find externally maintained test data.
#-----------------------------------------------------------------------
include(ANTSExternalData)

# Tell ExternalData commands to transform raw files to content links.
set(ExternalData_LINK_CONTENT SHA512)
set(ExternalData_SOURCE_ROOT ${${PROJECT_NAME}_SOURCE_DIR})
include(ExternalData)


# Set testing environment
set(ANTS_TEST_BIN_DIR ${CMAKE_BINARY_DIR}/Examples)

set(TEST_DATA_DIR ${CMAKE_SOURCE_DIR}/TestData)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/TEST_LOCATION)

## Test antsRegistration
##
file(GLOB_RECURSE INCS ${CMAKE_SOURCE_DIR}/Examples/"*.h")

MakeTestDriverFromANTSbinary(antsRegistrationTest antsRegistrationTest.cxx antsRegistration)
MakeTestDriverFromANTSbinary(antsApplyTransformsTest antsApplyTransformsTest.cxx antsApplyTransforms)


###
#  Perform testing
###
set(TestDataSHA512_DIR ${CMAKE_SOURCE_DIR}/TestData)
set(SMALL_DATA_DIR ${TestDataSHA512_DIR}/Data)

ExternalData_expand_arguments(${PROJECT_NAME}FetchData R16_IMAGE DATA{${SMALL_DATA_DIR}/r16slice.nii.gz})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData R64_IMAGE DATA{${SMALL_DATA_DIR}/r64slice.nii.gz})
###
# PSE sub-tests:  Check to see if .txt files and .vtk files also run correctly
###
ExternalData_expand_arguments(${PROJECT_NAME}FetchData DEVIL_IMAGE DATA{${SMALL_DATA_DIR}/Frown.nii})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData ANGEL_IMAGE DATA{${SMALL_DATA_DIR}/Smile.nii})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData DEVIL_IMAGE_TXT DATA{${SMALL_DATA_DIR}/Frown.txt})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData ANGEL_IMAGE_TXT DATA{${SMALL_DATA_DIR}/Smile.txt})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData DEVIL_IMAGE_VTK DATA{${SMALL_DATA_DIR}/Frown.vtk})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData ANGEL_IMAGE_VTK DATA{${SMALL_DATA_DIR}/Smile.vtk})
#ExternalData_expand_arguments(${PROJECT_NAME}FetchData SEG_IMAGE DATA{${SMALL_DATA_DIR}/nslice.nii.gz})

if(RUN_LONG_TESTS)

## A work around for the crazy ANTS command line that is not compatible
## with the testing framework of external tests

ExternalData_expand_arguments(${PROJECT_NAME}FetchData test_mask_VAR DATA{${TestDataSHA512_DIR}/test_mask.nii.gz})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData scale_test_mask_VAR DATA{${TestDataSHA512_DIR}/scale.test_mask.nii.gz})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData rotation_test_mask_VAR DATA{${TestDataSHA512_DIR}/rotation.test_mask.nii.gz})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData translation_test_mask_VAR DATA{${TestDataSHA512_DIR}/translation.test_mask.nii.gz})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData test_image_VAR DATA{${TestDataSHA512_DIR}/test.nii.gz})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData scale_test_image_VAR DATA{${TestDataSHA512_DIR}/scale.test.nii.gz})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData rotation_test_image_VAR DATA{${TestDataSHA512_DIR}/rotation.test.nii.gz})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData rotation_geom_test_image_VAR DATA{${TestDataSHA512_DIR}/rotation.geom.test.nii.gz})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData translation_test_image_VAR DATA{${TestDataSHA512_DIR}/translation.test.nii.gz})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData translation_rescale_test_image_VAR DATA{${TestDataSHA512_DIR}/translation.rescale.test.nii.gz})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData rotation_rescale_test_image_VAR DATA{${TestDataSHA512_DIR}/rotation.rescale.test.nii.gz})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData ANON0006_20_T1_dbg_splayed_image_VAR DATA{${TestDataSHA512_DIR}/ANON0006_20_T1_dbg_splayed.nii.gz})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData ANON0006_20_T1_sag_twisted_image_VAR DATA{${TestDataSHA512_DIR}/ANON0006_20_T1_sag_twisted.nii.gz})

set(antsRegistrationTestName antsRegistrationTest_AffineScaleMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000

  antsRegistrationTest
  -v --dimensionality 3
  --convergence 25x20x5 #numberOfIterations 2500
  --shrink-factors 1x1x1
  --smoothing-sigmas 0.1x0.1x0.1
  --metric Mattes[${test_image_VAR},${scale_test_image_VAR},1,200,Regular,0.10] #costMetric & fixedVolume & movingVolume & metricWeight? & numberOfHistogramBins & samplingStrategy? & samplingPercentage(numberOfSamples)
  --transform "Affine[0.001]" #transformType & minimumStepLength
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_0.05_${antsRegistrationTestName}.mat}
  --masks [${test_mask_VAR},${scale_test_mask_VAR}] #fixedBinaryVolume & movingBinaryVolume
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz,${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}_inverse.test.nii.gz]
)

set(antsRegistrationTestName antsRegistrationTest_AffineScaleNoMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000

  antsRegistrationTest
  -v --dimensionality 3
  --convergence 25x20x5
  --shrink-factors 1x1x1
  --smoothing-sigmas 0.1x0.1x0.1
  -u
  --metric Mattes[${test_image_VAR},${scale_test_image_VAR},1,200,Regular,0.10]
  --transform "Affine[0.001]" #transformType & minimumStepLength
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_0.05_antsRegistrationTest_AffineScaleMasks.mat}
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

set(antsRegistrationTestName antsRegistrationTest_AffineRotationMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 25
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 1000

  antsRegistrationTest
  -v --dimensionality 3
  --convergence 25x20x5
  --shrink-factors 1x1x1
  --smoothing-sigmas 0.1x0.1x0.1
  --metric Mattes[${test_image_VAR},${rotation_test_image_VAR},1,200,Regular,0.10]
  --transform "Affine[0.001]"
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_0.05_${antsRegistrationTestName}.mat}
  --masks [${test_mask_VAR},${rotation_test_mask_VAR}]
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

set(antsRegistrationTestName antsRegistrationTest_AffineRotationNoMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777

  antsRegistrationTest
  -v --dimensionality 3
  --convergence 25x20x5
  --shrink-factors 1x1x1
  --smoothing-sigmas 0.1x0.1x0.1
  --metric Mattes[${test_image_VAR},${rotation_test_image_VAR},1,200,Regular,0.10]
  --transform "Affine[0.001]"
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_0.05_${antsRegistrationTestName}.mat}
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

set(antsRegistrationTestName antsRegistrationTest_AffineTranslationMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777

  antsRegistrationTest
  -v --dimensionality 3
  --convergence 25x20x5
  --shrink-factors 1x1x1
  --smoothing-sigmas 0.1x0.1x0.1
  --metric MeanSquares[${test_image_VAR},${translation_test_image_VAR},1,radius=NA,Regular,0.10]
  --transform "Affine[0.001]"
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_${antsRegistrationTestName}.mat}
  --masks [${test_mask_VAR},${translation_test_mask_VAR}]
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

set(antsRegistrationTestName antsRegistrationTest_AffineTranslationNoMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 777

  antsRegistrationTest
  -v --dimensionality 3
  --convergence 25x20x5
  --shrink-factors 1x1x1
  --smoothing-sigmas 0.1x0.1x0.1
  --metric MeanSquares[${test_image_VAR},${translation_test_image_VAR},1,radius=NA,Regular,0.10]
  --transform "Affine[0.001]"
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_antsRegistrationTest_AffineTranslationMasks.mat}
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

set(antsRegistrationTestName antsRegistrationTest_RigidAnisotropicMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 1500

  antsRegistrationTest
  -v --dimensionality 3
  --convergence 25x20x5
  --shrink-factors 1x1x1
  --smoothing-sigmas 0.1x0.1x0.1
  --metric Mattes[${ANON0006_20_T1_dbg_splayed_image_VAR},${ANON0006_20_T1_sag_twisted_image_VAR},1,200,Regular,0.20]
  --transform "Rigid[0.001]"
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_antsRegistrationTest_RigidAnisotropicMasks.mat}
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

set(antsRegistrationTestName antsRegistrationTest_RigidRotationHeadMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 777

  antsRegistrationTest
  -v --dimensionality 3
  --convergence 25x20x5
  --shrink-factors 1x1x1
  --smoothing-sigmas 0.1x0.1x0.1
  --metric Mattes[${test_image_VAR},${rotation_test_image_VAR},1,200,Regular,0.20]
  --transform "Rigid[0.001]"
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_0.05_${antsRegistrationTestName}.mat}
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

set(antsRegistrationTestName antsRegistrationTest_RigidRotationMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 777

  antsRegistrationTest
  -v --dimensionality 3
  --convergence 25x20x5
  --shrink-factors 1x1x1
  --smoothing-sigmas 0.1x0.1x0.1
  --metric Mattes[${test_image_VAR},${rotation_test_image_VAR},1,200,Regular,0.20]
  --transform "Rigid[0.001]"
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_0.05_antsRegistrationTest_RigidRotationHeadMasks.mat}
  --masks [${test_mask_VAR},${rotation_test_mask_VAR}]
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

set(antsRegistrationTestName antsRegistrationTest_RigidRotationNoMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 777

  antsRegistrationTest
  -v --dimensionality 3
  --convergence 25x20x5
  --shrink-factors 1x1x1
  --smoothing-sigmas 0.1x0.1x0.1
  --metric Mattes[${test_image_VAR},${rotation_test_image_VAR},1,200,Regular,0.20]
  --transform "Rigid[0.001]"
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_0.05_${antsRegistrationTestName}.mat}
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

set(antsRegistrationTestName antsRegistrationTest_MSEAffineRotationMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 250

  antsRegistrationTest
  -v --dimensionality 3
  --convergence 25x20x5
  --shrink-factors 1x1x1
  --smoothing-sigmas 0.1x0.1x0.1
  --metric MeanSquares[${test_image_VAR},${rotation_test_image_VAR},1,200,Regular,0.10]
  --transform "Affine[0.001]"
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_0.05_antsRegistrationTest_AffineRotationMasks.mat}
  --masks [${test_mask_VAR},${rotation_test_mask_VAR}]
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

set(antsRegistrationTestName antsRegistrationTest_RigidRotGeomNoMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 777

  antsRegistrationTest
  -v --dimensionality 3
  --metric Mattes[${test_image_VAR},${rotation_geom_test_image_VAR},1,200,Regular,0.10]
  --transform "Rigid[0.05]"
  --convergence 25x20x5
  --shrink-factors 3x2x1
  --smoothing-sigmas 0x0x0
  --metric Mattes[${test_image_VAR},${rotation_geom_test_image_VAR},1,200,Regular,0.10]
  --transform "Rigid[0.005]"
  --convergence 25x20x5
  --shrink-factors 3x2x1
  --smoothing-sigmas 0x0x0
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_${antsRegistrationTestName}.mat}
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

set(antsRegistrationTestName antsRegistrationTest_RigidRotaRotaRotNoMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 1000

  antsRegistrationTest
  -v --dimensionality 3
  --metric Mattes[${test_image_VAR},${rotation_test_image_VAR},1,200,Regular,0.10]
  --transform "Rigid[0.05]"
  --convergence 25x20x5
  --shrink-factors 3x2x1
  --smoothing-sigmas 0x0x0
  --metric Mattes[${test_image_VAR},${rotation_test_image_VAR},1,200,Regular,0.10]
  --transform "Rigid[0.001]"
  --convergence 25x20x5
  --shrink-factors 3x2x1
  --smoothing-sigmas 0x0x0
  --metric Mattes[${test_image_VAR},${rotation_test_image_VAR},1,200,Regular,0.10]
  --transform "Rigid[0.00075]"
  --convergence 25x20x5
  --shrink-factors 3x2x1
  --smoothing-sigmas 0x0x0
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_0.05_antsRegistrationTest_RigidRotationNoMasks.mat}
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

set(antsRegistrationTestName antsRegistrationTest_RigidRotaRotaRotNoMasks_Float)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 1000

  antsRegistrationTest
  -v --dimensionality 3
  --float 1
  --metric Mattes[${test_image_VAR},${rotation_test_image_VAR},1,200,Regular,0.10]
  --transform "Rigid[0.05]"
  --convergence 25x20x5
  --shrink-factors 3x2x1
  --smoothing-sigmas 0x0x0
  --metric Mattes[${test_image_VAR},${rotation_test_image_VAR},1,200,Regular,0.10]
  --transform "Rigid[0.001]"
  --convergence 25x20x5
  --shrink-factors 3x2x1
  --smoothing-sigmas 0x0x0
  --metric Mattes[${test_image_VAR},${rotation_test_image_VAR},1,200,Regular,0.10]
  --transform "Rigid[0.00075]"
  --convergence 25x20x5
  --shrink-factors 3x2x1
  --smoothing-sigmas 0x0x0
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_0.05_antsRegistrationTest_RigidRotationNoMasks_Float.mat}
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

set(antsRegistrationTestName antsRegistrationTest_SimilarityRotationMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 11
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000

  antsRegistrationTest
  -v --dimensionality 3
  --convergence 25x20x5
  --shrink-factors 3x2x1
  --smoothing-sigmas 0x0x0
  --metric Mattes[${test_image_VAR},${rotation_test_image_VAR},1,200,Regular,0.20]
  --transform "Similarity[0.001]"
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_0.05_antsRegistrationTest_RigidRotationNoMasks.mat}
  --masks [${test_mask_VAR},${rotation_test_mask_VAR}]
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

# set(antsRegistrationTestName antsRegistrationTest_SimilarityRotationNoMasks)
# ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
#   COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
#   --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
#             ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
#   --compareIntensityTolerance 11
#   --compareRadiusTolerance 0
#   --compareNumberOfPixelsTolerance 1000
#
#   antsRegistrationTest
#   -v --dimensionality 3
#   --convergence 25x20x5
#   --shrink-factors 3x2x1
#   --smoothing-sigmas 0x0x0
#   --metric Mattes[${test_image_VAR},${rotation_test_image_VAR},1,200,Regular,0.20]
#   --transform "Similarity[0.001]"
#   --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_0.05_${antsRegistrationTestName}.mat}
#   --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
# )

# set(antsRegistrationTestName antsRegistrationTest_SimilarityScaleMasks)
# ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
#   COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
#   --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
#             ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
#   --compareIntensityTolerance 9
#   --compareRadiusTolerance 0
#   --compareNumberOfPixelsTolerance 1200
#
#   antsRegistrationTest
#   -v --dimensionality 3
#   --convergence 25x20x5
#   --shrink-factors 3x2x1
#   --smoothing-sigmas 0x0x0
#   --metric Mattes[${test_image_VAR},${scale_test_image_VAR},1,200,Regular,0.20]
#   --transform "Similarity[0.001]"
#   --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_0.05_antsRegistrationTest_AffineScaleMasks.mat}
#   --masks [${test_mask_VAR},${scale_test_mask_VAR}]
#   --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
# )

# set(antsRegistrationTestName antsRegistrationTest_SimilarityScaleNoMasks)
# ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
#   COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
#   --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
#             ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
#   --compareIntensityTolerance 9
#   --compareRadiusTolerance 0
#   --compareNumberOfPixelsTolerance 1000
#
#   antsRegistrationTest
#   -v --dimensionality 3
#   --convergence 25x20x5
#   --shrink-factors 3x2x1
#   --smoothing-sigmas 0x0x0
#   --metric Mattes[${test_image_VAR},${scale_test_image_VAR},1,200,Regular,0.20]
#   --transform "Similarity[0.001]"
#   --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_0.05_${antsRegistrationTestName}.mat}
#   --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
# )

# set(antsRegistrationTestName antsRegistrationTest_SimilarityTranslationRescaleNoMasks)
# ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
#   COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
#   --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
#             ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
#   --compareIntensityTolerance 7
#   --compareRadiusTolerance 0
#   --compareNumberOfPixelsTolerance 777
#
#   antsRegistrationTest
#   -v --dimensionality 3
#   --convergence 25x20x5
#   --shrink-factors 3x2x1
#   --smoothing-sigmas 0x0x0
#   --metric Mattes[${test_image_VAR},${translation_rescale_test_image_VAR},1,200,Regular,0.20]
#   --transform "Similarity[0.001]"
#   --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_antsRegistrationTest_TranslationRescaleHeadMasks.mat}
#   --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
# )

# set(antsRegistrationTestName antsRegistrationTest_SimilarityTranslationRescaleNoMasks_Float)
# ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
#   COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
#   --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
#             ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
#   --compareIntensityTolerance 7
#   --compareRadiusTolerance 0
#   --compareNumberOfPixelsTolerance 777
#
#   antsRegistrationTest
#   -v --dimensionality 3
#   --float 1
#   --convergence 25x20x5
#   --shrink-factors 3x2x1
#   --smoothing-sigmas 0x0x0
#   --metric Mattes[${test_image_VAR},${translation_rescale_test_image_VAR},1,200,Regular,0.20]
#   --transform "Similarity[0.001]"
#   --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_antsRegistrationTest_TranslationRescaleHeadMasks_Float.mat}
#   --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
# )

set(antsRegistrationTestName antsRegistrationTest_CCSimilarityRotationMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 11
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000

  antsRegistrationTest
  -v --dimensionality 3
  --convergence 25x20x5
  --shrink-factors 3x2x1
  --smoothing-sigmas 0x0x0
  --metric CC[${test_image_VAR},${rotation_test_image_VAR},1,5]
  --transform "Similarity[0.001]"
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_0.05_antsRegistrationTest_RigidRotationNoMasks.mat}
  --masks [${test_mask_VAR},${rotation_test_mask_VAR}]
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)


set(antsRegistrationTestName antsRegistrationTest_MSESimilarityRotationMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 11
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000

  antsRegistrationTest
  -v --dimensionality 3
  --convergence 25x20x5
  --shrink-factors 3x2x1
  --smoothing-sigmas 0x0x0
  --metric MeanSquares[${test_image_VAR},${rotation_test_image_VAR},1,NA]
  --transform "Similarity[0.001]"
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_0.05_antsRegistrationTest_RigidRotationNoMasks.mat}
  --masks [${test_mask_VAR},${rotation_test_mask_VAR}]
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

set(antsRegistrationTestName antsRegistrationTest_MSESimilarityRotationMasks_Float)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 11
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000

  antsRegistrationTest
  -v --dimensionality 3
  --float 1
  --convergence 25x20x5
  --shrink-factors 3x2x1
  --smoothing-sigmas 0x0x0
  --metric MeanSquares[${test_image_VAR},${rotation_test_image_VAR},1,NA]
  --transform "Similarity[0.001]"
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_0.05_antsRegistrationTest_RigidRotationNoMasks_Float.mat}
  --masks [${test_mask_VAR},${rotation_test_mask_VAR}]
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

set(antsRegistrationTestName antsRegistrationTest_ScaleTranslationRescaleHeadMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 7
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 777

  antsRegistrationTest
  -v --dimensionality 3
  --metric Mattes[${test_image_VAR},${translation_rescale_test_image_VAR},1,200]
  --transform "Rigid[0.001]"
  --convergence 25x20x5
  --shrink-factors 3x2x1
  --smoothing-sigmas 0x0x0
  --metric Mattes[${test_image_VAR},${translation_rescale_test_image_VAR},1,200]
  --transform "Similarity[0.0001]"
  --convergence 25x20x5
  --shrink-factors 3x2x1
  --smoothing-sigmas 0x0x0
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_antsRegistrationTest_TranslationRescaleHeadMasks.mat}
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

set(antsRegistrationTestName antsRegistrationTest_ScaleRotationRescaleHeadMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 1000

  antsRegistrationTest
  -v --dimensionality 3
  --metric Mattes[${test_image_VAR},${rotation_rescale_test_image_VAR},1,200]
  --transform "Rigid[0.01]"
  --convergence 25x20x5
  --shrink-factors 3x2x1
  --smoothing-sigmas 0x0x0
  --metric Mattes[${test_image_VAR},${rotation_rescale_test_image_VAR},1,200]
  --transform "Similarity[0.003]"
  --convergence 25x20x5
  --shrink-factors 3x2x1
  --smoothing-sigmas 0x0x0
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_${antsRegistrationTestName}.mat}
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

set(antsRegistrationTestName antsRegistrationTest_ScaleRotationRescaleHeadMasks_Float)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 1000

  antsRegistrationTest
  -v --dimensionality 3
  --float 1
  --metric Mattes[${test_image_VAR},${rotation_rescale_test_image_VAR},1,200]
  --transform "Rigid[0.01]"
  --convergence 25x20x5
  --shrink-factors 3x2x1
  --smoothing-sigmas 0x0x0
  --metric Mattes[${test_image_VAR},${rotation_rescale_test_image_VAR},1,200]
  --transform "Similarity[0.003]"
  --convergence 25x20x5
  --shrink-factors 3x2x1
  --smoothing-sigmas 0x0x0
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_${antsRegistrationTestName}.mat}
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

set(antsRegistrationTestName antsRegistrationTest_SyNScaleNoMasks)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 1000

  antsRegistrationTest
  -v --dimensionality 3
  --metric Mattes[${test_image_VAR},${scale_test_image_VAR},1,200]
  --transform "SyN[0.25,3.0,0.0]"
  --convergence 100x70x20
  --shrink-factors 3x2x1
  --smoothing-sigmas 2x1x0
  --use-histogram-matching
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_0.05_antsRegistrationTest_AffineScaleMasks.mat}
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

set(antsRegistrationTestName antsRegistrationTest_SyNScaleNoMasks_Float)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsRegistrationTestName}.result.nii.gz}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 1000

  antsRegistrationTest
  -v --dimensionality 3
  --float 1
  --metric Mattes[${test_image_VAR},${scale_test_image_VAR},1,200]
  --transform "SyN[0.25,3.0,0.0]"
  --convergence 100x70x20
  --shrink-factors 3x2x1
  --smoothing-sigmas 2x1x0
  --use-histogram-matching
  --initial-moving-transform DATA{${TestDataSHA512_DIR}/Initializer_0.05_antsRegistrationTest_AffineScaleMasks_Float.mat}
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)


## Test initializing stages at each iteration, ensure that only 2 stages result.
## This test is not very robust.  It ensures that the command line options
## can be accepted without failure, but no real testing of results is done.
## The output transform should have 2 levels (1 affine, 1 SyN)
## The next test uses "compareTwoTransforms" to compare the output transform to the baseline.
set(antsRegistrationTestName antsRegistrationTest_initialize_transforms_per_stage )
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  antsRegistrationTest
  -v --dimensionality 3
  --float 0
  --initialize-transforms-per-stage
  --write-composite-transform
  --collapse-output-transforms 0
  --metric Mattes[${test_image_VAR},${scale_test_image_VAR},1,200]
  --transform "Rigid[0.01]"
  --convergence 15
  --shrink-factors 3
  --smoothing-sigmas 2
  --metric Mattes[${test_image_VAR},${scale_test_image_VAR},1,200]
  --transform "Affine[0.001]"
  --convergence 15
  --shrink-factors 3
  --smoothing-sigmas 2
  --metric Mattes[${test_image_VAR},${scale_test_image_VAR},1,200]
  --transform "SyN[0.25,3.0,0.0]"
  --convergence 15
  --shrink-factors 3
  --smoothing-sigmas 2
  --metric Mattes[${test_image_VAR},${scale_test_image_VAR},1,200]
  --transform "SyN[0.25,3.0,0.0]"
  --convergence 11
  --shrink-factors 3
  --smoothing-sigmas 2
  --use-histogram-matching
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)

# Compare the output transform of the antsRegistrationTest_initialize_transforms_per_stage with the baseline file
set(antsRegistrationTestName antsRegistrationTest_initializePerStage_ComparisonTesting)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} compareTwoTransformsTestDriver

  compareTwoTransformsTest
  ${CMAKE_CURRENT_BINARY_DIR}/antsRegistrationTest_initialize_transforms_per_stageComposite.h5
  DATA{${TestDataSHA512_DIR}/antsRegistrationTest_initialize_transforms_per_stageComposite.result.h5}
  )
set_property(TEST ${antsRegistrationTestName} APPEND PROPERTY
  DEPENDS antsRegistrationTest_initialize_transforms_per_stage)

###############
#TEST Development SHOULD BE CONTINUED FROM HERE
###############


### checks that simpleSynRegistration runs without failure ####
set(simpleSynRegistrationTestName simpleSynRegistrationTesting)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${simpleSynRegistrationTestName}
  COMMAND ${LAUNCH_EXE} simpleSynRegistrationTestDriver

  simpleSynRegistrationTest
  ${test_image_VAR}
  ${scale_test_image_VAR}
  DATA{${TestDataSHA512_DIR}/Initializer_0.05_antsRegistrationTest_AffineScaleMasks.mat}
  ${CMAKE_CURRENT_BINARY_DIR}/${simpleSynRegistrationTestName}
  )

### antsApplyTransforms test #####
set(antsApplyTransformsTestName antsApplyTransformsTesting)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsApplyTransformsTestName}
  COMMAND ${LAUNCH_EXE} antsApplyTransformsTestDriver
  --compare DATA{${TestDataSHA512_DIR}/${antsApplyTransformsTestName}.result.nii.gz}
  ${CMAKE_CURRENT_BINARY_DIR}/${antsApplyTransformsTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 1000

  antsApplyTransformsTest
  -d 3
  -i ${scale_test_image_VAR} #{moving image}
  -r ${test_image_VAR} #{fix image}
  -n linear #{interpolation type}
  -t DATA{${TestDataSHA512_DIR}/${antsApplyTransformsTestName}_InputWarpTransform.nii.gz} #{Warp transform}
  -t DATA{${TestDataSHA512_DIR}/Initializer_0.05_antsRegistrationTest_AffineScaleMasks.mat} #{Affine Transform}
  -o ${CMAKE_CURRENT_BINARY_DIR}/${antsApplyTransformsTestName}.test.nii.gz #{output file name}
  )

### compare the outputs of the antsRegistration and simpleSynRegistration
# This needs four steps:
# 1- make the output SyN warp transform from antsRegistration
# 2- make the output Syn warp transform from simpleSynRegistration
# 3- apply the produced transforms to antsApplyTransforms
# 4- compare the results of step 3

# we use the results of the previous tests for steps 1 and 2.
# step 3 is done as follows:

set(antsApplyTransformsTestName antsApplyTransformsTestForAntsRegistration)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsApplyTransformsTestName}
  COMMAND ${LAUNCH_EXE} antsApplyTransformsTestDriver

  antsApplyTransformsTest
  -d 3
  -i ${scale_test_image_VAR} #{moving image}
  -r ${test_image_VAR} #{fix image}
  -n linear #{interpolation type}
  -t ${CMAKE_CURRENT_BINARY_DIR}/antsRegistrationTest_SyNScaleNoMasks1Warp.nii.gz #{Warp transform}
  -t DATA{${TestDataSHA512_DIR}/Initializer_0.05_antsRegistrationTest_AffineScaleMasks.mat} #{Affine Transform}
  -o ${CMAKE_CURRENT_BINARY_DIR}/${antsApplyTransformsTestName}.test.nii.gz #{output file name}
  )
set_property(TEST ${antsApplyTransformsTestName} APPEND PROPERTY DEPENDS antsRegistrationTest_SyNScaleNoMasks)

set(antsApplyTransformsTestName antsApplyTransformsTestForSimpleSynRegistration)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsApplyTransformsTestName}
  COMMAND ${LAUNCH_EXE} antsApplyTransformsTestDriver

  antsApplyTransformsTest
  -d 3
  -i ${scale_test_image_VAR} #{moving image}
  -r ${test_image_VAR} #{fix image}
  -n linear #{interpolation type}
  -t ${CMAKE_CURRENT_BINARY_DIR}/simpleSynRegistrationTestingWarp.nii.gz #{Warp transform}
  -t DATA{${TestDataSHA512_DIR}/Initializer_0.05_antsRegistrationTest_AffineScaleMasks.mat} #{Affine Transform}
  -o ${CMAKE_CURRENT_BINARY_DIR}/${antsApplyTransformsTestName}.test.nii.gz #{output file name}
  )
set_property(TEST ${antsApplyTransformsTestName} APPEND PROPERTY
  DEPENDS simpleSynRegistrationTesting)

# step 4: finally the results are compared together
set(antsApplyTransformsTestName simpleSynRegistration_Vs_antsRegistration)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsApplyTransformsTestName}
  COMMAND ${LAUNCH_EXE} antsApplyTransformsTestDriver
  --compare
  ${CMAKE_CURRENT_BINARY_DIR}/antsApplyTransformsTestForAntsRegistration.test.nii.gz
  ${CMAKE_CURRENT_BINARY_DIR}/antsApplyTransformsTestForSimpleSynRegistration.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 1
  --compareNumberOfPixelsTolerance 1000

  antsApplyTransformsTest -h
)
set_property(TEST simpleSynRegistration_Vs_antsRegistration APPEND PROPERTY
  DEPENDS antsApplyTransformsTestForSimpleSynRegistration antsApplyTransformsTestForAntsRegistration )
############ End of four step compare tests#################

### checks that CompositeTransformUtil runs without failure using --assemble flag####
# set(CompositeTransformUtilTestName CompositeTransformUtilTest)
# ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${CompositeTransformUtilTestName}
#  COMMAND ${LAUNCH_EXE} CompositeTransformUtilTestDriver

#  CompositeTransformUtilTest
#  --assemble ${CMAKE_CURRENT_BINARY_DIR}/${CompositeTransformUtilTestName}.test.h5
#  DATA{${TestDataSHA512_DIR}/${CompositeTransformUtilTestName}_RigidTransform.mat}
#  DATA{${TestDataSHA512_DIR}/${CompositeTransformUtilTestName}_AffineTransform.mat}
#  DATA{${TestDataSHA512_DIR}/${CompositeTransformUtilTestName}_SyNTransform.nii.gz}
#  )

# Compare the results of the CompositeTransformUtil with the baseline file
#set(CompositeTransformUtilTestName CompositeTransformUtil_ComparisonTesting)
# ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${CompositeTransformUtilTestName}
#  COMMAND ${LAUNCH_EXE} compareTwoTransformsTestDriver

#  compareTwoTransformsTest
#  ${CMAKE_CURRENT_BINARY_DIR}/CompositeTransformUtilTest.test.h5
#  DATA{${TestDataSHA512_DIR}/CompositeTransformUtilTest.result.h5}
#  )
# set_property(TEST ${CompositeTransformUtilTestName} APPEND PROPERTY
#  DEPENDS CompositeTransformUtilTest)

# Tests CompositeTransformUtil when it is run by --disassemble flag
# set(CompositeTransformUtilTestName CompositeTransformUtilTest_disassemble)
# ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${CompositeTransformUtilTestName}
#  COMMAND ${LAUNCH_EXE} CompositeTransformUtilTestDriver
#  --compare
#  DATA{${TestDataSHA512_DIR}/CompositeTransformUtilTest_SyNTransform.nii.gz}
#  ${CMAKE_CURRENT_BINARY_DIR}/02_${CompositeTransformUtilTestName}_DisplacementFieldTransform.nii.gz
#  --compareIntensityTolerance 9
#  --compareRadiusTolerance 1
#  --compareNumberOfPixelsTolerance 1000

#  CompositeTransformUtilTest
#  --disassemble DATA{${TestDataSHA512_DIR}/CompositeTransformUtilTest.result.h5}
#  ${CompositeTransformUtilTestName}
#  )

## Based on conversations with Nick and Brian, it is strange that
## these ToItself functions fail.  It is likely in the ITK scales estimators
## that the problem exists.  It will be addressed later.
if(0) #HACK -- REMOVING THESE TWO FAILING TEST
set(antsRegistrationTestName antsRegistrationTest_CCMetricBrainToItself)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare ${test_image_VAR}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000

  antsRegistrationTest
  -v --dimensionality 3
  --convergence 25x20x5
  --shrink-factors 3x2x1
  --smoothing-sigmas 0x0x0
  --metric CC[${test_image_VAR},${test_image_VAR},1,5]
  --transform "Similarity[0.001]"
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)
set(antsRegistrationTestName antsRegistrationTest_MSEMetricBrainToItself)
ExternalData_add_test( ${PROJECT_NAME}FetchData NAME ${antsRegistrationTestName}
  COMMAND ${LAUNCH_EXE} antsRegistrationTestDriver
  --compare ${test_image_VAR}
            ${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz
  --compareIntensityTolerance 9
  --compareRadiusTolerance 0
  --compareNumberOfPixelsTolerance 1000

  antsRegistrationTest
  -v --dimensionality 3
  --convergence 25x20x5
  --shrink-factors 3x2x1
  --smoothing-sigmas 0x0x0
  --metric MeanSquares[${test_image_VAR},${test_image_VAR},1,NA]
  --transform "Similarity[0.001]"
  --output [${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName},${CMAKE_CURRENT_BINARY_DIR}/${antsRegistrationTestName}.test.nii.gz]
)
endif()

############ End of CompositeTransformUtil tests#################

endif()

##############################################################
##############################################################
##############################################################
##############################################################
if(RUN_SHORT_TESTS)
###
#  ANTS metric testing
###

set(AllANTSPrograms
 ANTS
 antsApplyTransforms
 antsApplyTransformsToPoints
 ANTSIntegrateVectorField
 ANTSIntegrateVelocityField
 ANTSJacobian
 antsMotionCorr
 antsRegistration
 ANTSUseDeformationFieldToGetAffineTransform
 ANTSUseLandmarkImagesToGetAffineTransform
 ANTSUseLandmarkImagesToGetBSplineDisplacementField
 Atropos
 AverageAffineTransform
 AverageImages
 AverageTensorImages
 ClusterImageStatistics
 ComposeMultiTransform
 CompositeTransformUtil
 ConvertImagePixelType
 ConvertScalarImageToRGB
 ConvertToJpg
 CopyImageHeaderInformation
 CreateDisplacementField
 CreateDTICohort
 CreateImage
 CreateWarpedGridImage
 ExtractRegionFromImage
 ExtractRegionFromImageByMask
 ExtractSliceFromImage
 ImageCompare
 ImageMath
 ImageSetStatistics
 KellyKapowski
 KellySlater
 LabelClustersUniquely
 LabelGeometryMeasures
 LabelOverlapMeasures
 LaplacianThickness
 MeasureImageSimilarity
 MeasureMinMaxMean
 MemoryTest
 MultiplyImages
 N3BiasFieldCorrection
 N4BiasFieldCorrection
 PasteImageIntoImage
 PermuteFlipImageOrientationAxes
 PrintHeader
 RebaseTensorImage
 ReorientTensorImage
 ResampleImage
 ResampleImageBySpacing
 ResetDirection
 sccan
 SetDirectionByMatrix
 SetOrigin
 SetSpacing
 SmoothImage
 StackSlices
 SurfaceBasedSmoothing
 SurfaceCurvature
 ThresholdImage
 TileImages
 WarpImageMultiTransform
 WarpTensorImageMultiTransform
 WarpTimeSeriesImageMultiTransform
)

foreach(CurrProg ${AllANTSPrograms})
  set(HELP_FLAG "--help")
  add_test(NAME ${CurrProg}_HELP_LONG  COMMAND ${CurrProg} ${HELP_FLAG} ) ## Just print the help screen
  set(HELP_FLAG "-h" )
  add_test(NAME ${CurrProg}_HELP_SHORT COMMAND ${CurrProg} ${HELP_FLAG} ) ## Just print the help screen
  set_property(TEST ${CurrProg}_HELP_SHORT APPEND PROPERTY DEPENDS ${CurrProg}_HELP_LONG)
endforeach()

###
#  ANTS transform testing
###
# include(ANTS_CC_1_test.cmake)
# include(ANTS_CC_2_test.cmake)
# include(ANTS_CC_3_test.cmake)
# include(ANTS_MSQ_test.cmake)
# include(ANTS_MI_1_test.cmake)
# include(ANTS_MI_2_test.cmake)
# include(ANTS_ELASTIC_test.cmake)
# include(ANTS_GSYN_test.cmake)
# include(ANTS_EXP_test.cmake)
# include(ANTS_SYN_test.cmake)

###
#  ANTS labeled data testing
###
if(RUN_LONG_TESTS)
#   include(ANTS_PSE_MSQ_IMG_test.cmake)

###
#  ANTS images with non-trival rotation header test
###
#ExternalData_expand_arguments(${PROJECT_NAME}FetchData ROT_REF_IMAGE DATA{${SMALL_DATA_DIR}/ref2.nii.gz})
#ExternalData_expand_arguments(${PROJECT_NAME}FetchData ROT_MOV_IMAGE DATA{${SMALL_DATA_DIR}/mov2.nii.gz})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData ROT_REF_IMAGE DATA{${SMALL_DATA_DIR}/r16roth.nii.gz})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData ROT_MOV_IMAGE DATA{${SMALL_DATA_DIR}/r64roth.nii.gz})

  include(ANTS_ROT_GSYN_test.cmake)
  include(ANTS_ROT_EXP_test.cmake)
#add_test(NAME ANTS_ROT_EXP_CLEAN COMMAND rm ${ROT_WARP_FILES})
###
#  Test SyN with time
###
ExternalData_expand_arguments(${PROJECT_NAME}FetchData CHALF_IMAGE DATA{${SMALL_DATA_DIR}/chalf.nii.gz})
ExternalData_expand_arguments(${PROJECT_NAME}FetchData C_IMAGE DATA{${SMALL_DATA_DIR}/c.nii.gz})

include(ANTS_SYN_WITH_TIME_test.cmake)
include(APOC_OTSU_INIT_test.cmake)
endif()

# include(ANTS_PSE_MSQ_TXT_test.cmake)
# include(ANTS_PSE_MSQ_VTK_test.cmake)

endif()

ExternalData_add_target( ${PROJECT_NAME}FetchData )  # Name of data management target
