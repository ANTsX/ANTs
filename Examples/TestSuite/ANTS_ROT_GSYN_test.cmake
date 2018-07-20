#####################################################################################
#####################################################################################
set(THIS_TEST_NAME ANTS_ROT_GSYN)
set(OUTPUT_PREFIX ${CMAKE_BINARY_DIR}/TEST_${THIS_TEST_NAME} )
set(WARP ${OUTPUT_PREFIX}Warp.nii.gz ${OUTPUT_PREFIX}Affine.txt )
set(INVERSEWARP -i ${OUTPUT_PREFIX}Affine.txt ${OUTPUT_PREFIX}InverseWarp.nii.gz )
set(WARP_IMAGE ${CMAKE_BINARY_DIR}/${THIS_TEST_NAME}_warped.nii.gz)
set(INVERSEWARP_IMAGE ${CMAKE_BINARY_DIR}/${THIS_TEST_NAME}_inversewarped.nii.gz)

add_test(NAME ${THIS_TEST_NAME} COMMAND ANTS 3
 -m MSQ[${ROT_REF_IMAGE},${ROT_MOV_IMAGE},1,0]
 -t SyN[0.25] -i 50x5x0 -r Gauss[3,0.0,32]
 -o ${OUTPUT_PREFIX}.nii.gz)

add_test(NAME ${THIS_TEST_NAME}_WARP COMMAND WarpImageMultiTransform 3
 ${ROT_MOV_IMAGE} ${WARP_IMAGE} -R ${ROT_REF_IMAGE} ${WARP} )
set_property(TEST ${THIS_TEST_NAME}_WARP APPEND PROPERTY DEPENDS ${THIS_TEST_NAME})

# add_test(NAME ${THIS_TEST_NAME}_WARP_METRIC_0 COMMAND MeasureImageSimilarity 3 0
#  ${ROT_REF_IMAGE} ${WARP_IMAGE}
#  ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz
#  35.4473 0.5)
# set_property(TEST ${THIS_TEST_NAME}_WARP_METRIC_0 APPEND PROPERTY DEPENDS ${THIS_TEST_NAME}_WARP)

add_test(NAME ${THIS_TEST_NAME}_INVERSEWARP COMMAND WarpImageMultiTransform 3
 ${ROT_REF_IMAGE} ${INVERSEWARP_IMAGE} -R ${ROT_MOV_IMAGE} ${INVERSEWARP})
set_property(TEST ${THIS_TEST_NAME}_INVERSEWARP APPEND PROPERTY DEPENDS ${THIS_TEST_NAME})

# add_test(NAME ${THIS_TEST_NAME}_INVERSEWARP_METRIC_0 COMMAND MeasureImageSimilarity 3 0
#  ${ROT_MOV_IMAGE} ${INVERSEWARP_IMAGE}
#  ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz
#  35.5439 0.5)
# set_property(TEST ${THIS_TEST_NAME}_INVERSEWARP_METRIC_0 APPEND PROPERTY DEPENDS ${THIS_TEST_NAME}_INVERSEWARP)
