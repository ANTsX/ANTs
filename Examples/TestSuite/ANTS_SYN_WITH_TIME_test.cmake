#####################################################################################
#####################################################################################
set(THIS_TEST_NAME ANTS_SYN_WITH_TIME)
set(OUTPUT_PREFIX ${CMAKE_BINARY_DIR}/TEST_${THIS_TEST_NAME} )
set(WARP ${OUTPUT_PREFIX}Warp.nii.gz ${OUTPUT_PREFIX}Affine.txt )
set(INVERSEWARP -i ${OUTPUT_PREFIX}Affine.txt ${OUTPUT_PREFIX}InverseWarp.nii.gz )
set(WARP_IMAGE ${CMAKE_BINARY_DIR}/${THIS_TEST_NAME}_warped.nii.gz)
set(INVERSEWARP_IMAGE ${CMAKE_BINARY_DIR}/${THIS_TEST_NAME}_inversewarped.nii.gz)


add_test(NAME ${THIS_TEST_NAME} COMMAND ANTS 2
 -m MSQ[${CHALF_IMAGE},${C_IMAGE},1,0] -r Gauss[0.5,0.1]
 -t SyN[1,10,0.05] -i 150x100x2x2 -o ${OUTPUT_PREFIX}
 --geodesic 1 --number-of-affine-iterations 0)

add_test(NAME ${THIS_TEST_NAME}_WARP COMMAND WarpImageMultiTransform 2
 ${C_IMAGE} ${OUTPUT_PREFIX}.nii.gz ${OUTPUT_PREFIX}Warp.nii.gz -R ${CHALF_IMAGE} )
set_property(TEST ${THIS_TEST_NAME}_WARP APPEND PROPERTY DEPENDS ${THIS_TEST_NAME})

# add_test(NAME ${THIS_TEST_NAME}_WARP_METRIC_0 COMMAND MeasureImageSimilarity 2 0
#  ${CHALF_IMAGE} ${OUTPUT_PREFIX}.nii.gz
#  ${OUTPUT_PREFIX}log.txt ${OUTPUT_PREFIX}metric.nii.gz
#  0.0943736 0.1)
# set_property(TEST ${THIS_TEST_NAME}_WARP_METRIC_0 APPEND PROPERTY DEPENDS ${THIS_TEST_NAME}_WARP)
