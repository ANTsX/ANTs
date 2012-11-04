
##
# Apocrita tests
##
set(R16_MASK ${DATA_DIR}/r16mask.nii.gz)
set(R16_PRIORS ${DATA_DIR}/r16priors.nii.gz)

#add_test(NAME APOC_OTSU_INIT COMMAND $<TARGET_FILE:Apocrita> 2 -i otsu[${R16_IMAGE},3] -x ${R16_MASK} -n 10 -m [0.3,1x1,0.2,0.1] -o ${OUTPUT_PREFIX}.nii.gz )
#add_test(NAME APOC_OTSU_INIT_RADIUS_2x2 COMMAND $<TARGET_FILE:Apocrita> 2 -i otsu[${R16_IMAGE},3] -x ${R16_MASK} -n 10 -m [0.3,2,0.2,0.1] -o ${OUTPUT_PREFIX}.nii.gz )
#add_test(NAME APOC_KMEANS_INIT COMMAND $<TARGET_FILE:Apocrita> 2 -i kmeans[${R16_IMAGE},3] -x ${R16_MASK} -n 10 -m [0.3,1x1,0.2,0.1] -o [${OUTPUT_PREFIX}.nii.gz,${OUTPUT_PREFIX}_posteriors%d.nii.gz])
#add_test(NAME APOC_PRIORLABELIMAGE_INIT COMMAND $<TARGET_FILE:Apocrita> 2 -i priorlabelimage[${R16_IMAGE},5,${R16_PRIORS},0.5] -x ${R16_MASK} -n 10 -m [0.3,1x1,0.2,0.1] -o [${OUTPUT_PREFIX}.nii.gz,${OUTPUT_PREFIX}_posteriors%d.nii.gz] -l 1[1,0.75] -l 2[1,1.0] -l 3[0.5,0.5] -l 4[1,1])
