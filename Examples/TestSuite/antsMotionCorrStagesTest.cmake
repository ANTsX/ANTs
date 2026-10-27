# Run antsMotionCorr with three stages on a 2D+t series made from r16slice and r64slice.
# The stages name the moving files A, B, A, so the moving series is read again between stages.
#  1. A, A, A and A, Acopy, A must give byte-identical outputs: a series read again
#     gives the same result as a series kept from the previous stage.
#  2. A, C, A (C has other frames and another time spacing) must give other parameters
#     (stage 2 reads C) and the time spacing of A in the output (stage 3 reads A again).
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")

function(run_tool name)
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1 ${ARGN}
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(NOT "${result}" STREQUAL "0")
    message(FATAL_ERROR "${name} failed: ${result}\n${stdout}\n${stderr}")
  endif()
endfunction()

set(D "${TEST_OUTPUT_DIR}")
# TileImages stacks the 2D slices along the third axis, which antsMotionCorr -d 2 reads as time.
run_tool("stack A" "${TILE_IMAGES}" 3 "${D}/A.nii.gz" 1x1x3 "${R16_IMAGE}" "${R64_IMAGE}" "${R16_IMAGE}")
configure_file("${D}/A.nii.gz" "${D}/Acopy.nii.gz" COPYONLY)
run_tool("stack C" "${TILE_IMAGES}" 3 "${D}/C1.nii.gz" 1x1x3 "${R64_IMAGE}" "${R16_IMAGE}" "${R64_IMAGE}")
run_tool("time spacing C" "${IMAGE_MATH}" 3 "${D}/C.nii.gz" SetTimeSpacing "${D}/C1.nii.gz" 2)

function(moco name first second third)
  set(stage_options -t "Rigid[0.1]" -i 10x5 -s 1x0 -f 2x1)
  run_tool("antsMotionCorr ${name}" "${ANTS_MOTION_CORR}" -d 2 -u 1 -w 1 --random-seed 42
    -m "Demons[${R16_IMAGE},${D}/${first},1]" ${stage_options}
    -m "Demons[${R16_IMAGE},${D}/${second},1]" ${stage_options}
    -m "Demons[${R16_IMAGE},${D}/${third},1]" ${stage_options}
    -o "[${D}/${name}_,${D}/${name}_out.mha,${D}/${name}_avg.mha]")
endfunction()

moco(same A.nii.gz A.nii.gz A.nii.gz)
moco(copy A.nii.gz Acopy.nii.gz A.nii.gz)
moco(other A.nii.gz C.nii.gz A.nii.gz)

foreach(output out.mha avg.mha MOCOparams.csv Warp.nii.gz InverseWarp.nii.gz)
  file(SHA256 "${D}/same_${output}" same_hash)
  file(SHA256 "${D}/copy_${output}" copy_hash)
  if(NOT same_hash STREQUAL copy_hash)
    message(FATAL_ERROR "A,A,A and A,Acopy,A give different ${output}")
  endif()
endforeach()

file(SHA256 "${D}/same_MOCOparams.csv" same_hash)
file(SHA256 "${D}/other_MOCOparams.csv" other_hash)
if(same_hash STREQUAL other_hash)
  message(FATAL_ERROR "A,C,A gives the same parameters as A,A,A: stage 2 did not read C")
endif()

file(STRINGS "${D}/same_out.mha" same_spacing REGEX "^ElementSpacing = ")
file(STRINGS "${D}/other_out.mha" other_spacing REGEX "^ElementSpacing = ")
if(NOT same_spacing STREQUAL other_spacing)
  message(FATAL_ERROR "A,C,A output has the time spacing of C: stage 3 did not read A again\n"
    "${same_spacing}\n${other_spacing}")
endif()
