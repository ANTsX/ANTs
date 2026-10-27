# antsMotionCorr writes the linear stages of each time point to MOCOparams.csv as one composed
# transform. This test makes a 2D+t series from two slices, runs two linear stages, and applies
# the CSV row of each frame with antsApplyTransforms. The result must agree with the
# motion-corrected frame. A deformable stage between two linear stages must give an error.
file(REMOVE_RECURSE "${TEST_OUTPUT_DIR}")
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")
set(ENV{ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS} 1)

function(run_tool)
  execute_process(COMMAND ${ARGN} RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(NOT "${result}" STREQUAL "0")
    message(FATAL_ERROR "Command failed (${result}): ${ARGN}\n${stdout}\n${stderr}")
  endif()
  set(tool_output "${stdout}" PARENT_SCOPE)
endfunction()

# Frame 0 is the fixed image and frame 1 is the moving image. Frame 0 is also the reference.
set(series "${TEST_OUTPUT_DIR}/series.nii.gz")
run_tool("${STACK_SLICES}" "${series}" -1 -1 0 "${FIXED_IMAGE}" "${MOVING_IMAGE}")
foreach(frame RANGE 1)
  run_tool("${EXTRACT_SLICE}" 3 "${series}" "${TEST_OUTPUT_DIR}/series_${frame}.nii.gz" 2 ${frame})
endforeach()

set(metric "MI[${TEST_OUTPUT_DIR}/series_0.nii.gz,${series},1,32,Regular,0.5]")
set(stage_options -i 20x10 -s 1x0 -f 2x1)

# Run antsMotionCorr with the stages in ARGN and check each CSV row against the output frames.
# transform_type is the ITK transform that the CSV columns hold.
function(check_stages name transform_type number_of_parameters)
  set(prefix "${TEST_OUTPUT_DIR}/${name}_")
  set(corrected "${TEST_OUTPUT_DIR}/${name}.nii.gz")
  run_tool("${ANTS_MOTION_CORR}" -d 2 -u 1 --random-seed 42 -o "[${prefix},${corrected}]" ${ARGN})

  file(STRINGS "${prefix}MOCOparams.csv" rows)
  list(POP_FRONT rows header)
  string(REPLACE "," ";" header "${header}")
  list(LENGTH header number_of_columns)
  math(EXPR expected_columns "${number_of_parameters} + 2")
  if(NOT number_of_columns EQUAL expected_columns)
    message(FATAL_ERROR "${name}: expected ${expected_columns} CSV columns, got ${number_of_columns}")
  endif()

  foreach(frame RANGE 1)
    list(GET rows ${frame} row)
    string(REPLACE "," ";" row "${row}")
    list(SUBLIST row 2 -1 parameters)
    list(JOIN parameters " " parameters)
    set(transform "${TEST_OUTPUT_DIR}/${name}_${frame}.txt")
    file(WRITE "${transform}" "#Insight Transform File V1.0\n#Transform 0\nTransform: ${transform_type}\n"
                              "Parameters: ${parameters}\nFixedParameters: 0 0\n")
    set(moving "${TEST_OUTPUT_DIR}/series_${frame}.nii.gz")
    set(output_frame "${TEST_OUTPUT_DIR}/${name}_${frame}_output.nii.gz")
    set(csv_frame "${TEST_OUTPUT_DIR}/${name}_${frame}_csv.nii.gz")
    set(difference "${TEST_OUTPUT_DIR}/${name}_${frame}_difference.nii.gz")
    run_tool("${EXTRACT_SLICE}" 3 "${corrected}" "${output_frame}" 2 ${frame})
    run_tool("${ANTS_APPLY_TRANSFORMS}" -d 2 -i "${moving}" -r "${output_frame}" -t "${transform}" -n Linear
             -o "${csv_frame}")
    run_tool("${IMAGE_MATH}" 2 "${difference}" - "${csv_frame}" "${output_frame}")
    run_tool("${MEASURE_MIN_MAX_MEAN}" 2 "${difference}" "${TEST_OUTPUT_DIR}/${name}_difference.log" 1)
    if(NOT tool_output MATCHES "Max : \\[([^]]*)\\]")
      message(FATAL_ERROR "${name}: cannot read the difference of frame ${frame}: ${tool_output}")
    endif()
    set(max_difference "${CMAKE_MATCH_1}")
    if(NOT max_difference LESS 1e-4)
      message(FATAL_ERROR "${name}: frame ${frame} from the CSV differs from the output by ${max_difference}")
    endif()
  endforeach()
endfunction()

check_stages(RigidRigid Euler2DTransform_double_2_2 3
  -m "${metric}" -t Rigid[0.1] ${stage_options}
  -m "${metric}" -t Rigid[0.1] -i 10 -s 0 -f 1)
check_stages(RigidAffine AffineTransform_double_2_2 6
  -m "${metric}" -t Rigid[0.1] ${stage_options}
  -m "${metric}" -t Affine[0.1] -i 10 -s 0 -f 1)

execute_process(
  COMMAND "${ANTS_MOTION_CORR}" -d 2 -u 1 -o "[${TEST_OUTPUT_DIR}/RigidGDFRigid_]"
    -m "${metric}" -t Rigid[0.1] ${stage_options}
    -m "${metric}" -t GaussianDisplacementField[0.15,3,0.5] ${stage_options}
    -m "${metric}" -t Rigid[0.1] ${stage_options}
  RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
if("${result}" STREQUAL "0" OR NOT stderr MATCHES "A deformable stage is between linear stages 1 and 3")
  message(FATAL_ERROR "Rigid, GaussianDisplacementField, Rigid must fail (${result}):\n${stdout}\n${stderr}")
endif()
