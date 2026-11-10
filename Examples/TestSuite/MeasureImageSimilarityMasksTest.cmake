# Verify mask conversion through the CLI and its effect on the metric.
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")
set(fixed_values "")
set(moving_values "")
set(mask_values "")
set(binary_values "")
set(moving_mask_values "")
foreach(row RANGE 1 8)
  string(APPEND fixed_values "0 0 0 0 0 0 0 0 ")
  string(APPEND moving_values "1 2 3 4 1 2 3 4 ")
  string(APPEND mask_values "0.25 -0.25 256 0 0.25 -0.25 256 0 ")
  string(APPEND binary_values "1 1 1 0 1 1 1 0 ")
  string(APPEND moving_mask_values "-256 0.5 0 0 -256 0.5 0 0 ")
endforeach()
foreach(name fixed moving mask binary moving_mask)
  file(WRITE "${TEST_OUTPUT_DIR}/${name}.mha"
    "ObjectType = Image\nNDims = 2\nBinaryData = False\nDimSize = 8 8\nElementType = MET_FLOAT\nElementDataFile = LOCAL\n${${name}_values}\n")
endforeach()

function(check_metric name lower upper)
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=2
      "${MEASURE_IMAGE_SIMILARITY}" -d 2
      -m "MeanSquares[${TEST_OUTPUT_DIR}/fixed.mha,${TEST_OUTPUT_DIR}/moving.mha,1,0,None,1]"
      ${ARGN}
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  string(STRIP "${stdout}" value)
  if(NOT "${result}" STREQUAL "0" OR NOT value MATCHES "^[0-9]+([.][0-9]+)?$")
    message(FATAL_ERROR "${name} failed: ${result}\n${stdout}\n${stderr}")
  endif()
  if(value LESS lower OR value GREATER upper)
    message(FATAL_ERROR "${name}: expected metric in [${lower}, ${upper}], got ${value}")
  endif()
endfunction()

# Unmasked mean square = (1 + 4 + 9 + 16) / 4 = 7.5.
check_metric(unmasked 7.4999 7.5001)
# Fixed mask selects values 1, 2, 3: mean square = 14 / 3.
check_metric(binary_reference 4.6666 4.6668 -x "${TEST_OUTPUT_DIR}/binary.mha")
check_metric(fixed_mask 4.6666 4.6668 -x "${TEST_OUTPUT_DIR}/mask.mha")
check_metric(bracketed_fixed_mask 4.6666 4.6668 -x "[${TEST_OUTPUT_DIR}/mask.mha]")
check_metric(fixed_mask_null_moving 4.6666 4.6668 -x "[${TEST_OUTPUT_DIR}/mask.mha,NULL]")
# Moving mask limits the intersection to values 1 and 2: mean square = 2.5.
check_metric(both_masks 2.4999 2.5001
  -x "[${TEST_OUTPUT_DIR}/mask.mha,${TEST_OUTPUT_DIR}/moving_mask.mha]")
check_metric(moving_mask_only 2.4999 2.5001 -x "[NULL,${TEST_OUTPUT_DIR}/moving_mask.mha]")
