# Run DenoiseImage with a mask and a shrink factor of 2 on a 2D slice.
# The mask is an Otsu threshold of the input image. The test checks that:
# - the run succeeds,
# - the noise image is zero outside the mask and not zero inside the mask,
# - a missing mask file makes DenoiseImage fail.
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")
set(mask "${TEST_OUTPUT_DIR}/mask.nii.gz")
set(outside "${TEST_OUTPUT_DIR}/outside.nii.gz")
set(denoised "${TEST_OUTPUT_DIR}/denoised.nii.gz")
set(noise "${TEST_OUTPUT_DIR}/noise.nii.gz")

function(run)
  execute_process(COMMAND ${ARGN}
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(NOT "${result}" STREQUAL "0")
    message(FATAL_ERROR "Command failed (${result}): ${ARGN}\n${stdout}\n${stderr}")
  endif()
  set(stdout "${stdout}" PARENT_SCOPE)
endfunction()

# Get the maximum absolute value of the noise image in a region.
function(max_abs_noise region result_var)
  run("${MEASURE_MIN_MAX_MEAN}" 2 "${noise}" "${TEST_OUTPUT_DIR}/minmax.txt" 1 "${region}")
  if(NOT stdout MATCHES "Max : \\[?([-+0-9.eE]+)")
    message(FATAL_ERROR "Cannot read the maximum from: ${stdout}")
  endif()
  set(${result_var} "${CMAKE_MATCH_1}" PARENT_SCOPE)
endfunction()

run("${THRESHOLD_IMAGE}" 2 "${INPUT_IMAGE}" "${mask}" Otsu 1)
run("${THRESHOLD_IMAGE}" 2 "${mask}" "${outside}" 0 0 1 0)

run("${DENOISE_IMAGE}" -d 2 -i "${INPUT_IMAGE}" -x "${mask}" -s 2 -o "[${denoised},${noise}]")

max_abs_noise("${outside}" outside_max)
if(NOT outside_max EQUAL 0)
  message(FATAL_ERROR "The noise image is not zero outside the mask: max |noise| = ${outside_max}")
endif()

max_abs_noise("${mask}" inside_max)
if(NOT inside_max GREATER 0)
  message(FATAL_ERROR "The noise image is zero inside the mask: no voxel was denoised")
endif()

execute_process(
  COMMAND "${DENOISE_IMAGE}" -d 2 -i "${INPUT_IMAGE}" -x "${TEST_OUTPUT_DIR}/missing.nii.gz" -s 2
    -o "${TEST_OUTPUT_DIR}/missing_out.nii.gz"
  RESULT_VARIABLE result OUTPUT_QUIET ERROR_QUIET)
if("${result}" STREQUAL "0")
  message(FATAL_ERROR "DenoiseImage did not fail with a missing mask file")
endif()
