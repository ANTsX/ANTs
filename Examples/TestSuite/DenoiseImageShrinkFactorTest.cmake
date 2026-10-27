# DenoiseImage does not shrink the input at shrink factors 0 and 1. The test
# uses an oblique image, made here, so the test needs no data files. It checks:
# - at -s 1 the noise image is exactly input - denoised (no resample),
# - -s 0 and -s 1 give identical denoised and noise images,
# - a negative shrink factor makes DenoiseImage fail.
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")

set(input "${TEST_OUTPUT_DIR}/oblique.mha")
set(values "")
foreach(i RANGE 0 1727)
  math(EXPR value "100 + (${i} * 37) % 101")
  string(APPEND values "${value}.5 ")
endforeach()
# 12x12x12 voxels, rotated 30 degrees about z, anisotropic spacing. All values
# are in [100.5, 200.5], so each denoised value is within a factor of 2 of its
# input value. Then input - noise is exact in float, and the check of
# input - denoised against the noise image can require 0.
file(WRITE "${input}"
  "ObjectType = Image\nNDims = 3\nBinaryData = False\n"
  "TransformMatrix = 0.8660254 0.5 0 -0.5 0.8660254 0 0 0 1\n"
  "Offset = 3.7 -11.3 25.1\nElementSpacing = 1.1 0.9 1.3\n"
  "DimSize = 12 12 12\nElementType = MET_FLOAT\nElementDataFile = LOCAL\n${values}\n")

function(run)
  execute_process(COMMAND "${CMAKE_COMMAND}" -E env ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1 ${ARGN}
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(NOT "${result}" STREQUAL "0")
    message(FATAL_ERROR "Command failed (${result}): ${ARGN}\n${stdout}\n${stderr}")
  endif()
  set(stdout "${stdout}" PARENT_SCOPE)
endfunction()

function(denoise shrink)
  run("${DENOISE_IMAGE}" -d 3 -i "${input}" -s ${shrink} -p 1 -r 1
    -o "[${TEST_OUTPUT_DIR}/s${shrink}_denoised.mha,${TEST_OUTPUT_DIR}/s${shrink}_noise.mha]")
endfunction()

function(compare_outputs name)
  file(SHA256 "${TEST_OUTPUT_DIR}/s0_${name}.mha" first_hash)
  file(SHA256 "${TEST_OUTPUT_DIR}/s1_${name}.mha" second_hash)
  if(NOT first_hash STREQUAL second_hash)
    message(FATAL_ERROR "DenoiseImage -s 0 and -s 1 give different ${name} images")
  endif()
endfunction()

denoise(0)
denoise(1)

set(difference "${TEST_OUTPUT_DIR}/s1_input_minus_denoised.mha")
set(error "${TEST_OUTPUT_DIR}/s1_error.mha")
set(absolute_error "${TEST_OUTPUT_DIR}/s1_absolute_error.mha")
run("${IMAGE_MATH}" 3 "${difference}" - "${input}" "${TEST_OUTPUT_DIR}/s1_denoised.mha")
run("${IMAGE_MATH}" 3 "${error}" - "${difference}" "${TEST_OUTPUT_DIR}/s1_noise.mha")
run("${IMAGE_MATH}" 3 "${absolute_error}" abs "${error}")
run("${MEASURE_MIN_MAX_MEAN}" 3 "${absolute_error}")
if(NOT stdout MATCHES "Max : \\[?([-+0-9.eE]+)")
  message(FATAL_ERROR "Cannot read the maximum from: ${stdout}")
endif()
if(NOT CMAKE_MATCH_1 EQUAL 0)
  message(FATAL_ERROR "At -s 1 the noise image is not input - denoised: max difference = ${CMAKE_MATCH_1}")
endif()

compare_outputs(denoised)
compare_outputs(noise)

execute_process(
  COMMAND "${DENOISE_IMAGE}" -d 3 -i "${input}" -s -1 -o "${TEST_OUTPUT_DIR}/negative.mha"
  RESULT_VARIABLE result OUTPUT_QUIET ERROR_QUIET)
if("${result}" STREQUAL "0")
  message(FATAL_ERROR "DenoiseImage did not fail with shrink factor -1")
endif()
