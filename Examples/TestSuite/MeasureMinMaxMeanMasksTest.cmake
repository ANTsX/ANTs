file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")

function(write_image name components values offset)
  file(WRITE "${TEST_OUTPUT_DIR}/${name}.mha"
    "ObjectType = Image\nNDims = 2\nBinaryData = False\nTransformMatrix = -1 0 0 1\nOffset = ${offset}\nElementSpacing = 1.5 2.5\nDimSize = 3 2\nElementNumberOfChannels = ${components}\nElementType = MET_FLOAT\nElementDataFile = LOCAL\n${values}\n")
endfunction()

write_image(scalar 1 "10 5 0 100 100 100" "7 11")
write_image(vector 3 "6 8 0 0 3 4 0 0 0 100 100 100 100 100 100 100 100 100" "7 11")
write_image(mask 1 "0.25 -0.25 256 0 0 0" "7 11")
write_image(binary 1 "1 1 1 0 0 0" "7 11")
write_image(shifted 1 "1 1 1 0 0 0" "8 11")

function(measure image mask output)
  execute_process(
    COMMAND "${MEASURE_MIN_MAX_MEAN}" 2 "${TEST_OUTPUT_DIR}/${image}.mha"
      "${TEST_OUTPUT_DIR}/statistics.log" 0 "${TEST_OUTPUT_DIR}/${mask}.mha"
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(NOT "${result}" STREQUAL "0")
    message(FATAL_ERROR "${image}/${mask} failed: ${stdout}\n${stderr}")
  endif()
  set(${output} "${stdout}" PARENT_SCOPE)
endfunction()

foreach(image scalar vector)
  measure(${image} mask actual)
  measure(${image} binary expected)
  if(NOT actual STREQUAL expected)
    message(FATAL_ERROR "${image}: nonzero mask and binary mask results differ")
  endif()
  if(image STREQUAL "scalar")
    set(pattern "Max : \\[10\\] Min : \\[0\\] Mean : \\[5\\] Var : 16[.]6667")
  else()
    set(pattern "Max : \\[6, 8, 0\\] Min : \\[0, 0, 0\\] Mean : \\[2, 3[.]66667, 1[.]33333\\] Var : 22[.]4444")
  endif()
  if(NOT actual MATCHES "${pattern}")
    message(FATAL_ERROR "${image}: unexpected statistics: ${actual}")
  endif()
endforeach()

foreach(mask vector shifted)
  execute_process(
    COMMAND "${MEASURE_MIN_MAX_MEAN}" 2 "${TEST_OUTPUT_DIR}/vector.mha"
      "${TEST_OUTPUT_DIR}/statistics.log" 0 "${TEST_OUTPUT_DIR}/${mask}.mha"
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if("${result}" STREQUAL "0")
    message(FATAL_ERROR "Accepted invalid ${mask} mask")
  endif()
  if(mask STREQUAL "vector" AND NOT stderr MATCHES "Mask must be a scalar image")
    message(FATAL_ERROR "Missing scalar-mask diagnostic: ${stderr}")
  elseif(mask STREQUAL "shifted" AND NOT stderr MATCHES "Mask must have the same size, origin, spacing, and direction")
    message(FATAL_ERROR "Missing geometry diagnostic: ${stderr}")
  endif()
endforeach()
