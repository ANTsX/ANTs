file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")
set(input_values "")
set(mask_values "")
set(binary_values "")
set(full_values "")
set(nonzero_values 0.25 -0.25 256)
foreach(i RANGE 0 215)
  math(EXPR x "${i} % 6")
  math(EXPR y "(${i} / 6) % 6")
  math(EXPR z "${i} / 36")
  math(EXPR intensity "1 + (${i} * 17) % 97")
  string(APPEND input_values "${intensity} ")
  string(APPEND full_values "1 ")
  if(x GREATER 0 AND x LESS 5 AND y GREATER 0 AND y LESS 5 AND z GREATER 0 AND z LESS 5)
    math(EXPR index "${i} % 3")
    list(GET nonzero_values ${index} value)
    string(APPEND mask_values "${value} ")
    string(APPEND binary_values "1 ")
  else()
    string(APPEND mask_values "0 ")
    string(APPEND binary_values "0 ")
  endif()
endforeach()
foreach(name input mask binary full)
  file(WRITE "${TEST_OUTPUT_DIR}/${name}.mha"
    "ObjectType = Image\nNDims = 3\nBinaryData = False\nDimSize = 6 6 6\nElementType = MET_FLOAT\nElementDataFile = LOCAL\n${${name}_values}\n")
endforeach()

function(run)
  execute_process(COMMAND "${CMAKE_COMMAND}" -E env ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=2 ${ARGV}
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(NOT "${result}" STREQUAL "0")
    message(FATAL_ERROR "Command failed: ${ARGV}\n${stdout}\n${stderr}")
  endif()
endfunction()

function(compare first second)
  file(SHA256 "${TEST_OUTPUT_DIR}/${first}.mha" first_hash)
  file(SHA256 "${TEST_OUTPUT_DIR}/${second}.mha" second_hash)
  if(NOT first_hash STREQUAL second_hash)
    message(FATAL_ERROR "Different outputs: ${first} and ${second}")
  endif()
endfunction()

foreach(mask mask binary)
  run("${CONVERT_RGB}" 3 "${TEST_OUTPUT_DIR}/input.mha" "${TEST_OUTPUT_DIR}/rgb-${mask}.mha"
    "${TEST_OUTPUT_DIR}/${mask}.mha" red)
endforeach()
compare(rgb-mask rgb-binary)
run("${CONVERT_RGB}" 3 "${TEST_OUTPUT_DIR}/input.mha" "${TEST_OUTPUT_DIR}/rgb-full.mha" none red)
file(SHA256 "${TEST_OUTPUT_DIR}/rgb-full.mha" full_hash)
file(SHA256 "${TEST_OUTPUT_DIR}/rgb-mask.mha" masked_hash)
if(full_hash STREQUAL masked_hash)
  message(FATAL_ERROR "RGB mask did not affect the output")
endif()

foreach(mask mask binary)
  # Main RGB mask, including mask-derived cropping and slice selection.
  run("${CREATE_MOSAIC}" -i "${TEST_OUTPUT_DIR}/input.mha" -r "${TEST_OUTPUT_DIR}/rgb-full.mha"
    -x "${TEST_OUTPUT_DIR}/${mask}.mha" -o "${TEST_OUTPUT_DIR}/main-${mask}.mha"
    -d 2 -p mask+0 -s "[1,mask,mask]")
  # Functional overlay mask: exercise byte slice extraction, padding and flipping.
  run("${CREATE_MOSAIC}" -i "${TEST_OUTPUT_DIR}/input.mha"
    -e "[${TEST_OUTPUT_DIR}/rgb-full.mha,${TEST_OUTPUT_DIR}/${mask}.mha,1]"
    -o "${TEST_OUTPUT_DIR}/overlay-${mask}.mha" -d 2 -p 1 -f 1x0)
endforeach()
compare(main-mask main-binary)
compare(overlay-mask overlay-binary)

# Without an explicit mask, the RGB overlay gets an all-ones byte mask.
run("${CREATE_MOSAIC}" -i "${TEST_OUTPUT_DIR}/input.mha" -r "${TEST_OUTPUT_DIR}/rgb-full.mha"
  -o "${TEST_OUTPUT_DIR}/implicit-full.mha" -d 2)
run("${CREATE_MOSAIC}" -i "${TEST_OUTPUT_DIR}/input.mha" -r "${TEST_OUTPUT_DIR}/rgb-full.mha"
  -x "${TEST_OUTPUT_DIR}/full.mha" -o "${TEST_OUTPUT_DIR}/explicit-full.mha" -d 2)
compare(implicit-full explicit-full)
