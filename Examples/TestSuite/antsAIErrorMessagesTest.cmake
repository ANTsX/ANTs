# antsAI must say why it fails without -v. Several input errors printed their
# message only with -v, so antsAI failed with no output. The test checks each
# of them without -v. The image is made here, so the test needs no data files.
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")

set(values "")
foreach(y RANGE 15)
  foreach(x RANGE 15)
    math(EXPR value "(${x} * 7 + ${y} * 3) % 17")
    string(APPEND values "${value} ")
  endforeach()
endforeach()
set(image "${TEST_OUTPUT_DIR}/image.mha")
file(WRITE "${image}"
  "ObjectType = Image\nNDims = 2\nBinaryData = False\nElementSpacing = 1 1\n"
  "DimSize = 16 16\nElementType = MET_FLOAT\nElementDataFile = LOCAL\n${values}\n")

# Runs antsAI with the given arguments and checks that it fails with a message
# that matches pattern.
function(expect_error pattern)
  execute_process(COMMAND "${ANTS_AI}" ${ARGN} -o "${TEST_OUTPUT_DIR}/out.mat"
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if("${result}" STREQUAL "0")
    message(FATAL_ERROR "antsAI succeeded with ${ARGN}")
  endif()
  if(NOT "${stdout}${stderr}" MATCHES "${pattern}")
    message(FATAL_ERROR "antsAI failed without saying why (expected '${pattern}') with ${ARGN}:\n${stdout}\n${stderr}")
  endif()
endfunction()

expect_error("dimensionality" -m "Mattes[${image},${image}]" -t "Rigid[0.1]")
expect_error("images not specified" -d 2 -m "Mattes[${image}]" -t "Rigid[0.1]")
expect_error("sampling strategy" -d 2 -m "Mattes[${image},${image},32,Sparse]" -t "Rigid[0.1]")
expect_error("transform" -d 2 -m "Mattes[${image},${image}]" -t "Bogus[0.1]")
expect_error("metric" -d 2 -m "Bogus[${image},${image}]" -t "Rigid[0.1]" -c 1)
