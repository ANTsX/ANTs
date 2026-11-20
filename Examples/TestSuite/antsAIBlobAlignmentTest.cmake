# antsAI -b aligns two images by matching blobs. The test checks two cases:
# - an image aligned to itself gives the identity, for a rigid and an affine
#   transform. Blob pairs whose distances to the other blobs are the same in
#   both images become landmarks; the distances were measured wrongly, so no
#   landmark was kept and the translation was NaN.
# - a blank moving image has no blobs to match, so antsAI must fail and write
#   no transform.
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")

function(run)
  execute_process(COMMAND ${ARGN} RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(NOT "${result}" STREQUAL "0")
    message(FATAL_ERROR "Command failed (${result}): ${ARGN}\n${stdout}\n${stderr}")
  endif()
  set(stdout "${stdout}" PARENT_SCOPE)
endfunction()

set(number "([-+0-9.eE]+|-?nan|-?inf)")

foreach(transform Rigid Affine)
  set(output "${TEST_OUTPUT_DIR}/self_${transform}.mat")
  run("${ANTS_AI}" -d 2 -m "Mattes[${INPUT_IMAGE},${INPUT_IMAGE},32,None]" -t "${transform}[0.1]" -b 20 -c 0
    -o "${output}")
  run("${ANTS_TRANSFORM_INFO}" "${output}")
  if(NOT stdout MATCHES "Translation: \\[${number}, ${number}\\]")
    message(FATAL_ERROR "Cannot read the translation from:\n${stdout}")
  endif()
  set(translation "${CMAKE_MATCH_1}" "${CMAKE_MATCH_2}")
  if(NOT stdout MATCHES "Matrix \\(2x2\\)[ \n]+${number} ${number} *\n *${number} ${number}")
    message(FATAL_ERROR "Cannot read the matrix from:\n${stdout}")
  endif()
  set(matrix "${CMAKE_MATCH_1}" "${CMAKE_MATCH_2}" "${CMAKE_MATCH_3}" "${CMAKE_MATCH_4}")
  # The identity: translation 0 0, matrix 1 0 0 1.
  set(actual ${translation} ${matrix})
  set(expected 0 0 1 0 0 1)
  foreach(i RANGE 5)
    list(GET actual ${i} value)
    list(GET expected ${i} e)
    if(NOT value MATCHES "^[-+0-9.eE]+$")
      message(FATAL_ERROR "${transform} alignment of an image to itself is not finite: "
        "translation ${translation}, matrix ${matrix}")
    endif()
    math(EXPR low "${e} * 1000000 - 1")
    math(EXPR high "${e} * 1000000 + 1")
    if(value GREATER "${high}e-6" OR value LESS "${low}e-6")
      message(FATAL_ERROR "${transform} alignment of an image to itself is not the identity: "
        "translation ${translation}, matrix ${matrix}")
    endif()
  endforeach()
endforeach()

set(blank "${TEST_OUTPUT_DIR}/blank.nii.gz")
run("${IMAGE_MATH}" 2 "${blank}" m "${INPUT_IMAGE}" 0)
set(output "${TEST_OUTPUT_DIR}/blank.mat")
file(REMOVE "${output}")
execute_process(
  COMMAND "${ANTS_AI}" -d 2 -m "Mattes[${INPUT_IMAGE},${blank},32,None]" -t "Rigid[0.1]" -b 20 -c 0 -o "${output}"
  RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
if("${result}" STREQUAL "0")
  message(FATAL_ERROR "antsAI -b succeeded with a blank moving image:\n${stdout}\n${stderr}")
endif()
if(EXISTS "${output}")
  message(FATAL_ERROR "antsAI -b wrote a transform for a blank moving image")
endif()
