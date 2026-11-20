# antsAI with Random metric sampling must give the same transform in two runs
# with the same --random-seed. The random image iterator that picks the sample
# points has its own generator, which was not seeded, so the sample points and
# the transform changed from run to run. The images are made here, so the test
# needs no data files.
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")

foreach(name fixed moving)
  if(name STREQUAL "fixed")
    set(shift 0)
  else()
    set(shift 2)
  endif()
  set(values "")
  foreach(y RANGE 31)
    foreach(x RANGE 31)
      math(EXPR value "((${x} + ${shift}) * 7 + ${y} * 3) % 23 + ((${x} + ${shift}) / 8) * 10")
      string(APPEND values "${value} ")
    endforeach()
  endforeach()
  file(WRITE "${TEST_OUTPUT_DIR}/${name}.mha"
    "ObjectType = Image\nNDims = 2\nBinaryData = False\nElementSpacing = 1 1\n"
    "DimSize = 32 32\nElementType = MET_FLOAT\nElementDataFile = LOCAL\n${values}\n")
endforeach()

foreach(run 1 2)
  set(output "${TEST_OUTPUT_DIR}/run${run}.mat")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env --unset=ANTS_RANDOM_SEED
      "${ANTS_AI}" -d 2 -m "Mattes[${TEST_OUTPUT_DIR}/fixed.mha,${TEST_OUTPUT_DIR}/moving.mha,16,Random,0.3]"
      -t "Affine[0.1]" -s "[90,0.25]" -c 5 --random-seed 3 -o "${output}"
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(NOT "${result}" STREQUAL "0")
    message(FATAL_ERROR "antsAI failed (${result}):\n${stdout}\n${stderr}")
  endif()
  execute_process(COMMAND "${ANTS_TRANSFORM_INFO}" "${output}" OUTPUT_VARIABLE info RESULT_VARIABLE result)
  if(NOT "${result}" STREQUAL "0" OR NOT info MATCHES "Matrix \\(2x2\\)[^\n]*\n([^\n]*\n[^\n]*)\n  Offset: ([^\n]*)")
    message(FATAL_ERROR "Cannot read the transform from:\n${info}")
  endif()
  set(transform_${run} "${CMAKE_MATCH_1} ${CMAKE_MATCH_2}")
endforeach()

if(NOT transform_1 STREQUAL transform_2)
  message(FATAL_ERROR "Two antsAI runs with --random-seed 3 give different transforms:\n${transform_1}\n${transform_2}")
endif()
