# antsAI -t AlignCentersOfMass and -t AlignGeometricCenters only align the
# centers; there is no search. With -c greater than 0, antsAI wrote an empty
# transform file and still succeeded. The test checks that, for both
# transforms, -c 5 writes the same transform as -c 0. The images are made
# here, so the test needs no data files.
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")

foreach(name fixed moving)
  if(name STREQUAL "fixed")
    set(offset "0 0")
  else()
    set(offset "7 -3")
  endif()
  set(values "")
  foreach(y RANGE 15)
    foreach(x RANGE 15)
      math(EXPR value "(${x} * 7 + ${y} * 3) % 17")
      string(APPEND values "${value} ")
    endforeach()
  endforeach()
  file(WRITE "${TEST_OUTPUT_DIR}/${name}.mha"
    "ObjectType = Image\nNDims = 2\nBinaryData = False\nElementSpacing = 1 1\nOffset = ${offset}\n"
    "DimSize = 16 16\nElementType = MET_FLOAT\nElementDataFile = LOCAL\n${values}\n")
endforeach()

foreach(transform AlignCentersOfMass AlignGeometricCenters)
  foreach(iterations 0 5)
    set(output "${TEST_OUTPUT_DIR}/${transform}_${iterations}.mat")
    file(REMOVE "${output}")
    execute_process(
      COMMAND "${ANTS_AI}" -d 2 -m "Mattes[${TEST_OUTPUT_DIR}/fixed.mha,${TEST_OUTPUT_DIR}/moving.mha,32,None]"
        -t ${transform} -c ${iterations} -o "${output}"
      RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
    if(NOT "${result}" STREQUAL "0")
      message(FATAL_ERROR "antsAI -t ${transform} -c ${iterations} failed (${result}):\n${stdout}\n${stderr}")
    endif()
    file(SIZE "${output}" size)
    if(size EQUAL 0)
      message(FATAL_ERROR "antsAI -t ${transform} -c ${iterations} wrote an empty transform file")
    endif()
    execute_process(COMMAND "${ANTS_TRANSFORM_INFO}" "${output}" OUTPUT_VARIABLE info RESULT_VARIABLE result)
    if(NOT "${result}" STREQUAL "0" OR NOT info MATCHES "Offset: \\[([^]]*)\\]")
      message(FATAL_ERROR "Cannot read the transform written by antsAI -t ${transform} -c ${iterations}:\n${info}")
    endif()
    set(offset_${iterations} "${CMAKE_MATCH_1}")
  endforeach()
  if(NOT offset_0 STREQUAL offset_5)
    message(FATAL_ERROR "antsAI -t ${transform}: -c 0 gives offset ${offset_0}, -c 5 gives ${offset_5}")
  endif()
  if(NOT offset_0 MATCHES "^7, -3$")
    message(FATAL_ERROR "antsAI -t ${transform}: the offset is ${offset_0}, not 7, -3")
  endif()
endforeach()
