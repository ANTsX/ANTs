# antsAI must not optimize the same starting transform twice. The rotation grid
# repeats rotations: -180 and +180 degrees are the same turn, and in 3D the
# turns about axis 1, axis 2 and axis 1 again repeat rotations. With a search
# factor of 90 degrees over the full circle, 2D has 4 different rotations of
# the 5 angles, and 3D has the 24 rotations of a cube of the 125 angle triples.
# The test checks the numbers of starting points that antsAI -v reports. The
# images are made here, so the test needs no data files.
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")

set(values "")
foreach(y RANGE 15)
  foreach(x RANGE 15)
    math(EXPR value "(${x} * 7 + ${y} * 3) % 17 + ${x}")
    string(APPEND values "${value} ")
  endforeach()
endforeach()
set(image2d "${TEST_OUTPUT_DIR}/image2d.mha")
file(WRITE "${image2d}"
  "ObjectType = Image\nNDims = 2\nBinaryData = False\nElementSpacing = 1 1\n"
  "DimSize = 16 16\nElementType = MET_FLOAT\nElementDataFile = LOCAL\n${values}\n")
set(values "")
foreach(z RANGE 7)
  foreach(y RANGE 7)
    foreach(x RANGE 7)
      math(EXPR value "(${x} * 7 + ${y} * 3 + ${z} * 5) % 17 + ${x} + 2 * ${y}")
      string(APPEND values "${value} ")
    endforeach()
  endforeach()
endforeach()
set(image3d "${TEST_OUTPUT_DIR}/image3d.mha")
file(WRITE "${image3d}"
  "ObjectType = Image\nNDims = 3\nBinaryData = False\nElementSpacing = 1 1 1\n"
  "DimSize = 8 8 8\nElementType = MET_FLOAT\nElementDataFile = LOCAL\n${values}\n")

# Runs antsAI -v and checks the reported numbers of starting points.
function(check_starts ndims transform expected_unique expected_total)
  if(ndims EQUAL 2)
    set(image "${image2d}")
  else()
    set(image "${image3d}")
  endif()
  execute_process(
    COMMAND "${ANTS_AI}" -d ${ndims} -v 1 -m "Mattes[${image},${image},8,None]" -t "${transform}[0.1]" -s 90 -c 1
      -o "${TEST_OUTPUT_DIR}/out.mat" ${ARGN}
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(NOT "${result}" STREQUAL "0")
    message(FATAL_ERROR "antsAI failed (${result}):\n${stdout}\n${stderr}")
  endif()
  if(NOT stdout MATCHES "Starting optimizer with ([0-9]+) starting points \\(([0-9]+) before removing repeats\\)")
    message(FATAL_ERROR "Cannot read the number of starting points from:\n${stdout}")
  endif()
  if(NOT CMAKE_MATCH_1 EQUAL expected_unique OR NOT CMAKE_MATCH_2 EQUAL expected_total)
    message(FATAL_ERROR "${ndims}D ${transform} ${ARGN}: expected ${expected_unique} of ${expected_total} starting points, "
      "got ${CMAKE_MATCH_1} of ${CMAKE_MATCH_2}")
  endif()
endfunction()

foreach(transform Rigid Similarity Affine)
  check_starts(2 ${transform} 4 5)
  check_starts(2 ${transform} 12 15 -g "[5,5x0]")
  check_starts(3 ${transform} 24 125)
endforeach()
