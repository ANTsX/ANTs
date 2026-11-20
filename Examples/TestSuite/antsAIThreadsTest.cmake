# antsAI optimizes each start on one thread, with the starts spread over the
# threads, so its transform must not depend on the number of threads. The
# multi-start optimizer it replaced split each metric evaluation over the
# threads, and the transforms differed in the last digits. The test aligns two
# slices with 1, 2 and 4 threads and checks that the transform files, written
# with full precision, are the same.
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")

foreach(threads 1 2 4)
  set(output "${TEST_OUTPUT_DIR}/transform_${threads}.txt")
  file(REMOVE "${output}")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${threads}
      "${ANTS_AI}" -d 2 -m "Mattes[${FIXED_IMAGE},${MOVING_IMAGE},32,Regular,0.25]" -t "Affine[0.1]" -s "[30,0.5]"
      -c 5 --random-seed 1 -o "${output}"
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(NOT "${result}" STREQUAL "0" OR NOT EXISTS "${output}")
    message(FATAL_ERROR "antsAI failed with ${threads} threads (${result}):\n${stdout}\n${stderr}")
  endif()
  file(READ "${output}" transform_${threads})
endforeach()

if(NOT transform_1 STREQUAL transform_2 OR NOT transform_1 STREQUAL transform_4)
  message(FATAL_ERROR "antsAI gives different transforms with 1, 2 and 4 threads:\n"
    "${transform_1}\n${transform_2}\n${transform_4}")
endif()
