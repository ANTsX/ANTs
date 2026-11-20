# ImageMath BlobDetector must run with several threads. The blob detector keeps
# one heap per work unit, and the pool threader runs more work units than
# threads, so a heap per thread was indexed past its end. The test uses a 2D
# image of three bright squares, made here, so it needs no data files. It runs
# the detector with 1, 2 and 4 threads and checks that all three runs succeed
# and write the same blob image.
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")

set(input "${TEST_OUTPUT_DIR}/squares.mha")
set(values "")
foreach(y RANGE 0 31)
  foreach(x RANGE 0 31)
    set(value 0)
    if(x GREATER_EQUAL 6 AND x LESS_EQUAL 10 AND y GREATER_EQUAL 7 AND y LESS_EQUAL 11)
      set(value 100)
    elseif(x GREATER_EQUAL 19 AND x LESS_EQUAL 25 AND y GREATER_EQUAL 8 AND y LESS_EQUAL 14)
      set(value 80)
    elseif(x GREATER_EQUAL 13 AND x LESS_EQUAL 17 AND y GREATER_EQUAL 21 AND y LESS_EQUAL 27)
      set(value 120)
    endif()
    string(APPEND values "${value} ")
  endforeach()
endforeach()
file(WRITE "${input}"
  "ObjectType = Image\nNDims = 2\nBinaryData = False\nElementSpacing = 1 1\n"
  "DimSize = 32 32\nElementType = MET_FLOAT\nElementDataFile = LOCAL\n${values}\n")

foreach(threads 1 2 4)
  set(output "${TEST_OUTPUT_DIR}/blobs_${threads}.mha")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${threads}
      "${IMAGE_MATH}" 2 "${output}" BlobDetector "${input}" 3
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(NOT "${result}" STREQUAL "0")
    message(FATAL_ERROR "ImageMath BlobDetector failed with ${threads} threads (${result}):\n${stdout}\n${stderr}")
  endif()
  file(SHA256 "${output}" hash_${threads})
endforeach()

if(NOT hash_1 STREQUAL hash_2 OR NOT hash_1 STREQUAL hash_4)
  message(FATAL_ERROR "ImageMath BlobDetector writes different blob images with 1, 2 and 4 threads")
endif()
