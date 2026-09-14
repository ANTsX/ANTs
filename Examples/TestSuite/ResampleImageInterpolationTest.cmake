# Exercise the positional CLI without external image fixtures. In particular,
# a pixel type following a numeric interpolator must not become its parameter.
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")

function(resample dimension output)
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=2
      "${RESAMPLE_IMAGE}" "${dimension}" "${TEST_OUTPUT_DIR}/input${dimension}.mha"
      "${TEST_OUTPUT_DIR}/${output}.mha" "${output_size_${dimension}}" 1 ${ARGN}
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(NOT "${result}" STREQUAL "0")
    message(FATAL_ERROR "ResampleImage ${dimension} ${ARGN} failed: ${result}\n${stdout}\n${stderr}")
  endif()
endfunction()

function(compare_outputs first second)
  file(SHA256 "${TEST_OUTPUT_DIR}/${first}.mha" first_hash)
  file(SHA256 "${TEST_OUTPUT_DIR}/${second}.mha" second_hash)
  if(NOT first_hash STREQUAL second_hash)
    message(FATAL_ERROR "Different resampling results: ${first} and ${second}")
  endif()
endfunction()

set(interpolators Linear NearestNeighbor Gaussian WindowedSinc BSpline)
set(pixel_types MET_CHAR MET_UCHAR MET_SHORT MET_USHORT MET_INT MET_UINT MET_FLOAT MET_DOUBLE)
foreach(dimension RANGE 2 4)
  set(size "")
  set(spacing "")
  set(output_size_${dimension} "")
  set(count 1)
  foreach(axis RANGE 1 ${dimension})
    string(APPEND size "6 ")
    string(APPEND spacing "${axis} ")
    list(APPEND output_size_${dimension} 4)
    math(EXPR count "${count} * 6")
  endforeach()
  list(JOIN output_size_${dimension} x output_size_${dimension})
  set(values "")
  math(EXPR last "${count} - 1")
  foreach(i RANGE 0 ${last})
    math(EXPR value "(${i} * 17) % 97")
    string(APPEND values "${value}.25 ")
  endforeach()
  file(WRITE "${TEST_OUTPUT_DIR}/input${dimension}.mha"
    "ObjectType = Image\nNDims = ${dimension}\nBinaryData = False\nElementSpacing = ${spacing}\nDimSize = ${size}\nElementType = MET_FLOAT\nElementDataFile = LOCAL\n${values}\n")

  foreach(index RANGE 0 4)
    list(GET interpolators ${index} name)
    foreach(pixel_type RANGE 0 7)
      set(prefix "${dimension}-${index}-${pixel_type}")
      resample(${dimension} "${prefix}-numeric" ${index} ${pixel_type})
      resample(${dimension} "${prefix}-named" "${name}" ${pixel_type})
      compare_outputs("${prefix}-numeric" "${prefix}-named")
      list(GET pixel_types ${pixel_type} expected_type)
      file(STRINGS "${TEST_OUTPUT_DIR}/${prefix}-numeric.mha" type_line REGEX "^ElementType = ")
      if(NOT type_line STREQUAL "ElementType = ${expected_type}")
        message(FATAL_ERROR "${prefix}: expected ${expected_type}, got ${type_line}")
      endif()
    endforeach()
    resample(${dimension} "${dimension}-${index}-default-type" ${index})
    compare_outputs("${dimension}-${index}-default-type" "${dimension}-${index}-6-numeric")
  endforeach()
  resample(${dimension} "${dimension}-defaults")
  compare_outputs("${dimension}-defaults" "${dimension}-0-6-numeric")
endforeach()

# Named parameter lists must reach the interpolator, including vector sigma,
# optional Gaussian alpha, case-insensitive names, and window aliases.
resample(2 gaussian-default "Gaussian[spacing,1]")
compare_outputs(gaussian-default 2-2-6-numeric)
resample(2 gaussian-vector "Gaussian[1x2,1]")
compare_outputs(gaussian-vector gaussian-default)
resample(2 gaussian-scalar "gAuSsIaN[0.8,2]")
resample(2 gaussian-repeated "Gaussian[0.8x0.8,2]")
compare_outputs(gaussian-scalar gaussian-repeated)
resample(2 bspline-default "BSpline[3]")
compare_outputs(bspline-default 2-4-6-numeric)
resample(2 bspline-linear "BSpline[1]")
compare_outputs(bspline-linear 2-0-6-numeric)
foreach(window cosine welch blackman lanczos hamming)
  string(SUBSTRING "${window}" 0 1 alias)
  resample(2 "sinc-${window}" "WindowedSinc[${window}]")
  resample(2 "sinc-${alias}" "WindowedSinc[${alias}]")
  compare_outputs("sinc-${window}" "sinc-${alias}")
endforeach()
compare_outputs(sinc-hamming 2-3-6-numeric)

foreach(invalid "2[1]" "Linear[1]" "Gaussian[bad]" "Gaussian[1,bad]"
    "Gaussian[1x2x3]" "Gaussian[1,2,3]" "BSpline[bad]" "BSpline[-1]"
    "BSpline[6]" "WindowedSinc[bad]")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=2
      "${RESAMPLE_IMAGE}" 2 "${TEST_OUTPUT_DIR}/input2.mha"
      "${TEST_OUTPUT_DIR}/invalid.mha" 4x4 1 "${invalid}"
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(NOT "${result}" STREQUAL "1" OR NOT "${stderr}" MATCHES "interpolation")
    message(FATAL_ERROR "Expected a reported interpolation error for ${invalid}: ${result}\n${stdout}\n${stderr}")
  endif()
endforeach()
