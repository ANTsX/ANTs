# Exercise antsMotionCorr -w 1 on a small 2D+t series without external image
# fixtures. The Warp and InverseWarp fields must be written with float
# components. A fixed image that is smaller than the moving frames must give
# an error with -u 0, not an out-of-bounds copy.
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")

# Write a MetaImage of a Gaussian-like blob, one frame per center in ARGN.
function(write_blob_image file width)
  set(values "")
  foreach(center IN LISTS ARGN)
    string(REPLACE ":" ";" center "${center}")
    list(GET center 0 cx)
    list(GET center 1 cy)
    math(EXPR last "${width} - 1")
    foreach(y RANGE 0 ${last})
      foreach(x RANGE 0 ${last})
        math(EXPR value "400 - ((${x} - ${cx}) * (${x} - ${cx}) + (${y} - ${cy}) * (${y} - ${cy})) * 8")
        if(value LESS 0)
          set(value 0)
        endif()
        string(APPEND values "${value} ")
      endforeach()
    endforeach()
  endforeach()
  list(LENGTH ARGN frames)
  if(frames EQUAL 1)
    set(header "NDims = 2\nElementSpacing = 1 1\nDimSize = ${width} ${width}")
  else()
    set(header "NDims = 3\nElementSpacing = 1 1 1\nDimSize = ${width} ${width} ${frames}")
  endif()
  file(WRITE "${file}"
    "ObjectType = Image\n${header}\nBinaryData = False\nElementType = MET_FLOAT\nElementDataFile = LOCAL\n${values}\n")
endfunction()

function(run_moco prefix fixed use_fixed result_var output_var)
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
      "${ANTS_MOTION_CORR}" -d 2
      -o "[${TEST_OUTPUT_DIR}/${prefix},${TEST_OUTPUT_DIR}/${prefix}.nii.gz]"
      -m "Demons[${TEST_OUTPUT_DIR}/${fixed},${TEST_OUTPUT_DIR}/series.mha,1,1]"
      -t "Rigid[0.1]" -i 5 -u ${use_fixed} -e 1 -s 0 -f 1 -w 1 --random-seed 1
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  set(${result_var} "${result}" PARENT_SCOPE)
  set(${output_var} "${stdout}\n${stderr}" PARENT_SCOPE)
endfunction()

write_blob_image("${TEST_OUTPUT_DIR}/fixed.mha" 16 "8:8")
write_blob_image("${TEST_OUTPUT_DIR}/fixed_small.mha" 8 "4:4")
write_blob_image("${TEST_OUTPUT_DIR}/series.mha" 16 "8:8" "9:8" "8:9")

run_moco(moco fixed.mha 1 result output)
if(NOT "${result}" STREQUAL "0")
  message(FATAL_ERROR "antsMotionCorr -w 1 failed: ${result}\n${output}")
endif()
foreach(field Warp InverseWarp)
  set(file "${TEST_OUTPUT_DIR}/moco${field}.nii.gz")
  if(NOT EXISTS "${file}")
    message(FATAL_ERROR "antsMotionCorr -w 1 did not write ${file}\n${output}")
  endif()
  execute_process(COMMAND "${PRINT_HEADER}" "${file}"
    RESULT_VARIABLE result OUTPUT_VARIABLE header ERROR_VARIABLE stderr)
  if(NOT "${result}" STREQUAL "0")
    message(FATAL_ERROR "PrintHeader ${file} failed: ${result}\n${stderr}")
  endif()
  # NIfTI datatype 16 is float32 (float64 is 64).
  if(NOT header MATCHES "datatype = 16\n" OR NOT header MATCHES "bitpix = 32\n")
    message(FATAL_ERROR "${field} does not have float components:\n${header}")
  endif()
endforeach()

run_moco(mismatch fixed_small.mha 0 result output)
if(NOT "${result}" STREQUAL "1" OR NOT output MATCHES "are not inside the output image")
  message(FATAL_ERROR "antsMotionCorr -u 0 with a small fixed image must fail with an error: ${result}\n${output}")
endif()
