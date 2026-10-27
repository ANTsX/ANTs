# Resample a 3-D time series with -e 3, the same data as a multichannel image with -e 4 and as a
# five-dimensional image with -e 5. Each result must be equal, voxel for voxel, to the volumes of the
# series resampled one at a time with -e 0. The inputs are written here, so no external data is needed.
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")

set(nx 6)
set(ny 5)
set(nz 4)
set(nt 3)
math(EXPR voxels_per_volume "${nx} * ${ny} * ${nz}")
math(EXPR last_voxel "${voxels_per_volume} - 1")
math(EXPR last_volume "${nt} - 1")

function(write_meta_image file ndims dimsize values)
  file(WRITE "${file}"
    "ObjectType = Image\nNDims = ${ndims}\nBinaryData = False\nDimSize = ${dimsize}\n"
    "ElementType = MET_FLOAT\nElementDataFile = LOCAL\n${values}\n")
endfunction()

# Volume t of the series has the value (17 * i + 29 * t) % 97 + 0.25 at voxel i
set(series_values "")
foreach(t RANGE 0 ${last_volume})
  set(volume_values "")
  foreach(i RANGE 0 ${last_voxel})
    math(EXPR value "(17 * ${i} + 29 * ${t}) % 97")
    string(APPEND volume_values "${value}.25 ")
  endforeach()
  string(APPEND series_values "${volume_values}")
  write_meta_image("${TEST_OUTPUT_DIR}/volume${t}.mha" 3 "${nx} ${ny} ${nz}" "${volume_values}")
endforeach()
set(multichannel_values "")
foreach(i RANGE 0 ${last_voxel})
  foreach(t RANGE 0 ${last_volume})
    math(EXPR value "(17 * ${i} + 29 * ${t}) % 97")
    string(APPEND multichannel_values "${value}.25 ")
  endforeach()
endforeach()
write_meta_image("${TEST_OUTPUT_DIR}/series.mha" 4 "${nx} ${ny} ${nz} ${nt}" "${series_values}")
write_meta_image("${TEST_OUTPUT_DIR}/fivedimensional.mha" 5 "${nx} ${ny} ${nz} 1 ${nt}" "${series_values}")
# ITK does not read a multichannel ASCII MetaImage correctly, so use an ASCII NRRD file
file(WRITE "${TEST_OUTPUT_DIR}/multichannel.nrrd"
  "NRRD0004\ntype: float\ndimension: 4\nsizes: ${nt} ${nx} ${ny} ${nz}\nkinds: vector domain domain domain\n"
  "encoding: ascii\n\n${multichannel_values}\n")

# The reference grid differs from the input grid in size, spacing and origin
file(WRITE "${TEST_OUTPUT_DIR}/reference.mha"
  "ObjectType = Image\nNDims = 3\nBinaryData = False\nElementSpacing = 1.3 0.8 1.1\nOffset = 0.5 -0.3 0.2\n"
  "DimSize = 4 7 3\nElementType = MET_FLOAT\nElementDataFile = LOCAL\n")
foreach(i RANGE 1 84)
  file(APPEND "${TEST_OUTPUT_DIR}/reference.mha" "0 ")
endforeach()
file(WRITE "${TEST_OUTPUT_DIR}/affine.txt"
  "#Insight Transform File V1.0\n#Transform 0\nTransform: AffineTransform_double_3_3\n"
  "Parameters: 0.9 0.1 0 -0.05 1.1 0.05 0 0.02 0.95 0.3 -0.2 0.1\nFixedParameters: 2.5 2 1.5\n")

function(apply_transforms input_type input output)
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=2
      "${ANTS_APPLY_TRANSFORMS}" -d 3 -e ${input_type} -i "${TEST_OUTPUT_DIR}/${input}"
      -r "${TEST_OUTPUT_DIR}/reference.mha" -t "${TEST_OUTPUT_DIR}/affine.txt"
      -o "${TEST_OUTPUT_DIR}/${output}" ${ARGN}
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(NOT "${result}" STREQUAL "0")
    message(FATAL_ERROR "antsApplyTransforms -e ${input_type} ${ARGN} failed: ${result}\n${stdout}\n${stderr}")
  endif()
endfunction()

# The voxel data of an uncompressed NIfTI-1 file written by ITK starts at byte 352. ITK writes the
# channels of a vector image one after the other, so all three outputs have the layout of the series.
function(read_voxel_data file variable)
  file(READ "${TEST_OUTPUT_DIR}/${file}" data OFFSET 352 HEX)
  set(${variable} "${data}" PARENT_SCOPE)
endfunction()

set(option_sets "-n Linear -u float" "-n Gaussian -u float" "-n Linear -u int")
set(option_set_index 0)
foreach(option_set IN LISTS option_sets)
  separate_arguments(options UNIX_COMMAND "${option_set}")
  set(prefix "case${option_set_index}")
  file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}/${prefix}")

  set(expected "")
  foreach(t RANGE 0 ${last_volume})
    apply_transforms(0 "volume${t}.mha" "${prefix}/volume${t}.nii" ${options})
    read_voxel_data("${prefix}/volume${t}.nii" volume_data)
    string(APPEND expected "${volume_data}")
  endforeach()
  string(LENGTH "${expected}" expected_length)
  if(expected_length EQUAL 0)
    message(FATAL_ERROR "${option_set}: the per-volume outputs contain no voxel data")
  endif()

  foreach(case "3;series.mha" "4;multichannel.nrrd" "5;fivedimensional.mha")
    list(GET case 0 input_type)
    list(GET case 1 input)
    apply_transforms(${input_type} "${input}" "${prefix}/e${input_type}.nii" ${options})
    read_voxel_data("${prefix}/e${input_type}.nii" actual)
    if(NOT actual STREQUAL expected)
      message(FATAL_ERROR "${option_set}: -e ${input_type} output is not equal to the per-volume -e 0 outputs")
    endif()
  endforeach()

  math(EXPR option_set_index "${option_set_index} + 1")
endforeach()

# The -e 3 output has the reference grid in space and the input number of time points
file(READ "${TEST_OUTPUT_DIR}/case0/e3.nii" dims OFFSET 40 LIMIT 10 HEX)
if(NOT dims STREQUAL "04000400070003000300")
  message(FATAL_ERROR "-e 3 output has dim[0..4] ${dims}, expected 4, 4, 7, 3, ${nt}")
endif()
