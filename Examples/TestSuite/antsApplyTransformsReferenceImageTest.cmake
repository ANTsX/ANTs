# Exercise the reference image (-r) handling of antsApplyTransforms without
# external image fixtures. The tool reads only the geometry of the reference
# image. A reference image that cannot be read is an error when an image is
# written, but not when only a linear or composite transform is written.
file(REMOVE_RECURSE "${TEST_OUTPUT_DIR}")
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")

# Set the environment here and not with "cmake -E env", because "cmake -E env"
# reports a crash of the child process as exit status 1
set(ENV{ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS} 2)

function(apply_transforms expected_result)
  execute_process(
    COMMAND "${ANTS_APPLY_TRANSFORMS}" ${ARGN}
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(NOT "${result}" STREQUAL "${expected_result}")
    message(FATAL_ERROR "Expected exit status ${expected_result}, got ${result}: ${ARGN}\n${stdout}\n${stderr}")
  endif()
endfunction()

function(check_exists file)
  if(NOT EXISTS "${TEST_OUTPUT_DIR}/${file}")
    message(FATAL_ERROR "Expected output ${file} was not written")
  endif()
endfunction()

function(check_not_exists file)
  if(EXISTS "${TEST_OUTPUT_DIR}/${file}")
    message(FATAL_ERROR "Unexpected output ${file} was written")
  endif()
endfunction()

# Scalar 3D input and a reference image with a different geometry
set(values "")
foreach(i RANGE 0 63)
  math(EXPR value "(${i} * 17) % 97")
  string(APPEND values "${value} ")
endforeach()
file(WRITE "${TEST_OUTPUT_DIR}/input.mha"
  "ObjectType = Image\nNDims = 3\nBinaryData = False\nElementSpacing = 1 1 1\nDimSize = 4 4 4\nElementType = MET_FLOAT\nElementDataFile = LOCAL\n${values}\n")
set(values "")
foreach(i RANGE 0 26)
  string(APPEND values "0 ")
endforeach()
file(WRITE "${TEST_OUTPUT_DIR}/reference.mha"
  "ObjectType = Image\nNDims = 3\nBinaryData = False\nOffset = 0.5 0.5 0.5\nElementSpacing = 1.5 1.5 1.5\nDimSize = 3 3 3\nElementType = MET_FLOAT\nElementDataFile = LOCAL\n${values}\n")
file(WRITE "${TEST_OUTPUT_DIR}/not_an_image.txt" "This is not an image\n")
file(WRITE "${TEST_OUTPUT_DIR}/translation.txt"
  "#Insight Transform File V1.0\n#Transform 0\nTransform: AffineTransform_double_3_3\nParameters: 1 0 0 0 1 0 0 0 1 0.5 0 0\nFixedParameters: 0 0 0\n")

# A readable reference image gives an output image with the reference geometry
apply_transforms(0 -d 3 -i "${TEST_OUTPUT_DIR}/input.mha" -r "${TEST_OUTPUT_DIR}/reference.mha"
  -t "${TEST_OUTPUT_DIR}/translation.txt" -o "${TEST_OUTPUT_DIR}/scalar.mha")
check_exists(scalar.mha)
file(STRINGS "${TEST_OUTPUT_DIR}/scalar.mha" size_line REGEX "^DimSize = ")
if(NOT size_line STREQUAL "DimSize = 3 3 3")
  message(FATAL_ERROR "Output size does not match the reference image: ${size_line}")
endif()

# A missing or unreadable reference image is an error when an image is written
apply_transforms(1 -d 3 -i "${TEST_OUTPUT_DIR}/input.mha" -r "${TEST_OUTPUT_DIR}/missing.mha"
  -t "${TEST_OUTPUT_DIR}/translation.txt" -o "${TEST_OUTPUT_DIR}/missing_out.mha")
check_not_exists(missing_out.mha)
apply_transforms(1 -d 3 -i "${TEST_OUTPUT_DIR}/input.mha" -r "${TEST_OUTPUT_DIR}/not_an_image.txt"
  -t "${TEST_OUTPUT_DIR}/translation.txt" -o "${TEST_OUTPUT_DIR}/bad_out.mha")
check_not_exists(bad_out.mha)
apply_transforms(1 -d 3 -r "${TEST_OUTPUT_DIR}/missing.mha"
  -t "${TEST_OUTPUT_DIR}/translation.txt" -o "[${TEST_OUTPUT_DIR}/missing_field.nii.gz,1]")
check_not_exists(missing_field.nii.gz)

# The reference image is not necessary when only a transform is written
apply_transforms(0 -d 3 -r "${TEST_OUTPUT_DIR}/missing.mha"
  -t "${TEST_OUTPUT_DIR}/translation.txt" -o "Linear[${TEST_OUTPUT_DIR}/linear.mat]")
check_exists(linear.mat)
apply_transforms(0 -d 3 -r "${TEST_OUTPUT_DIR}/missing.mha"
  -t "${TEST_OUTPUT_DIR}/translation.txt" -o "CompositeTransform[${TEST_OUTPUT_DIR}/composite.h5]")
check_exists(composite.h5)

# Tensor input (-e 2) uses the reference geometry for the output tensor image
set(values "")
foreach(i RANGE 0 63)
  string(APPEND values "1 0 0 1 0 1 ")
endforeach()
file(WRITE "${TEST_OUTPUT_DIR}/tensor.mha"
  "ObjectType = Image\nNDims = 3\nBinaryData = False\nElementSpacing = 1 1 1\nDimSize = 4 4 4\nElementNumberOfChannels = 6\nElementType = MET_FLOAT\nElementDataFile = LOCAL\n${values}\n")
apply_transforms(0 -d 3 -e 2 -i "${TEST_OUTPUT_DIR}/tensor.mha" -r "${TEST_OUTPUT_DIR}/reference.mha"
  -t "${TEST_OUTPUT_DIR}/translation.txt" -o "${TEST_OUTPUT_DIR}/tensor_out.nii.gz")
check_exists(tensor_out.nii.gz)
apply_transforms(1 -d 3 -e 2 -i "${TEST_OUTPUT_DIR}/tensor.mha" -r "${TEST_OUTPUT_DIR}/missing.mha"
  -t "${TEST_OUTPUT_DIR}/translation.txt" -o "${TEST_OUTPUT_DIR}/tensor_missing.nii.gz")
check_not_exists(tensor_missing.nii.gz)
