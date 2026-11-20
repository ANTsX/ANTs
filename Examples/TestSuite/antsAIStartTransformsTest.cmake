# antsAI starts its search from a grid of rotations about the center of mass
# and of translations (-s, -g). With a learning rate of 1e-6, the result is the
# best of those starting transforms, so the test checks the starting
# transforms of the rigid, similarity and affine searches. The images are made
# here, so the test needs no data files. Two cases:
# - translation: the fixed image has objects A and B of the same mass; in the
#   moving image only A is shifted, by 6 mm along x. The centers of mass then
#   differ by 3 mm, so the initial translation is 3 mm and only the +3 mm start
#   of the grid -g [3,...] aligns A, which the fixed mask selects. The result
#   must be the 6 mm translation, in 2D and in 3D.
# - rotation (2D): the moving image is the fixed image turned 180 degrees about
#   the image center, and the image origin is far from 0. A start turned 180
#   degrees about the center of mass aligns the images; the result must be
#   that 180 degree rotation.
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")

function(run)
  execute_process(COMMAND "${CMAKE_COMMAND}" -E env ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1 ${ARGN}
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(NOT "${result}" STREQUAL "0")
    message(FATAL_ERROR "Command failed (${result}): ${ARGN}\n${stdout}\n${stderr}")
  endif()
  set(stdout "${stdout}" PARENT_SCOPE)
endfunction()

# Value of a voxel of object A or B: a ramp, so that the objects have texture.
function(object_value x0 y0 x y out)
  math(EXPR value "50 + 10 * (${x} - ${x0}) + 7 * (${y} - ${y0})")
  set(${out} ${value} PARENT_SCOPE)
endfunction()

# Writes a MetaImage of nx x ny (x nz) voxels with object A at ax.., object B at
# bx.., each 6 x 8 voxels (x all z). mask=1 writes the mask of A, grown by 1.
function(write_image file ndims nz ax bx mask)
  set(nx 32)
  set(ny 24)
  set(values "")
  math(EXPR zlast "${nz} - 1")
  foreach(z RANGE ${zlast})
    foreach(y RANGE 23)
      foreach(x RANGE 31)
        set(value 0)
        math(EXPR axl "${ax} + 5")
        math(EXPR bxl "${bx} + 5")
        if(mask)
          math(EXPR axm "${ax} - 1")
          math(EXPR axlm "${axl} + 1")
          if(x GREATER_EQUAL axm AND x LESS_EQUAL axlm AND y GREATER_EQUAL 7 AND y LESS_EQUAL 16)
            set(value 1)
          endif()
        elseif(y GREATER_EQUAL 8 AND y LESS_EQUAL 15)
          if(x GREATER_EQUAL ax AND x LESS_EQUAL axl)
            object_value(${ax} 8 ${x} ${y} value)
          elseif(x GREATER_EQUAL bx AND x LESS_EQUAL bxl)
            object_value(${bx} 8 ${x} ${y} value)
          endif()
        endif()
        string(APPEND values "${value} ")
      endforeach()
    endforeach()
  endforeach()
  if(ndims EQUAL 2)
    set(header "NDims = 2\nElementSpacing = 1 1\nOffset = 40 -25\nDimSize = ${nx} ${ny}\n")
  else()
    set(header "NDims = 3\nElementSpacing = 1 1 1\nOffset = 40 -25 60\nDimSize = ${nx} ${ny} ${nz}\n")
  endif()
  file(WRITE "${file}" "ObjectType = Image\n${header}ElementType = MET_FLOAT\nBinaryData = False\n"
    "ElementDataFile = LOCAL\n${values}\n")
endfunction()

set(number "([-+0-9.eE]+|-?nan|-?inf)")

# Reads the matrix and the translation of the transform in file.
function(read_transform file ndims matrix_out translation_out)
  run("${ANTS_TRANSFORM_INFO}" "${file}")
  if(ndims EQUAL 2)
    set(row_pattern "Matrix \\(2x2\\)[ \n]+${number} ${number} *\n *${number} ${number}")
    set(translation_pattern "Translation: \\[${number}, ${number}\\]")
  else()
    set(row_pattern "Matrix \\(3x3\\)[ \n]+${number} ${number} ${number} *\n *${number} ${number} ${number} *\n *${number} ${number} ${number}")
    set(translation_pattern "Translation: \\[${number}, ${number}, ${number}\\]")
  endif()
  if(NOT stdout MATCHES "${row_pattern}")
    message(FATAL_ERROR "Cannot read the matrix from:\n${stdout}")
  endif()
  set(matrix "")
  foreach(i RANGE 1 ${CMAKE_MATCH_COUNT})
    list(APPEND matrix "${CMAKE_MATCH_${i}}")
  endforeach()
  if(NOT stdout MATCHES "${translation_pattern}")
    message(FATAL_ERROR "Cannot read the translation from:\n${stdout}")
  endif()
  set(translation "")
  foreach(i RANGE 1 ${CMAKE_MATCH_COUNT})
    list(APPEND translation "${CMAKE_MATCH_${i}}")
  endforeach()
  set(${matrix_out} "${matrix}" PARENT_SCOPE)
  set(${translation_out} "${translation}" PARENT_SCOPE)
endfunction()

# Converts a number printed by antsTransformInfo, such as -0.738693 or
# 2.5804e-17, to an integer number of thousandths, rounded toward zero.
function(to_thousandths value out)
  if(value MATCHES "[eE]-")
    set(${out} 0 PARENT_SCOPE) # antsTransformInfo prints at most 6 digits, so e-4 and below are 0 here
    return()
  endif()
  if(NOT value MATCHES "^([-+]?)([0-9]*)(\\.([0-9]*))?$")
    message(FATAL_ERROR "Cannot read the number ${value}")
  endif()
  set(sign "${CMAKE_MATCH_1}")
  set(whole "${CMAKE_MATCH_2}")
  set(fraction "${CMAKE_MATCH_4}000")
  if(whole STREQUAL "")
    set(whole 0)
  endif()
  string(SUBSTRING "${fraction}" 0 3 fraction)
  string(REGEX REPLACE "^0+([0-9])" "\\1" fraction "${fraction}")
  math(EXPR result "${whole} * 1000 + ${fraction}")
  if(sign STREQUAL "-")
    math(EXPR result "-${result}")
  endif()
  set(${out} ${result} PARENT_SCOPE)
endfunction()

# Fails unless the two lists of numbers agree to within 0.01.
function(check_close actual expected what)
  list(LENGTH actual n)
  math(EXPR last "${n} - 1")
  foreach(i RANGE ${last})
    list(GET actual ${i} a)
    list(GET expected ${i} e)
    to_thousandths("${a}" a_thousandths)
    to_thousandths("${e}" e_thousandths)
    math(EXPR difference "${a_thousandths} - ${e_thousandths}")
    if(difference GREATER 10 OR difference LESS -10)
      message(FATAL_ERROR "${what}: expected ${expected}, got ${actual}")
    endif()
  endforeach()
endfunction()

# Checks the transform in file against the expected translation along x and
# matrix diagonal (1 or -1); all other matrix entries and translations must be 0.
function(check_transform file ndims expected_tx diagonal what)
  read_transform("${file}" ${ndims} matrix translation)
  if(ndims EQUAL 2)
    set(expected_matrix ${diagonal} 0 0 ${diagonal})
  else()
    set(expected_matrix ${diagonal} 0 0 0 ${diagonal} 0 0 0 ${diagonal})
  endif()
  if(diagonal EQUAL 1)
    set(expected_translation ${expected_tx} 0)
    if(ndims EQUAL 3)
      list(APPEND expected_translation 0)
    endif()
  else()
    set(expected_translation "") # a 180 degree turn: only the matrix is checked
  endif()
  set(actual ${matrix})
  set(expected ${expected_matrix})
  if(expected_translation)
    list(APPEND actual ${translation})
    list(APPEND expected ${expected_translation})
  endif()
  list(LENGTH actual n)
  math(EXPR last "${n} - 1")
  foreach(i RANGE ${last})
    list(GET actual ${i} value)
    list(GET expected ${i} e)
    if(NOT value MATCHES "^[-+0-9.eE]+$")
      message(FATAL_ERROR "${what}: transform is not finite: matrix ${matrix}, translation ${translation}")
    endif()
    # 0.01 tolerance: the 1e-6 learning rate moves the start by much less.
    math(EXPR low "${e} * 100 - 1")
    math(EXPR high "${e} * 100 + 1")
    if(value LESS "${low}e-2" OR value GREATER "${high}e-2")
      message(FATAL_ERROR "${what}: expected matrix ${expected_matrix} and translation ${expected_translation}, "
        "got matrix ${matrix} and translation ${translation}")
    endif()
  endforeach()
endfunction()

# Translation case, 2D and 3D.
foreach(ndims 2 3)
  if(ndims EQUAL 2)
    set(nz 1)
    set(grid "3x0")
  else()
    set(nz 4)
    set(grid "3x0x0")
  endif()
  set(fixed "${TEST_OUTPUT_DIR}/fixed${ndims}d.mha")
  set(moving "${TEST_OUTPUT_DIR}/moving${ndims}d.mha")
  set(mask "${TEST_OUTPUT_DIR}/mask${ndims}d.mha")
  write_image("${fixed}" ${ndims} ${nz} 4 20 0)
  write_image("${moving}" ${ndims} ${nz} 10 20 0)
  write_image("${mask}" ${ndims} ${nz} 4 20 1)
  foreach(transform Rigid Similarity Affine)
    set(output "${TEST_OUTPUT_DIR}/translation_${ndims}d_${transform}.mat")
    run("${ANTS_AI}" -d ${ndims} -m "Mattes[${fixed},${moving},32,None]" -x "${mask}" -t "${transform}[1e-6]"
      -s "[180,0]" -g "[3,${grid}]" -c 1 -o "${output}")
    check_transform("${output}" ${ndims} 6 1 "${ndims}D ${transform} translation search")
  endforeach()
endforeach()

# Rotation case, 2D: the moving image is the fixed image with its voxels in
# reverse order, which turns it 180 degrees about the image center.
set(fixed "${TEST_OUTPUT_DIR}/fixed2d.mha")
file(STRINGS "${fixed}" fixed_lines)
list(GET fixed_lines -1 fixed_values)
string(REGEX REPLACE " +$" "" fixed_values "${fixed_values}")
string(REPLACE " " ";" voxel_list "${fixed_values}")
list(REVERSE voxel_list)
string(REPLACE ";" " " turned_values "${voxel_list}")
set(turned "${TEST_OUTPUT_DIR}/turned2d.mha")
file(WRITE "${turned}"
  "ObjectType = Image\nNDims = 2\nElementSpacing = 1 1\nOffset = 40 -25\nDimSize = 32 24\n"
  "ElementType = MET_FLOAT\nBinaryData = False\nElementDataFile = LOCAL\n${turned_values}\n")
foreach(transform Rigid Similarity Affine)
  set(output "${TEST_OUTPUT_DIR}/rotation_2d_${transform}.mat")
  run("${ANTS_AI}" -d 2 -m "Mattes[${fixed},${turned},32,None]" -t "${transform}[1e-6]" -s "[180,1]" -c 1
    -o "${output}")
  check_transform("${output}" 2 0 -1 "2D ${transform} rotation search")
endforeach()

# Principal axes case (-p 1), 2D and 3D: the moving image is the fixed image
# with x and y swapped, so aligning the principal axes turns the image. A
# single start (-s [180,0], no translation grid) must be the initial transform
# that -c 0 writes.
foreach(ndims 2 3)
  if(ndims EQUAL 2)
    set(nz 1)
  else()
    set(nz 4)
  endif()
  set(fixed "${TEST_OUTPUT_DIR}/fixed${ndims}d.mha")
  file(STRINGS "${fixed}" fixed_lines)
  list(GET fixed_lines -1 fixed_values)
  string(REGEX REPLACE " +$" "" fixed_values "${fixed_values}")
  string(REPLACE " " ";" voxel_list "${fixed_values}")
  set(swapped_values "")
  math(EXPR zlast "${nz} - 1")
  foreach(z RANGE ${zlast})
    foreach(y RANGE 31)
      foreach(x RANGE 23)
        math(EXPR index "${z} * 768 + ${x} * 32 + ${y}")
        list(GET voxel_list ${index} value)
        string(APPEND swapped_values "${value} ")
      endforeach()
    endforeach()
  endforeach()
  set(swapped "${TEST_OUTPUT_DIR}/swapped${ndims}d.mha")
  if(ndims EQUAL 2)
    set(header "NDims = 2\nElementSpacing = 1 1\nOffset = 40 -25\nDimSize = 24 32\n")
  else()
    set(header "NDims = 3\nElementSpacing = 1 1 1\nOffset = 40 -25 60\nDimSize = 24 32 ${nz}\n")
  endif()
  file(WRITE "${swapped}" "ObjectType = Image\n${header}ElementType = MET_FLOAT\nBinaryData = False\n"
    "ElementDataFile = LOCAL\n${swapped_values}\n")
  foreach(transform Rigid Similarity Affine)
    set(initial "${TEST_OUTPUT_DIR}/axes_${ndims}d_${transform}_initial.mat")
    set(start "${TEST_OUTPUT_DIR}/axes_${ndims}d_${transform}_start.mat")
    run("${ANTS_AI}" -d ${ndims} -m "Mattes[${fixed},${swapped},32,None]" -t "${transform}[1e-6]" -p 1 -c 0
      -o "${initial}")
    run("${ANTS_AI}" -d ${ndims} -m "Mattes[${fixed},${swapped},32,None]" -t "${transform}[1e-6]" -p 1
      -s "[180,0]" -c 1 -o "${start}")
    read_transform("${initial}" ${ndims} initial_matrix initial_translation)
    read_transform("${start}" ${ndims} start_matrix start_translation)
    # The test needs the principal axes to turn the image.
    list(GET initial_matrix 0 m00)
    to_thousandths("${m00}" m00)
    if(m00 GREATER 900 OR m00 LESS -900)
      message(FATAL_ERROR "${ndims}D ${transform}: aligning the principal axes did not turn the image: ${initial_matrix}")
    endif()
    check_close("${start_matrix};${start_translation}" "${initial_matrix};${initial_translation}"
      "${ndims}D ${transform} start from the principal axes")
  endforeach()
endforeach()
