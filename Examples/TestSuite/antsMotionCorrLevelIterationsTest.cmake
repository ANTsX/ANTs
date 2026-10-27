# antsMotionCorr must use the number of iterations given for each level by -i.
# The optimizer can converge before it reaches that number, so the test checks
# upper bounds. Before the fix, every level used the first value of -i.
# The 2D+t input is stacked at test time from the r16 and r64 slices.
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")

function(run_tool)
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1 ${ARGN}
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(NOT "${result}" STREQUAL "0")
    message(FATAL_ERROR "Command failed (${result}): ${ARGN}\n${stdout}\n${stderr}")
  endif()
  set(stdout "${stdout}" PARENT_SCOPE)
endfunction()

set(series "${TEST_OUTPUT_DIR}/series.nii.gz")
run_tool("${STACK_SLICES}" "${series}" -1 -1 0 "${R16_IMAGE}" "${R64_IMAGE}")

# One registration for each of the two time points, with three levels each.
set(maximum_iterations 20 5 0)
foreach(transform Rigid Affine GaussianDisplacementField)
  if(transform STREQUAL "GaussianDisplacementField")
    set(metric "Demons[${R16_IMAGE},${series},1,1]")
    set(parameters "0.15,3,0.5")
  else()
    set(metric "MI[${R16_IMAGE},${series},1,32,Regular,0.5]")
    set(parameters "0.1")
  endif()
  run_tool("${ANTS_MOTION_CORR}" -d 2 -v 1 -m "${metric}"
    -t "${transform}[${parameters}]" -i 20x5x0 -f 4x2x1 -s 2x1x0 -u 1 -e 1
    --random-seed 42
    -o "[${TEST_OUTPUT_DIR}/${transform},${TEST_OUTPUT_DIR}/${transform}Warped.nii.gz]")

  string(REGEX MATCHALL "Level [0-9]+: iterations run = [0-9]+" runs "${stdout}")
  list(LENGTH runs number_of_runs)
  if(NOT number_of_runs EQUAL 6)
    message(FATAL_ERROR "${transform}: expected 6 levels (2 time points x 3 levels), got ${number_of_runs}\n${stdout}")
  endif()
  foreach(run IN LISTS runs)
    string(REGEX REPLACE "Level ([0-9]+): iterations run = ([0-9]+)" "\\1;\\2" level_and_count "${run}")
    list(GET level_and_count 0 level)
    list(GET level_and_count 1 count)
    list(GET maximum_iterations ${level} maximum)
    # Level 0 must run, and no level can run more than its -i value.
    if(count GREATER maximum OR (level EQUAL 0 AND count EQUAL 0))
      message(FATAL_ERROR "${transform}: level ${level} ran ${count} iterations, the maximum is ${maximum}\n  ${runs}\n${stdout}")
    endif()
  endforeach()
endforeach()
