# antsAI must reject search options that give no search, or one that never
# ends, before it starts. A search factor or step size that is not positive
# made the search loops run until memory ran out, a negative extent gave no
# starting point, and a translation grid with fewer values than dimensions was
# read past its end. The images are made here, so the test needs no data files.
# A step size given without brackets, and a grid with extra values (as
# antsBrainExtraction.sh passes for 2D images), must still work.
file(MAKE_DIRECTORY "${TEST_OUTPUT_DIR}")

set(values "")
foreach(y RANGE 15)
  foreach(x RANGE 15)
    math(EXPR value "(${x} * 7 + ${y} * 3) % 17")
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
      math(EXPR value "(${x} * 7 + ${y} * 3 + ${z} * 5) % 17")
      string(APPEND values "${value} ")
    endforeach()
  endforeach()
endforeach()
set(image3d "${TEST_OUTPUT_DIR}/image3d.mha")
file(WRITE "${image3d}"
  "ObjectType = Image\nNDims = 3\nBinaryData = False\nElementSpacing = 1 1 1\n"
  "DimSize = 8 8 8\nElementType = MET_FLOAT\nElementDataFile = LOCAL\n${values}\n")

# Runs antsAI with the given search options. expect is "pass" or "fail".
function(search ndims expect)
  if(ndims EQUAL 2)
    set(image "${image2d}")
  else()
    set(image "${image3d}")
  endif()
  set(output "${TEST_OUTPUT_DIR}/search.mat")
  file(REMOVE "${output}")
  # The time limit stops a search that never ends, which also grows memory.
  execute_process(
    COMMAND "${ANTS_AI}" -d ${ndims} -m "Mattes[${image},${image},8,None]" -t "Rigid[0.1]" -c 1 -o "${output}" ${ARGN}
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr TIMEOUT 60)
  if(expect STREQUAL "pass")
    if(NOT "${result}" STREQUAL "0" OR NOT EXISTS "${output}")
      message(FATAL_ERROR "antsAI failed (${result}) with ${ARGN}:\n${stdout}\n${stderr}")
    endif()
  else()
    if("${result}" STREQUAL "0")
      message(FATAL_ERROR "antsAI accepted ${ARGN}")
    endif()
    if(NOT stderr MATCHES "\\((-s|-g)\\)")
      message(FATAL_ERROR "antsAI did not say which option is wrong for ${ARGN} (${result}):\n${stdout}\n${stderr}")
    endif()
  endif()
endfunction()

search(2 fail -s 0)
search(2 fail -s "[-20,1]")
search(2 fail -s "[20,-0.5]")
search(2 fail -g "[0,10x10]")
search(2 fail -g "[-5,10x10]")
search(2 fail -g "[10,-10x10]")
search(3 fail -g "[10,10x10]")
search(2 pass -s "[90,0.5]" -g 10)
search(2 pass -s "[90,0.5]" -g "[10,10x10x10]")
search(3 pass -s "[90,0.5]" -g "[10,10x0x0]")
