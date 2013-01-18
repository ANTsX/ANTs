#
# Follow the boost suggestions
#
execute_process(COMMAND ./b2 install
  WORKING_DIRECTORY ${BUILD_DIR} RESULT_VARIABLE build_result)

return(${build_result})
