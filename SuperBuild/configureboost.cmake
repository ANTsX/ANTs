#
# Follow the boost suggestions
#
execute_process(COMMAND ./bootstrap.sh --prefix=${BOOST_INSTALL_DIR}
  WORKING_DIRECTORY ${BUILD_DIR} RESULT_VARIABLE boostrap_result)

return(${bootstrap_result})
