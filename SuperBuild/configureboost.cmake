#
# Follow the boost suggestions
#
execute_process(COMMAND ./bootstrap.sh --prefix=${BOOST_INSTALL_DIR}
  --with-libraries=
  WORKING_DIRECTORY ${BUILD_DIR} RESULT_VARIABLE boostrap_result)
### NOTE:  --with-libraries= is purposefull left blank to avoid building
###        any of the unnecessary boost libraries.  ANTS only needs
###        the header-only components of boost!

return(${bootstrap_result})
