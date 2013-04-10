#
# Follow the boost suggestions
#
execute_process(COMMAND ./bootstrap.sh --prefix=${BOOST_INSTALL_DIR}
  --without-libraries=atomic,chrono,context,date_time,exception,filesystem,graph,graph_parallel,iostreams,locale,math,mpi,program_options,python,random,regex,serialization,signals,system,test,thread,timer,wave
  WORKING_DIRECTORY ${BUILD_DIR} RESULT_VARIABLE boostrap_result)
### NOTE:  --with-libraries= is purposefull left blank to avoid building
###        any of the unnecessary boost libraries.  ANTS only needs
###        the header-only components of boost!

return(${bootstrap_result})
