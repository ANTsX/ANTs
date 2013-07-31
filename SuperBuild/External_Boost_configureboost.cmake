#
# Follow the boost suggestions
#

if(WIN32) # bootstrap.bat has no options, the options are given to ./b2 when BUILDING: see buildboost.cmake

  execute_process(COMMAND bootstrap.bat
    WORKING_DIRECTORY ${BUILD_DIR} RESULT_VARIABLE boostrap_result)

else(WIN32)

  execute_process(COMMAND ./bootstrap.sh --prefix=${BOOST_INSTALL_DIR}
    --without-libraries=atomic,chrono,context,date_time,exception,filesystem,graph,graph_parallel,iostreams,locale,log,math,mpi,program_options,python,random,regex,serialization,signals,test,timer,wave
## Needed for UKF system,thread,
    WORKING_DIRECTORY ${BUILD_DIR} RESULT_VARIABLE boostrap_result)
### NOTE:  --with-libraries= is purposefull left blank to avoid building
###        any of the unnecessary boost libraries.  ANTS only needs
###        the header-only components of boost!

endif(WIN32)

return(${bootstrap_result})
