set(Boost_Install_Dir ${CMAKE_CURRENT_BINARY_DIR}/boost-install)
set(Boost_Configure_Script ${CMAKE_CURRENT_LIST_DIR}/configureboost.cmake)
set(Boost_Build_Script ${CMAKE_CURRENT_LIST_DIR}/buildboost.cmake)



ExternalProject_add(Boost
# SVN is too slow SVN_REPOSITORY http://svn.boost.org/svn/boost/trunk
# SVN is too slow SVN_REVISION -r "82586"
  URL http://sourceforge.net/projects/boost/files/boost/1.53.0/boost_1_53_0.tar.gz
  URL_MD5 57a9e2047c0f511c4dfcf00eb5eb2fbb
  SOURCE_DIR Boost
  ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
  CONFIGURE_COMMAND ${CMAKE_COMMAND}
  -DBUILD_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/Boost
  -DBOOST_INSTALL_DIR:PATH=${Boost_Install_Dir}
  -P ${Boost_Configure_Script}
  INSTALL_COMMAND ""
  BUILD_IN_SOURCE 1
  BUILD_COMMAND ${CMAKE_COMMAND}
  -DBUILD_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/Boost
  -DBOOST_INSTALL_DIR:PATH=${Boost_Install_Dir} -P ${Boost_Build_Script}
)
set(BOOST_ROOT ${CMAKE_CURRENT_BINARY_DIR}/boost-install)

set(BOOST_INCLUDE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/boost-install/include)
