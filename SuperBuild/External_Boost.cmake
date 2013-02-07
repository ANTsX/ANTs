set(Boost_Install_Dir ${CMAKE_CURRENT_BINARY_DIR}/boost-install)
set(Boost_Configure_Script ${CMAKE_CURRENT_LIST_DIR}/configureboost.cmake)
set(Boost_Build_Script ${CMAKE_CURRENT_LIST_DIR}/buildboost.cmake)

ExternalProject_add(Boost
  SVN_REPOSITORY http://svn.boost.org/svn/boost/trunk
  SVN_REVISION -r "82586"
#  URL http://sourceforge.net/projects/boost/files/boost/1.49.0/boost_1_49_0.tar.gz
#  URL_MD5 e0defc8c818e4f1c5bbb29d0292b76ca
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
