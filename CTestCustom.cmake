#-- #NOTES from http: // www.cmake.org/Wiki/CMake_Testing_With_CTest
set(CTEST_CUSTOM_MEMCHECK_IGNORE
    $ {CTEST_CUSTOM_MEMCHECK_IGNORE}
    DummyExcludeMemcheckIgnoreTestSetGet
    )

#-- #set(CTEST_CUSTOM_WARNING_MATCH
#-- #${CTEST_CUSTOM_WARNING_MATCH}
#-- #"{standard input}:[0-9][0-9]*: Warning: "
#-- #)

#-- #if("@CMAKE_SYSTEM@" MATCHES "OSF")
set(CTEST_CUSTOM_WARNING_EXCEPTION
    $ {CTEST_CUSTOM_WARNING_EXCEPTION}
    "fltk"
    "xforms"
    "vtkKWApplication"
    "vtkKWObject"
    )
#-- #endif("@CMAKE_SYSTEM@" MATCHES "OSF")

#-- #The following are brains2 warnings that just need to be suppressed because they are caused
#-- #by third parties, and will never be fixed.
set(CTEST_CUSTOM_WARNING_EXCEPTION
  "Point.*may be used uninitialized"
  ".*was declared here.*"
    )
##Intel compiler does not like itkLegacyMacro warning #1170

#-- #Reset maximum number of warnings so that they all show up.
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 1000)

set(CTEST_CUSTOM_COVERAGE_EXCLUDE $ {CTEST_CUSTOM_COVERAGE_EXCLUDE}
    "./SlicerExecutionModel/"
    "./SlicerExecutionModel/.*"
    "./SlicerExecutionModel/.*/.*"
    ".*SlicerExecutionModel.*"
    "SlicerExecutionModel"
    )
