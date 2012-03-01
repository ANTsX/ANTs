MESSAGE("fixfile = ${fixfile}")
FILE(READ ${fixfile} cmakefile)
string(REPLACE "-mlong-branch" ""
  cmakefile "${cmakefile}")
file(WRITE ${fixfile} "${cmakefile}")
