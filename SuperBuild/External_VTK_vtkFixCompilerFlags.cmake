message("fixfile = ${fixfile}")
file(READ ${fixfile} cmakefile)
string(REPLACE "-mlong-branch" ""
  cmakefile "${cmakefile}")
file(WRITE ${fixfile} "${cmakefile}")
