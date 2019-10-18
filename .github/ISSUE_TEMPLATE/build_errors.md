---
name:  Build problems
about: Report a problem compiling or installing ANTs from source

---

Before opening an issue, please review the build documentation here

https://github.com/ANTsX/ANTs/wiki/Compiling-ANTs-on-Linux-and-Mac-OS

If you are building the latest code snapshot, please also check Travis
to see if the code is building successfully

https://travis-ci.org/ANTsX/ANTs

If you are using system ITK or VTK, please verify you have the correct
versions (see the Wiki link above). If you can, please also attempt a
Superbuild and let us know if that compiles successfully.

**When did the error occur?**
[ ] CMake configuration (cmake / ccmake)
[ ] Compilation (make)
[ ] Installation (make install)

**Build environment**
 - OS name and version: [e.g. Ubuntu 18.04]
 
**ANTs version**
Specify the release tag or commit hash of the code you are building
from. If you downloaded a snapshot as a ZIP file, please provide the
date of the download.

**Build configuration and logs**
Please attach the following files (relative to your build directory)

  - build.log (the terminal output from make)
  - CMakeCache.txt
  - CMakeFiles/CMakeError.log
  - CMakeFiles/CMakeOutput.log

**Additional context**
Add any other context about the problem here.
