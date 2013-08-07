set(vtkDetCFLAGS
  ${VTKSource}/CMake/vtkDetermineCompilerFlags.cmake)

file(READ ${vtkDetCFLAGS}
  code)

string(REPLACE
"SET(VTK_REQUIRED_C_FLAGS \"\${VTK_REQUIRED_C_FLAGS} -mlong-branch\")"
""
code "${code}")
string(REPLACE
"SET(VTK_REQUIRED_CXX_FLAGS \"\${VTK_REQUIRED_CXX_FLAGS} -mlong-branch\")"
""
code "${code}")

file(WRITE ${vtkDetCFLAGS}
  "${code}"
  )

set(ftglCMakeLists_txt ${VTKSource}/Utilities/ftgl/CMakeLists.txt)
file(READ ${ftglCMakeLists_txt}
  code)
string(REPLACE " -fpascal-strings" "" code "${code}")

file(WRITE ${ftglCMakeLists_txt} "${code}")

find_file(vtkVRMLImporter vtkVRMLImporter.cxx
  HINTS ${VTKSource}/Hybrid ${VTKSource}/IO/IMPORT
)

file(READ ${vtkVRMLImporter}
  code)

string(REPLACE
"#ifdef __GNUC__
#undef alloca
#define alloca __builtin_alloca
"
"#ifdef __GNUC__
#ifndef __clang__
#undef alloca
#define alloca __builtin_alloca
#endif
"
code "${code}")

file(WRITE ${vtkVRMLImporter}
  "${code}"
  )
