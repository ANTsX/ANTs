#include "ReadWriteData.h"

#include "itksys/SystemTools.hxx"

bool ANTSFileExists(const std::string & strFilename)
{
  return itksys::SystemTools::FileExists(strFilename);
}
