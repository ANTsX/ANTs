#include "ReadWriteData.h"

bool ANTSFileExists(const std::string & strFilename)
{

  // Attempt to get the file attributes
  struct stat stFileInfo;
  const int intStat = stat(strFilename.c_str(), &stFileInfo);
  bool blnReturn;
  if( intStat == 0 )
    {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
    }
  else
    {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
    }
  return blnReturn;
}
