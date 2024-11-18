#include "ReadWriteData.h"

#include "itksys/SystemTools.hxx"

bool
ANTSFileExists(const std::string & strFilename)
{
  // ITK checks file existence on all platforms, also read permissions on POSIX systems
  return itksys::SystemTools::FileExists(strFilename);
}

bool
ANTSFileIsImage(const std::string &filename)
{
  if (!ANTSFileExists(filename))
  {
    return false;
  }

  // Check if the file is recognized as a valid image
  itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
    filename.c_str(), itk::ImageIOFactory::IOFileModeEnum::ReadMode);
  if (!imageIO)
  {
    return false;
  }

  // File passed both checks
  return true;
}