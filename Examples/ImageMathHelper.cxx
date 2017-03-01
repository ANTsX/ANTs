// Non-templated functions used by ImageMath templates.
#include <string>
#include <map>
#include <iostream>
#include <fstream>

namespace ants
{

std::string ANTSGetFilePrefix(const char *str)
{
  const std::string      filename = str;
  const std::string::size_type pos = filename.rfind( "." );
  const std::string            filepre = std::string( filename, 0, pos );

#if 0 // HACK:  This does nothing useful
  if( pos != std::string::npos )
    {
    std::string extension = std::string( filename, pos, filename.length() - 1);
    if( extension == std::string(".gz") )
      {
      pos = filepre.rfind( "." );
      // extension = std::string( filepre, pos, filepre.length() - 1 );
      }
    }
#endif
  return filepre;
}

std::string ANTSOptionName(const char *str)
{
  std::string            filename = str;
  std::string::size_type pos = filename.rfind( "=" );
  std::string            name = std::string( filename, 0, pos );

  return name;
}

std::string ANTSOptionValue(const char *str)
{
  std::string            filename = str;
  std::string::size_type pos = filename.rfind( "=" );
  std::string            value = std::string( filename, pos + 1, filename.length() );

  return value;
}


// int is the key, string the return value
std::map<unsigned int, std::string> RoiList(std::string file)
{
  unsigned int wordindex = 0;
  std::string  tempstring = "";

  std::map<unsigned int, std::string> RoiList;
  //  RoiList[0]=std::string("Background");
  char         str[2000];
  std::fstream file_op(file.c_str(), std::ios::in);

  while( file_op >> str )
    {
    tempstring = std::string(str);
    RoiList[wordindex] = tempstring;
    wordindex++;
    }

  return RoiList; // returns the maximum index
}


} // namespace ants
