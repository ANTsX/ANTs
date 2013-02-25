
#include "antsUtilities.h"

#include <deque>
#include <string>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>

namespace ants
{
// We need to ensure that only one of these exists!
boost::iostreams::stream<ants_Sink> antscout( ( ants_Sink() ) );
}

TRAN_FILE_TYPE CheckFileType(const char * const str)
{
  std::string            filename = str;
  std::string::size_type pos = filename.rfind( "." );
  std::string            filepre = std::string( filename, 0, pos );

  if( pos != std::string::npos )
    {
    std::string extension = std::string( filename, pos, filename.length() - 1);
    if( extension == std::string(".gz") )
      {
      pos = filepre.rfind( "." );
      extension = std::string( filepre, pos, filepre.length() - 1 );
      }
    if( extension == ".txt" || extension == ".mat" || extension == ".hdf5" || extension == ".hdf" )
      {
      return AFFINE_FILE;
      }
    else
      {
      return DEFORMATION_FILE;
      }
    }
  else
    {
    return INVALID_FILE;
    }
  return AFFINE_FILE;
}

TRAN_FILE_TYPE CheckFileType(const std::string & str)
{
  return CheckFileType( str.c_str() );
}

void SetAffineInvFlag(TRAN_OPT & opt, bool & set_current_affine_inv)
{
  opt.do_affine_inv = set_current_affine_inv;
  if( set_current_affine_inv )
    {
    set_current_affine_inv = false;
    }
}

void FilePartsWithgz(const std::string & filename, std::string & path, std::string & name, std::string & ext)
{
  std::string            extension;
  std::string::size_type pos = filename.rfind( "." );
  std::string            filepre = std::string( filename, 0, pos );

  if( pos != std::string::npos )
    {
    extension = std::string( filename, pos, filename.length() - 1);
    if( extension == std::string(".gz") )
      {
      pos = filepre.rfind( "." );
      if( pos != std::string::npos )
        {
        extension = std::string( filepre, pos, filepre.length() - 1 ) + ".gz";
        filepre = std::string(filepre, 0, pos);
        }
      }
    }
  else
    {
    extension = std::string("");
    }

  ext = extension;

  pos = filepre.rfind('/');

  if( pos != std::string::npos )
    {
    path = std::string(filepre, 0, pos + 1);
    name = std::string(filepre, pos + 1, filepre.length() - 1);
    }
  else
    {
    path = std::string("");
    name = filepre;
    }
}

bool CheckFileExistence(const char * const str)
{
  std::ifstream myfile(str);
  bool          b = myfile.is_open();

  myfile.close();
  return b;
}

// adapted from http://stackoverflow.com/questions/194465/how-to-parse-a-string-to-an-int-in-c
bool get_a_double_number(const char * const str, double & v)
{
  errno = 0;
  char *end;
  v = strtod(str, &end);
  if( (errno == ERANGE && v == HUGE_VAL) )
    {
    return false;     // OVERFLOW
    }
  if( (errno == ERANGE && v == 0.0) )
    {
    return false;     // UNDERFLOW
    }
  if( *str == '\0' || *end != '\0' )
    {
    return false;     // INCONVERTIBLE;
    }
  return true;
}

void DisplayOptQueue(const TRAN_OPT_QUEUE & opt_queue)
{
  const int kQueueSize = opt_queue.size();

  for( int i = 0; i < kQueueSize; i++ )
    {
    ants::antscout << "[" << i << "/" << kQueueSize << "]: ";

    switch( opt_queue[i].file_type )
      {
      case AFFINE_FILE:
        {
        ants::antscout << "AFFINE";
        }
        break;
      case DEFORMATION_FILE:
        {
        ants::antscout << "FIELD";
        }
        break;
      case IDENTITY_TRANSFORM:
        {
        ants::antscout << "IDENTITY";
        }
        break;
      case IMAGE_AFFINE_HEADER:
        {
        ants::antscout << "HEADER";
        }
        break;
      default:
        {
        ants::antscout << "Invalid Format!!!";
        }
        break;
      }
    if( opt_queue[i].do_affine_inv )
      {
      ants::antscout << "-INV";
      }
    ants::antscout << ": " << opt_queue[i].filename << std::endl;
    }
}

void DisplayOpt(const TRAN_OPT & opt)
{
  switch( opt.file_type )
    {
    case AFFINE_FILE:
      {
      ants::antscout << "AFFINE";
      }
      break;
    case DEFORMATION_FILE:
      {
      ants::antscout << "FIELD";
      }
      break;
    case IDENTITY_TRANSFORM:
      {
      ants::antscout << "IDENTITY";
      }
      break;
    case IMAGE_AFFINE_HEADER:
      {
      ants::antscout << "HEADER";
      }
      break;
    default:
      {
      ants::antscout << "Invalid Format!!!";
      }
      break;
    }
  if( opt.do_affine_inv )
    {
    ants::antscout << "-INV";
    }
  ants::antscout << ": " << opt.filename << std::endl;
}

std::string GetPreferredTransformFileType(void)
{
  // return ".mat";
  return ".txt";
}

void ConvertToLowerCase( std::string& str )
{
  std::transform( str.begin(), str.end(), str.begin(), tolower );
  // You may need to cast the above line to (int(*)(int))
  // tolower - this works as is on VC 7.1 but may not work on
  // other compilers
}
