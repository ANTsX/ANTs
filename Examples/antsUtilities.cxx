
#include "antsUtilities.h"

#include <deque>
#include <string>
#include <vector>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>

namespace ants
{
// We need to ensure that only one of these exists!
// boost::iostreams::stream<ants_Sink> std::cout( ( ants_Sink() ) );
}

TRAN_FILE_TYPE
CheckFileType(const char * const str)
{
  std::string            filename = str;
  std::string::size_type pos = filename.rfind(".");
  std::string            filepre = std::string(filename, 0, pos);

  if (pos != std::string::npos)
  {
    std::string extension = std::string(filename, pos, filename.length() - 1);
    if (extension == std::string(".gz"))
    {
      pos = filepre.rfind(".");
      extension = std::string(filepre, pos, filepre.length() - 1);
    }
    if (extension == ".txt" || extension == ".mat" || extension == ".hdf5" || extension == ".hdf" ||
        extension == ".xfm")
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
}

TRAN_FILE_TYPE
CheckFileType(const std::string & str)
{
  return CheckFileType(str.c_str());
}

void
SetAffineInvFlag(TRAN_OPT & opt, bool & set_current_affine_inv)
{
  opt.do_affine_inv = set_current_affine_inv;
  if (set_current_affine_inv)
  {
    set_current_affine_inv = false;
  }
}

void
FilePartsWithgz(const std::string & filename, std::string & path, std::string & name, std::string & ext)
{
  std::string            extension;
  std::string::size_type pos = filename.rfind(".");
  std::string            filepre = std::string(filename, 0, pos);

  if (pos != std::string::npos)
  {
    extension = std::string(filename, pos, filename.length() - 1);
    if (extension == std::string(".gz"))
    {
      pos = filepre.rfind(".");
      if (pos != std::string::npos)
      {
        extension = std::string(filepre, pos, filepre.length() - 1) + ".gz";
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

  if (pos != std::string::npos)
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

bool
CheckFileExistence(const char * const str)
{
  std::ifstream myfile(str);
  bool          b = myfile.is_open();

  myfile.close();
  return b;
}

// adapted from http://stackoverflow.com/questions/194465/how-to-parse-a-string-to-an-int-in-c
bool
get_a_double_number(const char * const str, double & v)
{
  errno = 0;
  char * end;
  v = strtod(str, &end);
  if ((errno == ERANGE && v == HUGE_VAL))
  {
    return false; // OVERFLOW
  }
  if ((errno == ERANGE && itk::Math::FloatAlmostEqual(v, 0.0)))
  {
    return false; // UNDERFLOW
  }
  if (*str == '\0' || *end != '\0')
  {
    return false; // INCONVERTIBLE;
  }
  return true;
}

void
DisplayOptQueue(const TRAN_OPT_QUEUE & opt_queue)
{
  const itk::SizeValueType kQueueSize = opt_queue.size();

  for (itk::SizeValueType i = 0; i < kQueueSize; i++)
  {
    std::cout << "[" << i << "/" << kQueueSize << "]: ";

    switch (opt_queue[i].file_type)
    {
      case AFFINE_FILE:
      {
        std::cout << "AFFINE";
      }
      break;
      case DEFORMATION_FILE:
      {
        std::cout << "FIELD";
      }
      break;
      case IDENTITY_TRANSFORM:
      {
        std::cout << "IDENTITY";
      }
      break;
      case IMAGE_AFFINE_HEADER:
      {
        std::cout << "HEADER";
      }
      break;
      case INVALID_FILE:
      default:
      {
        std::cout << "Invalid Format!!!";
      }
      break;
    }
    if (opt_queue[i].do_affine_inv)
    {
      std::cout << "-INV";
    }
    std::cout << ": " << opt_queue[i].filename << std::endl;
  }
}

void
DisplayOpt(const TRAN_OPT & opt)
{
  switch (opt.file_type)
  {
    case AFFINE_FILE:
    {
      std::cout << "AFFINE";
    }
    break;
    case DEFORMATION_FILE:
    {
      std::cout << "FIELD";
    }
    break;
    case IDENTITY_TRANSFORM:
    {
      std::cout << "IDENTITY";
    }
    break;
    case IMAGE_AFFINE_HEADER:
    {
      std::cout << "HEADER";
    }
    break;
    case INVALID_FILE:
    default:
    {
      std::cout << "Invalid Format!!!";
    }
    break;
  }
  if (opt.do_affine_inv)
  {
    std::cout << "-INV";
  }
  std::cout << ": " << opt.filename << std::endl;
}

std::string
GetPreferredTransformFileType()
{
  // return ".mat";
  return ".txt";
}

void
ConvertToLowerCase(std::string & str)
{
  std::transform(str.begin(), str.end(), str.begin(), tolower);
  // You may need to cast the above line to (int(*)(int))
  // tolower - this works as is on VC 7.1 but may not work on
  // other compilers
}
