/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
  https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "antsUtilities.h"

#include "itkTrxFileReader.h"
#include "itkTrxGroup.h"
#include "itkTrxStreamlineData.h"

#include <iostream>
#include <string>
#include <vector>

namespace ants
{

namespace
{

const char *
CoordinateTypeToString(itk::TrxStreamlineData::CoordinateType ct)
{
  switch (ct)
  {
    case itk::TrxStreamlineData::CoordinateType::Float16:
      return "float16";
    case itk::TrxStreamlineData::CoordinateType::Float32:
      return "float32";
    case itk::TrxStreamlineData::CoordinateType::Float64:
      return "float64";
    default:
      return "unknown";
  }
}

void
PrintDpsSummary(const itk::TrxStreamlineData * data, bool verboseFields)
{
  const auto dpsNames = data->GetDpsFieldNames();
  std::cout << "  DPS Fields         : " << dpsNames.size() << "\n";
  const auto nStreamlines = static_cast<size_t>(data->GetNumberOfStreamlines());
  for (const auto & name : dpsNames)
  {
    if (verboseFields)
    {
      const auto values = data->GetDpsField(name);
      size_t     nCols = 0;
      if (nStreamlines > 0 && values.size() % nStreamlines == 0)
      {
        nCols = values.size() / nStreamlines;
      }
      std::cout << "    - " << name << ": values=" << values.size();
      if (nCols > 0)
      {
        std::cout << " (" << nCols << " col)";
      }
      std::cout << "\n";
    }
    else
    {
      std::cout << "    - " << name << "\n";
    }
  }
}

void
PrintDpvSummary(const itk::TrxStreamlineData * data, bool verboseFields)
{
  const auto dpvNames = data->GetDpvFieldNames();
  std::cout << "  DPV Fields         : " << dpvNames.size() << "\n";
  const auto nVertices = static_cast<size_t>(data->GetNumberOfVertices());
  for (const auto & name : dpvNames)
  {
    if (verboseFields)
    {
      const auto values = data->GetDpvField(name);
      size_t     nCols = 0;
      if (nVertices > 0 && values.size() % nVertices == 0)
      {
        nCols = values.size() / nVertices;
      }
      std::cout << "    - " << name << ": values=" << values.size();
      if (nCols > 0)
      {
        std::cout << " (" << nCols << " col)";
      }
      std::cout << "\n";
    }
    else
    {
      std::cout << "    - " << name << "\n";
    }
  }
}

void
PrintGroupSummary(const itk::TrxStreamlineData * data)
{
  const auto groupNames = data->GetGroupNames();
  std::cout << "  Groups             : " << groupNames.size() << "\n";
  for (const auto & name : groupNames)
  {
    const size_t count = static_cast<size_t>(data->GetGroupStreamlineCount(name));
    std::cout << "    - " << name << ": " << count << " streamlines\n";
  }
}

} // anonymous namespace

int
PrintTRXHeader(std::vector<std::string> args, std::ostream * /*out_stream*/)
{
  args.insert(args.begin(), "PrintTRXHeader");

  int    argc = static_cast<int>(args.size());
  char** argv = new char *[args.size() + 1];
  for (unsigned int i = 0; i < args.size(); ++i)
  {
    argv[i] = new char[args[i].length() + 1];
    std::strncpy(argv[i], args[i].c_str(), args[i].length());
    argv[i][args[i].length()] = '\0';
  }
  argv[argc] = nullptr;

  class Cleanup_argv
  {
  public:
    Cleanup_argv(char ** argv_, int argc_plus_one_)
      : argv(argv_)
      , argc_plus_one(argc_plus_one_)
    {}
    ~Cleanup_argv()
    {
      for (unsigned int i = 0; i < argc_plus_one; ++i)
      {
        delete[] argv[i];
      }
      delete[] argv;
    }

  private:
    char **      argv;
    unsigned int argc_plus_one;
  };
  Cleanup_argv cleanup_argv(argv, argc + 1);

  if (argc < 2 || (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)))
  {
    std::cout << "Usage: " << argv[0] << " tractogram.trx [--verbose-fields]\n";
    if (argc < 2)
    {
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }

  std::string inputPath;
  bool        verboseFields = false;

  for (int i = 1; i < argc; ++i)
  {
    const std::string arg = argv[i];
    if (arg == "--verbose-fields")
    {
      verboseFields = true;
    }
    else if (inputPath.empty() && arg[0] != '-')
    {
      inputPath = arg;
    }
    else
    {
      std::cerr << "Unknown argument: " << arg << "\n";
      return EXIT_FAILURE;
    }
  }

  if (inputPath.empty())
  {
    std::cerr << "Missing input TRX file.\n";
    std::cout << "Usage: " << argv[0] << " tractogram.trx [--verbose-fields]\n";
    return EXIT_FAILURE;
  }

  try
  {
    auto reader = itk::TrxFileReader::New();
    reader->SetFileName(inputPath);
    reader->Update();
    auto data = reader->GetOutput();
    if (!data)
    {
      std::cerr << "Failed to load TRX data from: " << inputPath << "\n";
      return EXIT_FAILURE;
    }

    std::cout << "\n";
    std::cout << "  File               : " << inputPath << "\n";
    std::cout << "  Streamlines        : " << data->GetNumberOfStreamlines() << "\n";
    std::cout << "  Vertices           : " << data->GetNumberOfVertices() << "\n";
    std::cout << "  Coord dtype        : " << CoordinateTypeToString(data->GetCoordinateType()) << "\n";
    std::cout << "  Dimensions         : " << (data->HasDimensions() ? "present" : "absent") << "\n";
    std::cout << "  VoxelToRAS         : " << (data->HasVoxelToRasMatrix() ? "present" : "absent") << "\n";

    std::cout << "\n";
    PrintDpsSummary(data, verboseFields);
    PrintDpvSummary(data, verboseFields);
    PrintGroupSummary(data);
    std::cout << "\n";
  }
  catch (const itk::ExceptionObject & e)
  {
    std::cerr << "ITK error: " << e << "\n";
    return EXIT_FAILURE;
  }
  catch (const std::exception & e)
  {
    std::cerr << "Error: " << e.what() << "\n";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

} // namespace ants
