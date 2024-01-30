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
#include "antsAllocImage.h"
#include <algorithm>

#include <iostream>
#include <fstream>
#include <cstdio>
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

namespace ants
{
template <unsigned int ImageDimension>
int
SetOrigin(int itkNotUsed(argc), char * argv[])
{
  using outPixelType = float;
  using inPixelType = float;
  using ImageType = itk::Image<inPixelType, ImageDimension>;
  using OutImageType = itk::Image<outPixelType, ImageDimension>;
  using readertype = itk::ImageFileReader<ImageType>;
  using writertype = itk::ImageFileWriter<OutImageType>;

  typename readertype::Pointer reader = readertype::New();
  reader->SetFileName(argv[1]);
  reader->Update();

  typename OutImageType::Pointer   outim = reader->GetOutput();
  typename OutImageType::PointType orig = outim->GetOrigin();

  std::cout << " Old Orig " << outim->GetOrigin();

  for (unsigned int d = 0; d < ImageDimension; d++)
  {
    orig[d] = atof(argv[d+3]);
  }
  std::cout << "  New Orig " << orig << std::endl;

  outim->SetOrigin(orig);

  using Iterator = itk::ImageRegionIteratorWithIndex<ImageType>;
  typename ImageType::Pointer varimage =
    AllocImage<ImageType>(outim->GetLargestPossibleRegion(), outim->GetSpacing(), orig, outim->GetDirection());

  Iterator vfIter2(varimage, varimage->GetLargestPossibleRegion());
  for (vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2)
  {
    vfIter2.Set(outim->GetPixel(vfIter2.GetIndex()));
  }

  typename writertype::Pointer writer = writertype::New();
  writer->SetFileName(argv[2]);
  writer->SetInput(varimage);
  writer->Update();
  writer->Write();

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
SetOrigin(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "SetOrigin");

  int     argc = args.size();
  char ** argv = new char *[args.size() + 1];
  for (unsigned int i = 0; i < args.size(); ++i)
  {
    // allocate space for the string plus a null character
    argv[i] = new char[args[i].length() + 1];
    std::strncpy(argv[i], args[i].c_str(), args[i].length());
    // place the null character in the end
    argv[i][args[i].length()] = '\0';
  }
  argv[argc] = nullptr;
  // class to automatically cleanup argv upon destruction
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

  // antscout->set_stream( out_stream );

  if (argc < 3)
  {
    std::cout << "Usage:   " << argv[0] << "  Dimension infile.hdr outfile.nii  OriginX OriginY {OriginZ} "
              << std::endl;
    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  // Get the image dimension
  switch (std::stoi(argv[1]))
  {
    case 2:
    {
      return SetOrigin<2>(argc - 1, argv + 1);
    }
    break;
    case 3:
    {
      return SetOrigin<3>(argc - 1, argv + 1);
    }
    break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
} // namespace ants
