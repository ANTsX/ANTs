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
#include <algorithm>

#include "itkImage.h"
#include "itkConstantPadImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkFlipImageFilter.h"

#include "ReadWriteData.h"

namespace ants
{
template <unsigned int Dimension>
int
PermuteFlipImageOrientationAxes(int argc, char * argv[])
{
  using InputPixelType = float;
  using OutputPixelType = float;

  using InputImageType = itk::Image<InputPixelType, Dimension>;
  using OutputImageType = itk::Image<OutputPixelType, Dimension>;

  typename InputImageType::Pointer inputImage = nullptr;
  ReadImage<InputImageType>(inputImage, argv[1]);

  // Create a filter
  using ShortImage = OutputImageType;
  typename itk::PermuteAxesImageFilter<ShortImage>::Pointer permute;
  permute = itk::PermuteAxesImageFilter<ShortImage>::New();
  permute->SetInput(inputImage);

  using FlipType = itk::FlipImageFilter<ShortImage>;
  typename FlipType::FlipAxesArrayType lowerFactors;

  typename itk::PermuteAxesImageFilter<ShortImage>::PermuteOrderArrayType upperFactors;

  for (unsigned int q = 0; q < Dimension; ++q)
  {
    upperFactors[q] = 0;
    lowerFactors[q] = 0;
  }

  bool flipaboutorigin = false;
  if (Dimension == 2)
  {
    if (argc > 3)
    {
      upperFactors[0] = std::stoi(argv[3]);
    }
    if (argc > 4)
    {
      upperFactors[1] = std::stoi(argv[4]);
    }
    if (argc > 5)
    {
      lowerFactors[0] = std::stoi(argv[5]);
    }
    if (argc > 6)
    {
      lowerFactors[1] = std::stoi(argv[6]);
    }
    if (argc > 7)
    {
      flipaboutorigin = std::stoi(argv[7]);
    }
  }
  else if (Dimension == 3)
  {
    if (argc > 3)
    {
      upperFactors[0] = std::stoi(argv[3]);
    }
    if (argc > 4)
    {
      upperFactors[1] = std::stoi(argv[4]);
    }
    if (argc > 5)
    {
      upperFactors[2] = std::stoi(argv[5]);
    }
    if (argc > 6)
    {
      lowerFactors[0] = std::stoi(argv[6]);
    }
    if (argc > 7)
    {
      lowerFactors[1] = std::stoi(argv[7]);
    }
    if (argc > 8)
    {
      lowerFactors[2] = std::stoi(argv[8]);
    }
    if (argc > 9)
    {
      flipaboutorigin = std::stoi(argv[9]);
    }
  }


  permute->SetOrder(upperFactors);
  permute->Update();

  typename FlipType::FlipAxesArrayType flip;
  for (unsigned int i = 0; i < Dimension; i++)
  {
    flip[i] = lowerFactors[i];
  }
  typename FlipType::Pointer flipper = FlipType::New();
  flipper->SetFlipAboutOrigin(flipaboutorigin);
  flipper->SetFlipAxes(flip);
  flipper->SetInput(permute->GetOutput());
  flipper->Update();

  typename InputImageType::Pointer image = flipper->GetOutput();
  ANTs::WriteImage<OutputImageType>(image, argv[2]);

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
PermuteFlipImageOrientationAxes(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "PermuteFlipImageOrientationAxes");

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
    std::cout << "Usage: " << std::endl;
    std::cout
      << argv[0]
      << " ImageDimension  inputImageFile  outputImageFile xperm yperm {zperm}  xflip yflip {zflip}  {FlipAboutOrigin}"
      << std::endl;
    std::cout << " for 3D:  " << argv[0]
              << " 3  in.nii out.nii   2 0 1  1 1 1  \n would map z=>x, x=>y, y=>z and flip each " << std::endl;
    std::cout << " for 2D:  " << argv[0] << " 2  in.nii out.nii   1 0  1 0  \n would map x=>y, y=>x and flip x  "
              << std::endl;
    std::cout << std::endl;
    std::cout << " 0 1 2 for permute factors gives no axis permutation " << std::endl;
    std::cout << " 1 2 0 maps y to x,  z to y and x to z " << std::endl;
    std::cout << " the flip values are boolean -  0 1 0 would flip the y-axis only " << std::endl;
    std::cout << std::endl
              << " The FlipAboutOrigin boolean lets you flip about the coordinate set in the origin " << std::endl;
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
      return PermuteFlipImageOrientationAxes<2>(argc - 1, argv + 1);
    }
    break;
    case 3:
    {
      return PermuteFlipImageOrientationAxes<3>(argc - 1, argv + 1);
    }
    break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
} // namespace ants
