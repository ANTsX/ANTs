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

#include "itkDiscreteGaussianImageFilter.h"

//  RecursiveAverageImages img1  img2 weightonimg2 outputname

// We divide the 2nd input image by its mean and add it to the first
// input image with weight 1/n.
// The output overwrites the 1st img with the sum.

#include <list>
#include <vector>
#include <fstream>
#include "vnl/vnl_vector.h"

#include "itkMinimumMaximumImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkCastImageFilter.h"
#include "ReadWriteData.h"

namespace ants
{
template <unsigned int ImageDimension>
int
LabelUniquely(int argc, char * argv[])
{
  using PixelType = float;

  using ImageType = itk::Image<PixelType, ImageDimension>;

  using ULPixelType = unsigned int;
  using labelimagetype = itk::Image<ULPixelType, ImageDimension>;
  using CastFilterType = itk::CastImageFilter<ImageType, labelimagetype>;
  using CastFilterType2 = itk::CastImageFilter<labelimagetype, ImageType>;

  using FilterType = itk::ConnectedComponentImageFilter<labelimagetype, labelimagetype, labelimagetype>;
  using RelabelType = itk::RelabelComponentImageFilter<labelimagetype, labelimagetype>;

  // want the average value in each cluster as defined by the mask and the value thresh and the clust thresh

  if (argc < 2)
  {
    std::cout << "missing 1st filename" << std::endl;
    throw;
  }
  if (argc < 3)
  {
    std::cout << "missing 2nd filename" << std::endl;
    throw;
  }
  if (argc < 4)
  {
    std::cout << "missing cluster thresholod" << std::endl;
    throw;
  }
  bool fullyConnected = false;
  if (argc > 5)
  {
    fullyConnected = static_cast<bool>(std::stoi(argv[4]));
  }
  std::string fn1 = std::string(argv[1]);
  float       clusterthresh = atof(argv[3]);

  typename ImageType::Pointer image1 = nullptr;

  ReadImage<ImageType>(image1, fn1.c_str());

  typename FilterType::Pointer  filter = FilterType::New();
  typename RelabelType::Pointer relabel = RelabelType::New();

  typename CastFilterType::Pointer castInput = CastFilterType::New();
  castInput->SetInput(image1);

  filter->SetInput(castInput->GetOutput());
  filter->SetFullyConnected(fullyConnected); // old default was false
  filter->SetBackgroundValue(0);
  filter->SetMaskImage(castInput->GetOutput());
  relabel->SetInput(filter->GetOutput());
  relabel->SetMinimumObjectSize((unsigned int)clusterthresh);

  try
  {
    relabel->Update();
  }
  catch (const itk::ExceptionObject & excep)
  {
    std::cout << "Relabel: exception caught !" << std::endl;
    std::cout << excep << std::endl;
  }

  //  float maximum=relabel->GetNumberOfObjects();
  typename CastFilterType2::Pointer castRegions = CastFilterType2::New();
  castRegions->SetInput(relabel->GetOutput());
  castRegions->Update();
  ANTs::WriteImage<ImageType>(castRegions->GetOutput(), argv[2]);

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
LabelClustersUniquely(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "LabelClustersUniquely");

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
    std::cout << "Usage:  " << std::endl;
    std::cout << argv[0]
              << " ImageDimension clustersin.hdr labeledclustersout.hdr   sizethresh optionalBoolFullyConnected"
              << std::endl;
    if (argc >= 2 && (std::string(argv[1]) == std::string("--help") || std::string(argv[1]) == std::string("-h")))
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }

  switch (std::stoi(argv[1]))
  {
    case 2:
    {
      return LabelUniquely<2>(argc, argv + 1);
    }
    break;
    case 3:
    {
      return LabelUniquely<3>(argc, argv + 1);
    }
    break;
    case 4:
    {
      return LabelUniquely<4>(argc, argv + 1);
    }
    break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
} // namespace ants
