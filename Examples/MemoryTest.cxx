/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>

#include "ReadWriteData.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkAvantsMutualInformationRegistrationFunction.h"
#include "itkProbabilisticRegistrationFunction.h"
#include "itkCrossCorrelationRegistrationFunction.h"

namespace ants
{
template <unsigned int ImageDimension>
int
MemoryTest(unsigned int argc, char * argv[])
{
  using PixelType = float;
  using VectorType = itk::Vector<float, ImageDimension>;
  using FieldType = itk::Image<VectorType, ImageDimension>;
  using ImageType = itk::Image<PixelType, ImageDimension>;
  using Iterator = itk::ImageRegionIteratorWithIndex<ImageType>;

  // get command line params
  unsigned int argct = 2;
  unsigned int whichmetric = std::stoi(argv[argct]);
  argct++;
  std::string fn1 = std::string(argv[argct]);
  argct++;
  std::string fn2 = std::string(argv[argct]);
  argct++;
  unsigned int numberoffields = 11;
  if (argc > argct)
  {
    numberoffields = std::stoi(argv[argct]);
  }
  argct++;

  typename ImageType::Pointer image1 = nullptr;
  ReadImage<ImageType>(image1, fn1.c_str());
  typename ImageType::Pointer image2 = nullptr;
  ReadImage<ImageType>(image2, fn2.c_str());

  std::vector<typename FieldType::Pointer> fieldvec;
  for (unsigned int i = 0; i < numberoffields; i++)
  {
    std::cout << " NFields " << i << " of " << numberoffields << std::endl;
    VectorType zero;
    zero.Fill(0);
    // ORIENTATION ALERT: Original code set spacing/origin without
    // also setting directions.
    typename FieldType::Pointer field = AllocImage<FieldType>(image1, zero);
    fieldvec.push_back(field);
  }

  typename ImageType::Pointer metricimg = AllocImage<ImageType>(image1, 0);
  Iterator                    iter(metricimg, metricimg->GetLargestPossibleRegion());

  using FixedImageType = ImageType;
  using MovingImageType = ImageType;
  using DisplacementFieldType = FieldType;

  // Choose the similarity metric
  using MIMetricType =
    itk::AvantsMutualInformationRegistrationFunction<FixedImageType, MovingImageType, DisplacementFieldType>;
  using CCMetricType =
    itk::CrossCorrelationRegistrationFunction<FixedImageType, MovingImageType, DisplacementFieldType>;
  // typedef itk::LandmarkCrossCorrelationRegistrationFunction<FixedImageType,MovingImageType,DisplacementFieldType>
  // MetricType;
  // typename
  typename MIMetricType::Pointer mimet = MIMetricType::New();
  typename CCMetricType::Pointer ccmet = CCMetricType::New();

  //  int nbins=32;

  typename CCMetricType::RadiusType ccradius;
  ccradius.Fill(4);
  typename MIMetricType::RadiusType miradius;
  miradius.Fill(0);

  //  mimet->SetDisplacementField(field);
  mimet->SetFixedImage(image1);
  mimet->SetMovingImage(image2);
  mimet->SetRadius(miradius);
  mimet->SetGradientStep(1.e2);
  mimet->SetNormalizeGradient(false);

  //  ccmet->SetDisplacementField(field);
  ccmet->SetFixedImage(image1);
  ccmet->SetMovingImage(image2);
  ccmet->SetRadius(ccradius);
  ccmet->SetGradientStep(1.e2);
  ccmet->SetNormalizeGradient(false);

  if (whichmetric == 1) // imagedifference
  {
    ccmet->InitializeIteration();
  }
  else if (whichmetric != 0)
  {
    mimet->InitializeIteration();
  }

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
MemoryTest(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "MemoryTest");

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
    std::cout << "Basic useage ex: " << std::endl;
    std::cout << argv[0] << " ImageDimension whichmetric image1.ext image2.ext NumberOfFieldsToAllocate " << std::endl;
    std::cout << "  outimage and logfile are optional  " << std::endl;
    std::cout << "  Metric 0 - MeanSquareDifference, 1 - Cross-Correlation, 2-Mutual Information  " << std::endl;
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
      return MemoryTest<2>(argc, argv);
    }
    break;
    case 3:
    {
      return MemoryTest<3>(argc, argv);
    }
    break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
} // namespace ants
