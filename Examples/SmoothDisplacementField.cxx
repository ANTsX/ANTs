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
#include "ReadWriteData.h"

#include "itkDisplacementFieldToBSplineImageFilter.h"
#include "itkGaussianOperator.h"
#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkTimeProbe.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkVectorNeighborhoodOperatorImageFilter.h"

namespace ants
{

template <unsigned int ImageDimension>
int
SmoothDisplacementField(int argc, char * argv[])
{

  using RealType = float;
  using RealImageType = itk::Image<RealType, ImageDimension>;
  using VectorType = itk::Vector<RealType, ImageDimension>;
  using DisplacementFieldType = itk::Image<VectorType, ImageDimension>;

  typename DisplacementFieldType::Pointer field = nullptr;
  ReadImage<DisplacementFieldType>(field, argv[2]);

  typename DisplacementFieldType::Pointer smoothField = DisplacementFieldType::New();

  float elapsedTime = 0.0;

  std::vector<float> var = ConvertVector<float>(std::string(argv[4]));
  if (var.size() == 1)
  {
    float variance = var[0];

    using GaussianSmoothingOperatorType = itk::GaussianOperator<float, ImageDimension>;
    using GaussianSmoothingSmootherType =
      itk::VectorNeighborhoodOperatorImageFilter<DisplacementFieldType, DisplacementFieldType>;

    GaussianSmoothingOperatorType gaussianSmoothingOperator;

    using DuplicatorType = itk::ImageDuplicator<DisplacementFieldType>;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(field);
    duplicator->Update();

    smoothField = duplicator->GetOutput();

    itk::TimeProbe timer;
    timer.Start();
    typename GaussianSmoothingSmootherType::Pointer smoother = GaussianSmoothingSmootherType::New();

    for (unsigned int dimension = 0; dimension < ImageDimension; ++dimension)
    {
      // smooth along this dimension
      gaussianSmoothingOperator.SetDirection(dimension);
      gaussianSmoothingOperator.SetVariance(variance);
      gaussianSmoothingOperator.SetMaximumError(0.001);
      gaussianSmoothingOperator.SetMaximumKernelWidth(smoothField->GetRequestedRegion().GetSize()[dimension]);
      gaussianSmoothingOperator.CreateDirectional();

      // todo: make sure we only smooth within the buffered region
      smoother->SetOperator(gaussianSmoothingOperator);
      smoother->SetInput(smoothField);
      smoother->Update();

      smoothField = smoother->GetOutput();
      smoothField->Update();
      smoothField->DisconnectPipeline();
    }

    const VectorType zeroVector(0.0);

    // make sure boundary does not move
    float weight1 = itk::NumericTraits<float>::OneValue();
    if (variance < static_cast<float>(0.5))
    {
      weight1 = itk::NumericTraits<float>::OneValue() -
                itk::NumericTraits<float>::OneValue() * (variance / static_cast<float>(0.5));
    }
    float weight2 = itk::NumericTraits<float>::OneValue() - weight1;

    const typename DisplacementFieldType::RegionType region = field->GetLargestPossibleRegion();
    const typename DisplacementFieldType::SizeType   size = region.GetSize();
    const typename DisplacementFieldType::IndexType  startIndex = region.GetIndex();

    itk::ImageRegionConstIteratorWithIndex<DisplacementFieldType> fieldIt(field, field->GetLargestPossibleRegion());
    itk::ImageRegionIteratorWithIndex<DisplacementFieldType>      smoothedFieldIt(smoothField,
                                                                             smoothField->GetLargestPossibleRegion());
    for (fieldIt.GoToBegin(), smoothedFieldIt.GoToBegin(); !fieldIt.IsAtEnd(); ++fieldIt, ++smoothedFieldIt)
    {
      typename DisplacementFieldType::IndexType index = fieldIt.GetIndex();
      bool                                      isOnBoundary = false;
      for (unsigned int dimension = 0; dimension < ImageDimension; ++dimension)
      {
        if (index[dimension] == startIndex[dimension] ||
            index[dimension] == static_cast<int>(size[dimension]) - startIndex[dimension] - 1)
        {
          isOnBoundary = true;
          break;
        }
      }
      if (isOnBoundary)
      {
        smoothedFieldIt.Set(zeroVector);
      }
      else
      {
        smoothedFieldIt.Set(smoothedFieldIt.Get() * weight1 + fieldIt.Get() * weight2);
      }
    }
    timer.Stop();
    elapsedTime = timer.GetMean();
  }
  else if (var.size() == ImageDimension)
  {
    using BSplineFilterType = itk::DisplacementFieldToBSplineImageFilter<DisplacementFieldType>;

    unsigned int numberOfLevels = 1;
    if (argc > 5)
    {
      numberOfLevels = std::stoi(argv[5]);
    }

    unsigned int splineOrder = 3;
    if (argc > 6)
    {
      splineOrder = std::stoi(argv[6]);
    }

    typename BSplineFilterType::ArrayType ncps;
    for (unsigned int d = 0; d < ImageDimension; d++)
    {
      ncps[d] = var[d] + splineOrder;
    }

    typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
    bspliner->SetDisplacementField(field);
    bspliner->SetBSplineDomainFromImage(field);
    bspliner->SetNumberOfControlPoints(ncps);
    bspliner->SetSplineOrder(splineOrder);
    bspliner->SetNumberOfFittingLevels(numberOfLevels);
    bspliner->SetEnforceStationaryBoundary(true);
    bspliner->SetEstimateInverse(false);

    if (argc > 7)
    {
      bspliner->SetEstimateInverse(static_cast<bool>(std::stoi(argv[7])));
    }

    if (argc > 8)
    {
      typename BSplineFilterType::RealImageType::Pointer confidenceImage = nullptr;
      ReadImage<typename BSplineFilterType::RealImageType>(confidenceImage, argv[8]);
      bspliner->SetConfidenceImage(confidenceImage);
    }

    itk::TimeProbe timer;
    timer.Start();
    bspliner->Update();
    timer.Stop();
    elapsedTime = timer.GetMean();

    smoothField = bspliner->GetOutput();
    smoothField->DisconnectPipeline();
  }
  else
  {
    std::cerr << "Error: unexpected variance format." << std::endl;
    return EXIT_FAILURE;
  }

  typename RealImageType::Pointer rmseImage = RealImageType::New();
  rmseImage->SetOrigin(field->GetOrigin());
  rmseImage->SetDirection(field->GetDirection());
  rmseImage->SetSpacing(field->GetSpacing());
  rmseImage->SetRegions(field->GetLargestPossibleRegion());
  rmseImage->AllocateInitialized();

  itk::ImageRegionConstIterator<DisplacementFieldType> fieldIt(field, field->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<DisplacementFieldType> smoothedFieldIt(smoothField,
                                                                       smoothField->GetLargestPossibleRegion());
  itk::ImageRegionIterator<RealImageType>              ItR(rmseImage, rmseImage->GetLargestPossibleRegion());

  float             rmse = 0.0;
  vnl_vector<float> rmse_comp(2);
  rmse_comp.fill(0.0);
  float N = 0.0;
  for (fieldIt.GoToBegin(), smoothedFieldIt.GoToBegin(), ItR.GoToBegin(); !fieldIt.IsAtEnd();
       ++fieldIt, ++smoothedFieldIt, ++ItR)
  {
    ItR.Set((fieldIt.Get() - smoothedFieldIt.Get()).GetSquaredNorm());
    for (unsigned int d = 0; d < ImageDimension; d++)
    {
      rmse_comp[d] += itk::Math::sqr(fieldIt.Get()[d] - smoothedFieldIt.Get()[d]);
    }
    rmse += static_cast<float>((fieldIt.Get() - smoothedFieldIt.Get()).GetSquaredNorm());
    N += itk::NumericTraits<float>::OneValue();
  }
  rmse = std::sqrt(rmse / N);

  std::cout << "Elapsed time: " << elapsedTime << std::endl;
  std::cout << "RMSE = " << rmse << std::endl;
  for (unsigned int d = 0; d < ImageDimension; d++)
  {
    std::cout << "  rmse[" << d << "] = " << std::sqrt(rmse_comp[d] / N) << std::endl;
  }

  ANTs::WriteImage<DisplacementFieldType>(smoothField, argv[3]);

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
SmoothDisplacementField(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "SmoothDisplacementField");

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

  if (argc < 5)
  {
    std::cout << argv[0] << " imageDimension inputField outputField variance_or_mesh_size_base_level "
              << "[numberOfevels=1] [splineOrder=3] [estimateInverse=0] [confidenceImage]" << std::endl;

    return EXIT_FAILURE;
  }

  switch (std::stoi(argv[1]))
  {
    case 2:
      return SmoothDisplacementField<2>(argc, argv);
      break;
    case 3:
      return SmoothDisplacementField<3>(argc, argv);
      break;
    case 4:
      return SmoothDisplacementField<4>(argc, argv);
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

} // namespace ants
