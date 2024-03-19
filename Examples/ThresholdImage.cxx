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

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>

#include "itkArray.h"
#include "itkExtractImageFilter.h"
#include "itkImage.h"
#include "ReadWriteData.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkListSample.h"
#include "itkMinimumDecisionRule.h"
#include "itkMultiplyImageFilter.h"
#include "itkNeighborhoodIterator.h"
#include "itkOtsuMultipleThresholdsImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkSampleClassifierFilter.h"
#include "itkWeightedCentroidKdTreeGenerator.h"

namespace ants
{
template <typename TImage>
typename TImage::Pointer
MultiplyImage(typename TImage::Pointer image1, typename TImage::Pointer image2)
{
  std::cout << " Multiply " << std::endl;
  // Begin Multiply Images
  using tImageType = TImage;
  //  output will be the speed image for FMM
  using MultFilterType = itk::MultiplyImageFilter<tImageType, tImageType, tImageType>;
  typename MultFilterType::Pointer filter = MultFilterType::New();
  filter->SetInput1(image1);
  filter->SetInput2(image2);
  filter->Update();
  return filter->GetOutput(); // this is the speed image

  // write a function to threshold the speedimage so
  // if the dist is g.t. D then speed = 1
}

template <typename TImage>
typename TImage::Pointer
BinaryThreshold_AltInsideOutside_threashold(typename TImage::PixelType low,
                                            typename TImage::PixelType high,
                                            typename TImage::PixelType insideval,
                                            typename TImage::PixelType outsideval,
                                            typename TImage::Pointer   input)
{

  using PixelType = typename TImage::PixelType;
  // Begin Threshold Image
  using InputThresholderType = itk::BinaryThresholdImageFilter<TImage, TImage>;
  typename InputThresholderType::Pointer inputThresholder = InputThresholderType::New();

  inputThresholder->SetInput(input);
  inputThresholder->SetInsideValue(insideval);
  inputThresholder->SetOutsideValue(outsideval);

  inputThresholder->SetLowerThreshold((PixelType)low);
  inputThresholder->SetUpperThreshold((PixelType)high);
  inputThresholder->Update();

  return inputThresholder->GetOutput();
}

template <typename TImage>
typename TImage::Pointer
LabelSurface(typename TImage::PixelType foreground, typename TImage::PixelType newval, typename TImage::Pointer input)
{
  std::cout << " Label Surf " << std::endl;
  using ImageType = TImage;
  enum
  {
    ImageDimension = ImageType::ImageDimension
  };
  // ORIENTATION ALERT: Original code set spacing & origin without
  // setting directions.
  typename ImageType::Pointer Image = AllocImage<ImageType>(input);

  using iteratorType = itk::NeighborhoodIterator<ImageType>;

  typename iteratorType::RadiusType rad;
  for (int j = 0; j < ImageDimension; j++)
  {
    rad[j] = 1;
  }
  iteratorType GHood(rad, input, input->GetLargestPossibleRegion());

  GHood.GoToBegin();

  //  std::cout << " foreg " << (int) foreground;
  while (!GHood.IsAtEnd())
  {
    typename TImage::PixelType p = GHood.GetCenterPixel();
    typename TImage::IndexType ind = GHood.GetIndex();
    typename TImage::IndexType ind2;
    if (p == foreground)
    {
      bool atedge = false;
      for (int i = 0; i < GHood.Size(); i++)
      {
        ind2 = GHood.GetIndex(i);
        float dist = 0.0;
        for (int j = 0; j < ImageDimension; j++)
        {
          dist += (float)(ind[j] - ind2[j]) * (float)(ind[j] - ind2[j]);
        }
        dist = sqrt(dist);
        if (GHood.GetPixel(i) != foreground && dist < 1.1f)
        {
          atedge = true;
        }
      }
      if (atedge && p == foreground)
      {
        Image->SetPixel(ind, newval);
      }
      else if (p == foreground)
      {
        Image->SetPixel(ind, 0);
      }
    }
    ++GHood;
  }

  return Image;
}

template <typename TImage, typename TMaskImage>
typename TImage::Pointer
OtsuThreshold(int NumberOfThresholds, typename TImage::Pointer input, typename TMaskImage::Pointer maskImage)
{
  // std::cout << " Otsu Thresh with " << NumberOfThresholds << " thresholds" << std::endl;

  if (maskImage.IsNull())
  {
    // Begin Threshold Image
    using InputThresholderType = itk::OtsuMultipleThresholdsImageFilter<TImage, TImage>;
    typename InputThresholderType::Pointer inputThresholder = InputThresholderType::New();

    inputThresholder->SetInput(input);
    /*
    inputThresholder->SetInsideValue(  replaceval );
    int outval=0;
    if ((float) replaceval == (float) -1) outval=1;
    inputThresholder->SetOutsideValue( outval );
    */
    inputThresholder->SetNumberOfThresholds(NumberOfThresholds);

    inputThresholder->Update();

    return inputThresholder->GetOutput();
  }
  else
  {
    using ImageType = TImage;
    using MaskImageType = TMaskImage;
    using PixelType = float;

    unsigned int numberOfBins = 200;
    int          maskLabel = 1;

    itk::ImageRegionIterator<ImageType>     ItI(input, input->GetLargestPossibleRegion());
    itk::ImageRegionIterator<MaskImageType> ItM(maskImage, maskImage->GetLargestPossibleRegion());
    PixelType                               maxValue = itk::NumericTraits<PixelType>::min();
    PixelType                               minValue = itk::NumericTraits<PixelType>::max();
    for (ItM.GoToBegin(), ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItM, ++ItI)
    {
      if (ItM.Get() == maskLabel)
      {
        if (ItI.Get() < minValue)
        {
          minValue = ItI.Get();
        }
        else if (ItI.Get() > maxValue)
        {
          maxValue = ItI.Get();
        }
      }
    }

    using StatsType = itk::LabelStatisticsImageFilter<ImageType, MaskImageType>;
    typename StatsType::Pointer stats = StatsType::New();
    stats->SetInput(input);
    stats->SetLabelInput(maskImage);
    stats->UseHistogramsOn();
    stats->SetHistogramParameters(numberOfBins, minValue, maxValue);
    stats->Update();

    using OtsuType = itk::OtsuMultipleThresholdsCalculator<typename StatsType::HistogramType>;
    typename OtsuType::Pointer otsu = OtsuType::New();
    otsu->SetInputHistogram(stats->GetHistogram(maskLabel));
    otsu->SetNumberOfThresholds(NumberOfThresholds);
    otsu->Compute();

    typename OtsuType::OutputType thresholds = otsu->GetOutput();

    typename ImageType::Pointer output = ImageType::New();
    output->CopyInformation(maskImage);
    output->SetRegions(maskImage->GetLargestPossibleRegion());
    output->AllocateInitialized();

    itk::ImageRegionIterator<ImageType> ItO(output, output->GetLargestPossibleRegion());
    for (unsigned int i = 0; i < thresholds.size(); i++)
    {

      ItI.GoToBegin();
      ItM.GoToBegin();
      ItO.GoToBegin();
      while (!ItM.IsAtEnd())
      {
        if (ItM.Get() == maskLabel)
        {
          if (itk::Math::FloatAlmostEqual(ItO.Get(), itk::NumericTraits<typename ImageType::PixelType>::ZeroValue()) &&
              ItI.Get() < static_cast<typename ImageType::PixelType>(thresholds[i]))
          {
            ItO.Set(i + 1);
          }
        }
        ++ItI;
        ++ItM;
        ++ItO;
      }
    }

    ItI.GoToBegin();
    ItM.GoToBegin();
    ItO.GoToBegin();
    while (!ItM.IsAtEnd())
    {
      if (itk::Math::FloatAlmostEqual(ItO.Get(), itk::NumericTraits<typename ImageType::PixelType>::ZeroValue()) &&
          ItM.Get() == maskLabel)
      {
        ItO.Set(thresholds.size() + 1);
      }
      ++ItI;
      ++ItM;
      ++ItO;
    }
    return output;
  }
}

template <typename TImage, typename TMaskImage>
typename TImage::Pointer
KmeansThreshold(int NumberOfThresholds, typename TImage::Pointer input, typename TMaskImage::Pointer maskImage)
{
  std::cout << " Kmeans with " << NumberOfThresholds << " thresholds" << std::endl;

  using ImageType = TImage;
  using LabelImageType = TImage;
  using MaskImageType = TMaskImage;
  using RealType = float;
  using LabelType = int;

  using MeasurementVectorType = itk::Array<RealType>;
  using SampleType = typename itk::Statistics::ListSample<MeasurementVectorType>;

  int maskLabel = 1;
  if (maskImage.IsNull())
  {
    maskImage = AllocImage<MaskImageType>(input, maskLabel);
  }

  unsigned int                     numberOfTissueClasses = NumberOfThresholds + 1;
  typename LabelImageType::Pointer output = AllocImage<LabelImageType>(input, 0);

  using StatsType = itk::LabelStatisticsImageFilter<ImageType, MaskImageType>;
  typename StatsType::Pointer stats = StatsType::New();
  stats->SetInput(input);
  stats->SetLabelInput(maskImage);
  stats->UseHistogramsOff();
  stats->Update();

  RealType minValue = stats->GetMinimum(maskLabel);
  RealType maxValue = stats->GetMaximum(maskLabel);

  // The code below can be replaced by itkListSampleToImageFilter when we
  // migrate over to the Statistics classes current in the Review/ directory.
  //
  typename SampleType::Pointer sample = SampleType::New();
  sample->SetMeasurementVectorSize(1);

  itk::ImageRegionConstIteratorWithIndex<ImageType> ItI(input, input->GetRequestedRegion());
  for (ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI)
  {
    if (!maskImage || maskImage->GetPixel(ItI.GetIndex()) == maskLabel)
    {
      typename SampleType::MeasurementVectorType measurement;
      measurement.SetSize(1);
      measurement[0] = ItI.Get();
      sample->PushBack(measurement);
    }
  }

  using TreeGeneratorType = itk::Statistics::WeightedCentroidKdTreeGenerator<SampleType>;
  typename TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
  treeGenerator->SetSample(sample);
  treeGenerator->SetBucketSize(16);
  treeGenerator->Update();

  using TreeType = typename TreeGeneratorType::KdTreeType;
  using EstimatorType = itk::Statistics::KdTreeBasedKmeansEstimator<TreeType>;
  typename EstimatorType::Pointer estimator = EstimatorType::New();
  estimator->SetKdTree(treeGenerator->GetOutput());
  estimator->SetMaximumIteration(200);
  estimator->SetCentroidPositionChangesThreshold(0.0);

  typename EstimatorType::ParametersType initialMeans(numberOfTissueClasses);

  for (unsigned int n = 0; n < numberOfTissueClasses; n++)
  {
    initialMeans[n] = static_cast<double>(minValue + (maxValue - minValue) * (static_cast<RealType>(n) + 0.5f) /
                                                       static_cast<RealType>(numberOfTissueClasses));
  }
  estimator->SetParameters(initialMeans);
  estimator->StartOptimization();

  //
  // Classify the samples
  //
  using DecisionRuleType = itk::Statistics::MinimumDecisionRule;
  typename DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();

  using ClassifierType = itk::Statistics::SampleClassifierFilter<SampleType>;
  typename ClassifierType::Pointer classifier = ClassifierType::New();
  classifier->SetDecisionRule(decisionRule);
  classifier->SetInput(sample);
  classifier->SetNumberOfClasses(numberOfTissueClasses);

  typename ClassifierType::ClassLabelVectorObjectType::Pointer classLabels =
    ClassifierType::ClassLabelVectorObjectType::New();
  classifier->SetClassLabels(classLabels);
  typename ClassifierType::ClassLabelVectorType & classLabelVector = classLabels->Get();

  //
  // Order the cluster means so that the lowest mean of the input image
  // corresponds to label '1', the second lowest to label '2', etc.
  //
  std::vector<RealType> estimatorParameters;
  for (unsigned int n = 0; n < numberOfTissueClasses; n++)
  {
    estimatorParameters.push_back(estimator->GetParameters()[n]);
  }
  std::sort(estimatorParameters.begin(), estimatorParameters.end());

  using MembershipFunctionType = itk::Statistics::DistanceToCentroidMembershipFunction<MeasurementVectorType>;
  typename ClassifierType::MembershipFunctionVectorObjectType::Pointer membershipFunctions =
    ClassifierType::MembershipFunctionVectorObjectType::New();
  typename ClassifierType::MembershipFunctionVectorType & membershipFunctionsVector = membershipFunctions->Get();

  classifier->SetMembershipFunctions(membershipFunctions);
  for (unsigned int n = 0; n < numberOfTissueClasses; n++)
  {
    typename MembershipFunctionType::Pointer membershipFunction = MembershipFunctionType::New();
    membershipFunction->SetMeasurementVectorSize(sample->GetMeasurementVectorSize());
    typename MembershipFunctionType::CentroidType centroid;
    itk::NumericTraits<typename MembershipFunctionType::CentroidType>::SetLength(centroid,
                                                                                 sample->GetMeasurementVectorSize());
    centroid[0] = estimatorParameters[n];
    membershipFunction->SetCentroid(centroid);
    membershipFunctionsVector.push_back(membershipFunction.GetPointer());

    classLabelVector.push_back(static_cast<typename ClassifierType::ClassLabelType>(n + 1));
  }
  classifier->Update();

  //
  // Classify the voxels
  //
  using ClassifierOutputType = typename ClassifierType::MembershipSampleType;
  using LabelIterator = typename ClassifierOutputType::ConstIterator;

  itk::ImageRegionIteratorWithIndex<LabelImageType> ItO(output, output->GetRequestedRegion());
  ItO.GoToBegin();
  LabelIterator it = classifier->GetOutput()->Begin();
  while (it != classifier->GetOutput()->End())
  {
    if (!maskImage || maskImage->GetPixel(ItO.GetIndex()) == maskLabel)
    {
      ItO.Set(it.GetClassLabel());
      ++it;
    }
    else
    {
      ItO.Set(itk::NumericTraits<LabelType>::ZeroValue());
    }
    ++ItO;
  }

  return output;
}

template <unsigned int InImageDimension>
int
ThresholdImage(int argc, char * argv[])
{
  //  const     unsigned int   InImageDimension = AvantsImageDimension;
  using PixelType = float;
  using FixedImageType = itk::Image<PixelType, InImageDimension>;
  using LabelType = int;
  using MaskImageType = itk::Image<LabelType, InImageDimension>;

  typename FixedImageType::Pointer fixed;
  ReadImage<FixedImageType>(fixed, argv[2]);

  typename MaskImageType::Pointer maskImage = nullptr;
  if (argc > 6)
  {
    ReadImage<MaskImageType>(maskImage, argv[6]);
  }
  // Label the surface of the image
  typename FixedImageType::Pointer thresh;
  std::string                      threshtype = std::string(argv[4]);
  if (strcmp(threshtype.c_str(), "Otsu") == 0)
  {
    thresh = OtsuThreshold<FixedImageType, MaskImageType>(std::stoi(argv[5]), fixed, maskImage);
  }
  else if (strcmp(threshtype.c_str(), "Kmeans") == 0)
  {
    thresh = KmeansThreshold<FixedImageType, MaskImageType>(std::stoi(argv[5]), fixed, maskImage);
  }
  else
  {
    PixelType insideValue = 1;
    PixelType outsideValue = 0;
    if (argc > 6)
    {
      insideValue = static_cast<PixelType>(atof(argv[6]));
    }
    if (argc > 7)
    {
      outsideValue = static_cast<PixelType>(atof(argv[7]));
    }
    thresh = BinaryThreshold_AltInsideOutside_threashold<FixedImageType>(
      atof(argv[4]), atof(argv[5]), insideValue, outsideValue, fixed);
  }

  ANTs::WriteImage<FixedImageType>(thresh, argv[3]);
  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
ThresholdImage(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "ThresholdImage");

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
    std::cout << "Usage: " << argv[0];
    std::cout << "   ImageDimension ImageIn.ext outImage.ext  threshlo threshhi <insideValue> <outsideValue>"
              << std::endl;
    std::cout << "   ImageDimension ImageIn.ext outImage.ext  Otsu NumberofThresholds <maskImage.ext>" << std::endl;
    std::cout << "   ImageDimension ImageIn.ext outImage.ext  Kmeans NumberofThresholds <maskImage.ext>" << std::endl;

    std::cout << " Inclusive thresholds " << std::endl;
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
      return ThresholdImage<2>(argc, argv);
    }
    break;
    case 3:
    {
      return ThresholdImage<3>(argc, argv);
    }
    break;
    case 4:
    {
      return ThresholdImage<4>(argc, argv);
    }
    break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
} // namespace ants
