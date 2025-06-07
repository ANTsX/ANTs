#include "antsCommandLineParser.h"
#include "antsUtilities.h"
#include "ReadWriteData.h"

#include "itkantsRegistrationHelper.h"

#include "itkImageRegionIterator.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkNumericTraits.h"

namespace ants
{

template <unsigned int ImageDimension>
int
MeasureImageSimilarity(itk::ants::CommandLineParser * parser)
{
  using TComputeType = double;

  using RegistrationHelperType = typename ants::RegistrationHelper<TComputeType, ImageDimension>;
  using ImageType = typename RegistrationHelperType::ImageType;
  using DisplacementFieldType = typename RegistrationHelperType::DisplacementFieldType;
  using DisplacementVectorType = typename DisplacementFieldType::PixelType;
  using DisplacementFieldTransformType = typename RegistrationHelperType::DisplacementFieldTransformType;

  using ImageMetricType = itk::ImageToImageMetricv4<ImageType, ImageType, ImageType, TComputeType>;
  using ImageMaskSpatialObjectType = itk::ImageMaskSpatialObject<ImageDimension>;
  using MaskImageType = typename ImageMaskSpatialObjectType::ImageType;

  using LabeledPointSetType = typename RegistrationHelperType::LabeledPointSetType;


  bool                verbose = false;
  OptionType::Pointer verboseOption = parser->GetOption("verbose");
  if (verboseOption && verboseOption->GetNumberOfFunctions())
  {
    verbose = parser->Convert<bool>(verboseOption->GetFunction(0)->GetName());
  }

  typename ImageMaskSpatialObjectType::Pointer fixedImageMask = nullptr;
  typename ImageMaskSpatialObjectType::Pointer movingImageMask = nullptr;

  typename ImageType::Pointer fixedImage = nullptr;
  typename ImageType::Pointer movingImage = nullptr;

  OptionType::Pointer maskOption = parser->GetOption("masks");
  if (maskOption && maskOption->GetNumberOfFunctions())
  {
    if (verbose)
    {
      std::cout << "  Reading mask(s)." << std::endl;
    }
    if (maskOption->GetFunction(0)->GetNumberOfParameters() > 0)
    {
      for (unsigned m = 0; m < maskOption->GetFunction(0)->GetNumberOfParameters(); m++)
      {
        std::string                     fname = maskOption->GetFunction(0)->GetParameter(m);
        typename MaskImageType::Pointer maskImage;
        ReadImage<MaskImageType>(maskImage, fname.c_str());
        if (m == 0)
        {
          if (maskImage.IsNotNull())
          {
            fixedImageMask = ImageMaskSpatialObjectType::New();
            fixedImageMask->SetImage(maskImage);
            if (verbose)
            {
              std::cout << "      Fixed mask = " << fname.c_str() << std::endl;
            }
          }
          else if (verbose)
          {
            std::cout << "      No fixed mask" << std::endl;
          }
        }
        else if (m == 1)
        {
          if (maskImage.IsNotNull())
          {
            movingImageMask = ImageMaskSpatialObjectType::New();
            movingImageMask->SetImage(maskImage);
            if (verbose)
            {
              std::cout << "      Moving mask = " << fname << std::endl;
            }
          }
          else if (verbose)
          {
            std::cout << "      No moving mask" << std::endl;
          }
        }
      }
    }
    else
    {
      std::string                     fname = maskOption->GetFunction(0)->GetName();
      typename MaskImageType::Pointer maskImage;
      ReadImage<MaskImageType>(maskImage, fname.c_str());
      if (maskImage.IsNotNull())
      {
        fixedImageMask = ImageMaskSpatialObjectType::New();
        fixedImageMask->SetImage(maskImage);
        if (verbose)
        {
          std::cout << "      Fixed mask = " << fname << std::endl;
        }
      }
      else if (verbose)
      {
        std::cout << "      No fixed mask" << std::endl;
      }
    }
  }

  OptionType::Pointer metricOption = parser->GetOption("metric");
  if (!metricOption || metricOption->GetNumberOfFunctions() == 0)
  {
    if (verbose)
    {
      std::cerr << "ERROR: the metric option ('-m') must be specified.  See help menu." << std::endl;
    }
    return EXIT_FAILURE;
  }

  typename RegistrationHelperType::Pointer regHelper = RegistrationHelperType::New();
  std::string                              whichMetric = metricOption->GetFunction(0)->GetName();
  ConvertToLowerCase(whichMetric);
  typename RegistrationHelperType::MetricEnumeration currentMetric = regHelper->StringToMetricType(whichMetric);

  // The metric weighting is irrelevant for this program.  We keep it to maintain consistency
  //   with other ANTs programs.
  float metricWeighting = 1.0;
  if (metricOption->GetFunction(0)->GetNumberOfParameters() > 2)
  {
    metricWeighting = parser->Convert<float>(metricOption->GetFunction(0)->GetParameter(2));
  }

  bool isImageMetric = true;
  if (currentMetric == RegistrationHelperType::ICP || currentMetric == RegistrationHelperType::PSE ||
      currentMetric == RegistrationHelperType::JHCT)
  {
    isImageMetric = false;
  }

  if (isImageMetric)
  {
    std::string fixedImageName = metricOption->GetFunction(0)->GetParameter(0);
    if (!ReadImage<ImageType>(fixedImage, fixedImageName.c_str()))
    {
      if (verbose)
      {
        std::cerr << "ERROR:  The fixed image file " << fixedImageName.c_str() << " does not exist." << std::endl;
      }
      return EXIT_FAILURE;
    }

    std::string movingImageName = metricOption->GetFunction(0)->GetParameter(1);
    if (!ReadImage<ImageType>(movingImage, movingImageName.c_str()))
    {
      if (verbose)
      {
        std::cerr << "ERROR:  The moving image file " << movingImageName.c_str() << " does not exist." << std::endl;
      }
      return EXIT_FAILURE;
    }

    const bool gradientfilter = false;

    TComputeType samplingPercentage = 1.0;
    if (metricOption->GetFunction(0)->GetNumberOfParameters() > 5)
    {
      samplingPercentage = parser->Convert<TComputeType>(metricOption->GetFunction(0)->GetParameter(5));
    }

    typename RegistrationHelperType::SamplingStrategy samplingStrategy = RegistrationHelperType::none;
    std::string                                       strategy = "none";
    if (metricOption->GetFunction(0)->GetNumberOfParameters() > 4)
    {
      strategy = metricOption->GetFunction(0)->GetParameter(4);
    }
    ConvertToLowerCase(strategy);

    if (strategy == "random")
    {
      samplingStrategy = RegistrationHelperType::random;
    }
    else if (strategy == "regular")
    {
      samplingStrategy = RegistrationHelperType::regular;
    }
    else if ((strategy == "none") || (strategy.empty()))
    {
      samplingStrategy = RegistrationHelperType::none;
    }

    typename ImageMetricType::Pointer imageMetric = nullptr;

    using LabeledPointSetMetricType =
      itk::LabeledPointSetToPointSetMetricv4<LabeledPointSetType, LabeledPointSetType, TComputeType>;
    typename LabeledPointSetMetricType::Pointer labeledPointSetMetric = nullptr;

    switch (currentMetric)
    {
      case RegistrationHelperType::CC:
      {
        const auto radiusOption = parser->Convert<unsigned int>(metricOption->GetFunction(0)->GetParameter(3));
        if (verbose)
        {
          std::cout << "  using the CC metric (radius = " << radiusOption << ", weight = " << metricWeighting << ")"
                    << std::endl;
        }
        using CorrelationMetricType =
          itk::ANTSNeighborhoodCorrelationImageToImageMetricv4<ImageType, ImageType, ImageType, TComputeType>;
        typename CorrelationMetricType::Pointer    correlationMetric = CorrelationMetricType::New();
        typename CorrelationMetricType::RadiusType radius;
        radius.Fill(radiusOption);
        correlationMetric->SetRadius(radius);
        correlationMetric->SetUseMovingImageGradientFilter(gradientfilter);
        correlationMetric->SetUseFixedImageGradientFilter(gradientfilter);

        imageMetric = correlationMetric;
      }
      break;
      case RegistrationHelperType::Mattes:
      {
        const auto binOption = parser->Convert<unsigned int>(metricOption->GetFunction(0)->GetParameter(3));
        if (verbose)
        {
          std::cout << "  using the Mattes MI metric (number of bins = " << binOption
                    << ", weight = " << metricWeighting << ")" << std::endl;
        }
        using MutualInformationMetricType =
          itk::MattesMutualInformationImageToImageMetricv4<ImageType, ImageType, ImageType, TComputeType>;
        typename MutualInformationMetricType::Pointer mutualInformationMetric = MutualInformationMetricType::New();
        // mutualInformationMetric = mutualInformationMetric;
        mutualInformationMetric->SetNumberOfHistogramBins(binOption);
        mutualInformationMetric->SetUseMovingImageGradientFilter(gradientfilter);
        mutualInformationMetric->SetUseFixedImageGradientFilter(gradientfilter);
        mutualInformationMetric->SetUseSampledPointSet(false);

        imageMetric = mutualInformationMetric;
      }
      break;
      case RegistrationHelperType::MI:
      {
        const auto binOption = parser->Convert<unsigned int>(metricOption->GetFunction(0)->GetParameter(3));
        if (verbose)
        {
          std::cout << "  using the joint histogram MI metric (number of bins = " << binOption
                    << ", weight = " << metricWeighting << ")" << std::endl;
        }
        using MutualInformationMetricType =
          itk::JointHistogramMutualInformationImageToImageMetricv4<ImageType, ImageType, ImageType, TComputeType>;
        typename MutualInformationMetricType::Pointer mutualInformationMetric = MutualInformationMetricType::New();
        // mutualInformationMetric = mutualInformationMetric;
        mutualInformationMetric->SetNumberOfHistogramBins(binOption);
        mutualInformationMetric->SetUseMovingImageGradientFilter(gradientfilter);
        mutualInformationMetric->SetUseFixedImageGradientFilter(gradientfilter);
        mutualInformationMetric->SetUseSampledPointSet(false);
        mutualInformationMetric->SetVarianceForJointPDFSmoothing(1.0);

        imageMetric = mutualInformationMetric;
      }
      break;
      case RegistrationHelperType::MeanSquares:
      {
        if (verbose)
        {
          std::cout << "  using the MeanSquares metric (weight = " << metricWeighting << ")" << std::endl;
        }

        using MeanSquaresMetricType =
          itk::MeanSquaresImageToImageMetricv4<ImageType, ImageType, ImageType, TComputeType>;
        typename MeanSquaresMetricType::Pointer meanSquaresMetric = MeanSquaresMetricType::New();
        // meanSquaresMetric = meanSquaresMetric;
        imageMetric = meanSquaresMetric;
      }
      break;
      case RegistrationHelperType::Demons:
      {
        if (verbose)
        {
          std::cout << "  using the Demons metric (weight = " << metricWeighting << ")" << std::endl;
        }

        using DemonsMetricType = itk::DemonsImageToImageMetricv4<ImageType, ImageType, ImageType, TComputeType>;
        typename DemonsMetricType::Pointer demonsMetric = DemonsMetricType::New();

        imageMetric = demonsMetric;
      }
      break;
      case RegistrationHelperType::GC:
      {
        if (verbose)
        {
          std::cout << "  using the global correlation metric (weight = " << metricWeighting << ")" << std::endl;
        }
        using corrMetricType = itk::CorrelationImageToImageMetricv4<ImageType, ImageType, ImageType, TComputeType>;
        typename corrMetricType::Pointer corrMetric = corrMetricType::New();

        imageMetric = corrMetric;
      }
      break;
      default:
        if (verbose)
        {
          std::cout << "ERROR: Unrecognized metric. " << std::endl;
        }
        return EXIT_FAILURE;
    }

    imageMetric->SetVirtualDomainFromImage(fixedImage);

    imageMetric->SetFixedImage(fixedImage);
    imageMetric->SetFixedImageMask(fixedImageMask);
    imageMetric->SetUseFixedImageGradientFilter(gradientfilter);

    imageMetric->SetMovingImage(movingImage);
    imageMetric->SetMovingImageMask(movingImageMask);
    imageMetric->SetUseMovingImageGradientFilter(gradientfilter);

    /** Sample the image domain **/

    if (samplingStrategy != RegistrationHelperType::none)
    {
      const typename ImageType::SpacingType oneThirdVirtualSpacing = fixedImage->GetSpacing() / 3.0;

      using MetricSamplePointSetType = typename ImageMetricType::FixedSampledPointSetType;
      typename MetricSamplePointSetType::Pointer samplePointSet = MetricSamplePointSetType::New();
      samplePointSet->Initialize();

      using SamplePointType = typename MetricSamplePointSetType::PointType;

      using RandomizerType = itk::Statistics::MersenneTwisterRandomVariateGenerator;
      typename RandomizerType::Pointer randomizer = RandomizerType::New();
      randomizer->SetSeed(1234);

      unsigned long index = 0;

      switch (samplingStrategy)
      {
        case RegistrationHelperType::regular:
        {
          const auto    sampleCount = static_cast<unsigned long>(std::ceil(1.0 / samplingPercentage));
          unsigned long count =
            sampleCount; // Start at sampleCount to keep behavior backwards identical, using first element.
          itk::ImageRegionConstIteratorWithIndex<ImageType> It(fixedImage, fixedImage->GetRequestedRegion());
          for (It.GoToBegin(); !It.IsAtEnd(); ++It)
          {
            if (count == sampleCount)
            {
              count = 0; // Reset counter
              SamplePointType point;
              fixedImage->TransformIndexToPhysicalPoint(It.GetIndex(), point);

              // randomly perturb the point within a voxel (approximately)
              for (unsigned int d = 0; d < ImageDimension; d++)
              {
                point[d] += static_cast<typename SamplePointType::CoordRepType>(randomizer->GetNormalVariate()) *
                            static_cast<typename SamplePointType::CoordRepType>(oneThirdVirtualSpacing[d]);
              }
              if (!fixedImageMask || fixedImageMask->IsInsideInWorldSpace(point))
              {
                samplePointSet->SetPoint(index, point);
                ++index;
              }
            }
            ++count;
          }
          break;
        }
        case RegistrationHelperType::random:
        {
          const unsigned long totalVirtualDomainVoxels = fixedImage->GetRequestedRegion().GetNumberOfPixels();
          const auto          sampleCount = static_cast<unsigned long>(static_cast<float>(totalVirtualDomainVoxels) *
                                                              static_cast<float>(samplingPercentage));
          itk::ImageRandomConstIteratorWithIndex<ImageType> ItR(fixedImage, fixedImage->GetRequestedRegion());
          ItR.SetNumberOfSamples(sampleCount);
          for (ItR.GoToBegin(); !ItR.IsAtEnd(); ++ItR)
          {
            SamplePointType point;
            fixedImage->TransformIndexToPhysicalPoint(ItR.GetIndex(), point);

            // randomly perturb the point within a voxel (approximately)
            for (unsigned int d = 0; d < ImageDimension; d++)
            {
              point[d] += static_cast<typename SamplePointType::CoordRepType>(randomizer->GetNormalVariate()) *
                          static_cast<typename SamplePointType::CoordRepType>(oneThirdVirtualSpacing[d]);
            }
            if (!fixedImageMask || fixedImageMask->IsInsideInWorldSpace(point))
            {
              samplePointSet->SetPoint(index, point);
              ++index;
            }
          }
          break;
        }
        case RegistrationHelperType::none:
          break;
        case RegistrationHelperType::invalid:
          break;
      }
      imageMetric->SetFixedSampledPointSet(samplePointSet);
      imageMetric->SetUseSampledPointSet(true);
    }

    imageMetric->Initialize();

    if (verbose)
    {
      imageMetric->Print(std::cout, 3);
    }

    std::cout << imageMetric->GetValue() << std::endl;

    OptionType::Pointer outputOption = parser->GetOption("output");
    if (outputOption && outputOption->GetNumberOfFunctions())
    {
      const DisplacementVectorType zeroVector(0.0);

      typename DisplacementFieldType::Pointer identityField = DisplacementFieldType::New();
      identityField->CopyInformation(fixedImage);
      identityField->SetRegions(fixedImage->GetLargestPossibleRegion());
      identityField->AllocateInitialized();

      typename DisplacementFieldTransformType::Pointer identityDisplacementFieldTransform =
        DisplacementFieldTransformType::New();
      identityDisplacementFieldTransform->SetDisplacementField(identityField);
      identityDisplacementFieldTransform->SetInverseDisplacementField(identityField);

      imageMetric->SetFixedTransform(identityDisplacementFieldTransform);
      imageMetric->SetMovingTransform(identityDisplacementFieldTransform);

      imageMetric->Initialize();

      using MetricDerivativeType = typename ImageMetricType::DerivativeType;
      const typename MetricDerivativeType::SizeValueType metricDerivativeSize =
        fixedImage->GetLargestPossibleRegion().GetNumberOfPixels() * ImageDimension;
      MetricDerivativeType metricDerivative(metricDerivativeSize);

      metricDerivative.Fill(itk::NumericTraits<typename MetricDerivativeType::ValueType>::ZeroValue());
      TComputeType value;
      imageMetric->GetValueAndDerivative(value, metricDerivative);

      typename DisplacementFieldType::Pointer gradientField = DisplacementFieldType::New();
      gradientField->CopyInformation(fixedImage);
      gradientField->SetRegions(fixedImage->GetLargestPossibleRegion());
      gradientField->AllocateInitialized();

      itk::ImageRegionIterator<DisplacementFieldType> ItG(gradientField, gradientField->GetRequestedRegion());

      itk::SizeValueType count = 0;
      for (ItG.GoToBegin(); !ItG.IsAtEnd(); ++ItG)
      {
        DisplacementVectorType displacement;
        for (itk::SizeValueType d = 0; d < ImageDimension; d++)
        {
          displacement[d] = metricDerivative[count++];
        }
        ItG.Set(displacement);
      }

      ANTs::WriteImage<DisplacementFieldType>(gradientField, (outputOption->GetFunction(0)->GetName()).c_str());
    }

    return EXIT_SUCCESS;
  }
  else
  {
    typename LabeledPointSetType::Pointer fixedLabeledPointSet = nullptr;
    typename LabeledPointSetType::Pointer movingLabeledPointSet = nullptr;

    std::string fixedPointSetFileName = metricOption->GetFunction(0)->GetParameter(0);
    if (!ReadLabeledPointSet<LabeledPointSetType>(fixedLabeledPointSet, fixedPointSetFileName.c_str(), false, 1.0))
    {
      if (verbose)
      {
        std::cerr << "ERROR:  The fixed point set file " << fixedPointSetFileName.c_str() << " does not exist."
                  << std::endl;
      }
      return EXIT_FAILURE;
    }

    std::string movingPointSetFileName = metricOption->GetFunction(0)->GetParameter(1);
    if (!ReadLabeledPointSet<LabeledPointSetType>(movingLabeledPointSet, movingPointSetFileName.c_str(), false, 1.0))
    {
      if (verbose)
      {
        std::cerr << "ERROR:  The moving point set file " << movingPointSetFileName.c_str() << " does not exist."
                  << std::endl;
      }
      return EXIT_FAILURE;
    }

    using LabeledPointSetMetricType =
      itk::LabeledPointSetToPointSetMetricv4<LabeledPointSetType, LabeledPointSetType, TComputeType>;
    typename LabeledPointSetMetricType::Pointer labeledPointSetMetric = LabeledPointSetMetricType::New();

    switch (currentMetric)
    {
      case RegistrationHelperType::ICP:
      {
        if (verbose)
        {
          std::cout << "  using the ICP metric (weight = " << metricWeighting << ")" << std::endl;
        }
        using IcpPointSetMetricType =
          itk::EuclideanDistancePointSetToPointSetMetricv4<LabeledPointSetType, LabeledPointSetType, TComputeType>;
        typename IcpPointSetMetricType::Pointer icpMetric = IcpPointSetMetricType::New();

        labeledPointSetMetric->SetPointSetMetric(icpMetric.GetPointer());
      }
      break;
      case RegistrationHelperType::PSE:
      {
        TComputeType pointSetSigma = 1.0;
        if (metricOption->GetFunction(0)->GetNumberOfParameters() > 3)
        {
          pointSetSigma = parser->Convert<TComputeType>(metricOption->GetFunction(0)->GetParameter(3));
        }
        unsigned int kNeighborhood = 50;
        if (metricOption->GetFunction(0)->GetNumberOfParameters() > 4)
        {
          kNeighborhood = parser->Convert<unsigned int>(metricOption->GetFunction(0)->GetParameter(4));
        }

        if (verbose)
        {
          std::cout << "  using the PSE metric (weight = " << metricWeighting << ")" << std::endl;
        }
        using PsePointSetMetricType =
          itk::ExpectationBasedPointSetToPointSetMetricv4<LabeledPointSetType, LabeledPointSetType, TComputeType>;
        typename PsePointSetMetricType::Pointer pseMetric = PsePointSetMetricType::New();
        pseMetric->SetPointSetSigma(pointSetSigma);
        pseMetric->SetEvaluationKNeighborhood(kNeighborhood);

        labeledPointSetMetric->SetPointSetMetric(pseMetric.GetPointer());
      }
      break;
      case RegistrationHelperType::JHCT:
      {
        TComputeType pointSetSigma = 1.0;
        if (metricOption->GetFunction(0)->GetNumberOfParameters() > 3)
        {
          pointSetSigma = parser->Convert<TComputeType>(metricOption->GetFunction(0)->GetParameter(3));
        }
        unsigned int kNeighborhood = 50;
        if (metricOption->GetFunction(0)->GetNumberOfParameters() > 4)
        {
          kNeighborhood = parser->Convert<unsigned int>(metricOption->GetFunction(0)->GetParameter(4));
        }
        TComputeType alpha = 1.1;
        if (metricOption->GetFunction(0)->GetNumberOfParameters() > 5)
        {
          alpha = parser->Convert<TComputeType>(metricOption->GetFunction(0)->GetParameter(5));
        }
        bool useAnisotropicCovariances = true;
        if (metricOption->GetFunction(0)->GetNumberOfParameters() > 6)
        {
          useAnisotropicCovariances = parser->Convert<bool>(metricOption->GetFunction(0)->GetParameter(6));
        }

        if (verbose)
        {
          std::cout << "  using the JHCT metric (weight = " << metricWeighting << ")" << std::endl;
        }
        using JhctPointSetMetricType =
          itk::JensenHavrdaCharvatTsallisPointSetToPointSetMetricv4<LabeledPointSetType, TComputeType>;
        typename JhctPointSetMetricType::Pointer jhctMetric = JhctPointSetMetricType::New();
        jhctMetric->SetPointSetSigma(pointSetSigma);
        jhctMetric->SetKernelSigma(10.0);
        jhctMetric->SetUseAnisotropicCovariances(useAnisotropicCovariances);
        jhctMetric->SetCovarianceKNeighborhood(5);
        jhctMetric->SetEvaluationKNeighborhood(kNeighborhood);
        jhctMetric->SetAlpha(alpha);

        labeledPointSetMetric->SetPointSetMetric(jhctMetric.GetPointer());
      }
      break;

      default:
        if (verbose)
        {
          std::cout << "ERROR: Unrecognized metric. " << std::endl;
        }
        return EXIT_FAILURE;
    }
    labeledPointSetMetric->SetFixedPointSet(fixedLabeledPointSet);
    labeledPointSetMetric->SetMovingPointSet(movingLabeledPointSet);

    labeledPointSetMetric->Initialize();

    if (verbose)
    {
      labeledPointSetMetric->Print(std::cout, 3);
    }
    std::cout << labeledPointSetMetric->GetValue() << std::endl;

    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE; // should not ever get here
}

void
InitializeCommandLineOptions(itk::ants::CommandLineParser * parser)
{
  {
    std::string description = std::string("Dimensionality of the fixed/moving image pair.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("dimensionality");
    option->SetShortName('d');
    option->SetUsageOption(0, "2/3/4");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("These image metrics are available--- ") +
      std::string("CC:  ANTS neighborhood cross correlation, MI:  Mutual information, ") +
      std::string("Demons: (Thirion), MeanSquares, and GC: Global Correlation. ") +
      std::string("The \"metricWeight\" variable is not used.  ") +
      std::string("The metrics can also employ a sampling strategy defined by a ") +
      std::string("sampling percentage. The sampling strategy defaults to \'None\' (aka a dense sampling of ") +
      std::string("one sample per voxel), otherwise it defines a point set over which to optimize the metric. ") +
      std::string("In addition, three point set metrics are available: Euclidean (ICP), Point-set ") +
      std::string("expectation (PSE), and Jensen-Havrda-Charvet-Tsallis (JHCT). ");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("metric");
    option->SetShortName('m');
    option->SetUsageOption(0,
                           "CC[fixedImage,movingImage,metricWeight,radius,<samplingStrategy={None,Regular,Random}>,<"
                           "samplingPercentage=[0,1]>]");
    option->SetUsageOption(1,
                           "MI[fixedImage,movingImage,metricWeight,numberOfBins,<samplingStrategy={None,Regular,Random}"
                           ">,<samplingPercentage=[0,1]>]");
    option->SetUsageOption(2,
                           "Mattes[fixedImage,movingImage,metricWeight,numberOfBins,<samplingStrategy={None,Regular,"
                           "Random}>,<samplingPercentage=[0,1]>]");
    option->SetUsageOption(3,
                           "MeanSquares[fixedImage,movingImage,metricWeight,radius=NA,<samplingStrategy={None,Regular,"
                           "Random}>,<samplingPercentage=[0,1]>]");
    option->SetUsageOption(4,
                           "Demons[fixedImage,movingImage,metricWeight,radius=NA,<samplingStrategy={None,Regular,"
                           "Random}>,<samplingPercentage=[0,1]>]");
    option->SetUsageOption(5,
                           "GC[fixedImage,movingImage,metricWeight,radius=NA,<samplingStrategy={None,Regular,Random}>,<"
                           "samplingPercentage=[0,1]>]");
    option->SetUsageOption(6, "ICP[fixedPointSet,movingPointSet,metricWeight]");
    option->SetUsageOption(7, "PSE[fixedPointSet,movingPointSet,metricWeight,<pointSetSigma=1>,<kNeighborhood=50>]");
    option->SetUsageOption(8,
                           "JHCT[fixedPointSet,movingPointSet,metricWeight,<pointSetSigma=1>,<kNeighborhood=50>,<alpha="
                           "1.1>,<useAnisotropicCovariances=1>]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("Image masks to limit voxels considered by the metric. ");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("masks");
    option->SetShortName('x');
    option->SetUsageOption(0, "fixedImageMask");
    option->SetUsageOption(1, "[fixedImageMask,movingImageMask]");
    option->SetUsageOption(2, "[fixedImageMask,NULL]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Verbose output.");

    OptionType::Pointer option = OptionType::New();
    option->SetShortName('v');
    option->SetLongName("verbose");
    option->SetUsageOption(0, "(0)/1");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Output the metric gradient image (optional).");

    OptionType::Pointer option = OptionType::New();
    option->SetShortName('o');
    option->SetLongName("output");
    option->SetUsageOption(0, "gradientImage");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Print the help menu (short version).");

    OptionType::Pointer option = OptionType::New();
    option->SetShortName('h');
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Print the help menu.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("help");
    option->SetDescription(description);
    parser->AddOption(option);
  }
}

int
MeasureImageSimilarity(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{

  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle

  args.insert(args.begin(), "MeasureImageSimilarity");

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

  ParserType::Pointer parser = ParserType::New();
  parser->SetCommand(argv[0]);

  std::string commandDescription =
    std::string("Program to calculate the ") + std::string("similarity between two images using various metrics.");

  parser->SetCommandDescription(commandDescription);
  InitializeCommandLineOptions(parser);

  if (parser->Parse(argc, argv) == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  }

  bool                verbose = false;
  OptionType::Pointer verboseOption = parser->GetOption("verbose");
  if (verboseOption && verboseOption->GetNumberOfFunctions())
  {
    verbose = parser->Convert<bool>(verboseOption->GetFunction(0)->GetName());
  }

  if (argc == 1)
  {
    parser->PrintMenu(std::cout, 5, false);
    return EXIT_FAILURE;
  }
  else if (parser->GetOption("help")->GetFunction() &&
           parser->Convert<bool>(parser->GetOption("help")->GetFunction()->GetName()))
  {
    parser->PrintMenu(std::cout, 5, false);
    return EXIT_SUCCESS;
  }
  else if (parser->GetOption('h')->GetFunction() &&
           parser->Convert<bool>(parser->GetOption('h')->GetFunction()->GetName()))
  {
    parser->PrintMenu(std::cout, 5, true);
    return EXIT_SUCCESS;
  }

  unsigned int dimension = 3;

  OptionType::Pointer dimOption = parser->GetOption("dimensionality");
  if (dimOption && dimOption->GetNumberOfFunctions())
  {
    dimension = parser->Convert<unsigned int>(dimOption->GetFunction(0)->GetName());
  }
  else
  {
    if (verbose)
    {
      std::cerr << "Image dimensionality not specified.  See command line option --dimensionality" << std::endl;
    }
    return EXIT_FAILURE;
  }

  switch (dimension)
  {
    case 2:
    {
      return MeasureImageSimilarity<2>(parser);
      break;
    }
    case 3:
    {
      return MeasureImageSimilarity<3>(parser);
      break;
    }
    case 4:
    {
      return MeasureImageSimilarity<4>(parser);
      break;
    }
    default:
    {
      std::cerr << "Unrecognized dimensionality.  Please see help menu." << std::endl;
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}
} // namespace ants
