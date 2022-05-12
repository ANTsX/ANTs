#include "antsUtilities.h"
#include "antsAllocImage.h"
#include "antsCommandLineParser.h"

#include "ReadWriteData.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkExpImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkN3MRIBiasFieldCorrectionImageFilter.h"
#include "itkN3BiasFieldCorrectionImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkTimeProbe.h"
#include "itkShrinkImageFilter.h"

#include <string>
#include <algorithm>
#include <vector>

#include "ANTsVersion.h"

namespace ants
{
template <typename TFilter>
class CommandIterationUpdate final : public itk::Command
{
public:
  using Self = CommandIterationUpdate<TFilter>;
  using Superclass = itk::Command;
  using Pointer = itk::SmartPointer<Self>;
  itkNewMacro(Self);

protected:
  CommandIterationUpdate() = default;

public:
  void
  Execute(itk::Object * caller, const itk::EventObject & event) override
  {
    Execute((const itk::Object *)caller, event);
  }

  void
  Execute(const itk::Object * object, const itk::EventObject & event) override
  {
    const auto * filter = dynamic_cast<const TFilter *>(object);

    if (typeid(event) != typeid(itk::IterationEvent))
    {
      return;
    }

    std::cout << "Iteration " << filter->GetElapsedIterations() << " (of " << filter->GetMaximumNumberOfIterations()
              << ").  ";
    std::cout << " Current convergence value = " << filter->GetCurrentConvergenceMeasurement()
              << " (threshold = " << filter->GetConvergenceThreshold() << ")" << std::endl;
  }
};

template <unsigned int ImageDimension>
int
N3BiasFieldCorrection(int argc, char * argv[])
{
  using RealType = float;

  using ImageType = itk::Image<RealType, ImageDimension>;
  using MaskImageType = itk::Image<unsigned char, ImageDimension>;
  typename ImageType::Pointer image;
  ReadImage<ImageType>(image, argv[2]);

  using ShrinkerType = itk::ShrinkImageFilter<ImageType, ImageType>;
  typename ShrinkerType::Pointer shrinker = ShrinkerType::New();
  shrinker->SetInput(image);
  shrinker->SetShrinkFactors(1);

  typename MaskImageType::Pointer maskImage = nullptr;

  if (argc > 5)
  {
    ReadImage<MaskImageType>(maskImage, argv[5]);
  }
  if (!maskImage)
  {
    using ThresholderType = itk::OtsuThresholdImageFilter<ImageType, MaskImageType>;
    typename ThresholderType::Pointer otsu = ThresholderType::New();
    otsu->SetInput(image);
    otsu->SetNumberOfHistogramBins(200);
    otsu->SetInsideValue(0);
    otsu->SetOutsideValue(1);
    otsu->Update();

    maskImage = otsu->GetOutput();
  }
  using MaskShrinkerType = itk::ShrinkImageFilter<MaskImageType, MaskImageType>;
  typename MaskShrinkerType::Pointer maskshrinker = MaskShrinkerType::New();
  maskshrinker->SetInput(maskImage);
  maskshrinker->SetShrinkFactors(1);

  // 4 is [shrinkFactor]
  // 5 is [maskImage]
  if (argc > 4)
  {
    shrinker->SetShrinkFactors(std::stoi(argv[4]));
    maskshrinker->SetShrinkFactors(std::stoi(argv[4]));
  }
  shrinker->Update();
  maskshrinker->Update();

  using CorrecterType = itk::N3MRIBiasFieldCorrectionImageFilter<ImageType, MaskImageType, ImageType>;
  typename CorrecterType::Pointer correcter = CorrecterType::New();
  correcter->SetInput(shrinker->GetOutput());
  correcter->SetMaskImage(maskshrinker->GetOutput());

  // 6 is [numberOfIterations]
  if (argc > 6)
  {
    correcter->SetMaximumNumberOfIterations(std::stoi(argv[6]));
  }
  // 7 is [numberOfFittingLevels]
  if (argc > 7)
  {
    correcter->SetNumberOfFittingLevels(std::stoi(argv[7]));
  }

  // 8 is [outputBiasField]
  // 9 is [verbose]
  bool verbose = false;
  if (argc > 9)
  {
    verbose = static_cast<bool>(std::stoi(argv[9]));
  }

  using CommandType = CommandIterationUpdate<CorrecterType>;
  typename CommandType::Pointer observer = CommandType::New();
  if (verbose)
  {
    correcter->AddObserver(itk::IterationEvent(), observer);
  }
  try
  {
    correcter->Update();
  }
  catch (...)
  {
    std::cout << "Exception caught." << std::endl;
    return EXIT_FAILURE;
  }

  //  correcter->Print( std::cout, 3 );

  /**
   * Reconstruct the bias field at full image resolution.  Divide
   * the original input image by the bias field to get the final
   * corrected image.
   */
  using BSplinerType = itk::BSplineControlPointImageFilter<typename CorrecterType::BiasFieldControlPointLatticeType,
                                                           typename CorrecterType::ScalarImageType>;
  typename BSplinerType::Pointer bspliner = BSplinerType::New();
  bspliner->SetInput(correcter->GetLogBiasFieldControlPointLattice());
  bspliner->SetSplineOrder(correcter->GetSplineOrder());
  bspliner->SetSize(image->GetLargestPossibleRegion().GetSize());
  bspliner->SetOrigin(image->GetOrigin());
  bspliner->SetDirection(image->GetDirection());
  bspliner->SetSpacing(image->GetSpacing());
  bspliner->Update();

  typename ImageType::Pointer logField = AllocImage<ImageType>(bspliner->GetOutput());

  itk::ImageRegionIterator<typename CorrecterType::ScalarImageType> ItB(
    bspliner->GetOutput(), bspliner->GetOutput()->GetLargestPossibleRegion());
  itk::ImageRegionIterator<ImageType> ItF(logField, logField->GetLargestPossibleRegion());
  for (ItB.GoToBegin(), ItF.GoToBegin(); !ItB.IsAtEnd(); ++ItB, ++ItF)
  {
    ItF.Set(ItB.Get()[0]);
  }

  using ExpFilterType = itk::ExpImageFilter<ImageType, ImageType>;
  typename ExpFilterType::Pointer expFilter = ExpFilterType::New();
  expFilter->SetInput(logField);
  expFilter->Update();

  using DividerType = itk::DivideImageFilter<ImageType, ImageType, ImageType>;
  typename DividerType::Pointer divider = DividerType::New();
  divider->SetInput1(image);
  divider->SetInput2(expFilter->GetOutput());
  divider->Update();

  using WriterType = itk::ImageFileWriter<ImageType>;
  typename WriterType::Pointer writer = WriterType::New();
  if (argc < 4)
  {
    std::cout << "missing divider image filename" << std::endl;
    throw;
  }
  ANTs::WriteImage<ImageType>(divider->GetOutput(), argv[3]);

  if (argc > 8)
  {
    ANTs::WriteImage<ImageType>(expFilter->GetOutput(), argv[8]);
  }

  return EXIT_SUCCESS;
}

template <unsigned int ImageDimension>
int
N3(itk::ants::CommandLineParser * parser)
{
  using RealType = float;

  using ImageType = itk::Image<RealType, ImageDimension>;
  typename ImageType::Pointer inputImage = nullptr;

  using MaskImageType = itk::Image<RealType, ImageDimension>;
  typename MaskImageType::Pointer maskImage = nullptr;

  bool                                                       verbose = false;
  typename itk::ants::CommandLineParser::OptionType::Pointer verboseOption = parser->GetOption("verbose");
  if (verboseOption && verboseOption->GetNumberOfFunctions())
  {
    verbose = parser->Convert<bool>(verboseOption->GetFunction(0)->GetName());
  }
  if (verbose)
  {
    std::cout << std::endl << "Running N3 for " << ImageDimension << "-dimensional images." << std::endl << std::endl;
  }

  using CorrecterType = itk::N3BiasFieldCorrectionImageFilter<ImageType, MaskImageType, ImageType>;
  typename CorrecterType::Pointer correcter = CorrecterType::New();

  typename itk::ants::CommandLineParser::OptionType::Pointer inputImageOption = parser->GetOption("input-image");
  if (inputImageOption && inputImageOption->GetNumberOfFunctions())
  {
    std::string inputFile = inputImageOption->GetFunction(0)->GetName();
    ReadImage<ImageType>(inputImage, inputFile.c_str());
  }
  else
  {
    if (verbose)
    {
      std::cerr << "Input image not specified." << std::endl;
    }
    return EXIT_FAILURE;
  }

  /**
   * handle the mask image
   */

  bool isMaskImageSpecified = false;

  typename itk::ants::CommandLineParser::OptionType::Pointer maskImageOption = parser->GetOption("mask-image");
  if (maskImageOption && maskImageOption->GetNumberOfFunctions())
  {
    std::string inputMaskFile = maskImageOption->GetFunction(0)->GetName();
    ReadImage<MaskImageType>(maskImage, inputMaskFile.c_str());

    isMaskImageSpecified = true;
  }
  if (maskImage.IsNull())
  {
    if (verbose)
    {
      std::cout << "Mask not read.  Using the entire image as the mask." << std::endl << std::endl;
    }
    maskImage = MaskImageType::New();
    maskImage->CopyInformation(inputImage);
    maskImage->SetRegions(inputImage->GetRequestedRegion());
    maskImage->Allocate();
    maskImage->FillBuffer(itk::NumericTraits<typename MaskImageType::PixelType>::OneValue());
  }

  /**
   * check for negative values in the masked region
   */
  using ThresholderType = itk::BinaryThresholdImageFilter<MaskImageType, MaskImageType>;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInsideValue(itk::NumericTraits<typename MaskImageType::PixelType>::ZeroValue());
  thresholder->SetOutsideValue(itk::NumericTraits<typename MaskImageType::PixelType>::OneValue());
  thresholder->SetLowerThreshold(itk::NumericTraits<typename MaskImageType::PixelType>::ZeroValue());
  thresholder->SetUpperThreshold(itk::NumericTraits<typename MaskImageType::PixelType>::ZeroValue());
  thresholder->SetInput(maskImage);

  using StatsType = itk::LabelStatisticsImageFilter<ImageType, MaskImageType>;
  typename StatsType::Pointer statsOriginal = StatsType::New();
  statsOriginal->SetInput(inputImage);
  statsOriginal->SetLabelInput(thresholder->GetOutput());
  statsOriginal->UseHistogramsOff();
  statsOriginal->Update();

  using StatsLabelType = typename StatsType::LabelPixelType;
  StatsLabelType maskLabel = itk::NumericTraits<StatsLabelType>::OneValue();

  RealType minOriginal = statsOriginal->GetMinimum(maskLabel);
  RealType maxOriginal = statsOriginal->GetMaximum(maskLabel);

  if (verbose)
  {
    std::cout << "Original intensity range:  [" << minOriginal << ", " << maxOriginal << "]" << std::endl;
  }

  if (minOriginal <= 0)
  {
    if (verbose)
    {
      std::cout << std::endl;
      std::cout << "***********************************************************" << std::endl;
      std::cout << "Warning:  Your input image contains nonpositive values" << std::endl;
      std::cout << "which could cause failure or problematic results.  A" << std::endl;
      std::cout << "possible workaround would be to:" << std::endl;
      std::cout << "   1. rescale your image to positive values e.g., [10,100]." << std::endl;
      std::cout << "   2. run N3 on your rescaled image." << std::endl;
      std::cout << "   3. (optional) rescale the N3 output to the original" << std::endl;
      std::cout << "      intensity range." << std::endl;
      std::cout << "***********************************************************" << std::endl;
      std::cout << std::endl;
    }
  }

  /**
   * handle the weight image
   */

  typename ImageType::Pointer weightImage = nullptr;

  typename itk::ants::CommandLineParser::OptionType::Pointer weightImageOption = parser->GetOption("weight-image");
  if (weightImageOption && weightImageOption->GetNumberOfFunctions())
  {
    std::string inputFile = weightImageOption->GetFunction(0)->GetName();
    ReadImage<ImageType>(weightImage, inputFile.c_str());
  }

  /**
   * convergence options
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer convergenceOption = parser->GetOption("convergence");
  if (convergenceOption && convergenceOption->GetNumberOfFunctions())
  {
    if (convergenceOption->GetFunction(0)->GetNumberOfParameters() > 0)
    {
      unsigned int maximumNumberOfIterations =
        parser->Convert<unsigned int>(convergenceOption->GetFunction(0)->GetParameter(0));
      correcter->SetMaximumNumberOfIterations(maximumNumberOfIterations);
      correcter->SetConvergenceThreshold(0.0);
    }
    if (convergenceOption->GetFunction(0)->GetNumberOfParameters() > 1)
    {
      correcter->SetConvergenceThreshold(parser->Convert<float>(convergenceOption->GetFunction(0)->GetParameter(1)));
    }
  }
  else // set default values
  {
    correcter->SetMaximumNumberOfIterations(50);
    correcter->SetConvergenceThreshold(0.0);
  }

  /**
   * B-spline options -- we place this here to take care of the case where
   * the user wants to specify things in terms of the spline distance.
   */

  typename ImageType::IndexType inputImageIndex = inputImage->GetLargestPossibleRegion().GetIndex();
  typename ImageType::SizeType  inputImageSize = inputImage->GetLargestPossibleRegion().GetSize();

  typename ImageType::PointType newOrigin = inputImage->GetOrigin();

  typename itk::ants::CommandLineParser::OptionType::Pointer bsplineOption = parser->GetOption("bspline-fitting");
  if (bsplineOption && bsplineOption->GetNumberOfFunctions())
  {
    if (bsplineOption->GetFunction(0)->GetNumberOfParameters() > 2)
    {
      correcter->SetSplineOrder(parser->Convert<unsigned int>(bsplineOption->GetFunction(0)->GetParameter(2)));
    }
    if (bsplineOption->GetFunction(0)->GetNumberOfParameters() > 1)
    {
      correcter->SetNumberOfFittingLevels(
        parser->Convert<unsigned int>(bsplineOption->GetFunction(0)->GetParameter(1)));
    }
    if (bsplineOption->GetFunction(0)->GetNumberOfParameters() > 0)
    {
      std::vector<float> array = parser->ConvertVector<float>(bsplineOption->GetFunction(0)->GetParameter(0));
      typename CorrecterType::ArrayType numberOfControlPoints;
      if (array.size() == 1)
      {
        // the user wants to specify things in terms of spline distance.
        //  1. need to pad the images to get as close to possible to the
        //     requested domain size.
        float splineDistance = array[0];

        typename ImageType::SizeType originalImageSize = inputImage->GetLargestPossibleRegion().GetSize();

        itk::Size<ImageDimension> lowerBound;
        itk::Size<ImageDimension> upperBound;
        for (unsigned int d = 0; d < ImageDimension; d++)
        {
          float domain = static_cast<float>(originalImageSize[d] - 1) * static_cast<float>(inputImage->GetSpacing()[d]);
          auto  numberOfSpans = static_cast<unsigned int>(std::ceil(domain / splineDistance));
          auto  extraPadding = static_cast<unsigned long>((numberOfSpans * splineDistance - domain) /
                                                           static_cast<float>(inputImage->GetSpacing()[d]) +
                                                         static_cast<float>(0.5));
          lowerBound[d] = static_cast<unsigned long>(0.5 * extraPadding);
          upperBound[d] = extraPadding - lowerBound[d];
          newOrigin[d] -= (static_cast<double>(lowerBound[d]) * inputImage->GetSpacing()[d]);
          numberOfControlPoints[d] = numberOfSpans + correcter->GetSplineOrder();
        }

        using PadderType = itk::ConstantPadImageFilter<ImageType, ImageType>;
        typename PadderType::Pointer padder = PadderType::New();
        padder->SetInput(inputImage);
        padder->SetPadLowerBound(lowerBound);
        padder->SetPadUpperBound(upperBound);
        padder->SetConstant(itk::NumericTraits<typename ImageType::PixelType>::ZeroValue());
        padder->Update();

        inputImage = padder->GetOutput();
        inputImage->DisconnectPipeline();

        using MaskPadderType = itk::ConstantPadImageFilter<MaskImageType, MaskImageType>;
        typename MaskPadderType::Pointer maskPadder = MaskPadderType::New();
        maskPadder->SetInput(maskImage);
        maskPadder->SetPadLowerBound(lowerBound);
        maskPadder->SetPadUpperBound(upperBound);
        maskPadder->SetConstant(0);
        maskPadder->Update();

        maskImage = maskPadder->GetOutput();
        maskImage->DisconnectPipeline();

        if (weightImage)
        {
          typename PadderType::Pointer weightPadder = PadderType::New();
          weightPadder->SetInput(weightImage);
          weightPadder->SetPadLowerBound(lowerBound);
          weightPadder->SetPadUpperBound(upperBound);
          weightPadder->SetConstant(0);
          weightPadder->Update();

          weightImage = weightPadder->GetOutput();
          weightImage->DisconnectPipeline();
        }

        if (verbose)
        {
          std::cout << "Specified spline distance: " << splineDistance << "mm" << std::endl;
          std::cout << "  original image size:  " << originalImageSize << std::endl;
          std::cout << "  padded image size:  " << inputImage->GetLargestPossibleRegion().GetSize() << std::endl;
          std::cout << "  number of control points:  " << numberOfControlPoints << std::endl;
          std::cout << std::endl;
        }
      }
      else if (array.size() == ImageDimension)
      {
        for (unsigned int d = 0; d < ImageDimension; d++)
        {
          numberOfControlPoints[d] = static_cast<unsigned int>(array[d]) + correcter->GetSplineOrder();
        }
      }
      else
      {
        if (verbose)
        {
          std::cerr << "Incorrect mesh resolution" << std::endl;
        }
        return EXIT_FAILURE;
      }
      correcter->SetNumberOfControlPoints(numberOfControlPoints);
    }
  }

  using ShrinkerType = itk::ShrinkImageFilter<ImageType, ImageType>;
  typename ShrinkerType::Pointer shrinker = ShrinkerType::New();
  shrinker->SetInput(inputImage);
  shrinker->SetShrinkFactors(1);

  using MaskShrinkerType = itk::ShrinkImageFilter<MaskImageType, MaskImageType>;
  typename MaskShrinkerType::Pointer maskshrinker = MaskShrinkerType::New();
  maskshrinker->SetInput(maskImage);
  maskshrinker->SetShrinkFactors(1);

  typename itk::ants::CommandLineParser::OptionType::Pointer shrinkFactorOption = parser->GetOption("shrink-factor");
  int                                                        shrinkFactor = 4;
  if (shrinkFactorOption && shrinkFactorOption->GetNumberOfFunctions())
  {
    shrinkFactor = parser->Convert<int>(shrinkFactorOption->GetFunction(0)->GetName());
  }
  shrinker->SetShrinkFactors(shrinkFactor);
  maskshrinker->SetShrinkFactors(shrinkFactor);
  if (ImageDimension == 4)
  {
    shrinker->SetShrinkFactor(3, 1);
    maskshrinker->SetShrinkFactor(3, 1);
  }
  shrinker->Update();
  maskshrinker->Update();

  itk::TimeProbe timer;
  timer.Start();

  correcter->SetInput(shrinker->GetOutput());
  correcter->SetMaskImage(maskshrinker->GetOutput());

  using WeightShrinkerType = itk::ShrinkImageFilter<ImageType, ImageType>;
  typename WeightShrinkerType::Pointer weightshrinker = WeightShrinkerType::New();
  if (weightImage)
  {
    weightshrinker->SetInput(weightImage);
    weightshrinker->SetShrinkFactors(shrinker->GetShrinkFactors());
    weightshrinker->Update();

    correcter->SetConfidenceImage(weightshrinker->GetOutput());
  }

  if (verbose)
  {
    using CommandType = CommandIterationUpdate<CorrecterType>;
    typename CommandType::Pointer observer = CommandType::New();
    correcter->AddObserver(itk::IterationEvent(), observer);
  }

  /**
   * histogram sharpening options
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer histOption = parser->GetOption("histogram-sharpening");
  if (histOption && histOption->GetNumberOfFunctions())
  {
    if (histOption->GetFunction(0)->GetNumberOfParameters() > 0)
    {
      correcter->SetBiasFieldFullWidthAtHalfMaximum(
        parser->Convert<float>(histOption->GetFunction(0)->GetParameter(0)));
    }
    if (histOption->GetFunction(0)->GetNumberOfParameters() > 1)
    {
      correcter->SetWienerFilterNoise(parser->Convert<float>(histOption->GetFunction(0)->GetParameter(1)));
    }
    if (histOption->GetFunction(0)->GetNumberOfParameters() > 2)
    {
      correcter->SetNumberOfHistogramBins(parser->Convert<unsigned int>(histOption->GetFunction(0)->GetParameter(2)));
    }
  }

  try
  {
    // correcter->DebugOn();
    correcter->Update();
  }
  catch (const itk::ExceptionObject & e)
  {
    if (verbose)
    {
      std::cerr << "Exception caught: " << e << std::endl;
    }
    return EXIT_FAILURE;
  }

  if (verbose)
  {
    correcter->Print(std::cout, 3);
  }

  timer.Stop();
  if (verbose)
  {
    std::cout << "Elapsed time: " << timer.GetMean() << std::endl;
  }

  /**
   * output
   */

  typename itk::ants::CommandLineParser::OptionType::Pointer outputOption = parser->GetOption("output");
  if (outputOption && outputOption->GetNumberOfFunctions())
  {
    /**
     * Reconstruct the bias field at full image resolution.  Divide
     * the original input image by the bias field to get the final
     * corrected image.
     */
    using BSplinerType = itk::BSplineControlPointImageFilter<typename CorrecterType::BiasFieldControlPointLatticeType,
                                                             typename CorrecterType::ScalarImageType>;
    typename BSplinerType::Pointer bspliner = BSplinerType::New();
    bspliner->SetInput(correcter->GetLogBiasFieldControlPointLattice());
    bspliner->SetSplineOrder(correcter->GetSplineOrder());
    bspliner->SetSize(inputImage->GetLargestPossibleRegion().GetSize());
    bspliner->SetOrigin(newOrigin);
    bspliner->SetDirection(inputImage->GetDirection());
    bspliner->SetSpacing(inputImage->GetSpacing());
    bspliner->Update();

    typename ImageType::Pointer logField = AllocImage<ImageType>(inputImage);

    itk::ImageRegionIterator<typename CorrecterType::ScalarImageType> ItB(
      bspliner->GetOutput(), bspliner->GetOutput()->GetLargestPossibleRegion());
    itk::ImageRegionIterator<ImageType> ItF(logField, logField->GetLargestPossibleRegion());
    for (ItB.GoToBegin(), ItF.GoToBegin(); !ItB.IsAtEnd(); ++ItB, ++ItF)
    {
      ItF.Set(ItB.Get()[0]);
    }

    using ExpFilterType = itk::ExpImageFilter<ImageType, ImageType>;
    typename ExpFilterType::Pointer expFilter = ExpFilterType::New();
    expFilter->SetInput(logField);
    expFilter->Update();

    using DividerType = itk::DivideImageFilter<ImageType, ImageType, ImageType>;
    typename DividerType::Pointer divider = DividerType::New();
    divider->SetInput1(inputImage);
    divider->SetInput2(expFilter->GetOutput());

    typename ImageType::Pointer dividedImage = divider->GetOutput();
    dividedImage->Update();
    dividedImage->DisconnectPipeline();

    if (maskImageOption && maskImageOption->GetNumberOfFunctions() > 0)
    {
      itk::ImageRegionIteratorWithIndex<ImageType> ItD(dividedImage, dividedImage->GetLargestPossibleRegion());
      itk::ImageRegionIterator<ImageType>          ItI(inputImage, inputImage->GetLargestPossibleRegion());
      for (ItD.GoToBegin(), ItI.GoToBegin(); !ItD.IsAtEnd(); ++ItD, ++ItI)
      {
        if (itk::Math::FloatAlmostEqual(maskImage->GetPixel(ItD.GetIndex()),
                                        itk::NumericTraits<typename MaskImageType::PixelType>::ZeroValue()))
        {
          ItD.Set(ItI.Get());
        }
      }
    }

    bool doRescale = true;

    typename itk::ants::CommandLineParser::OptionType::Pointer rescaleOption = parser->GetOption("rescale-intensities");
    if (!isMaskImageSpecified || (rescaleOption && rescaleOption->GetNumberOfFunctions() &&
                                  !parser->Convert<bool>(rescaleOption->GetFunction()->GetName())))
    {
      doRescale = false;
    }

    if (doRescale)
    {
      typename ThresholderType::Pointer thresholder2 = ThresholderType::New();
      thresholder2->SetInsideValue(itk::NumericTraits<typename MaskImageType::PixelType>::ZeroValue());
      thresholder2->SetOutsideValue(itk::NumericTraits<typename MaskImageType::PixelType>::OneValue());
      thresholder2->SetLowerThreshold(itk::NumericTraits<typename MaskImageType::PixelType>::ZeroValue());
      thresholder2->SetUpperThreshold(itk::NumericTraits<typename MaskImageType::PixelType>::ZeroValue());
      thresholder2->SetInput(maskImage);

      typename StatsType::Pointer statsBiasCorrected = StatsType::New();
      statsBiasCorrected->SetInput(dividedImage);
      statsBiasCorrected->SetLabelInput(thresholder2->GetOutput());
      statsBiasCorrected->UseHistogramsOff();
      statsBiasCorrected->Update();

      RealType minBiasCorrected = statsBiasCorrected->GetMinimum(maskLabel);
      RealType maxBiasCorrected = statsBiasCorrected->GetMaximum(maskLabel);

      RealType slope = (maxOriginal - minOriginal) / (maxBiasCorrected - minBiasCorrected);

      itk::ImageRegionIteratorWithIndex<ImageType> ItD(dividedImage, dividedImage->GetLargestPossibleRegion());
      for (ItD.GoToBegin(); !ItD.IsAtEnd(); ++ItD)
      {
        if (itk::Math::FloatAlmostEqual(maskImage->GetPixel(ItD.GetIndex()), static_cast<RealType>(maskLabel)))
        {
          RealType originalIntensity = ItD.Get();
          RealType rescaledIntensity = maxOriginal - slope * (maxBiasCorrected - originalIntensity);
          ItD.Set(rescaledIntensity);
        }
      }
    }

    typename ImageType::RegionType inputRegion;
    inputRegion.SetIndex(inputImageIndex);
    inputRegion.SetSize(inputImageSize);

    using CropperType = itk::ExtractImageFilter<ImageType, ImageType>;
    typename CropperType::Pointer cropper = CropperType::New();
    cropper->SetInput(dividedImage);
    cropper->SetExtractionRegion(inputRegion);
    cropper->SetDirectionCollapseToSubmatrix();
    cropper->Update();

    typename CropperType::Pointer biasFieldCropper = CropperType::New();
    biasFieldCropper->SetInput(expFilter->GetOutput());
    biasFieldCropper->SetExtractionRegion(inputRegion);
    biasFieldCropper->SetDirectionCollapseToSubmatrix();
    biasFieldCropper->Update();

    if (outputOption->GetFunction(0)->GetNumberOfParameters() == 0)
    {
      ANTs::WriteImage<ImageType>(cropper->GetOutput(), (outputOption->GetFunction(0)->GetName()).c_str());
    }
    else if (outputOption->GetFunction(0)->GetNumberOfParameters() > 0)
    {
      ANTs::WriteImage<ImageType>(cropper->GetOutput(), (outputOption->GetFunction(0)->GetParameter(0)).c_str());
      if (outputOption->GetFunction(0)->GetNumberOfParameters() > 1)
      {
        ANTs::WriteImage<ImageType>(biasFieldCropper->GetOutput(),
                                    (outputOption->GetFunction(0)->GetParameter(1)).c_str());
      }
    }
  }

  return EXIT_SUCCESS;
}

void
N3InitializeCommandLineOptions(itk::ants::CommandLineParser * parser)
{
  using OptionType = itk::ants::CommandLineParser::OptionType;

  {
    std::string description = std::string("This option forces the image to be treated as a specified-") +
                              std::string("dimensional image.  If not specified, N3 tries to ") +
                              std::string("infer the dimensionality from the input image.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("image-dimensionality");
    option->SetShortName('d');
    option->SetUsageOption(0, "2/3/4");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("A scalar image is expected as input for bias correction.  ") +
                              std::string("Since N3 log transforms the intensities, negative values ") +
                              std::string("or values close to zero should be processed prior to ") +
                              std::string("correction.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("input-image");
    option->SetShortName('i');
    option->SetUsageOption(0, "inputImageFilename");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("If a mask image is specified, the final bias correction is ") +
      std::string("only performed in the mask region.  If a weight image is not ") +
      std::string("specified, only intensity values inside the masked region are ") +
      std::string("used during the execution of the algorithm.  If a weight ") +
      std::string("image is specified, only the non-zero weights are used in the ") +
      std::string("execution of the algorithm although the mask region defines ") +
      std::string("where bias correction is performed in the final output. ") +
      std::string("Otherwise bias correction occurs over the entire image domain. ") +
      std::string("See also the option description for the weight image. ") +
      std::string("If a mask image is *not* specified then the entire image region ") +
      std::string("will be used as the mask region.  Note that this is different than ") +
      std::string("the N3 implementation which uses the results of Otsu thresholding ") +
      std::string("to define a mask.  However, this leads to unknown anatomical regions being ") +
      std::string("included and excluded during the bias correction.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("mask-image");
    option->SetShortName('x');
    option->SetUsageOption(0, "maskImageFilename");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("At each iteration, a new intensity mapping is calculated ") +
                              std::string("and applied but there is nothing which constrains the ") +
                              std::string("new intensity range to be within certain values.  The ") +
                              std::string("result is that the range can \"drift\" from the original ") +
                              std::string("at each iteration.  This option rescales to the [min,max] ") +
                              std::string("range of the original image intensities within the user-specified mask.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("rescale-intensities");
    option->SetShortName('r');
    option->SetUsageOption(0, "0/(1)");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("The weight image allows the user to perform a relative ") +
                              std::string("weighting of specific voxels during the B-spline fitting. ") +
                              std::string("For example, some studies have shown that N3 performed on ") +
                              std::string("white matter segmentations improves performance.  If one ") +
                              std::string("has a spatial probability map of the white matter, one can ") +
                              std::string("use this map to weight the b-spline fitting towards those ") +
                              std::string("voxels which are more probabilistically classified as white ") +
                              std::string("matter.  See also the option description for the mask image.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("weight-image");
    option->SetUsageOption(0, "weightImageFilename");
    option->SetShortName('w');
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Running N3 on large images can be time consuming. ") +
                              std::string("To lessen computation time, the input image can be resampled. ") +
                              std::string("The shrink factor, specified as a single integer, describes ") +
                              std::string("this resampling.  Shrink factors <= 4 are commonly used.") +
                              std::string("Note that the shrink factor is only applied to the first two or ") +
                              std::string("three dimensions which we assume are spatial.  ");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("shrink-factor");
    option->SetShortName('s');
    option->SetUsageOption(0, "1/2/3/(4)/...");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Convergence is determined by calculating the sigma of the coefficient of ") +
                              std::string("variation between subsequent iterations. When this value ") +
                              std::string("is less than the specified threshold ") +
                              std::string("from the previous iteration or the maximum number of ") +
                              std::string("iterations is exceeded the program terminates.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("convergence");
    option->SetShortName('c');
    option->SetUsageOption(0, "[<numberOfIterations=50>,<convergenceThreshold=0.0>]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("These options describe the b-spline fitting parameters. ") + std::string("The b-spline mesh is ") +
      std::string("specified either as the number of elements in each dimension, ") +
      std::string("e.g. 2x2x3 for 3-D images, or it can be specified as a ") +
      std::string("single scalar parameter which describes the isotropic sizing ") +
      std::string("of the mesh elements. The latter option is typically preferred. ") +
      std::string("The default setting ") +
      std::string("is to employ a single mesh element over the entire domain, i.e., ") + std::string("-b [1x1x1,4,3].");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("bspline-fitting");
    option->SetShortName('b');
    option->SetUsageOption(0, "[splineDistance,<numberOfFittingLevels=4>,<splineOrder=3>]");
    option->SetUsageOption(1, "[meshResolution,<numberOfFittingLevels=4>,<splineOrder=3>]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("These options describe the histogram sharpening parameters, ") +
                              std::string("i.e. the deconvolution step parameters described in the ") +
                              std::string("original N3 algorithm.  The default values have been shown ") +
                              std::string("to work fairly well.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("histogram-sharpening");
    option->SetShortName('t');
    option->SetUsageOption(0, "[<FWHM=0.15>,<wienerNoise=0.01>,<numberOfHistogramBins=200>]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("The output consists of the bias corrected version of the ") +
                              std::string("input image.  Optionally, one can also output the estimated ") +
                              std::string("bias field.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("output");
    option->SetShortName('o');
    option->SetUsageOption(0, "correctedImage");
    option->SetUsageOption(1, "[correctedImage,<biasField>]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("Get Version Information.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("version");
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

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
N3BiasFieldCorrection(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{

  if (args.size() >= 1 && std::string(args[0]) != std::string("2") && std::string(args[0]) != std::string("3") &&
      std::string(args[0]) != std::string("4"))
  {
    // put the arguments coming in as 'args' into standard (argc,argv) format;
    // 'args' doesn't have the command name as first, argument, so add it manually;
    // 'args' may have adjacent arguments concatenated into one argument,
    // which the parser should handle
    args.insert(args.begin(), "N3BiasFieldCorrection");

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

    itk::ants::CommandLineParser::Pointer parser = itk::ants::CommandLineParser::New();

    parser->SetCommand(argv[0]);

    std::string commandDescription = std::string("This N3 is a variant of the popular N3 (nonparameteric nonuniform ") +
                                     std::string("normalization) retrospective bias correction algorithm. Based ") +
                                     std::string("on the assumption that the corruption of the low frequency bias ") +
                                     std::string("field can be modeled as a convolution of the intensity histogram ") +
                                     std::string("by a Gaussian, the basic algorithmic protocol is to iterate ") +
                                     std::string("between deconvolving the intensity histogram by a Gaussian, ") +
                                     std::string("remapping the intensities, and then spatially smoothing this ") +
                                     std::string("result by a B-spline modeling of the bias field itself. ");

    parser->SetCommandDescription(commandDescription);
    N3InitializeCommandLineOptions(parser);

    if (parser->Parse(argc, argv) == EXIT_FAILURE)
    {
      return EXIT_FAILURE;
    }

    if (argc == 1)
    {
      parser->PrintMenu(std::cerr, 5, false);
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
    // Show automatic version
    itk::ants::CommandLineParser::OptionType::Pointer versionOption = parser->GetOption("version");
    if (versionOption && versionOption->GetNumberOfFunctions())
    {
      std::string versionFunction = versionOption->GetFunction(0)->GetName();
      ConvertToLowerCase(versionFunction);
      if (versionFunction.compare("1") == 0 || versionFunction.compare("true") == 0)
      {
        // Print Version Information
        std::cout << ANTs::Version::ExtendedVersionString() << std::endl;
        return EXIT_SUCCESS;
      }
    }
    // Get dimensionality
    unsigned int dimension = 3;

    itk::ants::CommandLineParser::OptionType::Pointer dimOption = parser->GetOption("image-dimensionality");
    if (dimOption && dimOption->GetNumberOfFunctions())
    {
      dimension = parser->Convert<unsigned int>(dimOption->GetFunction(0)->GetName());
    }
    else
    {
      // Read in the first intensity image to get the image dimension.
      std::string filename;

      itk::ants::CommandLineParser::OptionType::Pointer imageOption = parser->GetOption("input-image");
      if (imageOption && imageOption->GetNumberOfFunctions() > 0)
      {
        if (imageOption->GetFunction(0)->GetNumberOfParameters() > 0)
        {
          filename = imageOption->GetFunction(0)->GetParameter(0);
        }
        else
        {
          filename = imageOption->GetFunction(0)->GetName();
        }
      }
      else
      {
        std::cerr << "No input images were specified.  Specify an input image"
                  << " with the -i option" << std::endl;
        return EXIT_FAILURE;
      }
      itk::ImageIOBase::Pointer imageIO =
        itk::ImageIOFactory::CreateImageIO(filename.c_str(), itk::IOFileModeEnum::ReadMode);
      dimension = imageIO->GetNumberOfDimensions();
    }

    int returnValue = EXIT_FAILURE;

    switch (dimension)
    {
      case 2:
      {
        returnValue = N3<2>(parser);
      }
      break;
      case 3:
      {
        returnValue = N3<3>(parser);
      }
      break;
      case 4:
      {
        returnValue = N3<4>(parser);
      }
      break;
      default:
        std::cout << "Unsupported dimension" << std::endl;
        return EXIT_FAILURE;
    }
    return returnValue;
  }
  else
  {
    // put the arguments coming in as 'args' into standard (argc,argv) format;
    // 'args' doesn't have the command name as first, argument, so add it manually;
    // 'args' may have adjacent arguments concatenated into one argument,
    // which the parser should handle
    args.insert(args.begin(), "N3BiasFieldCorrection");

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

    if (argc < 4)
    {
      std::cout << "Usage: " << argv[0] << " imageDimension inputImage "
                << "outputImage [shrinkFactor] [maskImage] [numberOfIterations] "
                << "[numberOfFittingLevels] [outputBiasField] [verbose]" << std::endl;
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
        return N3BiasFieldCorrection<2>(argc, argv);
      }
      break;
      case 3:
      {
        return N3BiasFieldCorrection<3>(argc, argv);
      }
      break;
      case 4:
      {
        return N3BiasFieldCorrection<4>(argc, argv);
      }
      break;
      default:
        std::cout << "Unsupported dimension" << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }
}
} // namespace ants
