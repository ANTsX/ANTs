/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include "antsUtilities.h"
#include "antsAllocImage.h"
#include "ReadWriteData.h"
#include "antsCommandLineParser.h"
#include "itkCSVNumericObjectFileWriter.h"
#include "itkImageRegistrationMethodv4.h"
#include "itkSyNImageRegistrationMethod.h"
#include "itkDisplacementFieldTransform.h"
#include "itkANTSNeighborhoodCorrelationImageToImageMetricv4.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkCorrelationImageToImageMetricv4.h"
#include "itkImageToImageMetricv4.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkImageToHistogramFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkTransformToDisplacementFieldFilter.h"
#include "itkAffineTransform.h"
#include "itkBSplineTransform.h"
#include "itkBSplineSmoothingOnUpdateDisplacementFieldTransform.h"
#include "itkCompositeTransform.h"
#include "itkGaussianSmoothingOnUpdateDisplacementFieldTransform.h"
#include "itkIdentityTransform.h"
#include "itkEuler2DTransform.h"
#include "itkEuler3DTransform.h"
#include "itkSimilarity2DTransform.h"
#include "itkTransform.h"
#include "itkExtractImageFilter.h"
#include "itkBSplineTransformParametersAdaptor.h"
#include "itkBSplineSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor.h"
#include "itkGaussianSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor.h"
#include "itkTimeVaryingVelocityFieldTransformParametersAdaptor.h"
#include "itkGradientDescentOptimizerv4.h"
#include "itkConjugateGradientLineSearchOptimizerv4.h"
#include "itkQuasiNewtonOptimizerv4.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMacro.h"
#include "itkRegistrationParameterScalesFromPhysicalShift.h"
#include "itkResampleImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkTimeProbe.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkSimilarity2DTransform.h"
#include "itkSimilarity3DTransform.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkGaussianInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkLabelImageGaussianInterpolateImageFunction.h"
#include "itkLabelImageGenericInterpolateImageFunction.h"
#include <sstream>

namespace ants
{
/** \class antsRegistrationCommandIterationUpdate
 *  \brief change parameters between iterations of registration
 */
template <typename TFilter>
class antsSliceRegularizedRegistrationCommandIterationUpdate : public itk::Command
{
public:
  using Self = antsSliceRegularizedRegistrationCommandIterationUpdate<TFilter>;
  using Superclass = itk::Command;
  using Pointer = itk::SmartPointer<Self>;
  itkNewMacro(Self);

protected:
  antsSliceRegularizedRegistrationCommandIterationUpdate() { this->m_LogStream = &std::cout; }

public:
  void
  Execute(itk::Object * caller, const itk::EventObject & event) override
  {
    Execute((const itk::Object *)caller, event);
  }

  void
  Execute(const itk::Object * object, const itk::EventObject & event) override
  {
    auto * filter = const_cast<TFilter *>(dynamic_cast<const TFilter *>(object));

    unsigned int currentLevel = 0;

    if (typeid(event) == typeid(itk::IterationEvent))
    {
      currentLevel = filter->GetCurrentLevel() + 1;
    }
    if (currentLevel < this->m_NumberOfIterations.size())
    {
      typename TFilter::ShrinkFactorsPerDimensionContainerType shrinkFactors =
        filter->GetShrinkFactorsPerDimension(currentLevel);
      typename TFilter::SmoothingSigmasArrayType                 smoothingSigmas = filter->GetSmoothingSigmasPerLevel();
      typename TFilter::TransformParametersAdaptorsContainerType adaptors =
        filter->GetTransformParametersAdaptorsPerLevel();

      this->Logger() << "  Current level = " << currentLevel << std::endl;
      this->Logger() << "    number of iterations = " << this->m_NumberOfIterations[currentLevel] << std::endl;
      this->Logger() << "    shrink factors = " << shrinkFactors << std::endl;
      this->Logger() << "    smoothing sigmas = " << smoothingSigmas[currentLevel] << std::endl;
      this->Logger() << "    required fixed parameters = " << adaptors[currentLevel]->GetRequiredFixedParameters()
                     << std::endl;

      using GradientDescentOptimizerType = itk::ConjugateGradientLineSearchOptimizerv4;
      auto * optimizer = reinterpret_cast<GradientDescentOptimizerType *>(filter->GetModifiableOptimizer());
      optimizer->SetNumberOfIterations(this->m_NumberOfIterations[currentLevel]);
      optimizer->SetMinimumConvergenceValue(1.e-7);
      optimizer->SetConvergenceWindowSize(10);
      optimizer->SetLowerLimit(0);
      optimizer->SetUpperLimit(2);
      optimizer->SetEpsilon(0.1);
    }
  }

  void
  SetNumberOfIterations(const std::vector<unsigned int> & iterations)
  {
    this->m_NumberOfIterations = iterations;
  }

  void
  SetLogStream(std::ostream & logStream)
  {
    this->m_LogStream = &logStream;
  }

private:
  std::ostream &
  Logger() const
  {
    return *m_LogStream;
  }

  std::vector<unsigned int> m_NumberOfIterations;
  std::ostream *            m_LogStream;
};

template <typename T>
inline std::string
ants_slice_regularized_to_string(const T & t)
{
  std::stringstream ss;

  ss << t;
  return ss.str();
}

template <typename ImageType>
typename ImageType::Pointer
sliceRegularizedPreprocessImage(ImageType *                   inputImage,
                                typename ImageType::PixelType lowerScaleFunction,
                                typename ImageType::PixelType upperScaleFunction,
                                float                         winsorizeLowerQuantile,
                                float                         winsorizeUpperQuantile,
                                ImageType *                   histogramMatchSourceImage = nullptr)
{
  using HistogramFilterType = itk::Statistics::ImageToHistogramFilter<ImageType>;
  using InputBooleanObjectType = typename HistogramFilterType::InputBooleanObjectType;
  using HistogramSizeType = typename HistogramFilterType::HistogramSizeType;

  HistogramSizeType histogramSize(1);
  histogramSize[0] = 256;

  typename InputBooleanObjectType::Pointer autoMinMaxInputObject = InputBooleanObjectType::New();
  autoMinMaxInputObject->Set(true);

  typename HistogramFilterType::Pointer histogramFilter = HistogramFilterType::New();
  histogramFilter->SetInput(inputImage);
  histogramFilter->SetAutoMinimumMaximumInput(autoMinMaxInputObject);
  histogramFilter->SetHistogramSize(histogramSize);
  histogramFilter->SetMarginalScale(10.0);
  histogramFilter->Update();

  float lowerFunction = histogramFilter->GetOutput()->Quantile(0, winsorizeLowerQuantile);
  float upperFunction = histogramFilter->GetOutput()->Quantile(0, winsorizeUpperQuantile);
  using IntensityWindowingImageFilterType = itk::IntensityWindowingImageFilter<ImageType, ImageType>;

  typename IntensityWindowingImageFilterType::Pointer windowingFilter = IntensityWindowingImageFilterType::New();
  windowingFilter->SetInput(inputImage);
  windowingFilter->SetWindowMinimum(lowerFunction);
  windowingFilter->SetWindowMaximum(upperFunction);
  windowingFilter->SetOutputMinimum(lowerScaleFunction);
  windowingFilter->SetOutputMaximum(upperScaleFunction);
  windowingFilter->Update();

  typename ImageType::Pointer outputImage = nullptr;
  if (histogramMatchSourceImage)
  {
    using HistogramMatchingFilterType = itk::HistogramMatchingImageFilter<ImageType, ImageType>;
    typename HistogramMatchingFilterType::Pointer matchingFilter = HistogramMatchingFilterType::New();
    matchingFilter->SetSourceImage(windowingFilter->GetOutput());
    matchingFilter->SetReferenceImage(histogramMatchSourceImage);
    matchingFilter->SetNumberOfHistogramLevels(64);
    matchingFilter->SetNumberOfMatchPoints(12);
    matchingFilter->ThresholdAtMeanIntensityOn();
    matchingFilter->Update();

    outputImage = matchingFilter->GetOutput();
    outputImage->Update();
    outputImage->DisconnectPipeline();

    using CalculatorType = itk::MinimumMaximumImageCalculator<ImageType>;
    typename CalculatorType::Pointer calc = CalculatorType::New();
    calc->SetImage(inputImage);
    calc->ComputeMaximum();
    calc->ComputeMinimum();
    if (itk::Math::abs(calc->GetMaximum() - calc->GetMinimum()) < static_cast<float>(1.e-9))
    {
      std::cout << "Warning: bad time point - too little intensity variation" << calc->GetMinimum() << " "
                << calc->GetMaximum() << std::endl;
      return nullptr;
    }
  }
  else
  {
    outputImage = windowingFilter->GetOutput();
    outputImage->Update();
    outputImage->DisconnectPipeline();
  }
  return outputImage;
}

template <typename T>
struct ants_slice_regularized_index_cmp
{
  ants_slice_regularized_index_cmp(const T _arr)
    : arr(_arr)
  {}

  bool
  operator()(const size_t a, const size_t b) const
  {
    return arr[a] < arr[b];
  }

  const T arr;
};

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
    auto * filter = const_cast<TFilter *>(dynamic_cast<const TFilter *>(object));

    if (typeid(event) != typeid(itk::IterationEvent))
    {
      return;
    }

    unsigned int                                               currentLevel = filter->GetCurrentLevel();
    typename TFilter::TransformParametersAdaptorsContainerType adaptors =
      filter->GetTransformParametersAdaptorsPerLevel();

    using OptimizerType = itk::ConjugateGradientLineSearchOptimizerv4;
    auto * optimizer = reinterpret_cast<OptimizerType *>(filter->GetModifiableOptimizer());
    optimizer->SetNumberOfIterations(this->m_NumberOfIterations[currentLevel]);
    optimizer->SetMinimumConvergenceValue(1.e-7);
    optimizer->SetConvergenceWindowSize(10);
    optimizer->SetLowerLimit(0);
    optimizer->SetUpperLimit(2);
    optimizer->SetEpsilon(0.1);
  }

  void
  SetNumberOfIterations(std::vector<unsigned int> iterations)
  {
    this->m_NumberOfIterations = iterations;
  }

private:
  std::vector<unsigned int> m_NumberOfIterations;
};


void
ants_slice_poly_regularize(vnl_matrix<double>   A,
                           unsigned int         timedims,
                           vnl_vector<double> & solnx,
                           double &             interceptx,
                           vnl_matrix<double>   param_values,
                           unsigned int         whichCol)
{
  using RealType = double;
  using vVector = vnl_vector<RealType>;
  for (unsigned int z = 0; z < timedims; z++)
  {
    auto zz = static_cast<RealType>(z + 1);
    A(z, 0) = zz;
    for (unsigned int lcol = 1; lcol < A.cols(); lcol++)
    {
      A(z, lcol) = std::pow(zz, static_cast<RealType>(lcol + 1));
    }
  }
  for (unsigned int lcol = 0; lcol < A.cols(); lcol++)
  {
    vVector  acol = A.get_column(lcol);
    RealType acolsd = (acol - acol.mean()).rms();
    A.set_column(lcol, (acol - acol.mean()) / acolsd);
  }
  vnl_svd<double> svd(A);
  vVector         ob = param_values.get_column(whichCol);
  vVector         polyx = svd.solve(ob);
  interceptx = param_values.get_column(whichCol).mean();
  for (unsigned int Acol = 0; Acol < A.cols(); Acol++)
  {
    interceptx -= A.get_column(Acol).mean() * polyx(Acol);
  }
  solnx = A * polyx + interceptx;
}


template <unsigned int ImageDimension, typename TXType>
int
ants_slice_regularized_registration(itk::ants::CommandLineParser * parser)
{
  // We infer the number of stages by the number of transformations
  // specified by the user which should match the number of metrics.
  unsigned numberOfStages = 0;

  using PixelType = float;
  using RealType = double;
  using FixedIOImageType = itk::Image<PixelType, ImageDimension>;
  using FixedImageType = itk::Image<PixelType, ImageDimension - 1>;
  using MovingIOImageType = itk::Image<PixelType, ImageDimension>;
  using MovingImageType = itk::Image<PixelType, ImageDimension - 1>;

  using VectorIOType = itk::Vector<RealType, ImageDimension>;
  using DisplacementIOFieldType = itk::Image<VectorIOType, ImageDimension>;
  using VectorType = itk::Vector<RealType, ImageDimension - 1>;
  using DisplacementFieldType = itk::Image<VectorType, ImageDimension - 1>;

  using vMatrix = vnl_matrix<RealType>;
  using vVector = vnl_vector<RealType>;
  vMatrix param_values;
  using OptionType = typename itk::ants::CommandLineParser::OptionType;

  typename OptionType::Pointer transformOption = parser->GetOption("transform");
  if (transformOption && transformOption->GetNumberOfFunctions())
  {
    numberOfStages = transformOption->GetNumberOfFunctions();
  }
  else
  {
    std::cerr << "No transformations are specified." << std::endl;
    return EXIT_FAILURE;
  }
  if (numberOfStages != 1)
  {
    std::cerr << "Must specify one and only one stage." << std::endl;
    return EXIT_FAILURE;
  }
  typename OptionType::Pointer metricOption = parser->GetOption("metric");
  if (!metricOption || metricOption->GetNumberOfFunctions() != numberOfStages)
  {
    std::cerr << "The number of metrics specified does not match the number of stages." << std::endl;
    return EXIT_FAILURE;
  }

  std::string                                                whichInterpolator("linear");
  typename itk::ants::CommandLineParser::OptionType::Pointer interpolationOption = parser->GetOption("interpolation");
  if (interpolationOption && interpolationOption->GetNumberOfFunctions())
  {
    whichInterpolator = interpolationOption->GetFunction(0)->GetName();
    ConvertToLowerCase(whichInterpolator);
  }

  using ImageType = MovingImageType;
  typename ImageType::SpacingType cache_spacing_for_smoothing_sigmas(
    itk::NumericTraits<typename ImageType::SpacingType::ValueType>::ZeroValue());
  const unsigned int VImageDimension = ImageDimension - 1;
#include "make_interpolator_snip.tmpl"

  typename OptionType::Pointer iterationsOption = parser->GetOption("iterations");
  if (!iterationsOption || iterationsOption->GetNumberOfFunctions() != numberOfStages)
  {
    std::cerr << "The number of iteration sets specified does not match the number of stages." << std::endl;
    return EXIT_FAILURE;
  }

  typename OptionType::Pointer shrinkFactorsOption = parser->GetOption("shrinkFactors");
  if (!shrinkFactorsOption || shrinkFactorsOption->GetNumberOfFunctions() != numberOfStages)
  {
    std::cerr << "The number of shrinkFactor sets specified does not match the number of stages." << std::endl;
    return EXIT_FAILURE;
  }

  typename OptionType::Pointer smoothingSigmasOption = parser->GetOption("smoothingSigmas");
  if (!smoothingSigmasOption || smoothingSigmasOption->GetNumberOfFunctions() != numberOfStages)
  {
    std::cerr << "The number of smoothing sigma sets specified does not match the number of stages." << std::endl;
    return EXIT_FAILURE;
  }

  typename OptionType::Pointer outputOption = parser->GetOption("output");
  if (!outputOption)
  {
    std::cerr << "Output option not specified." << std::endl;
    return EXIT_FAILURE;
  }
  std::string outputPrefix = outputOption->GetFunction(0)->GetParameter(0);
  if (outputPrefix.length() < 3)
  {
    outputPrefix = outputOption->GetFunction(0)->GetName();
  }

  std::string                  maskfn = std::string("");
  typename OptionType::Pointer maskOption = parser->GetOption("mask");
  if (maskOption->GetNumberOfFunctions() > 0)
  {
    maskfn = maskOption->GetFunction(0)->GetName();
  }


  bool doEstimateLearningRateOnce(true);

  itk::TimeProbe totalTimer;
  totalTimer.Start();
  using TranslationRegistrationType = itk::ImageRegistrationMethodv4<FixedImageType, MovingImageType, TXType>;
  // We iterate backwards because the command line options are stored as a stack (first in last out)
  using SingleTransformItemType = typename TXType::Pointer;
  std::vector<SingleTransformItemType>          transformList;
  std::vector<SingleTransformItemType>          transformUList;
  std::vector<typename FixedImageType::Pointer> fixedSliceList;
  std::vector<typename FixedImageType::Pointer> movingSliceList;
  typename FixedIOImageType::Pointer            maskImage;
  using ImageMaskType = itk::Image<unsigned char, ImageDimension - 1>;
  typename ImageMaskType::Pointer mask_time_slice = nullptr;
  if (maskfn.length() > 3)
    ReadImage<FixedIOImageType>(maskImage, maskfn.c_str());

  for (int currentStage = numberOfStages - 1; currentStage >= 0; currentStage--)
  {
    std::stringstream currentStageString;
    currentStageString << currentStage;

    // Get the fixed and moving images
    std::string                        fixedImageFileName = metricOption->GetFunction(currentStage)->GetParameter(0);
    std::string                        movingImageFileName = metricOption->GetFunction(currentStage)->GetParameter(1);
    typename FixedImageType::Pointer   fixed_time_slice = nullptr;
    typename FixedImageType::Pointer   moving_time_slice = nullptr;
    typename FixedIOImageType::Pointer fixedImage;
    ReadImage<FixedIOImageType>(fixedImage, fixedImageFileName.c_str());
    unsigned int timedims = fixedImage->GetLargestPossibleRegion().GetSize()[ImageDimension - 1];
    unsigned int nparams = TXType::New()->GetNumberOfParameters();
    param_values.set_size(timedims, nparams);
    param_values.fill(0);

    typename MovingIOImageType::Pointer movingImage;
    ReadImage<MovingIOImageType>(movingImage, movingImageFileName.c_str());
    unsigned int moving_timedims = movingImage->GetLargestPossibleRegion().GetSize()[ImageDimension - 1];
    if (timedims != moving_timedims)
    {
      std::cerr << "We require that the n^th dimensions of the fixed and moving image are equal" << std::endl;
      return EXIT_FAILURE;
    }

    typename FixedIOImageType::Pointer outputImage;
    ReadImage<FixedIOImageType>(outputImage, fixedImageFileName.c_str());
    outputImage->FillBuffer(0);

    // Get the number of iterations and use that information to specify the number of levels
    std::vector<unsigned int> iterations =
      parser->ConvertVector<unsigned int>(iterationsOption->GetFunction(currentStage)->GetName());
    unsigned int numberOfLevels = iterations.size();

    // Get shrink factors
    std::vector<unsigned int> factors =
      parser->ConvertVector<unsigned int>(shrinkFactorsOption->GetFunction(currentStage)->GetName());
    typename TranslationRegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
    shrinkFactorsPerLevel.SetSize(factors.size());

    if (factors.size() != numberOfLevels)
    {
      std::cerr << "ERROR:  The number of shrink factors does not match the number of levels." << std::endl;
      return EXIT_FAILURE;
    }
    else
    {
      for (unsigned int n = 0; n < shrinkFactorsPerLevel.Size(); n++)
      {
        shrinkFactorsPerLevel[n] = factors[n];
      }
    }

    // Get smoothing sigmas

    std::vector<float> sigmas =
      parser->ConvertVector<float>(smoothingSigmasOption->GetFunction(currentStage)->GetName());
    typename TranslationRegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
    smoothingSigmasPerLevel.SetSize(sigmas.size());

    if (sigmas.size() != numberOfLevels)
    {
      std::cerr << "ERROR:  The number of smoothing sigmas does not match the number of levels." << std::endl;
      return EXIT_FAILURE;
    }
    else
    {
      for (unsigned int n = 0; n < smoothingSigmasPerLevel.Size(); n++)
      {
        smoothingSigmasPerLevel[n] = sigmas[n];
      }
    }

    bool verbose = false;
    if (parser->Convert<bool>(parser->GetOption('v')->GetFunction()->GetName()))
    {
      verbose = true;
    }


    std::vector<unsigned int>                         polydegree;
    itk::ants::CommandLineParser::OptionType::Pointer polyOption = parser->GetOption("polydegree");
    if (polyOption && polyOption->GetNumberOfFunctions())
    {
      polydegree = parser->ConvertVector<unsigned int>(polyOption->GetFunction(0)->GetName());
    }
    if ((polydegree.size() != nparams))
    {
      polydegree.resize(nparams, polydegree[0]);
    }
    for (unsigned int & pind : polydegree)
    {
      if (pind > (timedims - 2))
        pind = timedims - 2;
    }

    // the fixed image slice is a reference image in 2D while the moving is a 2D slice image
    // loop over every time point and register image_i_moving to image_i_fixed
    //
    // Set up the image metric and scales estimator
    for (unsigned int timedim = 0; timedim < timedims; timedim++)
    {
      using IdentityTransformType = itk::IdentityTransform<RealType, ImageDimension - 1>;
      typename IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();
      using ExtractFilterType = itk::ExtractImageFilter<FixedIOImageType, FixedImageType>;
      typename FixedIOImageType::RegionType extractRegion = movingImage->GetLargestPossibleRegion();
      extractRegion.SetSize(ImageDimension - 1, 0);
      extractRegion.SetIndex(ImageDimension - 1, timedim);

      typename ExtractFilterType::Pointer extractFilterF = ExtractFilterType::New();
      typename ExtractFilterType::Pointer extractFilterM = ExtractFilterType::New();
      extractFilterF->SetInput(fixedImage);
      extractFilterM->SetInput(movingImage);
      extractFilterF->SetDirectionCollapseToSubmatrix();
      extractFilterM->SetDirectionCollapseToSubmatrix();
      bool toidentity = false;
      if (toidentity)
        extractFilterF->SetDirectionCollapseToIdentity();
      if (toidentity)
        extractFilterM->SetDirectionCollapseToIdentity();
      extractFilterF->SetExtractionRegion(extractRegion);
      extractFilterM->SetExtractionRegion(extractRegion);
      extractFilterF->Update();
      extractFilterM->Update();
      fixed_time_slice = extractFilterF->GetOutput();
      moving_time_slice = extractFilterM->GetOutput();
      fixedSliceList.push_back(fixed_time_slice);
      movingSliceList.push_back(moving_time_slice);

      if (maskfn.length() > 3)
      {
        using ExtractFilterTypeX = itk::ExtractImageFilter<FixedIOImageType, ImageMaskType>;
        typename ExtractFilterTypeX::Pointer extractFilterX = ExtractFilterTypeX::New();
        extractFilterX->SetInput(maskImage);
        extractFilterX->SetDirectionCollapseToSubmatrix();
        if (toidentity)
          extractFilterX->SetDirectionCollapseToIdentity();
        extractFilterX->SetExtractionRegion(extractRegion);
        extractFilterX->Update();
        mask_time_slice = extractFilterX->GetOutput();
      }

      // set up initial transform parameters
      typename TXType::Pointer translationTransform = TXType::New();
      translationTransform->SetIdentity();
      transformList.push_back(translationTransform);

      // set up update transform parameters
      typename TXType::Pointer translationTransformU = TXType::New();
      translationTransformU->SetIdentity();
      transformUList.push_back(translationTransformU);
    }

    // implement a gradient descent on the polynomial parameters by looping over registration results
    using MetricType = itk::ImageToImageMetricv4<FixedImageType, FixedImageType>;
    typename MetricType::Pointer metric;
    unsigned int                 maxloop = 2;
    for (unsigned int loop = 0; loop < maxloop; loop++)
    {
      RealType metricval = 0;
      bool     skipThisTimePoint = false;
      for (unsigned int timedim = 0; timedim < timedims; timedim++)
      {
        typename FixedImageType::Pointer preprocessFixedImage =
          sliceRegularizedPreprocessImage<FixedImageType>(fixedSliceList[timedim], 0, 1, 0.005, 0.995, nullptr);

        typename FixedImageType::Pointer preprocessMovingImage = sliceRegularizedPreprocessImage<FixedImageType>(
          movingSliceList[timedim], 0, 1, 0.005, 0.995, preprocessFixedImage);
        if (preprocessFixedImage.IsNull() || preprocessMovingImage.IsNull())
        {
          preprocessFixedImage = fixedSliceList[timedim];
          preprocessMovingImage = movingSliceList[timedim];
          skipThisTimePoint = true;
        }

        std::string whichMetric = metricOption->GetFunction(currentStage)->GetName();
        ConvertToLowerCase(whichMetric);

        float samplingPercentage = 1.0;
        if (metricOption->GetFunction(0)->GetNumberOfParameters() > 5)
        {
          samplingPercentage = parser->Convert<float>(metricOption->GetFunction(currentStage)->GetParameter(5));
        }

        std::string samplingStrategy = "";
        if (metricOption->GetFunction(0)->GetNumberOfParameters() > 4)
        {
          samplingStrategy = metricOption->GetFunction(currentStage)->GetParameter(4);
        }
        ConvertToLowerCase(samplingStrategy);
        itk::ImageRegistrationMethodv4Enums::MetricSamplingStrategy metricSamplingStrategy =
		  itk::ImageRegistrationMethodv4Enums::MetricSamplingStrategy::NONE;
        if (std::strcmp(samplingStrategy.c_str(), "random") == 0)
        {
          metricSamplingStrategy = itk::ImageRegistrationMethodv4Enums::MetricSamplingStrategy::RANDOM;
        }
        if (std::strcmp(samplingStrategy.c_str(), "regular") == 0)
        {
          metricSamplingStrategy = itk::ImageRegistrationMethodv4Enums::MetricSamplingStrategy::REGULAR;
        }

        if (std::strcmp(whichMetric.c_str(), "cc") == 0)
        {
          auto radiusOption = parser->Convert<unsigned int>(metricOption->GetFunction(currentStage)->GetParameter(3));

          using CorrelationMetricType =
            itk::ANTSNeighborhoodCorrelationImageToImageMetricv4<FixedImageType, FixedImageType>;
          typename CorrelationMetricType::Pointer    correlationMetric = CorrelationMetricType::New();
          typename CorrelationMetricType::RadiusType radius;
          radius.Fill(radiusOption);
          correlationMetric->SetRadius(radius);
          correlationMetric->SetUseMovingImageGradientFilter(false);
          correlationMetric->SetUseFixedImageGradientFilter(false);
          metric = correlationMetric;
        }
        else if (std::strcmp(whichMetric.c_str(), "mi") == 0)
        {
          auto binOption = parser->Convert<unsigned int>(metricOption->GetFunction(currentStage)->GetParameter(3));
          using MutualInformationMetricType =
            itk::MattesMutualInformationImageToImageMetricv4<FixedImageType, FixedImageType>;
          typename MutualInformationMetricType::Pointer mutualInformationMetric = MutualInformationMetricType::New();
          // mutualInformationMetric = mutualInformationMetric;
          mutualInformationMetric->SetNumberOfHistogramBins(binOption);
          mutualInformationMetric->SetUseMovingImageGradientFilter(false);
          mutualInformationMetric->SetUseFixedImageGradientFilter(false);
          metric = mutualInformationMetric;
        }
        else if (std::strcmp(whichMetric.c_str(), "meansquares") == 0)
        {
          using MSQMetricType = itk::MeanSquaresImageToImageMetricv4<FixedImageType, FixedImageType>;
          typename MSQMetricType::Pointer demonsMetric = MSQMetricType::New();
          // demonsMetric = demonsMetric;
          metric = demonsMetric;
        }
        else if (std::strcmp(whichMetric.c_str(), "gc") == 0)
        {
          using corrMetricType = itk::CorrelationImageToImageMetricv4<FixedImageType, FixedImageType>;
          typename corrMetricType::Pointer corrMetric = corrMetricType::New();
          metric = corrMetric;
        }
        else
        {
          std::cerr << "ERROR: Unrecognized image metric: " << whichMetric << std::endl;
          return EXIT_FAILURE;
        }
        metric->SetVirtualDomainFromImage(fixedSliceList[timedim]);
        if (maskfn.length() > 3)
        {
          using spMaskType = itk::ImageMaskSpatialObject<ImageDimension - 1>;
          typename spMaskType::Pointer spatialObjectMask = spMaskType::New();
          spatialObjectMask->SetImage(mask_time_slice);
          metric->SetFixedImageMask(spatialObjectMask);
          if ((verbose) && (loop == 0) && (timedim == 0))
            std::cout << " setting mask " << maskfn << std::endl;
        }
        using ScalesEstimatorType = itk::RegistrationParameterScalesFromPhysicalShift<MetricType>;
        typename ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();
        scalesEstimator->SetMetric(metric);
        scalesEstimator->SetTransformForward(true);
        auto learningRate = parser->Convert<float>(transformOption->GetFunction(currentStage)->GetParameter(0));

        using OptimizerType = itk::ConjugateGradientLineSearchOptimizerv4;
        typename OptimizerType::Pointer optimizer = OptimizerType::New();
        optimizer->SetNumberOfIterations(iterations[0]);
        optimizer->SetMinimumConvergenceValue(1.e-7);
        optimizer->SetConvergenceWindowSize(8);
        optimizer->SetLowerLimit(0);
        optimizer->SetUpperLimit(2);
        optimizer->SetEpsilon(0.1);
        optimizer->SetScalesEstimator(scalesEstimator);
        optimizer->SetMaximumStepSizeInPhysicalUnits(learningRate);
        optimizer->SetDoEstimateLearningRateOnce(doEstimateLearningRateOnce);
        optimizer->SetDoEstimateLearningRateAtEachIteration(!doEstimateLearningRateOnce);

        // Set up the image registration methods along with the transforms
        std::string whichTransform = transformOption->GetFunction(currentStage)->GetName();
        ConvertToLowerCase(whichTransform);
        typename TranslationRegistrationType::Pointer translationRegistration = TranslationRegistrationType::New();
        if ((std::strcmp(whichTransform.c_str(), "translation") == 0) ||
            (std::strcmp(whichTransform.c_str(), "rigid") == 0) ||
            (std::strcmp(whichTransform.c_str(), "similarity") == 0))
        {
          transformList[timedim]->GetNumberOfParameters();
          metric->SetFixedImage(preprocessFixedImage);
          metric->SetVirtualDomainFromImage(preprocessFixedImage);
          metric->SetMovingImage(preprocessMovingImage);
          metric->SetMovingTransform(transformList[timedim]);
          typename ScalesEstimatorType::ScalesType scales(transformList[timedim]->GetNumberOfParameters());
          typename MetricType::ParametersType      newparams(transformList[timedim]->GetParameters());
          metric->SetParameters(newparams);
          metric->Initialize();
          scalesEstimator->SetMetric(metric);
          scalesEstimator->EstimateScales(scales);
          optimizer->SetScales(scales);
          translationRegistration->SetFixedImage(preprocessFixedImage);
          translationRegistration->SetMovingImage(preprocessMovingImage);
          translationRegistration->SetNumberOfLevels(numberOfLevels);
          translationRegistration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel);
          translationRegistration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);
          translationRegistration->SetMetricSamplingStrategy(metricSamplingStrategy);
          translationRegistration->SetMetricSamplingPercentage(samplingPercentage);
          translationRegistration->SetMetric(metric);
          translationRegistration->SetOptimizer(optimizer);

          using TranslationCommandType = CommandIterationUpdate<TranslationRegistrationType>;
          typename TranslationCommandType::Pointer translationObserver = TranslationCommandType::New();
          translationObserver->SetNumberOfIterations(iterations);
          translationRegistration->AddObserver(itk::IterationEvent(), translationObserver);
          if (!skipThisTimePoint)
            try
            {
              translationRegistration->Update();
            }
            catch (const itk::ExceptionObject & e)
            {
              std::cerr << "Exception caught: " << e << std::endl;
              return EXIT_FAILURE;
            }
        }
        else
        {
          std::cerr << "ERROR:  Unrecognized transform option - " << whichTransform << std::endl;
          return EXIT_FAILURE;
        }
        transformUList[timedim] = translationRegistration->GetModifiableTransform();
        metricval += metric->GetValue();
      }

      for (unsigned int i = 0; i < transformList.size(); i++)
      {
        typename TXType::ParametersType pu = transformUList[i]->GetParameters();
        for (unsigned int kk = 0; kk < pu.Size(); kk++)
          param_values(i, kk) = pu[kk];
      }

      // project updated solution back to polynomial space
      vVector  solnx;
      RealType interceptx = 0;
      if (polydegree[0] > 0)
      {
        // set up polynomial system of equations
        vMatrix A(timedims, polydegree[0], 0.0);
        ants_slice_poly_regularize(A, timedims, solnx, interceptx, param_values, 0);
      }

      // now y-dimension
      vVector  solny;
      RealType intercepty = 0;
      if (polydegree[1] > 0)
      {
        vMatrix A(timedims, polydegree[1], 0.0);
        ants_slice_poly_regularize(A, timedims, solny, intercepty, param_values, 1);
      }

      // now rotation
      vVector  solnr;
      RealType interceptr = 0;
      if (nparams >= 3)
        if (polydegree[2] > 0)
        {
          vMatrix A(timedims, polydegree[2], 0.0);
          ants_slice_poly_regularize(A, timedims, solnr, interceptr, param_values, 2);
        }

      // now scaling-dimension
      vVector  solns;
      RealType intercepts = 0;
      if (nparams >= 4)
        if (polydegree[3] > 0)
        {
          vMatrix A(timedims, polydegree[3], 0.0);
          ants_slice_poly_regularize(A, timedims, solns, intercepts, param_values, 3);
        }

      // now look at delta and do projection
      if ((solnx.size() != transformList.size()) && (polydegree[0] > 0))
        std::cerr << "solnx.size() != transformList.size()" << std::endl;
      if ((solny.size() != transformList.size()) && (polydegree[1] > 0))
        std::cerr << "solny.size() != transformList.size()" << std::endl;

      RealType err = 0;
      for (unsigned int i = 0; i < transformList.size(); i++)
      {
        /** FIXME - this should be vectorized: DRY */
        typename TXType::ParametersType pOld = transformList[i]->GetParameters();
        typename TXType::ParametersType p = transformUList[i]->GetParameters();
        if (polydegree[0] > 0)
          p[0] = solnx[i];
        if (polydegree[1] > 0)
          p[1] = solny[i];
        param_values(i, 0) = p[0];
        param_values(i, 1) = p[1];
        if (nparams >= 3)
        {
          if (polydegree[2] > 0)
            p[2] = solnr[i];
          param_values(i, 2) = p[2];
        }
        if (nparams >= 4)
        {
          if (polydegree[3] > 0)
            p[3] = solns[i];
          param_values(i, 3) = p[3];
        }
        transformList[i]->SetParameters(p);
        transformUList[i]->SetParameters(p);
        err += (p - pOld).rms();
      }
      err = err / static_cast<double>(transformList.size());
      if (verbose)
      {
        std::cout << "Loop" << loop << " polyerr: " << err << " image-metric " << metricval << std::endl;
      }
      //  transformList = transformUList;
      for (unsigned int i = 0; i < transformList.size(); i++)
      {
        typename TXType::ParametersType p = transformList[i]->GetParameters();
        if (polydegree[0] == 0)
          param_values(i, 0) = p[0];
        if (polydegree[1] == 0)
          param_values(i, 1) = p[1];
        if (nparams >= 3)
          if (polydegree[2] == 0)
            param_values(i, 2) = p[2];
        if (nparams >= 4)
          if (polydegree[3] == 0)
            param_values(i, 3) = p[3];
      }
    } // done with optimization, now move on to writing data ... closes for ( unsigned int loop = 0; loop < maxloop;
      // loop++ )

    // write polynomial predicted data
    {
      std::vector<std::string> ColumnHeaders;
      std::string              colname;
      colname = std::string("Tx");
      ColumnHeaders.push_back(colname);
      colname = std::string("Ty");
      ColumnHeaders.push_back(colname);
      if (nparams >= 3)
      {
        colname = std::string("Tr");
        ColumnHeaders.push_back(colname);
      }
      if (nparams >= 4)
      {
        colname = std::string("Ts");
        ColumnHeaders.push_back(colname);
      }
      using WriterType = itk::CSVNumericObjectFileWriter<double, 1, 1>;
      WriterType::Pointer writer = WriterType::New();
      std::string         fnmp;
      fnmp = outputPrefix + std::string("TxTy_poly.csv");
      writer->SetFileName(fnmp.c_str());
      writer->SetColumnHeaders(ColumnHeaders);
      writer->SetInput(&param_values);
      writer->Write();
    }
    /** Handle all output: images and displacement fields */
    using IdentityIOTransformType = itk::IdentityTransform<RealType, ImageDimension>;
    typename IdentityIOTransformType::Pointer identityIOTransform = IdentityIOTransformType::New();
    using ConverterType = typename itk::TransformToDisplacementFieldFilter<DisplacementIOFieldType, RealType>;
    typename ConverterType::Pointer idconverter = ConverterType::New();
    idconverter->SetOutputOrigin(outputImage->GetOrigin());
    idconverter->SetOutputStartIndex(outputImage->GetBufferedRegion().GetIndex());
    idconverter->SetSize(outputImage->GetBufferedRegion().GetSize());
    idconverter->SetOutputSpacing(outputImage->GetSpacing());
    idconverter->SetOutputDirection(outputImage->GetDirection());
    idconverter->SetTransform(identityIOTransform);
    idconverter->Update();
    typename DisplacementIOFieldType::Pointer displacementout = idconverter->GetOutput();
    typename DisplacementIOFieldType::Pointer displacementinv = DisplacementIOFieldType::New();
    displacementinv->CopyInformation(displacementout);
    displacementinv->SetRegions(displacementout->GetRequestedRegion());
    displacementinv->AllocateInitialized();
    typename DisplacementIOFieldType::IndexType dind;
    dind.Fill(0);
    displacementinv->FillBuffer(displacementout->GetPixel(dind));
    for (unsigned int timedim = 0; timedim < timedims; timedim++)
    {
      using _ConverterType = typename itk::TransformToDisplacementFieldFilter<DisplacementFieldType, RealType>;
      typename _ConverterType::Pointer converter = _ConverterType::New();
      converter->SetOutputOrigin(fixedSliceList[timedim]->GetOrigin());
      converter->SetOutputStartIndex(fixedSliceList[timedim]->GetBufferedRegion().GetIndex());
      converter->SetSize(fixedSliceList[timedim]->GetBufferedRegion().GetSize());
      converter->SetOutputSpacing(fixedSliceList[timedim]->GetSpacing());
      converter->SetOutputDirection(fixedSliceList[timedim]->GetDirection());
      converter->SetTransform(transformList[timedim]);
      converter->Update();

      // resample the moving image and then put it in its place
      interpolator->SetInputImage(movingSliceList[timedim]);
      using ResampleFilterType = itk::ResampleImageFilter<FixedImageType, FixedImageType>;
      typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
      resampler->SetTransform(transformList[timedim]);
      resampler->SetInterpolator(interpolator);
      resampler->SetInput(movingSliceList[timedim]);
      resampler->SetOutputParametersFromImage(fixedSliceList[timedim]);
      resampler->SetDefaultPixelValue(0);
      resampler->Update();

      /** Here, we put the resampled 2D image into the 3D volume */
      using Iterator = itk::ImageRegionIteratorWithIndex<FixedImageType>;
      Iterator vfIter2(resampler->GetOutput(), resampler->GetOutput()->GetLargestPossibleRegion());
      for (vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2)
      {
        typename FixedImageType::PixelType fval = vfIter2.Get();
        VectorType                         vec = converter->GetOutput()->GetPixel(vfIter2.GetIndex());
        VectorIOType                       vecout;
        vecout.Fill(0);
        typename MovingIOImageType::IndexType ind;
        for (itk::SizeValueType xx = 0; xx < ImageDimension - 1; xx++)
        {
          ind[xx] = vfIter2.GetIndex()[xx];
          vecout[xx] = vec[xx];
        }
        unsigned int tdim = timedim;
        if (tdim > (timedims - 1))
        {
          tdim = timedims - 1;
        }
        ind[MovingIOImageType::ImageDimension - 1] = tdim;
        outputImage->SetPixel(ind, fval);
        displacementout->SetPixel(ind, vecout);
      }
    }

    // now apply to the inverse map
    unsigned int timedimX = 0;
    for (timedimX = 0; timedimX < timedims; timedimX++)
    {
      using _ConverterType = typename itk::TransformToDisplacementFieldFilter<DisplacementFieldType, RealType>;
      typename _ConverterType::Pointer converter = _ConverterType::New();
      converter->SetOutputOrigin(movingSliceList[timedimX]->GetOrigin());
      converter->SetOutputStartIndex(movingSliceList[timedimX]->GetBufferedRegion().GetIndex());
      converter->SetSize(movingSliceList[timedimX]->GetBufferedRegion().GetSize());
      converter->SetOutputSpacing(movingSliceList[timedimX]->GetSpacing());
      converter->SetOutputDirection(movingSliceList[timedimX]->GetDirection());
      typename TXType::Pointer invtx = TXType::New();
      invtx->SetIdentity();
      transformList[timedimX]->GetInverse(invtx);
      converter->SetTransform(invtx);
      converter->Update();

      // resample the moving image and then put it in its place
      using ResampleFilterType = itk::ResampleImageFilter<FixedImageType, FixedImageType>;
      typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
      resampler->SetTransform(invtx);
      interpolator->SetInputImage(fixedSliceList[timedimX]);
      resampler->SetInterpolator(interpolator);
      resampler->SetInput(fixedSliceList[timedimX]);
      resampler->SetOutputParametersFromImage(movingSliceList[timedimX]);
      resampler->SetDefaultPixelValue(0);
      resampler->Update();

      /** Here, we put the resampled 2D image into the 3D volume */
      using Iterator = itk::ImageRegionIteratorWithIndex<FixedImageType>;
      Iterator vfIter2(resampler->GetOutput(), resampler->GetOutput()->GetLargestPossibleRegion());
      for (vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2)
      {
        VectorType   vec = converter->GetOutput()->GetPixel(vfIter2.GetIndex());
        VectorIOType vecout;
        vecout.Fill(0);
        typename MovingIOImageType::IndexType ind;
        for (itk::SizeValueType xx = 0; xx < ImageDimension - 1; xx++)
        {
          ind[xx] = vfIter2.GetIndex()[xx];
          vecout[xx] = vec[xx];
        }
        unsigned int tdim = timedimX;
        if (tdim > (timedims - 1))
        {
          tdim = timedims - 1;
        }
        ind[MovingIOImageType::ImageDimension - 1] = tdim;
        displacementinv->SetPixel(ind, vecout);
      }
    }

    if (outputOption && outputOption->GetFunction(0)->GetNumberOfParameters() > 1 && currentStage == 0)
    {
      std::string fileName = outputOption->GetFunction(0)->GetParameter(1);
      if (outputPrefix.length() < 3)
      {
        outputPrefix = outputOption->GetFunction(0)->GetName();
      }
      ANTs::WriteImage<MovingIOImageType>(outputImage, fileName.c_str());
    }
    if (outputOption && outputOption->GetFunction(0)->GetNumberOfParameters() > 2 && outputImage && currentStage == 0)
    {
      std::string fileName = outputOption->GetFunction(0)->GetParameter(2);
    }
    {
      std::string dispfn = outputPrefix + std::string("Warp.nii.gz");
      using DisplacementFieldWriterType = itk::ImageFileWriter<DisplacementIOFieldType>;
      typename DisplacementFieldWriterType::Pointer displacementFieldWriter = DisplacementFieldWriterType::New();
      displacementFieldWriter->SetInput(displacementout);
      displacementFieldWriter->SetFileName(dispfn.c_str());
      displacementFieldWriter->Update();
    }
    {
      std::string dispfn = outputPrefix + std::string("InverseWarp.nii.gz");
      using DisplacementFieldWriterType = itk::ImageFileWriter<DisplacementIOFieldType>;
      typename DisplacementFieldWriterType::Pointer displacementFieldWriter = DisplacementFieldWriterType::New();
      displacementFieldWriter->SetInput(displacementinv);
      displacementFieldWriter->SetFileName(dispfn.c_str());
      displacementFieldWriter->Update();
    }
  }
  totalTimer.Stop();
  //  std::cout << std::endl << "Total elapsed time: " << totalTimer.GetMean() << std::endl;
  return EXIT_SUCCESS;
}

void
antsSliceRegularizedRegistrationInitializeCommandLineOptions(itk::ants::CommandLineParser * parser)
{
  using OptionType = itk::ants::CommandLineParser::OptionType;

  {
    std::string description =
      std::string("Four image metrics are available--- ") +
      std::string("GC : global correlation, CC:  ANTS neighborhood cross correlation, MI:  Mutual information, and ") +
      std::string("MeanSquares:  mean-squares intensity difference. ") +
      std::string("Note that the metricWeight is currently not used.  ") +
      std::string("Rather, it is a temporary place holder until multivariate metrics ") +
      std::string("are available for a single stage.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("metric");
    option->SetShortName('m');
    option->SetUsageOption(0,
                           "CC[fixedImage,movingImage,metricWeight,radius,<samplingStrategy={Regular,Random,None}>,<"
                           "samplingPercentage=[0,1]>]");
    option->SetUsageOption(1,
                           "MI[fixedImage,movingImage,metricWeight,numberOfBins,<samplingStrategy={Regular,Random,None}"
                           ">,<samplingPercentage=[0,1]>]");
    option->SetUsageOption(2,

                           "MeanSquares[fixedImage,movingImage,metricWeight,radius,<samplingStrategy={Regular,Random,"
                           "None}>,<samplingPercentage=[0,1]>]");
    option->SetUsageOption(3,
                           "GC[fixedImage,movingImage,metricWeight,radius,<samplingStrategy={Regular,Random,None}>,<"
                           "samplingPercentage=[0,1]>]");
    option->SetDescription(description);
    parser->AddOption(option);
  }


  {
    std::string         description = "Fixed image mask to limit voxels considered by the metric.";
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("mask");
    option->SetShortName('x');
    option->SetUsageOption(0, "mask-in-fixed-image-space.nii.gz");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Several interpolation options are available in ITK. ") +
                              std::string("These have all been made available.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("interpolation");
    option->SetShortName('n');
    option->SetUsageOption(0, "Linear");
    option->SetUsageOption(1, "NearestNeighbor");
    option->SetUsageOption(2, "MultiLabel[<sigma=imageSpacing>,<alpha=4.0>]");
    option->SetUsageOption(3, "Gaussian[<sigma=imageSpacing>,<alpha=1.0>]");
    option->SetUsageOption(4, "BSpline[<order=3>]");
    option->SetUsageOption(5, "CosineWindowedSinc");
    option->SetUsageOption(6, "WelchWindowedSinc");
    option->SetUsageOption(7, "HammingWindowedSinc");
    option->SetUsageOption(8, "LanczosWindowedSinc");
    option->SetUsageOption(9, "GenericLabel[<interpolator=Linear>]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("Several transform options are available.  The gradientStep or") +
      std::string("learningRate characterizes the gradient descent optimization and is scaled appropriately ") +
      std::string("for each transform using the shift scales estimator.  Subsequent parameters are ") +
      std::string("transform-specific and can be determined from the usage. ");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("transform");
    option->SetShortName('t');
    option->SetUsageOption(0, "Translation[gradientStep]");
    option->SetUsageOption(1, "Rigid[gradientStep]");
    option->SetUsageOption(2, "Similarity[gradientStep]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Specify the number of iterations at each level.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("iterations");
    option->SetShortName('i');
    option->SetUsageOption(0, "MxNx0...");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description = std::string("Specify the amount of smoothing at each level.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("smoothingSigmas");
    option->SetShortName('s');
    option->SetUsageOption(0, "MxNx0...");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("Specify the shrink factor for the virtual domain (typically the fixed image) at each level.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("shrinkFactors");
    option->SetShortName('f');
    option->SetUsageOption(0, "MxNx0...");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("Specify the output transform prefix (output format is .nii.gz ).") +
      std::string("Optionally, one can choose to warp the moving image to the fixed space and, if the ") +
      std::string("inverse transform exists, one can also output the warped fixed image.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName("output");
    option->SetShortName('o');
    option->SetUsageOption(0, "[outputTransformPrefix,<outputWarpedImage>,<outputAverageImage>]");
    option->SetDescription(description);
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("Print the help menu (short version).");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("help");
    option->SetShortName('h');
    option->SetDescription(description);
    option->AddFunction(std::string("0"));
    parser->AddOption(option);
  }

  {
    std::string         description = std::string("verbose option");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("verbose");
    option->SetShortName('v');
    option->SetDescription(description);
    option->AddFunction(std::string("0"));
    parser->AddOption(option);
  }

  {
    std::string description =
      std::string("degree of polynomial - up to zDimension-2. ") +
      std::string("Controls the polynomial degree. 0 means no regularization. This may be a vector ") +
      std::string("denoted by 2x2x1 for a 3-parameter transform ( e.g. rigid ).  This would ") +
      std::string("regularize the translation by 2nd degree polynomial and the rotation by a ") +
      std::string("linear function.");
    OptionType::Pointer option = OptionType::New();
    option->SetLongName("polydegree");
    option->SetShortName('p');
    option->SetDescription(description);
    parser->AddOption(option);
  }
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int
antsSliceRegularizedRegistration(std::vector<std::string> args, std::ostream * /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert(args.begin(), "antsSliceRegularizedRegistration");

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

  std::string commandDescription =
    std::string("antsSliceRegularizedRegistration ") +
    std::string("This program is a user-level application for slice-by-slice translation registration. ") +
    std::string(
      "Results are regularized in z using polynomial regression.  The program is targeted at spinal cord MRI. ") +
    std::string("Only one stage is supported where a stage consists of a transform; an image metric; ") +
    std::string("and iterations, shrink factors, and smoothing sigmas for each level. ") +
    std::string("Specialized for 3D data: fixed image is 3D, moving image is 3D. ") +
    std::string("Registration is performed slice-by-slice then regularized in z. ") +
    std::string(" The parameter -p controls the polynomial degree. -p 0 means no regularization.") +
    std::string("Implemented by B. Avants and conceived by Julien Cohen-Adad.\n") + std::string("Outputs: \n\n") +
    std::string(" OutputPrefixTxTy_poly.csv: polynomial fit to Tx & Ty \n") +
    std::string(" OutputPrefix.nii.gz: transformed image \n") + std::string("Example call: \n\n") +
    std::string(" antsSliceRegularizedRegistration -p 4 --output [OutputPrefix,OutputPrefix.nii.gz]   ") +
    std::string("--transform Translation[0.1] --metric MI[ fixed.nii.gz, moving.nii.gz , 1 , 16 , Regular , 0.2 ] ") +
    std::string("--iterations 20 --shrinkFactors 1 --smoothingSigmas 0 \n\n");
  parser->SetCommandDescription(commandDescription);
  antsSliceRegularizedRegistrationInitializeCommandLineOptions(parser);

  if (parser->Parse(argc, argv) == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  }

  if (argc < 2 || parser->Convert<bool>(parser->GetOption("help")->GetFunction()->GetName()))
  {
    parser->PrintMenu(std::cout, 5, false);
    if (argc < 2)
    {
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }
  else if (parser->Convert<bool>(parser->GetOption('h')->GetFunction()->GetName()))
  {
    parser->PrintMenu(std::cout, 5, true);
    return EXIT_SUCCESS;
  }

  // Get dimensionality
  unsigned int dimension = 3;
  using RealType = double;
  using TranslationTransformType = itk::TranslationTransform<RealType, 2>;
  using EulerTransformType = itk::Euler2DTransform<RealType>;
  using SimilarityTransformType = itk::Similarity2DTransform<RealType>;

  itk::ants::CommandLineParser::OptionType::Pointer transformOption = parser->GetOption("transform");
  if (!(transformOption && transformOption->GetNumberOfFunctions()))
  {
    std::cerr << "No transformations are specified." << std::endl;
    return EXIT_FAILURE;
  }
  std::string whichTransform = transformOption->GetFunction(0)->GetName();
  ConvertToLowerCase(whichTransform);
  if (std::strcmp(whichTransform.c_str(), "translation") == 0)
  {
    switch (dimension)
    {
      case 3:
      {
        return ants_slice_regularized_registration<3, TranslationTransformType>(parser);
      }
      break;
      default:
        std::cerr << "Unsupported dimension" << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }
  else if (std::strcmp(whichTransform.c_str(), "rigid") == 0)
  {
    switch (dimension)
    {
      case 3:
      {
        return ants_slice_regularized_registration<3, EulerTransformType>(parser);
      }
      break;
      default:
        std::cerr << "Unsupported dimension" << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }
  else if (std::strcmp(whichTransform.c_str(), "similarity") == 0)
  {
    switch (dimension)
    {
      case 3:
      {
        return ants_slice_regularized_registration<3, SimilarityTransformType>(parser);
      }
      break;
      default:
        std::cerr << "Unsupported dimension" << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }

  return EXIT_SUCCESS;
}

} // namespace ants
