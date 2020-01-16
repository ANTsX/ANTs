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
#include "itkImageMomentsCalculator.h"
#include "itkImageToHistogramFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkTransformToDisplacementFieldFilter.h"
#include "itkIdentityTransform.h"

#include "itkAffineTransform.h"
#include "itkBSplineTransform.h"
#include "itkBSplineSmoothingOnUpdateDisplacementFieldTransform.h"
#include "itkCompositeTransform.h"
#include "itkGaussianSmoothingOnUpdateDisplacementFieldTransform.h"
#include "itkIdentityTransform.h"
#include "itkEuler2DTransform.h"
#include "itkEuler3DTransform.h"
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

// Headers for interpolating functions (to support the --interpolation choice)
#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"

#include <sstream>

namespace ants
{
/** \class antsRegistrationCommandIterationUpdate
 *  \brief change parameters between iterations of registration
 */
template <typename TFilter>
class antsRegistrationCommandIterationUpdate : public itk::Command
{
public:
  typedef antsRegistrationCommandIterationUpdate Self;
  typedef itk::Command                           Superclass;
  typedef itk::SmartPointer<Self>                Pointer;
  itkNewMacro( Self );
protected:
  antsRegistrationCommandIterationUpdate()
  {
    this->m_LogStream = &std::cout;
  }

public:

  void Execute(itk::Object *caller, const itk::EventObject & event) override
  {
    Execute( (const itk::Object *) caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event) override
  {
    auto * filter = const_cast<TFilter *>( dynamic_cast<const TFilter *>( object ) );

    unsigned int currentLevel = 0;

    if( typeid( event ) == typeid( itk::IterationEvent ) )
      {
      currentLevel = filter->GetCurrentLevel() + 1;
      }
    if( currentLevel < this->m_NumberOfIterations.size() )
      {
      typename TFilter::ShrinkFactorsPerDimensionContainerType shrinkFactors = filter->GetShrinkFactorsPerDimension(
          currentLevel );
      typename TFilter::SmoothingSigmasArrayType smoothingSigmas = filter->GetSmoothingSigmasPerLevel();
      typename TFilter::TransformParametersAdaptorsContainerType adaptors =
        filter->GetTransformParametersAdaptorsPerLevel();

      this->Logger() << "  Current level = " << currentLevel << std::endl;
      this->Logger() << "    number of iterations = " << this->m_NumberOfIterations[currentLevel] << std::endl;
      this->Logger() << "    shrink factors = " << shrinkFactors << std::endl;
      this->Logger() << "    smoothing sigmas = " << smoothingSigmas[currentLevel] << std::endl;
      this->Logger() << "    required fixed parameters = " << adaptors[currentLevel]->GetRequiredFixedParameters()
                     << std::endl;

      typedef itk::ConjugateGradientLineSearchOptimizerv4 GradientDescentOptimizerType;
      auto * optimizer = reinterpret_cast<GradientDescentOptimizerType *>( filter->GetModifiableOptimizer() );
      optimizer->SetNumberOfIterations( this->m_NumberOfIterations[currentLevel] );
      optimizer->SetMinimumConvergenceValue( 1.e-7 );
      optimizer->SetConvergenceWindowSize( 10 );
      optimizer->SetLowerLimit( 0 );
      optimizer->SetUpperLimit( 2 );
      optimizer->SetEpsilon( 0.1 );
      }
  }

  void SetNumberOfIterations( const std::vector<unsigned int> & iterations )
  {
    this->m_NumberOfIterations = iterations;
  }

  void SetLogStream(std::ostream & logStream)
  {
    this->m_LogStream = &logStream;
  }

private:
  std::ostream & Logger() const
  {
    return *m_LogStream;
  }

  std::vector<unsigned int> m_NumberOfIterations;
  std::ostream *            m_LogStream;
};

template <typename T>
inline std::string ants_moco_to_string(const T& t)
{
  std::stringstream ss;

  ss << t;
  return ss.str();
}

template <typename ImageType>
typename ImageType::Pointer PreprocessImage( ImageType * inputImage,
                                             typename ImageType::PixelType lowerScaleFunction,
                                             typename ImageType::PixelType upperScaleFunction,
                                             float winsorizeLowerQuantile, float winsorizeUpperQuantile,
                                             ImageType *histogramMatchSourceImage = nullptr )
{
  bool verbose = false;
  typedef itk::Statistics::ImageToHistogramFilter<ImageType>   HistogramFilterType;
  typedef typename HistogramFilterType::InputBooleanObjectType InputBooleanObjectType;
  typedef typename HistogramFilterType::HistogramSizeType      HistogramSizeType;

  HistogramSizeType histogramSize( 1 );
  histogramSize[0] = 256;

  typename InputBooleanObjectType::Pointer autoMinMaxInputObject = InputBooleanObjectType::New();
  autoMinMaxInputObject->Set( true );

  typename HistogramFilterType::Pointer histogramFilter = HistogramFilterType::New();
  histogramFilter->SetInput( inputImage );
  histogramFilter->SetAutoMinimumMaximumInput( autoMinMaxInputObject );
  histogramFilter->SetHistogramSize( histogramSize );
  histogramFilter->SetMarginalScale( 10.0 );
  histogramFilter->Update();

  float lowerFunction = histogramFilter->GetOutput()->Quantile( 0, winsorizeLowerQuantile );
  float upperFunction = histogramFilter->GetOutput()->Quantile( 0, winsorizeUpperQuantile );
  typedef itk::IntensityWindowingImageFilter<ImageType, ImageType> IntensityWindowingImageFilterType;

  typename IntensityWindowingImageFilterType::Pointer windowingFilter = IntensityWindowingImageFilterType::New();
  windowingFilter->SetInput( inputImage );
  windowingFilter->SetWindowMinimum( lowerFunction );
  windowingFilter->SetWindowMaximum( upperFunction );
  windowingFilter->SetOutputMinimum( lowerScaleFunction );
  windowingFilter->SetOutputMaximum( upperScaleFunction );
  windowingFilter->Update();

  typename ImageType::Pointer outputImage = nullptr;
  if( histogramMatchSourceImage )
    {
    typedef itk::HistogramMatchingImageFilter<ImageType, ImageType> HistogramMatchingFilterType;
    typename HistogramMatchingFilterType::Pointer matchingFilter = HistogramMatchingFilterType::New();
    matchingFilter->SetInput( windowingFilter->GetOutput() );
    matchingFilter->SetReferenceImage( histogramMatchSourceImage );
    matchingFilter->SetNumberOfHistogramLevels( 256 );
    matchingFilter->SetNumberOfMatchPoints( 12 );
    matchingFilter->ThresholdAtMeanIntensityOn();
    matchingFilter->Update();

    outputImage = matchingFilter->GetOutput();
    outputImage->Update();
    outputImage->DisconnectPipeline();

    typedef itk::MinimumMaximumImageCalculator<ImageType> CalculatorType;
    typename CalculatorType::Pointer calc = CalculatorType::New();
    calc->SetImage( inputImage );
    calc->ComputeMaximum();
    calc->ComputeMinimum();
    if( itk::Math::abs( calc->GetMaximum() - calc->GetMinimum() ) < static_cast<typename ImageType::PixelType>( 1.e-9 ) )
      {
      if ( verbose ) std::cout << "Warning: bad time point - too little intensity variation" << std::endl;
      return histogramMatchSourceImage;
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
struct ants_moco_index_cmp
  {
  ants_moco_index_cmp(const T _arr) : arr(_arr)
  {
  }

  bool operator()(const size_t a, const size_t b) const
  {
    return arr[a] < arr[b];
  }

  const T arr;
  };

template <typename TFilter>
class CommandIterationUpdate final : public itk::Command
{
public:
  typedef CommandIterationUpdate  Self;
  typedef itk::Command            Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() = default;
public:

  void Execute(itk::Object *caller, const itk::EventObject & event) override
  {
    Execute( (const itk::Object *) caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event) override
  {
    bool verbose = false;
    auto * filter = const_cast<TFilter *>( dynamic_cast<const TFilter *>( object ) );

    if( typeid( event ) != typeid( itk::IterationEvent ) )
      {
      return;
      }

    unsigned int currentLevel = filter->GetCurrentLevel();
    typename TFilter::ShrinkFactorsPerDimensionContainerType shrinkFactors = filter->GetShrinkFactorsPerDimension(
        currentLevel );
    typename TFilter::SmoothingSigmasArrayType smoothingSigmas = filter->GetSmoothingSigmasPerLevel();
    typename TFilter::TransformParametersAdaptorsContainerType adaptors =
      filter->GetTransformParametersAdaptorsPerLevel();

    if ( verbose ) std::cout << "  Current level = " << currentLevel << std::endl;
    if ( verbose ) std::cout << "    number of iterations = " << this->m_NumberOfIterations[currentLevel] << std::endl;
    if ( verbose ) std::cout << "    shrink factor = " << shrinkFactors[currentLevel] << std::endl;
    if ( verbose ) std::cout << "    smoothing sigma = " << smoothingSigmas[currentLevel] << std::endl;
    if ( verbose ) std::cout << "    required fixed parameters = " << adaptors[currentLevel]->GetRequiredFixedParameters()
              << std::endl;

    typedef itk::ConjugateGradientLineSearchOptimizerv4 OptimizerType;
    auto * optimizer = reinterpret_cast<OptimizerType *>( filter->GetModifiableOptimizer() );
    optimizer->SetNumberOfIterations( this->m_NumberOfIterations[currentLevel] );
    optimizer->SetMinimumConvergenceValue( 1.e-7 );
    optimizer->SetConvergenceWindowSize( 10 );
    optimizer->SetLowerLimit( 0 );
    optimizer->SetUpperLimit( 2 );
    optimizer->SetEpsilon( 0.1 );
  }

  void SetNumberOfIterations( std::vector<unsigned int> iterations )
  {
    this->m_NumberOfIterations = iterations;
  }

private:

  std::vector<unsigned int> m_NumberOfIterations;
};

// Transform traits to generalize the rigid transform
//
template <unsigned int ImageDimension>
class RigidTransformTraits
{
// Don't worry about the fact that the default option is the
// affine Transform, that one will not actually be instantiated.
public:
  typedef itk::AffineTransform<double, ImageDimension> TransformType;
};

template <>
class RigidTransformTraits<2>
{
public:
  typedef itk::Euler2DTransform<double> TransformType;
};

template <>
class RigidTransformTraits<3>
{
public:
  // typedef itk::VersorRigid3DTransform<double> TransformType;
  // typedef itk::QuaternionRigidTransform<double>  TransformType;
  typedef itk::Euler3DTransform<double> TransformType;
};

template <unsigned int ImageDimension>
class SimilarityTransformTraits
{
// Don't worry about the fact that the default option is the
// affine Transform, that one will not actually be instantiated.
public:
  typedef itk::AffineTransform<double, ImageDimension> TransformType;
};

template <>
class SimilarityTransformTraits<2>
{
public:
  typedef itk::Similarity2DTransform<double> TransformType;
};

template <>
class SimilarityTransformTraits<3>
{
public:
  typedef itk::Similarity3DTransform<double> TransformType;
};

/*
template <unsigned int ImageDimension>
class CompositeAffineTransformTraits
{
// Don't worry about the fact that the default option is the
// affine Transform, that one will not actually be instantiated.
public:
  typedef itk::AffineTransform<double, ImageDimension> TransformType;
};
template <>
class CompositeAffineTransformTraits<2>
{
public:
  typedef itk::ANTSCenteredAffine2DTransform<double> TransformType;
};
template <>
class CompositeAffineTransformTraits<3>
{
public:
  typedef itk::ANTSAffine3DTransform<double> TransformType;
};
*/

template <typename TImageIn, typename TImageOut>
void
AverageTimeImages( typename TImageIn::Pointer image_in,  typename TImageOut::Pointer image_avg,
                   std::vector<unsigned int> timelist )
{
  bool verbose = false;
  typedef TImageIn  ImageType;
  typedef TImageOut OutImageType;
  enum { ImageDimension = ImageType::ImageDimension };
  typedef itk::ImageRegionIteratorWithIndex<OutImageType> Iterator;
  image_avg->FillBuffer(0);
  unsigned int timedims = image_in->GetLargestPossibleRegion().GetSize()[ImageDimension - 1];
  if( timelist.empty() )
    {
    for( unsigned int timedim = 0; timedim < timedims; timedim++ )
      {
      timelist.push_back(timedim);
      }
    }
  if ( verbose ) std::cout << " averaging with " << timelist.size() << " images of " <<  timedims <<  " timedims " << std::endl;
  Iterator vfIter2(  image_avg, image_avg->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    typename OutImageType::PixelType  fval = 0;
    typename ImageType::IndexType ind;
    typename OutImageType::IndexType spind = vfIter2.GetIndex();
    for(unsigned int & xx : timelist)
      {
      for( unsigned int yy = 0; yy < ImageDimension - 1; yy++ )
        {
        ind[yy] = spind[yy];
        }
      ind[ImageDimension - 1] = xx;
      fval += image_in->GetPixel(ind);
      }
    fval /= static_cast<typename OutImageType::PixelType>( timelist.size() );
    image_avg->SetPixel(spind, fval);
    }
  if ( verbose ) std::cout << " averaging images done " << std::endl;
  return;
}

template <unsigned int ImageDimension>
int ants_motion( itk::ants::CommandLineParser *parser )
{
  unsigned int verbose = 0;
  itk::ants::CommandLineParser::OptionType::Pointer vOption =
    parser->GetOption( "verbose" );
  if( vOption && vOption->GetNumberOfFunctions() )
      {
      verbose = parser->Convert<unsigned int>( vOption->GetFunction( 0 )->GetName() );
      }
  if ( verbose ) std::cout << " verbose " << std::endl;
  // We infer the number of stages by the number of transformations
  // specified by the user which should match the number of metrics.
  unsigned numberOfStages = 0;

  typedef float                                     PixelType;
  typedef double                                    RealType;
  typedef itk::Image<PixelType, ImageDimension>     FixedIOImageType;
  typedef itk::Image<PixelType, ImageDimension>     FixedImageType;
  typedef itk::Image<PixelType, ImageDimension + 1> MovingIOImageType;
  typedef itk::Image<PixelType, ImageDimension + 1> MovingImageType;
  typedef itk::Vector<RealType, ImageDimension+1>     VectorIOType;
  typedef itk::Image<VectorIOType, ImageDimension+1>  DisplacementIOFieldType;
  typedef itk::Vector<RealType, ImageDimension>     VectorType;
  typedef itk::Image<VectorType, ImageDimension>    DisplacementFieldType;
  typedef vnl_matrix<RealType>                      vMatrix;
  vMatrix param_values;
  typedef itk::CompositeTransform<RealType, ImageDimension> CompositeTransformType;
  std::vector<typename CompositeTransformType::Pointer> CompositeTransformVector;
  typedef typename itk::ants::CommandLineParser ParserType;
  typedef typename ParserType::OptionType       OptionType;

  typename OptionType::Pointer averageOption = parser->GetOption( "average-image" );
  if( averageOption && averageOption->GetNumberOfFunctions() )
    {
    typename OptionType::Pointer outputOption = parser->GetOption( "output" );
    if( !outputOption )
      {
      std::cerr << "Output option not specified.  Should be the output average image name." << std::endl;
      return EXIT_FAILURE;
      }
    std::string outputPrefix = outputOption->GetFunction( 0 )->GetParameter( 0 );
    if( outputPrefix.length() < 3 )
      {
      outputPrefix = outputOption->GetFunction( 0 )->GetName();
      }
    std::string fn = averageOption->GetFunction( 0 )->GetName();
    typename MovingIOImageType::Pointer movingImage;
    ReadImage<MovingIOImageType>( movingImage, fn.c_str()  );
    typename FixedIOImageType::Pointer avgImage;
    typedef itk::ExtractImageFilter<MovingIOImageType, FixedIOImageType> ExtractFilterType;
    typename MovingIOImageType::RegionType extractRegion = movingImage->GetLargestPossibleRegion();
    extractRegion.SetSize(ImageDimension, 0);
    typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
    extractFilter->SetInput( movingImage );
    extractFilter->SetDirectionCollapseToSubmatrix();
    if( ImageDimension == 2 )
      {
      extractFilter->SetDirectionCollapseToIdentity();
      }
    unsigned int td = 0;
    extractRegion.SetIndex(ImageDimension, td );
    extractFilter->SetExtractionRegion( extractRegion );
    extractFilter->Update();
    avgImage = extractFilter->GetOutput();
    std::vector<unsigned int> timelist;
    AverageTimeImages<MovingIOImageType, FixedIOImageType>( movingImage, avgImage, timelist );
    if ( verbose ) std::cout << "average out " << outputPrefix <<  std::endl;
    WriteImage<FixedIOImageType>( avgImage, outputPrefix.c_str() );
    return EXIT_SUCCESS;
    }

  typename OptionType::Pointer transformOption = parser->GetOption( "transform" );
  if( transformOption && transformOption->GetNumberOfFunctions() )
    {
    numberOfStages = transformOption->GetNumberOfFunctions();
    }
  else
    {
    std::cerr << "No transformations are specified." << std::endl;
    return EXIT_FAILURE;
    }

  if ( verbose ) std::cout << "Registration using " << numberOfStages << " total stages." << std::endl;

  // Get the interpolator and possible parameters
  std::string whichInterpolator( "linear" );
  typename OptionType::Pointer interpolationOption = parser->GetOption( "interpolation" );
  if( interpolationOption && interpolationOption->GetNumberOfFunctions() )
    {
    whichInterpolator = interpolationOption->GetFunction( 0 )->GetName();
    ConvertToLowerCase( whichInterpolator );
    }

  typedef itk::Image<PixelType, ImageDimension>           ImageType; // Used only for templating interp functions
  typedef itk::InterpolateImageFunction<ImageType, RealType> InterpolatorType;
  typename InterpolatorType::Pointer interpolator = nullptr;

  if( !std::strcmp( whichInterpolator.c_str(), "linear" ) )
    {
    typedef itk::LinearInterpolateImageFunction<ImageType, RealType> LinearInterpolatorType;
    typename LinearInterpolatorType::Pointer linearInterpolator = LinearInterpolatorType::New();
    interpolator = linearInterpolator;
    }
  else if( !std::strcmp( whichInterpolator.c_str(), "nearestneighbor" ) )
    {
    typedef itk::NearestNeighborInterpolateImageFunction<ImageType, RealType> NearestNeighborInterpolatorType;
    typename NearestNeighborInterpolatorType::Pointer nearestNeighborInterpolator = NearestNeighborInterpolatorType::New();
    interpolator = nearestNeighborInterpolator;
    }
  else if( !std::strcmp( whichInterpolator.c_str(), "bspline" ) )
    {
    typedef itk::BSplineInterpolateImageFunction<ImageType, RealType> BSplineInterpolatorType;
    typename BSplineInterpolatorType::Pointer bSplineInterpolator = BSplineInterpolatorType::New();
    if( interpolationOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      unsigned int bsplineOrder = parser->Convert<unsigned int>( interpolationOption->GetFunction( 0 )->GetParameter( 0 ) );
      bSplineInterpolator->SetSplineOrder( bsplineOrder );
      }
    interpolator = bSplineInterpolator;
    }
  else if( !std::strcmp( whichInterpolator.c_str(), "CosineWindowedSinc" ) )
    {
    typedef itk::WindowedSincInterpolateImageFunction
                 <ImageType, 3, itk::Function::CosineWindowFunction<3, RealType, RealType>, itk::ConstantBoundaryCondition< ImageType >, RealType> CosineInterpolatorType;
    typename CosineInterpolatorType::Pointer cosineInterpolator = CosineInterpolatorType::New();
    interpolator = cosineInterpolator;
    }
  else if( !std::strcmp( whichInterpolator.c_str(), "hammingwindowedsinc" ) )
    {
    typedef itk::WindowedSincInterpolateImageFunction
                 <ImageType, 3, itk::Function::HammingWindowFunction<3, RealType, RealType >, itk::ConstantBoundaryCondition< ImageType >, RealType> HammingInterpolatorType;
    typename HammingInterpolatorType::Pointer hammingInterpolator = HammingInterpolatorType::New();
    interpolator = hammingInterpolator;
    }
  else if( !std::strcmp( whichInterpolator.c_str(), "lanczoswindowedsinc" ) )
    {
    typedef itk::WindowedSincInterpolateImageFunction
                 <ImageType, 3, itk::Function::LanczosWindowFunction<3, RealType, RealType>, itk::ConstantBoundaryCondition< ImageType >, RealType > LanczosInterpolatorType;
    typename LanczosInterpolatorType::Pointer lanczosInterpolator = LanczosInterpolatorType::New();
    interpolator = lanczosInterpolator;
    }
  else if( !std::strcmp( whichInterpolator.c_str(), "blackmanwindowedsinc" ) )
    {
    typedef itk::WindowedSincInterpolateImageFunction
                 <ImageType, 3, itk::Function::BlackmanWindowFunction<3, RealType, RealType>, itk::ConstantBoundaryCondition< ImageType >, RealType > BlackmanInterpolatorType;
    typename BlackmanInterpolatorType::Pointer blackmanInterpolator = BlackmanInterpolatorType::New();
    interpolator = blackmanInterpolator;
    }
  else if( !std::strcmp( whichInterpolator.c_str(), "welchwindowedsinc" ) )
    {
    typedef itk::WindowedSincInterpolateImageFunction
                 <ImageType, 3, itk::Function::WelchWindowFunction<3, RealType, RealType>, itk::ConstantBoundaryCondition< ImageType >, RealType > WelchInterpolatorType;
    typename WelchInterpolatorType::Pointer welchInterpolator = WelchInterpolatorType::New();
    interpolator = welchInterpolator;
    }

  typename OptionType::Pointer metricOption = parser->GetOption( "metric" );
  if( !metricOption || metricOption->GetNumberOfFunctions() != numberOfStages  )
    {
    std::cerr << "The number of metrics specified does not match the number of stages." << std::endl;
    return EXIT_FAILURE;
    }

  typename OptionType::Pointer iterationsOption = parser->GetOption( "iterations" );
  if( !iterationsOption || iterationsOption->GetNumberOfFunctions() != numberOfStages  )
    {
    std::cerr << "The number of iteration sets specified does not match the number of stages." << std::endl;
    return EXIT_FAILURE;
    }

  typename OptionType::Pointer shrinkFactorsOption = parser->GetOption( "shrinkFactors" );
  if( !shrinkFactorsOption || shrinkFactorsOption->GetNumberOfFunctions() != numberOfStages  )
    {
    std::cerr << "The number of shrinkFactor sets specified does not match the number of stages." << std::endl;
    return EXIT_FAILURE;
    }

  typename OptionType::Pointer smoothingSigmasOption = parser->GetOption( "smoothingSigmas" );
  if( !smoothingSigmasOption || smoothingSigmasOption->GetNumberOfFunctions() != numberOfStages  )
    {
    std::cerr << "The number of smoothing sigma sets specified does not match the number of stages." << std::endl;
    return EXIT_FAILURE;
    }

  typename OptionType::Pointer outputOption = parser->GetOption( "output" );
  if( !outputOption )
    {
    std::cerr << "Output option not specified." << std::endl;
    return EXIT_FAILURE;
    }
  std::string outputPrefix = outputOption->GetFunction( 0 )->GetParameter( 0 );
  if( outputPrefix.length() < 3 )
    {
    outputPrefix = outputOption->GetFunction( 0 )->GetName();
    }

  unsigned int                                      nimagestoavg = 0;
  itk::ants::CommandLineParser::OptionType::Pointer navgOption = parser->GetOption( "n-images" );
  if( navgOption && navgOption->GetNumberOfFunctions() )
    {
    nimagestoavg = parser->Convert<unsigned int>( navgOption->GetFunction( 0 )->GetName() );
    if ( verbose ) std::cout << " nimagestoavg " << nimagestoavg << std::endl;
    }

  unsigned int writeDisplacementField = 0;
  itk::ants::CommandLineParser::OptionType::Pointer wdopt = parser->GetOption( "write-displacement" );
  if( wdopt && wdopt->GetNumberOfFunctions() )
    {
    writeDisplacementField = parser->Convert<unsigned int>( wdopt->GetFunction( 0 )->GetName() );
    }

  bool                doEstimateLearningRateOnce(false);
  OptionType::Pointer rateOption = parser->GetOption( "use-estimate-learning-rate-once" );
  if( rateOption && rateOption->GetNumberOfFunctions() )
    {
    std::string rateFunction = rateOption->GetFunction( 0 )->GetName();
    ConvertToLowerCase( rateFunction );
    if( rateFunction.compare( "1" ) == 0 || rateFunction.compare( "true" ) == 0 )
      {
      doEstimateLearningRateOnce = true;
      }
    }

  bool                doHistogramMatch( true );
  OptionType::Pointer histogramMatchOption = parser->GetOption( "use-histogram-matching" );
  if( histogramMatchOption && histogramMatchOption->GetNumberOfFunctions() )
    {
    std::string histogramMatchFunction = histogramMatchOption->GetFunction( 0 )->GetName();
    ConvertToLowerCase( histogramMatchFunction );
    if( histogramMatchFunction.compare( "0" ) == 0 || histogramMatchFunction.compare( "false" ) == 0 )
      {
      doHistogramMatch = false;
      }
    }

  // Zero seed means use default behavior: registration randomizer seeds from system time
  // and does not re-seed iterator
  int antsRandomSeed = 0;

  itk::ants::CommandLineParser::OptionType::Pointer randomSeedOption = parser->GetOption( "random-seed" );
  if( randomSeedOption && randomSeedOption->GetNumberOfFunctions() )
    {
    antsRandomSeed = parser->Convert<int>( randomSeedOption->GetFunction(0)->GetName() );
    }
  else
    {
    char* envSeed = getenv( "ANTS_RANDOM_SEED" );

    if ( envSeed != nullptr )
      {
      antsRandomSeed = std::stoi( envSeed );
      }
    }

  unsigned int   nparams = 2;
  itk::TimeProbe totalTimer;
  totalTimer.Start();
  double metricmean = 0;

  typedef itk::AffineTransform<RealType, ImageDimension>                                      AffineTransformType;
  typedef itk::ImageRegistrationMethodv4<FixedImageType, FixedImageType, AffineTransformType> AffineRegistrationType;
  // We iterate backwards because the command line options are stored as a stack (first in last out)
  typename DisplacementIOFieldType::Pointer displacementout = nullptr;
  typename DisplacementIOFieldType::Pointer displacementinv = nullptr;

  for( int currentStage = numberOfStages - 1; currentStage >= 0; currentStage-- )
    {
    if ( verbose ) std::cout << std::endl << "Stage " << numberOfStages - currentStage << std::endl;
    std::stringstream currentStageString;
    currentStageString << currentStage;

    // Get the fixed and moving images

    std::string fixedImageFileName = metricOption->GetFunction( currentStage )->GetParameter(  0 );
    std::string movingImageFileName = metricOption->GetFunction( currentStage )->GetParameter(  1 );
    if ( verbose ) std::cout << "  fixed image: " << fixedImageFileName << std::endl;
    if ( verbose ) std::cout << "  moving image: " << movingImageFileName << std::endl;
    typename FixedImageType::Pointer fixed_time_slice = nullptr;
    typename FixedImageType::Pointer moving_time_slice = nullptr;
    typename FixedIOImageType::Pointer fixedInImage;
    ReadImage<FixedIOImageType>( fixedInImage, fixedImageFileName.c_str() );
    fixedInImage->Update();
    fixedInImage->DisconnectPipeline();
    typename FixedImageType::Pointer fixedImage;
    fixedImage = arCastImage<FixedIOImageType, FixedImageType>( fixedInImage );

    typename MovingIOImageType::Pointer movingInImage;
    typename MovingImageType::Pointer movingImage;
    ReadImage<MovingIOImageType>( movingInImage, movingImageFileName.c_str()  );
    movingInImage->Update();
    movingInImage->DisconnectPipeline();
    movingImage = arCastImage<MovingIOImageType, MovingImageType>( movingInImage );
    unsigned int              timedims = movingImage->GetLargestPossibleRegion().GetSize()[ImageDimension];

    typename MovingIOImageType::Pointer outputImage = MovingIOImageType::New();
    typename MovingIOImageType::RegionType outRegion;
    typename MovingIOImageType::SizeType outSize;
    typename MovingIOImageType::SpacingType outSpacing;
    typename MovingIOImageType::PointType outOrigin;
    typename MovingIOImageType::DirectionType outDirection;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      outSize[d] = fixedImage->GetLargestPossibleRegion().GetSize()[d];
      outSpacing[d] = fixedImage->GetSpacing()[d];
      outOrigin[d] = fixedImage->GetOrigin()[d];
      for( unsigned int e = 0; e < ImageDimension; e++ )
      	{
	      outDirection(e, d) = fixedImage->GetDirection() (e, d);
	      }
      }
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      outDirection(d, ImageDimension) = 0;
      outDirection(ImageDimension, d) = 0;
      }
    outDirection(ImageDimension, ImageDimension) = 1.0;

    outSize[ImageDimension] = timedims;
    outSpacing[ImageDimension] = movingImage->GetSpacing()[ImageDimension];
    outOrigin[ImageDimension] = movingImage->GetOrigin()[ImageDimension];

    outRegion.SetSize( outSize );
    outputImage->SetRegions( outRegion );
    outputImage->SetSpacing( outSpacing );
    outputImage->SetOrigin( outOrigin );
    outputImage->SetDirection( outDirection );
    outputImage->Allocate();
    outputImage->FillBuffer( 0 );


    if ( writeDisplacementField > 0 )
      {
      /** Handle all output: images and displacement fields */
      typedef itk::IdentityTransform<RealType, ImageDimension+1> IdentityIOTransformType;
      typename IdentityIOTransformType::Pointer identityIOTransform = IdentityIOTransformType::New();
      typedef typename itk::TransformToDisplacementFieldFilter<DisplacementIOFieldType, RealType> ConverterType;
      typename ConverterType::Pointer idconverter = ConverterType::New();
      idconverter->SetOutputOrigin( outputImage->GetOrigin() );
      idconverter->SetOutputStartIndex( outputImage->GetBufferedRegion().GetIndex() );
      idconverter->SetSize( outputImage->GetBufferedRegion().GetSize() );
      idconverter->SetOutputSpacing( outputImage->GetSpacing() );
      idconverter->SetOutputDirection( outputImage->GetDirection() );
      idconverter->SetTransform( identityIOTransform );
      idconverter->Update();
      displacementout = idconverter->GetOutput();


      typename ConverterType::Pointer invconverter = ConverterType::New();
      invconverter->SetOutputOrigin( movingInImage->GetOrigin() );
      invconverter->SetOutputStartIndex(
        movingInImage->GetBufferedRegion().GetIndex() );
      invconverter->SetSize( movingInImage->GetBufferedRegion().GetSize() );
      invconverter->SetOutputSpacing( movingInImage->GetSpacing() );
      invconverter->SetOutputDirection( movingInImage->GetDirection() );
      invconverter->SetTransform( identityIOTransform );
      invconverter->Update();
      displacementinv = invconverter->GetOutput();
      }


    // Get the number of iterations and use that information to specify the number of levels

    std::vector<unsigned int> iterations =
      parser->ConvertVector<unsigned int>( iterationsOption->GetFunction( currentStage )->GetName()  );
    unsigned int numberOfLevels = iterations.size();
    if ( verbose ) std::cout << "  number of levels = " << numberOfLevels << std::endl;

    // Get shrink factors

    std::vector<unsigned int> factors =
      parser->ConvertVector<unsigned int>( shrinkFactorsOption->GetFunction( currentStage )->GetName()  );
    typename AffineRegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
    shrinkFactorsPerLevel.SetSize( factors.size() );

    if( factors.size() != numberOfLevels )
      {
      std::cerr << "ERROR:  The number of shrink factors does not match the number of levels." << std::endl;
      return EXIT_FAILURE;
      }
    else
      {
      for( unsigned int n = 0; n < shrinkFactorsPerLevel.Size(); n++ )
        {
        shrinkFactorsPerLevel[n] = factors[n];
        }
      if ( verbose ) std::cout << "  shrink factors per level: " << shrinkFactorsPerLevel << std::endl;
      }

    // Get smoothing sigmas

    std::vector<float> sigmas = parser->ConvertVector<float>( smoothingSigmasOption->GetFunction(
                                                                currentStage )->GetName()  );
    typename AffineRegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
    smoothingSigmasPerLevel.SetSize( sigmas.size() );

    if( sigmas.size() != numberOfLevels )
      {
      std::cerr << "ERROR:  The number of smoothing sigmas does not match the number of levels." << std::endl;
      return EXIT_FAILURE;
      }
    else
      {
      for( unsigned int n = 0; n < smoothingSigmasPerLevel.Size(); n++ )
        {
        smoothingSigmasPerLevel[n] = sigmas[n];
        }
      if ( verbose ) std::cout << "  smoothing sigmas per level: " << smoothingSigmasPerLevel << std::endl;
      }

    // the fixed image is a reference image in 3D while the moving is a 4D image
    // loop over every time point and register image_i+1 to image_i
    //
    // Set up the image metric and scales estimator
    std::vector<unsigned int> timelist;
    std::vector<double>       metriclist;
    for( unsigned int timedim = 0; timedim < timedims; timedim++ )
      {
      timelist.push_back(timedim);
      }
    for( unsigned int timedim = 0; timedim < timedims; timedim++ )
      {
      typename CompositeTransformType::Pointer compositeTransform = nullptr;
      if( currentStage == static_cast<int>(numberOfStages) - 1 )
        {
        compositeTransform = CompositeTransformType::New();
        CompositeTransformVector.push_back(compositeTransform);
        }
      else if( CompositeTransformVector.size() == timedims && !CompositeTransformVector[timedim].IsNull() )
        {
        compositeTransform = CompositeTransformVector[timedim];
        if( timedim == 0 )
          {
          if ( verbose ) std::cout << " use existing transform " << compositeTransform->GetParameters() << std::endl;
          }
        }
      typedef itk::IdentityTransform<RealType, ImageDimension> IdentityTransformType;
      typename IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();
      //
      typedef itk::ExtractImageFilter<MovingImageType, FixedImageType> ExtractFilterType;
      typename MovingImageType::RegionType extractRegion = movingImage->GetLargestPossibleRegion();
      extractRegion.SetSize(ImageDimension, 0);
      bool maptoneighbor = true;
      typename OptionType::Pointer fixedOption =
        parser->GetOption( "useFixedReferenceImage" );
      if( fixedOption && fixedOption->GetNumberOfFunctions() )
        {
        std::string fixedFunction = fixedOption->GetFunction( 0 )->GetName();
        ConvertToLowerCase( fixedFunction );
        if( fixedFunction.compare( "1" ) == 0 || fixedFunction.compare( "true" ) == 0 )
          {
          if( timedim == 0 )
            {
            if ( verbose ) std::cout << "  using fixed reference image for all frames " << std::endl;
            }
          fixed_time_slice = fixedImage;
          extractRegion.SetIndex(ImageDimension, timedim );
          typename ExtractFilterType::Pointer extractFilter2 = ExtractFilterType::New();
          extractFilter2->SetInput( movingImage );
          extractFilter2->SetDirectionCollapseToSubmatrix();
          if( ImageDimension == 2 )
            {
            extractFilter2->SetDirectionCollapseToIdentity();
            }
          extractFilter2->SetExtractionRegion( extractRegion );
          extractFilter2->Update();
          moving_time_slice = extractFilter2->GetOutput();
          maptoneighbor = false;
          }
        }

      if( maptoneighbor )
        {
        extractRegion.SetIndex(ImageDimension, timedim );
        typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
        extractFilter->SetInput( movingImage );
        extractFilter->SetDirectionCollapseToSubmatrix();
        if( ImageDimension == 2 )
          {
          extractFilter->SetDirectionCollapseToIdentity();
          }
        extractFilter->SetExtractionRegion( extractRegion );
        extractFilter->Update();
        fixed_time_slice = extractFilter->GetOutput();
        unsigned int td = timedim + 1;
        if( td > timedims - 1 )
          {
          td = timedims - 1;
          }
        extractRegion.SetIndex(ImageDimension, td );
        typename ExtractFilterType::Pointer extractFilter2 = ExtractFilterType::New();
        extractFilter2->SetInput( movingImage );
        extractFilter2->SetDirectionCollapseToSubmatrix();
        if( ImageDimension == 2 )
          {
          extractFilter->SetDirectionCollapseToIdentity();
          }
        extractFilter2->SetExtractionRegion( extractRegion );
        extractFilter2->Update();
        moving_time_slice = extractFilter2->GetOutput();
        }

      typename FixedImageType::Pointer preprocessFixedImage =
        PreprocessImage<FixedImageType>( fixed_time_slice, 0,
                                         1, 0.001, 0.999,
                                         nullptr );

      if ( verbose ) std::cout << "  use histogram matching " << doHistogramMatch << std::endl;

      typename FixedImageType::Pointer preprocessMovingImage =
        PreprocessImage<FixedImageType>( moving_time_slice,
                                         0, 1,
                                         0.001, 0.999,
                                         preprocessFixedImage );

      typedef itk::ImageToImageMetricv4<FixedImageType, FixedImageType> MetricType;
      typename MetricType::Pointer metric;

      std::string whichMetric = metricOption->GetFunction( currentStage )->GetName();
      ConvertToLowerCase( whichMetric );

      float samplingPercentage = 1.0;
      if( metricOption->GetFunction( 0 )->GetNumberOfParameters() > 5 )
        {
        samplingPercentage = parser->Convert<float>( metricOption->GetFunction( currentStage )->GetParameter(  5 ) );
        }

      std::string samplingStrategy = "";
      if( metricOption->GetFunction( 0 )->GetNumberOfParameters() > 4 )
        {
        samplingStrategy = metricOption->GetFunction( currentStage )->GetParameter(  4 );
        }
      ConvertToLowerCase( samplingStrategy );
      typename AffineRegistrationType::MetricSamplingStrategyType metricSamplingStrategy = AffineRegistrationType::NONE;
      if( std::strcmp( samplingStrategy.c_str(), "random" ) == 0 )
        {
        if( timedim == 0 )
          {
          if ( verbose ) std::cout << "  random sampling (percentage = " << samplingPercentage << ")" << std::endl;
          }
        metricSamplingStrategy = AffineRegistrationType::RANDOM;
        }
      if( std::strcmp( samplingStrategy.c_str(), "regular" ) == 0 )
        {
        if( timedim == 0 )
          {
          if ( verbose ) std::cout << "  regular sampling (percentage = " << samplingPercentage << ")" << std::endl;
          }
        metricSamplingStrategy = AffineRegistrationType::REGULAR;
        }

      if( std::strcmp( whichMetric.c_str(), "cc" ) == 0 )
        {
        auto radiusOption = parser->Convert<unsigned int>( metricOption->GetFunction( currentStage )->GetParameter(  3 ) );

        if( timedim == 0 )
          {
          if ( verbose ) std::cout << "  using the CC metric (radius = " << radiusOption << ")." << std::endl;
          }
        typedef itk::ANTSNeighborhoodCorrelationImageToImageMetricv4<FixedImageType,
          FixedImageType> CorrelationMetricType;
        typename CorrelationMetricType::Pointer correlationMetric = CorrelationMetricType::New();
        typename CorrelationMetricType::RadiusType radius;
        radius.Fill( radiusOption );
        correlationMetric->SetRadius( radius );
        correlationMetric->SetUseMovingImageGradientFilter( false );
        correlationMetric->SetUseFixedImageGradientFilter( false );

        metric = correlationMetric;
        }
      else if( std::strcmp( whichMetric.c_str(), "mi" ) == 0 )
        {
        auto binOption =
          parser->Convert<unsigned int>( metricOption->GetFunction( currentStage )->GetParameter(  3 ) );

        if( timedim == 0 )
          {
          if ( verbose ) std::cout << "  using the Mattes MI metric." << std::endl;
          }
        typedef itk::MattesMutualInformationImageToImageMetricv4<FixedImageType,
                                                                 FixedImageType> MutualInformationMetricType;
        typename MutualInformationMetricType::Pointer mutualInformationMetric = MutualInformationMetricType::New();
        //mutualInformationMetric = mutualInformationMetric;
        mutualInformationMetric->SetNumberOfHistogramBins( binOption );
        mutualInformationMetric->SetUseMovingImageGradientFilter( false );
        mutualInformationMetric->SetUseFixedImageGradientFilter( false );
        metric = mutualInformationMetric;
        }
      else if( std::strcmp( whichMetric.c_str(), "demons" ) == 0 )
        {
        if( timedim == 0 )
          {
          if ( verbose ) std::cout << "  using the Demons metric." << std::endl;
          }
        typedef itk::MeanSquaresImageToImageMetricv4<FixedImageType, FixedImageType> DemonsMetricType;
        typename DemonsMetricType::Pointer demonsMetric = DemonsMetricType::New();
        //demonsMetric = demonsMetric;
        metric = demonsMetric;
        }
      else if( std::strcmp( whichMetric.c_str(), "gc" ) == 0 )
        {
        if( timedim == 0 )
          {
          if ( verbose ) std::cout << "  using the global correlation metric." << std::endl;
          }
        typedef itk::CorrelationImageToImageMetricv4<FixedImageType, FixedImageType> corrMetricType;
        typename corrMetricType::Pointer corrMetric = corrMetricType::New();
        metric = corrMetric;
        if ( verbose ) std::cout << "  global corr metric set " << std::endl;
        }
      else
        {
        std::cerr << "ERROR: Unrecognized image metric: " << whichMetric << std::endl;
        return EXIT_FAILURE;
        }
      metric->SetVirtualDomainFromImage(  fixed_time_slice );

      typedef itk::RegistrationParameterScalesFromPhysicalShift<MetricType> ScalesEstimatorType;
      typename ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();
      scalesEstimator->SetMetric( metric );
      scalesEstimator->SetTransformForward( true );

      auto learningRate = parser->Convert<float>( transformOption->GetFunction( currentStage )->GetParameter(  0 ) );

      typedef itk::ConjugateGradientLineSearchOptimizerv4 OptimizerType;
      OptimizerType::Pointer optimizer = OptimizerType::New();
      optimizer->SetNumberOfIterations( iterations[0] );
      optimizer->SetMinimumConvergenceValue( 1.e-7 );
      optimizer->SetConvergenceWindowSize( 10 );
      optimizer->SetLowerLimit( 0 );
      optimizer->SetUpperLimit( 2 );
      optimizer->SetEpsilon( 0.1 );

      typename OptionType::Pointer scalesOption = parser->GetOption( "useScalesEstimator" );
      if( scalesOption && scalesOption->GetNumberOfFunctions() )
        {
        std::string scalesFunction = scalesOption->GetFunction( 0 )->GetName();
        ConvertToLowerCase( scalesFunction );
        if( scalesFunction.compare( "1" ) == 0 || scalesFunction.compare( "true" ) == 0 )
          {
          if( timedim == 0 )
            {
            if ( verbose ) std::cout << "  employing scales estimator " << std::endl;
            }
          optimizer->SetScalesEstimator( scalesEstimator );
          }
        else
          {
          if( timedim == 0 )
            {
            if ( verbose ) std::cout << "  not employing scales estimator " << scalesFunction << std::endl;
            }
          }
        }
      optimizer->SetMaximumStepSizeInPhysicalUnits( learningRate );
      optimizer->SetDoEstimateLearningRateOnce( doEstimateLearningRateOnce );
      optimizer->SetDoEstimateLearningRateAtEachIteration( !doEstimateLearningRateOnce );
      //    optimizer->SetMaximumNewtonStepSizeInPhysicalUnits(sqrt(small_step)*learningR);

      // Set up the image registration methods along with the transforms
      std::string whichTransform = transformOption->GetFunction( currentStage )->GetName();
      ConvertToLowerCase( whichTransform );

      // initialize with moments
      typedef typename itk::ImageMomentsCalculator<FixedImageType> ImageCalculatorType;
      typename ImageCalculatorType::Pointer calculator1 =
        ImageCalculatorType::New();
      typename ImageCalculatorType::Pointer calculator2 =
        ImageCalculatorType::New();
      calculator1->SetImage(  fixed_time_slice );
      calculator2->SetImage(  moving_time_slice );
      typename ImageCalculatorType::VectorType fixed_center;
      fixed_center.Fill(0);
      typename ImageCalculatorType::VectorType moving_center;
      moving_center.Fill(0);
      try
        {
        calculator1->Compute();
        fixed_center = calculator1->GetCenterOfGravity();
        try
          {
          calculator2->Compute();
          moving_center = calculator2->GetCenterOfGravity();
          }
        catch( ... )
          {
          fixed_center.Fill(0);
          }
        }
      catch( ... )
        {
        // Rcpp::Rcerr << " zero image1 error ";
        }
      typename AffineTransformType::OffsetType trans;
      itk::Point<RealType, ImageDimension> trans2;
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        trans[i] = moving_center[i] - fixed_center[i];
        trans2[i] =  fixed_center[i];
        }
      if( std::strcmp( whichTransform.c_str(), "affine" ) == 0 )
        {
        typename AffineRegistrationType::Pointer affineRegistration = AffineRegistrationType::New();
        if ( antsRandomSeed != 0 )
          {
          affineRegistration->MetricSamplingReinitializeSeed( antsRandomSeed );
          }
        typename AffineTransformType::Pointer affineTransform = AffineTransformType::New();
        affineTransform->SetIdentity();
        affineTransform->SetOffset( trans );
        affineTransform->SetCenter( trans2 );
        nparams = affineTransform->GetNumberOfParameters() + 2;
        metric->SetFixedImage( preprocessFixedImage );
        metric->SetVirtualDomainFromImage( preprocessFixedImage );
        metric->SetMovingImage( preprocessMovingImage );
        metric->SetMovingTransform( affineTransform );
        typename ScalesEstimatorType::ScalesType scales(affineTransform->GetNumberOfParameters() );
        typename MetricType::ParametersType      newparams(  affineTransform->GetParameters() );
        metric->SetParameters( newparams );
        metric->Initialize();
        scalesEstimator->SetMetric(metric);
        scalesEstimator->EstimateScales(scales);
        optimizer->SetScales(scales);
        if( compositeTransform->GetNumberOfTransforms() > 0 )
          {
          affineRegistration->SetMovingInitialTransform( compositeTransform );
          }
        affineRegistration->SetFixedImage( preprocessFixedImage );
        affineRegistration->SetMovingImage( preprocessMovingImage );
        affineRegistration->SetNumberOfLevels( numberOfLevels );
        affineRegistration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
        affineRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
        affineRegistration->SetMetricSamplingStrategy( metricSamplingStrategy );
        affineRegistration->SetMetricSamplingPercentage( samplingPercentage );
        affineRegistration->SetMetric( metric );
        affineRegistration->SetOptimizer( optimizer );

        typedef CommandIterationUpdate<AffineRegistrationType> AffineCommandType;
        typename AffineCommandType::Pointer affineObserver = AffineCommandType::New();
        affineObserver->SetNumberOfIterations( iterations );

        affineRegistration->AddObserver( itk::IterationEvent(), affineObserver );

        try
          {
          if ( verbose ) std::cout << std::endl << "*** Running affine registration ***" << timedim << std::endl << std::endl;
          affineRegistration->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cerr << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }
        compositeTransform->AddTransform( affineRegistration->GetModifiableTransform() );
        // Write out the affine transform
        std::string filename = outputPrefix + std::string("TimeSlice") + ants_moco_to_string<unsigned int>(timedim)
          + std::string( "Affine.txt" );
        typedef itk::TransformFileWriter TransformWriterType;
        typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
        transformWriter->SetInput( affineRegistration->GetOutput()->Get() );
        transformWriter->SetFileName( filename.c_str() );
#if ITK_VERSION_MAJOR >= 5
        transformWriter->SetUseCompression(true);
#endif
        //      transformWriter->Update();
        if( timedim == 0 )
          {
          param_values.set_size(timedims, nparams);
          param_values.fill(0);
          }
        for( unsigned int i = 0; i < nparams - 2; i++ )
          {
          param_values(timedim, i + 2) = affineRegistration->GetOutput()->Get()->GetParameters()[i];
          }
        }
      else if( std::strcmp( whichTransform.c_str(), "rigid" ) == 0 )
        {
        typedef typename RigidTransformTraits<ImageDimension>::TransformType RigidTransformType;
        typename RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
        rigidTransform->SetOffset( trans );
        rigidTransform->SetCenter( trans2 );
        nparams = rigidTransform->GetNumberOfParameters() + 2;
        typedef itk::ImageRegistrationMethodv4<FixedImageType, FixedImageType,
                                               RigidTransformType> RigidRegistrationType;
        typename RigidRegistrationType::Pointer rigidRegistration = RigidRegistrationType::New();
        if ( antsRandomSeed != 0 )
          {
          rigidRegistration->MetricSamplingReinitializeSeed( antsRandomSeed );
          }
        metric->SetFixedImage( preprocessFixedImage );
        metric->SetVirtualDomainFromImage( preprocessFixedImage );
        metric->SetMovingImage( preprocessMovingImage );
        metric->SetMovingTransform( rigidTransform );
        typename ScalesEstimatorType::ScalesType
          scales(  rigidTransform->GetNumberOfParameters() );
        typename MetricType::ParametersType
          newparams(  rigidTransform->GetParameters() );
        metric->SetParameters( newparams );
        metric->Initialize();
        scalesEstimator->SetMetric(metric);
        scalesEstimator->EstimateScales(scales);
        optimizer->SetScales(scales);
        rigidRegistration->SetFixedImage( preprocessFixedImage );
        rigidRegistration->SetMovingImage( preprocessMovingImage );
        rigidRegistration->SetNumberOfLevels( numberOfLevels );
        rigidRegistration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
        rigidRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
        rigidRegistration->SetMetric( metric );
        rigidRegistration->SetMetricSamplingStrategy(
          static_cast<typename RigidRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
        rigidRegistration->SetMetricSamplingPercentage( samplingPercentage );
        rigidRegistration->SetOptimizer( optimizer );
        if( compositeTransform->GetNumberOfTransforms() > 0 )
          {
          rigidRegistration->SetMovingInitialTransform( compositeTransform );
          }

        typedef CommandIterationUpdate<RigidRegistrationType> RigidCommandType;
        typename RigidCommandType::Pointer rigidObserver = RigidCommandType::New();
        rigidObserver->SetNumberOfIterations( iterations );
        rigidRegistration->AddObserver( itk::IterationEvent(), rigidObserver );
        try
          {
          if ( verbose ) std::cout << std::endl << "*** Running rigid registration ***" << timedim  << std::endl << std::endl;
          rigidRegistration->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cerr << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }
        compositeTransform->AddTransform( rigidRegistration->GetModifiableTransform() );
        // Write out the rigid transform
        std::string filename = outputPrefix + std::string("TimeSlice") + ants_moco_to_string<unsigned int>(timedim)
          + std::string( "Rigid.txt" );
        typedef itk::TransformFileWriter TransformWriterType;
        typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
        transformWriter->SetInput( rigidRegistration->GetOutput()->Get() );
        transformWriter->SetFileName( filename.c_str() );
#if ITK_VERSION_MAJOR >= 5
        transformWriter->SetUseCompression(true);
#endif
        //      transformWriter->Update();
        if( timedim == 0 )
          {
          param_values.set_size(timedims, nparams);
          param_values.fill(0);
          }
        for( unsigned int i = 0; i < nparams - 2; i++ )
          {
          param_values(timedim, i + 2) = rigidRegistration->GetOutput()->Get()->GetParameters()[i];
          }
        }
      else if( std::strcmp( whichTransform.c_str(),
                            "gaussiandisplacementfield" ) == 0 ||  std::strcmp( whichTransform.c_str(), "gdf" ) == 0 )
        {
        RealType sigmaForUpdateField = parser->Convert<float>( transformOption->GetFunction(
                                                                 currentStage )->GetParameter(  1 ) );
        RealType sigmaForTotalField = parser->Convert<float>( transformOption->GetFunction(
                                                                currentStage )->GetParameter(  2 ) );
        const unsigned int VImageDimension = ImageDimension;
        VectorType zeroVector( 0.0 );
        // ORIENTATION ALERT: Original code set image size to
        // fixedImage buffered region, & if fixedImage BufferedRegion
        // != LargestPossibleRegion, this code would be wrong.
        typename DisplacementFieldType::Pointer displacementField = AllocImage<DisplacementFieldType>(
            preprocessFixedImage, zeroVector );
        typedef itk::GaussianSmoothingOnUpdateDisplacementFieldTransform<RealType,
                                                                         VImageDimension>
          GaussianDisplacementFieldTransformType;

        typedef itk::ImageRegistrationMethodv4<FixedImageType, FixedImageType,
                                               GaussianDisplacementFieldTransformType>
          DisplacementFieldRegistrationType;
        typename DisplacementFieldRegistrationType::Pointer displacementFieldRegistration =
          DisplacementFieldRegistrationType::New();

        typename GaussianDisplacementFieldTransformType::Pointer outputDisplacementFieldTransform =
                                                                      displacementFieldRegistration->GetModifiableTransform();

        // Create the transform adaptors

        typedef itk::GaussianSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor<GaussianDisplacementFieldTransformType> DisplacementFieldTransformAdaptorType;
        typename DisplacementFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;

        // Extract parameters
        outputDisplacementFieldTransform->SetGaussianSmoothingVarianceForTheUpdateField( sigmaForUpdateField );
        outputDisplacementFieldTransform->SetGaussianSmoothingVarianceForTheTotalField( sigmaForTotalField );
        outputDisplacementFieldTransform->SetDisplacementField( displacementField );
        for( unsigned int level = 0; level < numberOfLevels; level++ )
          {
          typedef itk::ShrinkImageFilter<DisplacementFieldType, DisplacementFieldType> ShrinkFilterType;
          typename ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
          shrinkFilter->SetShrinkFactors( shrinkFactorsPerLevel[level] );
          shrinkFilter->SetInput( displacementField );
          shrinkFilter->Update();
          typename DisplacementFieldTransformAdaptorType::Pointer fieldTransformAdaptor =
            DisplacementFieldTransformAdaptorType::New();
          fieldTransformAdaptor->SetRequiredSpacing( shrinkFilter->GetOutput()->GetSpacing() );
          fieldTransformAdaptor->SetRequiredSize( shrinkFilter->GetOutput()->GetBufferedRegion().GetSize() );
          fieldTransformAdaptor->SetRequiredDirection( shrinkFilter->GetOutput()->GetDirection() );
          fieldTransformAdaptor->SetRequiredOrigin( shrinkFilter->GetOutput()->GetOrigin() );
          fieldTransformAdaptor->SetTransform( outputDisplacementFieldTransform );
          adaptors.push_back( fieldTransformAdaptor.GetPointer() );
          }
        displacementFieldRegistration->SetFixedImage( 0, preprocessFixedImage );
        displacementFieldRegistration->SetMovingImage( 0, preprocessMovingImage );
        displacementFieldRegistration->SetMetric( metric );
        displacementFieldRegistration->SetNumberOfLevels( numberOfLevels );
        displacementFieldRegistration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
        displacementFieldRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
        displacementFieldRegistration->SetSmoothingSigmasAreSpecifiedInPhysicalUnits( false );
        displacementFieldRegistration->SetMetricSamplingStrategy(
          static_cast<typename DisplacementFieldRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
        displacementFieldRegistration->SetMetricSamplingPercentage( samplingPercentage );
        displacementFieldRegistration->SetOptimizer( optimizer );
        displacementFieldRegistration->SetTransformParametersAdaptorsPerLevel( adaptors );
        if( compositeTransform->GetNumberOfTransforms() > 0 )
          {
          displacementFieldRegistration->SetMovingInitialTransform( compositeTransform );
          }
        try
          {
          displacementFieldRegistration->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cerr << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }
        compositeTransform->AddTransform( outputDisplacementFieldTransform );
        if( timedim == 0 )
          {
          param_values.set_size(timedims, nparams);
          param_values.fill(0);
          }
        }
      else if( std::strcmp( whichTransform.c_str(),
                            "SyN" ) == 0 ||  std::strcmp( whichTransform.c_str(), "syn" ) == 0 )
        {
        RealType sigmaForUpdateField = parser->Convert<float>( transformOption->GetFunction(
                                                                 currentStage )->GetParameter(  1 ) );
        RealType sigmaForTotalField = parser->Convert<float>( transformOption->GetFunction(
                                                                currentStage )->GetParameter(  2 ) );
        const unsigned int VImageDimension = ImageDimension;
        VectorType zeroVector( 0.0 );

        typename DisplacementFieldType::Pointer displacementField = AllocImage<DisplacementFieldType>(
            preprocessFixedImage, zeroVector );

        typename DisplacementFieldType::Pointer inverseDisplacementField = AllocImage<DisplacementFieldType>(
            preprocessFixedImage, zeroVector );

        typedef itk::DisplacementFieldTransform<RealType, VImageDimension> DisplacementFieldTransformType;
        typedef itk::SyNImageRegistrationMethod<FixedImageType, FixedImageType,
                                                DisplacementFieldTransformType> DisplacementFieldRegistrationType;
        typename DisplacementFieldRegistrationType::Pointer displacementFieldRegistration =
          DisplacementFieldRegistrationType::New();

        typename DisplacementFieldTransformType::Pointer outputDisplacementFieldTransform =
                                                                  displacementFieldRegistration->GetModifiableTransform();

        // Create the transform adaptors

        typedef itk::DisplacementFieldTransformParametersAdaptor<DisplacementFieldTransformType>
          DisplacementFieldTransformAdaptorType;
        typename DisplacementFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;
        // Create the transform adaptors
        // For the gaussian displacement field, the specified variances are in image spacing terms
        // and, in normal practice, we typically don't change these values at each level.  However,
        // if the user wishes to add that option, they can use the class
        // GaussianSmoothingOnUpdateDisplacementFieldTransformAdaptor
        for( unsigned int level = 0; level < numberOfLevels; level++ )
          {
          // TODO:
          // We use the shrink image filter to calculate the fixed parameters of the virtual
          // domain at each level.  To speed up calculation and avoid unnecessary memory
          // usage, we could calculate these fixed parameters directly.
          typedef itk::ShrinkImageFilter<DisplacementFieldType, DisplacementFieldType> ShrinkFilterType;
          typename ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
          shrinkFilter->SetShrinkFactors( shrinkFactorsPerLevel[level] );
          shrinkFilter->SetInput( displacementField );
          shrinkFilter->Update();

          typename DisplacementFieldTransformAdaptorType::Pointer fieldTransformAdaptor =
            DisplacementFieldTransformAdaptorType::New();
          fieldTransformAdaptor->SetRequiredSpacing( shrinkFilter->GetOutput()->GetSpacing() );
          fieldTransformAdaptor->SetRequiredSize( shrinkFilter->GetOutput()->GetBufferedRegion().GetSize() );
          fieldTransformAdaptor->SetRequiredDirection( shrinkFilter->GetOutput()->GetDirection() );
          fieldTransformAdaptor->SetRequiredOrigin( shrinkFilter->GetOutput()->GetOrigin() );
          fieldTransformAdaptor->SetTransform( outputDisplacementFieldTransform );

          adaptors.push_back( fieldTransformAdaptor.GetPointer() );
          }

        // Extract parameters
        typename DisplacementFieldRegistrationType::NumberOfIterationsArrayType numberOfIterationsPerLevel;
        numberOfIterationsPerLevel.SetSize( numberOfLevels );
        if( timedim == 0 )
          {
          if ( verbose ) std::cout << "SyN iterations:";
          }
        for( unsigned int d = 0; d < numberOfLevels; d++ )
          {
          numberOfIterationsPerLevel[d] = iterations[d]; // currentStageIterations[d];
          if( timedim == 0 )
            {
            if ( verbose ) std::cout << numberOfIterationsPerLevel[d] << " ";
            }
          }
        if( timedim == 0 )
          {
          if ( verbose ) std::cout << std::endl;
          }

        const RealType varianceForUpdateField = sigmaForUpdateField;
        const RealType varianceForTotalField = sigmaForTotalField;
        displacementFieldRegistration->SetFixedImage( 0, preprocessFixedImage );
        displacementFieldRegistration->SetMovingImage( 0, preprocessMovingImage );
        displacementFieldRegistration->SetMetric( metric );

        if( compositeTransform->GetNumberOfTransforms() > 0 )
          {
          displacementFieldRegistration->SetMovingInitialTransform( compositeTransform );
          }
        displacementFieldRegistration->SetDownsampleImagesForMetricDerivatives( true );
        displacementFieldRegistration->SetAverageMidPointGradients( false );
        displacementFieldRegistration->SetNumberOfLevels( numberOfLevels );
        displacementFieldRegistration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
        displacementFieldRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
        displacementFieldRegistration->SetSmoothingSigmasAreSpecifiedInPhysicalUnits( false );
        displacementFieldRegistration->SetLearningRate( learningRate );
        displacementFieldRegistration->SetConvergenceThreshold( 1.e-8 );
        displacementFieldRegistration->SetConvergenceWindowSize( 10 );
        displacementFieldRegistration->SetNumberOfIterationsPerLevel( numberOfIterationsPerLevel );
        displacementFieldRegistration->SetTransformParametersAdaptorsPerLevel( adaptors );
        displacementFieldRegistration->SetGaussianSmoothingVarianceForTheUpdateField( varianceForUpdateField );
        displacementFieldRegistration->SetGaussianSmoothingVarianceForTheTotalField( varianceForTotalField );
        outputDisplacementFieldTransform->SetDisplacementField( displacementField );
        outputDisplacementFieldTransform->SetInverseDisplacementField( inverseDisplacementField );
        try
          {
          displacementFieldRegistration->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cerr << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }
        // Add calculated transform to the composite transform
        compositeTransform->AddTransform( outputDisplacementFieldTransform );
        if( timedim == 0 )
          {
          param_values.set_size(timedims, nparams);
          param_values.fill(0);
          }
        }
      else
        {
        std::cerr << "ERROR:  Unrecognized transform option - " << whichTransform << std::endl;
        return EXIT_FAILURE;
        }
      if( currentStage == static_cast<int>(numberOfStages) - 1 )
        {
        param_values(timedim, 1) = metric->GetValue();
        }
      metriclist.push_back( param_values(timedim, 1) );
      metricmean +=  param_values(timedim, 1) / ( double ) timedims;
      // resample the moving image and then put it in its place
      typedef itk::ResampleImageFilter<FixedImageType, FixedImageType> ResampleFilterType;
      typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
      resampler->SetTransform( compositeTransform );
      resampler->SetInput( moving_time_slice );
      resampler->SetOutputParametersFromImage( fixed_time_slice );
      resampler->SetDefaultPixelValue( 0 );
      resampler->SetInterpolator( interpolator );
      resampler->Update();
      if ( verbose ) std::cout << " done resampling timepoint : " << timedim << std::endl;

      /** Here, we put the resampled 3D image into the 4D volume */
      typedef itk::ImageRegionIteratorWithIndex<FixedImageType> Iterator;
      Iterator vfIter2(  resampler->GetOutput(), resampler->GetOutput()->GetLargestPossibleRegion() );
      for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
        {
        typename FixedImageType::PixelType  fval = vfIter2.Get();
        typename MovingImageType::IndexType ind;
        for( unsigned int xx = 0; xx < ImageDimension; xx++ )
          {
          ind[xx] = vfIter2.GetIndex()[xx];
          }
        unsigned int tdim = timedim;
        if( tdim > ( timedims - 1 ) )
          {
          tdim = timedims - 1;
          }
        ind[ImageDimension] = tdim;
        outputImage->SetPixel(ind, fval);
        }
      if ( writeDisplacementField > 0 )
        {
        typedef typename
        itk::TransformToDisplacementFieldFilter<DisplacementFieldType, RealType>
          ConverterType;
        typename ConverterType::Pointer converter = ConverterType::New();
        converter->SetOutputOrigin( fixed_time_slice->GetOrigin() );
        converter->SetOutputStartIndex(
          fixed_time_slice->GetBufferedRegion().GetIndex() );
        converter->SetSize( fixed_time_slice->GetBufferedRegion().GetSize() );
        converter->SetOutputSpacing( fixed_time_slice->GetSpacing() );
        converter->SetOutputDirection( fixed_time_slice->GetDirection() );
        converter->SetTransform( compositeTransform );
        converter->Update();
        /** Here, we put the 3d tx into a 4d displacement field */
        for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
          {
          VectorType vec =
            converter->GetOutput()->GetPixel( vfIter2.GetIndex() );
          VectorIOType vecout;
          vecout.Fill( 0 );
          typename MovingIOImageType::IndexType ind;
          for( unsigned int xx = 0; xx < ImageDimension; xx++ )
            {
            ind[xx] = vfIter2.GetIndex()[xx];
            vecout[xx] = vec[xx];
            }
          unsigned int tdim = timedim;
          if( tdim > ( timedims - 1 ) )
            {
            tdim = timedims - 1;
            }
          ind[ImageDimension] = tdim;
          displacementout->SetPixel( ind, vecout );
          }
#
        typename ConverterType::Pointer converter2 = ConverterType::New();
        converter2->SetOutputOrigin( moving_time_slice->GetOrigin() );
        converter2->SetOutputStartIndex(
          moving_time_slice->GetBufferedRegion().GetIndex() );
        converter2->SetSize( moving_time_slice->GetBufferedRegion().GetSize() );
        converter2->SetOutputSpacing( moving_time_slice->GetSpacing() );
        converter2->SetOutputDirection( moving_time_slice->GetDirection() );
        converter2->SetTransform( compositeTransform->GetInverseTransform() );
        converter2->Update();
        /** Here, we put the 3d tx into a 4d displacement field */
        Iterator vfIterInv(  moving_time_slice,
          moving_time_slice->GetLargestPossibleRegion() );
        for(  vfIterInv.GoToBegin(); !vfIterInv.IsAtEnd(); ++vfIterInv )
          {
          VectorType vec =
            converter2->GetOutput()->GetPixel( vfIterInv.GetIndex() );
          VectorIOType vecout;
          vecout.Fill( 0 );
          typename MovingIOImageType::IndexType ind;
          for( unsigned int xx = 0; xx < ImageDimension; xx++ )
            {
            ind[xx] = vfIterInv.GetIndex()[xx];
            vecout[xx] = vec[xx];
            }
          unsigned int tdim = timedim;
          if( tdim > ( timedims - 1 ) )
            {
            tdim = timedims - 1;
            }
          ind[ImageDimension] = tdim;
          displacementinv->SetPixel( ind, vecout );
          }
        }
      }
    if( outputOption && outputOption->GetFunction( 0 )->GetNumberOfParameters() > 1  && currentStage == 0 )
      {
      std::string fileName = outputOption->GetFunction( 0 )->GetParameter( 1 );
      if( outputPrefix.length() < 3 )
        {
        outputPrefix = outputOption->GetFunction( 0 )->GetName();
        }
      if ( verbose ) std::cout << "motion corrected out " << fileName <<  std::endl;
      WriteImage<MovingIOImageType>( outputImage, fileName.c_str()  );
      }
    if( outputOption && outputOption->GetFunction( 0 )->GetNumberOfParameters() > 2 && outputImage && currentStage ==
        0 )
      {
      std::string fileName = outputOption->GetFunction( 0 )->GetParameter( 2 );
      typename FixedIOImageType::Pointer avgImage;
      typedef itk::ExtractImageFilter<MovingImageType, FixedIOImageType> ExtractFilterType;
      typename MovingImageType::RegionType extractRegion = movingImage->GetLargestPossibleRegion();
      extractRegion.SetSize(ImageDimension, 0);
      typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
      extractFilter->SetInput( movingImage );
      extractFilter->SetDirectionCollapseToSubmatrix();
      if( ImageDimension == 2 )
        {
        extractFilter->SetDirectionCollapseToIdentity();
        }
      unsigned int td = 0;
      extractRegion.SetIndex(ImageDimension, td );
      extractFilter->SetExtractionRegion( extractRegion );
      extractFilter->Update();
      avgImage = extractFilter->GetOutput();
      std::sort(timelist.begin(), timelist.end(), ants_moco_index_cmp<std::vector<double> &>(metriclist) );
      if( nimagestoavg == 0 )
        {
        nimagestoavg = timelist.size();
        }
      std::vector<unsigned int> timelistsort;
      for( unsigned int i = 0; i < nimagestoavg; i++ )
        {
        if( i < timelist.size() )
          {
          timelistsort.push_back(timelist[i]);
          }
        if ( verbose ) std::cout << " i^th value " << i << "  is " << metriclist[timelist[i]] << std::endl;
        }
      AverageTimeImages<MovingIOImageType, FixedIOImageType>( outputImage, fixed_time_slice, timelistsort );
      if ( verbose ) std::cout << " write average post " << fileName << std::endl;
      WriteImage<FixedIOImageType>( fixed_time_slice, fileName.c_str() );
      }
    }
  if ( writeDisplacementField > 0 )
    {
    std::string dfn = outputPrefix + std::string("Warp.nii.gz");
    WriteImage<DisplacementIOFieldType>( displacementout, dfn.c_str()  );
    dfn = outputPrefix + std::string("InverseWarp.nii.gz");
    WriteImage<DisplacementIOFieldType>( displacementinv, dfn.c_str()  );
    }

  totalTimer.Stop();
  if ( verbose ) std::cout << std::endl << "Total elapsed time: " << totalTimer.GetMean() << " averagemetric " << metricmean
            << std::endl;
    {
    std::vector<std::string> ColumnHeaders;
    std::string              colname;
    colname = std::string("MetricPre");
    ColumnHeaders.push_back( colname );
    colname = std::string("MetricPost");
    ColumnHeaders.push_back( colname );
    for( unsigned int nv = 2; nv < nparams; nv++ )
      {
      std::string _colname = std::string("MOCOparam") + ants_moco_to_string<unsigned int>(nv - 2);
      ColumnHeaders.push_back( _colname );
      }
    typedef itk::CSVNumericObjectFileWriter<double, 1, 1> WriterType;
    WriterType::Pointer writer = WriterType::New();
    std::string         fnmp;
    if ( verbose ) std::cout << " get motion corr params " << outputPrefix << std::endl;
    if( outputPrefix[0] == '0' && outputPrefix[1] == 'x' )
      {
      void* ptr;
      std::sscanf(outputPrefix.c_str(), "%p", (void * *)&ptr);
      //      std::stringstream strstream;
      //      strstream << outputPrefix;
      //      void* ptr;
      //      strstream >> ptr;
      ( static_cast<std::pair<std::vector<std::string>, vnl_matrix<float> > *>( ptr ) )->first = ColumnHeaders;
      ( static_cast<std::pair<std::vector<std::string>, vnl_matrix<double> > *>( ptr ) )->second = param_values;
      if ( verbose ) std::cout << "motion-correction params written" << std::endl;
      }
    else
      {
      fnmp = outputPrefix + std::string("MOCOparams.csv");
      if ( verbose ) std::cout << " write " << fnmp << std::endl;
      writer->SetFileName( fnmp.c_str() );
      writer->SetColumnHeaders(ColumnHeaders);
      writer->SetInput( &param_values );
      writer->Write();
      }
    }

  return EXIT_SUCCESS;
}

void antsMotionCorrInitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

  {
  std::string description =
    std::string( "This option forces the image to be treated as a specified-" )
    + std::string( "dimensional image.  If not specified, the program tries to " )
    + std::string( "infer the dimensionality from the input image." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "dimensionality" );
  option->SetShortName( 'd' );
  option->SetUsageOption( 0, "2/3" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string(
      "turn on the option that lets you estimate the learning rate step size only at the beginning of each level.  * useful as a second stage of fine-scale registration." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "use-estimate-learning-rate-once" );
  option->SetShortName( 'l' );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "This option sets the number of images to use to construct the template image.");
  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "n-images" );
  option->SetShortName( 'n' );
  option->SetUsageOption( 0, "10" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Four image metrics are available--- " )
    + std::string( "GC : global correlation, CC:  ANTS neighborhood cross correlation, MI:  Mutual information, and " )
    + std::string( "Demons:  Thirion's Demons (modified mean-squares). " )
    + std::string( "Note that the metricWeight is currently not used.  " )
    + std::string( "Rather, it is a temporary place holder until multivariate metrics " )
    + std::string( "are available for a single stage. " )
    + std::string( "The fixed image should be a single time point (eg the average of the time series). " )
    + std::string( "By default, this image is not used, the fixed image for correction of each volume is the preceding volume " )
    + std::string( "in the time series. See below for the option to use a fixed reference image for all volumes. ");

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "metric" );
  option->SetShortName( 'm' );
  option->SetUsageOption(
    0,
    "CC[fixedImage,movingImage,metricWeight,radius,<samplingStrategy={Regular,Random}>,<samplingPercentage=[0,1]>]" );
  option->SetUsageOption(
    1,
    "MI[fixedImage,movingImage,metricWeight,numberOfBins,<samplingStrategy={Regular,Random}>,<samplingPercentage=[0,1]>]" );
  option->SetUsageOption(
    2,

    "Demons[fixedImage,movingImage,metricWeight,radius,<samplingStrategy={Regular,Random}>,<samplingPercentage=[0,1]>]" );
  option->SetUsageOption(
    3,
    "GC[fixedImage,movingImage,metricWeight,radius,<samplingStrategy={Regular,Random}>,<samplingPercentage=[0,1]>]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string(
      "use a fixed reference image to correct all volumes, instead of correcting each image ")
      + std::string( "to the prior volume in the time series." );
  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "useFixedReferenceImage" );
  option->SetShortName( 'u' );
  option->SetUsageOption( 0, "(0)/1" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string         description = std::string( "use the scale estimator to control optimization." );
  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "useScalesEstimator" );
  option->SetShortName( 'e' );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Several transform options are available.  The gradientStep or" )
    + std::string( "learningRate characterizes the gradient descent optimization and is scaled appropriately " )
    + std::string( "for each transform using the shift scales estimator.  Subsequent parameters are " )
    + std::string( "transform-specific and can be determined from the usage. " );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "transform" );
  option->SetShortName( 't' );
  option->SetUsageOption( 0, "Affine[gradientStep]" );
  option->SetUsageOption( 1, "Rigid[gradientStep]" );
  option->SetUsageOption(
    2, "GaussianDisplacementField[gradientStep,updateFieldSigmaInPhysicalSpace,totalFieldSigmaInPhysicalSpace]" );
  option->SetUsageOption(
    3, "SyN[gradientStep,updateFieldSigmaInPhysicalSpace,totalFieldSigmaInPhysicalSpace]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Specify the number of iterations at each level." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "iterations" );
  option->SetShortName( 'i' );
  option->SetUsageOption( 0, "MxNx0..." );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Specify the amount of smoothing at each level." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "smoothingSigmas" );
  option->SetShortName( 's' );
  option->SetUsageOption( 0, "MxNx0..." );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string(
      "Specify the shrink factor for the virtual domain (typically the fixed image) at each level." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "shrinkFactors" );
  option->SetShortName( 'f' );
  option->SetUsageOption( 0, "MxNx0..." );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Specify the output transform prefix (output format is .nii.gz )." )
    + std::string( "Optionally, one can choose to warp the moving image to the fixed space and, if the " )
    + std::string( "inverse transform exists, one can also output the warped fixed image." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "output" );
  option->SetShortName( 'o' );
  option->SetUsageOption( 0, "[outputTransformPrefix,<outputWarpedImage>,<outputAverageImage>]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string         description = std::string( "Average the input time series image." );
  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "average-image" );
  option->SetShortName( 'a' );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string         description = std::string( "Write the low-dimensional 3D transforms to a 4D displacement field." );
  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "write-displacement" );
  option->SetShortName( 'w' );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string         description = std::string( "Histogram match the moving images to the reference image." );
  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "use-histogram-matching" );
  option->SetUsageOption( 0, "0/(1)" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Use a fixed seed for random number generation. " )
    + std::string( "By default, the system clock is used to initialize the seeding. " )
    + std::string( "The fixed seed can be any nonzero int value." );
  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "random-seed" );
  option->SetUsageOption( 0, "seedValue" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Several interpolation options are available in ITK. " )
    + std::string( "The above are available (default Linear)." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "interpolation" );
  // n is already in use by --n-images. Unfortunately flag shortname is inconsistent with antsApplyTransforms.
  option->SetShortName( 'p' );
  option->SetUsageOption( 0, "Linear" );
  option->SetUsageOption( 1, "NearestNeighbor" );
  option->SetUsageOption( 2, "BSpline[<order=3>]" );
  option->SetUsageOption( 3, "BlackmanWindowedSinc" );
  option->SetUsageOption( 4, "CosineWindowedSinc" );
  option->SetUsageOption( 5, "WelchWindowedSinc" );
  option->SetUsageOption( 6, "HammingWindowedSinc" );
  option->SetUsageOption( 7, "LanczosWindowedSinc" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Verbose output." );

  OptionType::Pointer option = OptionType::New();
  option->SetShortName( 'v' );
  option->SetLongName( "verbose" );
  option->SetUsageOption( 0, "(0)/1" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Print the help menu (short version)." );

  OptionType::Pointer option = OptionType::New();
  option->SetShortName( 'h' );
  option->SetDescription( description );
  option->AddFunction( std::string( "0" ) );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Print the help menu." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "help" );
  option->SetDescription( description );
  option->AddFunction( std::string( "0" ) );
  parser->AddOption( option );
  }

}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int antsMotionCorr( std::vector<std::string> args, std::ostream * /*out_stream = nullptr */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "antsMotionCorr" );

  int     argc = args.size();
  char* * argv = new char *[args.size() + 1];
  for( unsigned int i = 0; i < args.size(); ++i )
    {
    // allocate space for the string plus a null character
    argv[i] = new char[args[i].length() + 1];
    std::strncpy( argv[i], args[i].c_str(), args[i].length() );
    // place the null character in the end
    argv[i][args[i].length()] = '\0';
    }
  argv[argc] = nullptr;
  // class to automatically cleanup argv upon destruction
  class Cleanup_argv
  {
public:
    Cleanup_argv( char* * argv_, int argc_plus_one_ ) : argv( argv_ ), argc_plus_one( argc_plus_one_ )
    {
    }

    ~Cleanup_argv()
    {
      for( unsigned int i = 0; i < argc_plus_one; ++i )
        {
        delete[] argv[i];
        }
      delete[] argv;
    }

private:
    char* *      argv;
    unsigned int argc_plus_one;
  };
  Cleanup_argv cleanup_argv( argv, argc + 1 );

  // antscout->set_stream( out_stream );

  itk::ants::CommandLineParser::Pointer parser = itk::ants::CommandLineParser::New();

  parser->SetCommand( argv[0] );

  std::string commandDescription = std::string( "antsMotionCorr = motion correction.  This program is a user-level " )
    + std::string( "registration application meant to utilize classes in ITK v4.0 or greater. The user can specify " )
    + std::string( "any number of \"stages\" where a stage consists of a transform; an image metric; " )
    + std::string( " and iterations, shrink factors, and smoothing sigmas for each level. " )
    + std::string(
      " Specialized for 4D time series data: fixed image is 3D, moving image should be the 4D time series. ")
    + std::string( " Fixed image is a reference space or time slice.");
  parser->SetCommandDescription( commandDescription );
  antsMotionCorrInitializeCommandLineOptions( parser );

  if( parser->Parse( argc, argv ) == EXIT_FAILURE )
    {
    return EXIT_FAILURE;
    }

  if( argc < 2 || parser->Convert<bool>( parser->GetOption( "help" )->GetFunction()->GetName() ) )
    {
    parser->PrintMenu( std::cout, 5, false );
    if( argc < 2 )
      {
      return EXIT_FAILURE;
      }
    return EXIT_SUCCESS;
    }
  else if( parser->Convert<bool>( parser->GetOption( 'h' )->GetFunction()->GetName() ) )
    {
    parser->PrintMenu( std::cout, 5, true );
    return EXIT_SUCCESS;
    }

  // Get dimensionality
  unsigned int dimension = 3;

  itk::ants::CommandLineParser::OptionType::Pointer dimOption = parser->GetOption( "dimensionality" );
  if( dimOption && dimOption->GetNumberOfFunctions() )
    {
    dimension = parser->Convert<unsigned int>( dimOption->GetFunction( 0 )->GetName() );
    }
  else
    {
    std::cerr << "Image dimensionality not specified.  See command line option --dimensionality" << std::endl;
    return EXIT_FAILURE;
    }

  switch( dimension )
    {
    case 2:
      {
      return ants_motion<2>( parser );
      }
      break;
    case 3:
      {
      return ants_motion<3>( parser );
      }
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}

} // namespace ants
