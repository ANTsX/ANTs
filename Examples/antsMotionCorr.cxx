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
#include "ReadWriteImage.h"
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

#include <sstream>

namespace ants
{
/** \class antsRegistrationCommandIterationUpdate
 *  \brief change parameters between iterations of registration
 */
template <class TFilter>
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

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *) caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    TFilter * filter = const_cast<TFilter *>( dynamic_cast<const TFilter *>( object ) );

    unsigned int currentLevel = 0;

    if( typeid( event ) == typeid( itk::IterationEvent ) )
      {
      currentLevel = filter->GetCurrentLevel() + 1;
      }
    if( currentLevel < this->m_NumberOfIterations.size() )
      {
      typename TFilter::ShrinkFactorsPerDimensionContainerType shrinkFactors = filter->GetShrinkFactorsPerDimension( currentLevel );
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
      GradientDescentOptimizerType * optimizer = reinterpret_cast<GradientDescentOptimizerType *>(
          const_cast<typename TFilter::OptimizerType *>( filter->GetOptimizer() ) );
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

template <class T>
inline std::string ants_moco_to_string(const T& t)
{
  std::stringstream ss;

  ss << t;
  return ss.str();
}

template <class ImageType>
typename ImageType::Pointer PreprocessImage( ImageType * inputImage,
                                             typename ImageType::PixelType lowerScaleFunction,
                                             typename ImageType::PixelType upperScaleFunction,
                                             float winsorizeLowerQuantile, float winsorizeUpperQuantile,
                                             ImageType *histogramMatchSourceImage = NULL )
{
  typedef itk::Statistics::ImageToHistogramFilter<ImageType>   HistogramFilterType;
  typedef typename HistogramFilterType::InputBooleanObjectType InputBooleanObjectType;
  typedef typename HistogramFilterType::HistogramSizeType      HistogramSizeType;
  typedef typename HistogramFilterType::HistogramType          HistogramType;

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

  typename ImageType::Pointer outputImage = NULL;
  if( histogramMatchSourceImage )
    {
    typedef itk::HistogramMatchingImageFilter<ImageType, ImageType> HistogramMatchingFilterType;
    typename HistogramMatchingFilterType::Pointer matchingFilter = HistogramMatchingFilterType::New();
    matchingFilter->SetSourceImage( windowingFilter->GetOutput() );
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
    if ( vnl_math_abs( calc->GetMaximum() - calc->GetMinimum() ) < 1.e-9 )
      {
      std::cout <<"Warning: bad time point - too little intensity variation" << std::endl;
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

template <class T>
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

template <class TFilter>
class CommandIterationUpdate : public itk::Command
{
public:
  typedef CommandIterationUpdate  Self;
  typedef itk::Command            Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate()
  {
  };
public:

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *) caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    TFilter * filter = const_cast<TFilter *>( dynamic_cast<const TFilter *>( object ) );

    if( typeid( event ) != typeid( itk::IterationEvent ) )
      {
      return;
      }

    unsigned int currentLevel = filter->GetCurrentLevel();
    typename TFilter::ShrinkFactorsPerDimensionContainerType shrinkFactors = filter->GetShrinkFactorsPerDimension( currentLevel );
    typename TFilter::SmoothingSigmasArrayType smoothingSigmas = filter->GetSmoothingSigmasPerLevel();
    typename TFilter::TransformParametersAdaptorsContainerType adaptors =
      filter->GetTransformParametersAdaptorsPerLevel();

    std::cout << "  Current level = " << currentLevel << std::endl;
    std::cout << "    number of iterations = " << this->m_NumberOfIterations[currentLevel] << std::endl;
    std::cout << "    shrink factor = " << shrinkFactors[currentLevel] << std::endl;
    std::cout << "    smoothing sigma = " << smoothingSigmas[currentLevel] << std::endl;
    std::cout << "    required fixed parameters = " << adaptors[currentLevel]->GetRequiredFixedParameters()
             << std::endl;

    typedef itk::ConjugateGradientLineSearchOptimizerv4 OptimizerType;
    OptimizerType * optimizer = reinterpret_cast<OptimizerType *>(
        const_cast<typename TFilter::OptimizerType *>( filter->GetOptimizer() ) );
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

template <class TImageIn, class TImageOut>
void
AverageTimeImages( typename TImageIn::Pointer image_in,  typename TImageOut::Pointer image_avg,
                   std::vector<unsigned int> timelist )
{
  typedef TImageIn  ImageType;
  typedef TImageOut OutImageType;
  enum { ImageDimension = ImageType::ImageDimension };
  typedef typename TImageIn::PixelType                    PixelType;
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
  std::cout << " averaging with " << timelist.size() << " images of " <<  timedims <<  " timedims " << std::endl;
  Iterator vfIter2(  image_avg, image_avg->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    typename OutImageType::PixelType  fval = 0;
    typename ImageType::IndexType ind;
    typename OutImageType::IndexType spind = vfIter2.GetIndex();
    for( unsigned int xx = 0; xx < timelist.size(); xx++ )
      {
      for( unsigned int yy = 0; yy < ImageDimension - 1; yy++ )
        {
        ind[yy] = spind[yy];
        }
      ind[ImageDimension - 1] = timelist[xx];
      fval += image_in->GetPixel(ind);
      }
    fval /= (double)timelist.size();
    image_avg->SetPixel(spind, fval);
    }
  std::cout << " averaging images done " << std::endl;
  return;
}

template <unsigned int ImageDimension>
int ants_motion( itk::ants::CommandLineParser *parser )
{
  // We infer the number of stages by the number of transformations
  // specified by the user which should match the number of metrics.
  unsigned numberOfStages = 0;

  typedef float                                     PixelType;
  typedef double                                    RealType;
  typedef itk::Image<PixelType, ImageDimension>     FixedIOImageType;
  typedef itk::Image<RealType , ImageDimension>     FixedImageType;
  typedef itk::Image<PixelType, ImageDimension + 1> MovingIOImageType;
  typedef itk::Image<RealType, ImageDimension + 1>  MovingImageType;
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
      std::cout << "Output option not specified.  Should be the output average image name." << std::endl;
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
    if ( ImageDimension == 2 ) extractFilter->SetDirectionCollapseToIdentity();
    unsigned int td = 0;
    extractRegion.SetIndex(ImageDimension, td );
    extractFilter->SetExtractionRegion( extractRegion );
    extractFilter->Update();
    avgImage = extractFilter->GetOutput();
    std::vector<unsigned int> timelist;
    AverageTimeImages<MovingIOImageType, FixedIOImageType>( movingImage, avgImage, timelist );
    std::cout << "average out " << outputPrefix <<  std::endl;
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
    std::cout << "No transformations are specified." << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Registration using " << numberOfStages << " total stages." << std::endl;

  typename OptionType::Pointer metricOption = parser->GetOption( "metric" );
  if( !metricOption || metricOption->GetNumberOfFunctions() != numberOfStages  )
    {
    std::cout << "The number of metrics specified does not match the number of stages." << std::endl;
    return EXIT_FAILURE;
    }

  typename OptionType::Pointer iterationsOption = parser->GetOption( "iterations" );
  if( !iterationsOption || iterationsOption->GetNumberOfFunctions() != numberOfStages  )
    {
    std::cout << "The number of iteration sets specified does not match the number of stages." << std::endl;
    return EXIT_FAILURE;
    }

  typename OptionType::Pointer shrinkFactorsOption = parser->GetOption( "shrinkFactors" );
  if( !shrinkFactorsOption || shrinkFactorsOption->GetNumberOfFunctions() != numberOfStages  )
    {
    std::cout << "The number of shrinkFactor sets specified does not match the number of stages." << std::endl;
    return EXIT_FAILURE;
    }

  typename OptionType::Pointer smoothingSigmasOption = parser->GetOption( "smoothingSigmas" );
  if( !smoothingSigmasOption || smoothingSigmasOption->GetNumberOfFunctions() != numberOfStages  )
    {
    std::cout << "The number of smoothing sigma sets specified does not match the number of stages." << std::endl;
    return EXIT_FAILURE;
    }

  typename OptionType::Pointer outputOption = parser->GetOption( "output" );
  if( !outputOption )
    {
    std::cout << "Output option not specified." << std::endl;
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
    std::cout << " nimagestoavg " << nimagestoavg << std::endl;
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

  unsigned int   nparams = 2;
  itk::TimeProbe totalTimer;
  totalTimer.Start();
  double metricmean = 0;

  typedef itk::AffineTransform<RealType, ImageDimension>                                      AffineTransformType;
  typedef itk::ImageRegistrationMethodv4<FixedImageType, FixedImageType, AffineTransformType> AffineRegistrationType;
  // We iterate backwards because the command line options are stored as a stack (first in last out)
  for( int currentStage = numberOfStages - 1; currentStage >= 0; currentStage-- )
    {
    std::cout << std::endl << "Stage " << numberOfStages - currentStage << std::endl;
    std::stringstream currentStageString;
    currentStageString << currentStage;

    // Get the fixed and moving images

    std::string fixedImageFileName = metricOption->GetFunction( currentStage )->GetParameter(  0 );
    std::string movingImageFileName = metricOption->GetFunction( currentStage )->GetParameter(  1 );
    std::cout << "  fixed image: " << fixedImageFileName << std::endl;
    std::cout << "  moving image: " << movingImageFileName << std::endl;
    typename FixedImageType::Pointer fixed_time_slice = NULL;
    typename FixedImageType::Pointer moving_time_slice = NULL;
    typename FixedIOImageType::Pointer fixedInImage;
    ReadImage<FixedIOImageType>( fixedInImage, fixedImageFileName.c_str() );
    fixedInImage->Update();
    fixedInImage->DisconnectPipeline();
    typename FixedImageType::Pointer fixedImage;
    fixedImage = arCastImage< FixedIOImageType, FixedImageType >( fixedInImage );

    typename MovingIOImageType::Pointer movingInImage;
    typename MovingImageType::Pointer movingImage;
    ReadImage<MovingIOImageType>( movingInImage, movingImageFileName.c_str()  );
    movingInImage->Update();
    movingInImage->DisconnectPipeline();
    movingImage = arCastImage<MovingIOImageType,MovingImageType>( movingInImage );

    typename MovingIOImageType::Pointer outputImage;
    ReadImage<MovingIOImageType>( outputImage, movingImageFileName.c_str() );
    outputImage->Update();
    outputImage->DisconnectPipeline();

    // Get the number of iterations and use that information to specify the number of levels

    std::vector<unsigned int> iterations =
      parser->ConvertVector<unsigned int>( iterationsOption->GetFunction( currentStage )->GetName()  );
    unsigned int numberOfLevels = iterations.size();
    std::cout << "  number of levels = " << numberOfLevels << std::endl;

    // Get shrink factors

    std::vector<unsigned int> factors =
      parser->ConvertVector<unsigned int>( shrinkFactorsOption->GetFunction( currentStage )->GetName()  );
    typename AffineRegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
    shrinkFactorsPerLevel.SetSize( factors.size() );

    if( factors.size() != numberOfLevels )
      {
      std::cout << "ERROR:  The number of shrink factors does not match the number of levels." << std::endl;
      return EXIT_FAILURE;
      }
    else
      {
      for( unsigned int n = 0; n < shrinkFactorsPerLevel.Size(); n++ )
        {
        shrinkFactorsPerLevel[n] = factors[n];
        }
      std::cout << "  shrink factors per level: " << shrinkFactorsPerLevel << std::endl;
      }

    // Get smoothing sigmas

    std::vector<float> sigmas = parser->ConvertVector<float>( smoothingSigmasOption->GetFunction(
                                                                currentStage )->GetName()  );
    typename AffineRegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
    smoothingSigmasPerLevel.SetSize( sigmas.size() );

    if( sigmas.size() != numberOfLevels )
      {
      std::cout << "ERROR:  The number of smoothing sigmas does not match the number of levels." << std::endl;
      return EXIT_FAILURE;
      }
    else
      {
      for( unsigned int n = 0; n < smoothingSigmasPerLevel.Size(); n++ )
        {
        smoothingSigmasPerLevel[n] = sigmas[n];
        }
      std::cout << "  smoothing sigmas per level: " << smoothingSigmasPerLevel << std::endl;
      }

    // the fixed image is a reference image in 3D while the moving is a 4D image
    // loop over every time point and register image_i+1 to image_i
    //
    // Set up the image metric and scales estimator
    unsigned int              timedims = movingImage->GetLargestPossibleRegion().GetSize()[ImageDimension];
    std::vector<unsigned int> timelist;
    std::vector<double>       metriclist;
    for( unsigned int timedim = 0; timedim < timedims; timedim++ )
      {
      timelist.push_back(timedim);
      }
    for( unsigned int timelistindex = 0;  timelistindex < timelist.size();  timelistindex++ )
      {
      unsigned int timedim = timelist[timelistindex];
      typename CompositeTransformType::Pointer compositeTransform = NULL;
      if( currentStage == static_cast<int>(numberOfStages) - 1 )
        {
        compositeTransform = CompositeTransformType::New();
        CompositeTransformVector.push_back(compositeTransform);
        }
      else if( CompositeTransformVector.size() == timedims && !CompositeTransformVector[timedim].IsNull() )
        {
        compositeTransform = CompositeTransformVector[timedim];
        if ( timedim == 0 ) std::cout << " use existing transform " << compositeTransform->GetParameters() << std::endl;
        }
      typedef itk::IdentityTransform<RealType, ImageDimension> IdentityTransformType;
      typename IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();
      //
      typedef itk::ExtractImageFilter<MovingImageType, FixedImageType> ExtractFilterType;
      typename MovingImageType::RegionType extractRegion = movingImage->GetLargestPossibleRegion();
      extractRegion.SetSize(ImageDimension, 0);
      bool maptoneighbor = true;
      typename OptionType::Pointer fixedOption = parser->GetOption( "useFixedReferenceImage" );
      if( fixedOption && fixedOption->GetNumberOfFunctions() )
        {
        std::string fixedFunction = fixedOption->GetFunction( 0 )->GetName();
        ConvertToLowerCase( fixedFunction );
        if( fixedFunction.compare( "1" ) == 0 || fixedFunction.compare( "true" ) == 0 )
          {
          if( timedim == 0 )
            {
            std::cout << "using fixed reference image for all frames " << std::endl;
            }
          fixed_time_slice = fixedImage;
          extractRegion.SetIndex(ImageDimension, timedim );
          typename ExtractFilterType::Pointer extractFilter2 = ExtractFilterType::New();
          extractFilter2->SetInput( movingImage );
          extractFilter2->SetDirectionCollapseToSubmatrix();
          if ( ImageDimension == 2 ) extractFilter2->SetDirectionCollapseToIdentity();
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
        if ( ImageDimension == 2 ) extractFilter->SetDirectionCollapseToIdentity();
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
        if ( ImageDimension == 2 ) extractFilter->SetDirectionCollapseToIdentity();
        extractFilter2->SetExtractionRegion( extractRegion );
        extractFilter2->Update();
        moving_time_slice = extractFilter2->GetOutput();
        }

      bool directionmatricesok = true;
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        for( unsigned int j = 0; j < ImageDimension; j++ )
          {
          if( fabs( moving_time_slice->GetDirection()[i][j] - fixed_time_slice->GetDirection()[i][j] ) > 1.e-6 )
            {
            directionmatricesok = false;
            }
          }
        }

      if ( ( !directionmatricesok ) && ( timedim == 0 )   )
        {
        std::cout << " WARNING!" << std::endl;
        std::cout << " fixed and moving DirectionMatrices not the same " << std::endl;
        std::cout << " Fixed Dir " << fixed_time_slice->GetDirection()  << std::endl;
        std::cout << " Moving Dir " << moving_time_slice->GetDirection()  << std::endl;
        std::cout << " setting moving direction matrix to equal fixed matrix " << std::endl;
        std::cout << " WARNING END!" << std::endl;
        std::cout <<  std::endl;
        moving_time_slice->SetDirection(  fixed_time_slice->GetDirection()  );
        }

      typename FixedImageType::Pointer preprocessFixedImage =
        PreprocessImage<FixedImageType>( fixed_time_slice, 0,
                                         1, 0.001, 0.999,
                                         NULL );

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
        if ( timedim == 0 ) std::cout << "  random sampling (percentage = " << samplingPercentage << ")" << std::endl;
        metricSamplingStrategy = AffineRegistrationType::RANDOM;
        }
      if( std::strcmp( samplingStrategy.c_str(), "regular" ) == 0 )
        {
        if ( timedim == 0 ) std::cout << "  regular sampling (percentage = " << samplingPercentage << ")" << std::endl;
        metricSamplingStrategy = AffineRegistrationType::REGULAR;
        }

      if( std::strcmp( whichMetric.c_str(), "cc" ) == 0 )
        {
        unsigned int radiusOption = parser->Convert<unsigned int>( metricOption->GetFunction(
                                                                     currentStage )->GetParameter(  3 ) );

        if ( timedim == 0 ) std::cout << "  using the CC metric (radius = " << radiusOption << ")." << std::endl;
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
        unsigned int binOption =
          parser->Convert<unsigned int>( metricOption->GetFunction( currentStage )->GetParameter(  3 ) );
        typedef itk::MattesMutualInformationImageToImageMetricv4<FixedImageType,
                                                                 FixedImageType> MutualInformationMetricType;
        typename MutualInformationMetricType::Pointer mutualInformationMetric = MutualInformationMetricType::New();
        mutualInformationMetric = mutualInformationMetric;
        mutualInformationMetric->SetNumberOfHistogramBins( binOption );
        mutualInformationMetric->SetUseMovingImageGradientFilter( false );
        mutualInformationMetric->SetUseFixedImageGradientFilter( false );
        metric = mutualInformationMetric;
        }
      else if( std::strcmp( whichMetric.c_str(), "demons" ) == 0 )
        {
        if ( timedim == 0 ) std::cout << "  using the Demons metric." << std::endl;
        typedef itk::MeanSquaresImageToImageMetricv4<FixedImageType, FixedImageType> DemonsMetricType;
        typename DemonsMetricType::Pointer demonsMetric = DemonsMetricType::New();
        demonsMetric = demonsMetric;
        metric = demonsMetric;
        }
      else if( std::strcmp( whichMetric.c_str(), "gc" ) == 0 )
        {
        if ( timedim == 0 ) std::cout << "  using the global correlation metric." << std::endl;
        typedef itk::CorrelationImageToImageMetricv4<FixedImageType, FixedImageType> corrMetricType;
        typename corrMetricType::Pointer corrMetric = corrMetricType::New();
        metric = corrMetric;
        std::cout << " global corr metric set " << std::endl;
        }
      else
        {
        std::cout << "ERROR: Unrecognized image metric: " << whichMetric << std::endl;
        }
      metric->SetVirtualDomainFromImage(  fixed_time_slice );
      // Set up the optimizer.  To change the iteration number for each level we rely
      // on the command observer.
      //    typedef itk::JointHistogramMutualInformationImageToImageMetricv4<FixedImageType, FixedImageType>
      // MutualInformationMetricType;

      typedef itk::RegistrationParameterScalesFromPhysicalShift<MetricType> ScalesEstimatorType;
      typename ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();
      scalesEstimator->SetMetric( metric );
      scalesEstimator->SetTransformForward( true );

      float learningRate = parser->Convert<float>( transformOption->GetFunction( currentStage )->GetParameter(  0 ) );

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
	  if ( timedim == 0 ) std::cout << " employing scales estimator " << std::endl;
          optimizer->SetScalesEstimator( scalesEstimator );
          }
        else
          {
          if ( timedim == 0 ) std::cout << " not employing scales estimator " << scalesFunction << std::endl;
          }
        }
      optimizer->SetMaximumStepSizeInPhysicalUnits( learningRate );
      optimizer->SetDoEstimateLearningRateOnce( doEstimateLearningRateOnce );
      optimizer->SetDoEstimateLearningRateAtEachIteration( !doEstimateLearningRateOnce );
      //    optimizer->SetMaximumNewtonStepSizeInPhysicalUnits(sqrt(small_step)*learningR);

      // Set up the image registration methods along with the transforms
      std::string whichTransform = transformOption->GetFunction( currentStage )->GetName();
      ConvertToLowerCase( whichTransform );
      if( std::strcmp( whichTransform.c_str(), "affine" ) == 0 )
        {
        typename AffineRegistrationType::Pointer affineRegistration = AffineRegistrationType::New();
        typedef itk::AffineTransform<double, ImageDimension> AffineTransformType;
        typename AffineTransformType::Pointer affineTransform = AffineTransformType::New();
        affineTransform->SetIdentity();
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
          std::cout << std::endl << "*** Running affine registration ***" << timedim << std::endl << std::endl;
          affineRegistration->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cout << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }
        compositeTransform->AddTransform( const_cast<AffineTransformType *>( affineRegistration->GetOutput()->Get() ) );
        // Write out the affine transform
        std::string filename = outputPrefix + std::string("TimeSlice") + ants_moco_to_string<unsigned int>(timedim)
          + std::string( "Affine.txt" );
        typedef itk::TransformFileWriter TransformWriterType;
        typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
        transformWriter->SetInput( affineRegistration->GetOutput()->Get() );
        transformWriter->SetFileName( filename.c_str() );
        //      transformWriter->Update();
        if( timedim == 0 )
          {
          param_values.set_size(timedims, nparams);
          param_values.fill(0);
          }
        for ( unsigned int i = 0; i < nparams - 2; i++ )
          {
          param_values(timedim, i + 2) = affineRegistration->GetOutput()->Get()->GetParameters()[i];
          }
        }
      else if( std::strcmp( whichTransform.c_str(), "rigid" ) == 0 )
        {
        typedef typename RigidTransformTraits<ImageDimension>::TransformType RigidTransformType;
        typename RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
        nparams = rigidTransform->GetNumberOfParameters() + 2;
        typedef itk::ImageRegistrationMethodv4<FixedImageType, FixedImageType,
                                               RigidTransformType> RigidRegistrationType;
        typename RigidRegistrationType::Pointer rigidRegistration = RigidRegistrationType::New();
        metric->SetFixedImage( preprocessFixedImage );
        metric->SetVirtualDomainFromImage( preprocessFixedImage );
        metric->SetMovingImage( preprocessMovingImage );
        metric->SetMovingTransform( rigidTransform );
        typename ScalesEstimatorType::ScalesType scales(rigidTransform->GetNumberOfParameters() );
        typename MetricType::ParametersType      newparams(  rigidTransform->GetParameters() );
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
          std::cout << std::endl << "*** Running rigid registration ***" << timedim  << std::endl << std::endl;
          rigidRegistration->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cout << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }
        compositeTransform->AddTransform( const_cast<RigidTransformType *>( rigidRegistration->GetOutput()->Get() ) );
        // Write out the rigid transform
        std::string filename = outputPrefix + std::string("TimeSlice") + ants_moco_to_string<unsigned int>(timedim)
          + std::string( "Rigid.txt" );
        typedef itk::TransformFileWriter TransformWriterType;
        typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
        transformWriter->SetInput( rigidRegistration->GetOutput()->Get() );
        transformWriter->SetFileName( filename.c_str() );
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
        typedef itk::Vector<RealType, VImageDimension> VectorType;
        VectorType zeroVector( 0.0 );
        typedef itk::Image<VectorType, VImageDimension> DisplacementFieldType;
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
          const_cast<GaussianDisplacementFieldTransformType *>( displacementFieldRegistration->GetOutput()->Get() );

        // Create the transform adaptors

        typedef itk::GaussianSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor<
            GaussianDisplacementFieldTransformType> DisplacementFieldTransformAdaptorType;
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
          std::cout << "Exception caught: " << e << std::endl;
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
        typedef itk::Vector<RealType, VImageDimension> VectorType;
        VectorType zeroVector( 0.0 );
        typedef itk::Image<VectorType, VImageDimension> DisplacementFieldType;

        typename DisplacementFieldType::Pointer displacementField = AllocImage<DisplacementFieldType>( preprocessFixedImage, zeroVector );

        typename DisplacementFieldType::Pointer inverseDisplacementField = AllocImage<DisplacementFieldType>(
            preprocessFixedImage, zeroVector );

	typedef itk::DisplacementFieldTransform<RealType, VImageDimension>         DisplacementFieldTransformType;
        typedef itk::SyNImageRegistrationMethod<FixedImageType, FixedImageType,
                                                DisplacementFieldTransformType> DisplacementFieldRegistrationType;
        typename DisplacementFieldRegistrationType::Pointer displacementFieldRegistration =
          DisplacementFieldRegistrationType::New();

        typename DisplacementFieldTransformType::Pointer outputDisplacementFieldTransform =
          const_cast<DisplacementFieldTransformType *>( displacementFieldRegistration->GetOutput()->Get() );

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
        if( timedim == 0 ) std::cout << "SyN iterations:";
        for( unsigned int d = 0; d < numberOfLevels; d++ )
          {
	  numberOfIterationsPerLevel[d] = iterations[d]; // currentStageIterations[d];
	  if( timedim == 0 ) std::cout << numberOfIterationsPerLevel[d] << " ";
          }
	if( timedim == 0 ) std::cout << std::endl;

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
          std::cout << "Exception caught: " << e << std::endl;
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
        std::cout << "ERROR:  Unrecognized transform option - " << whichTransform << std::endl;
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
      resampler->Update();
      std::cout << " done resampling timepoint : " << timedim << std::endl;

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
        unsigned int tdim = timedim ;
        if( tdim > ( timedims - 1 ) )
          {
          tdim = timedims - 1;
          }
        ind[ImageDimension] = tdim;
        outputImage->SetPixel(ind, fval);
        }
      }
    if( outputOption && outputOption->GetFunction( 0 )->GetNumberOfParameters() > 1  && currentStage == 0 )
      {
      std::string fileName = outputOption->GetFunction( 0 )->GetParameter( 1 );
      if( outputPrefix.length() < 3 )
        {
        outputPrefix = outputOption->GetFunction( 0 )->GetName();
        }
      std::cout << "motion corrected out " << fileName <<  std::endl;
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
      if ( ImageDimension == 2 ) extractFilter->SetDirectionCollapseToIdentity();
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
	if ( i < timelist.size() ) timelistsort.push_back(timelist[i]);
        std::cout << " i^th value " << i << "  is " << metriclist[timelist[i]] << std::endl;
        }
      AverageTimeImages<MovingIOImageType, FixedIOImageType>( outputImage, avgImage, timelistsort );
      std::cout << " write average post " << fileName << std::endl;
      WriteImage<FixedIOImageType>( avgImage, fileName.c_str() );
      }
    }
  totalTimer.Stop();
  std::cout << std::endl << "Total elapsed time: " << totalTimer.GetMean() << " averagemetric " << metricmean
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
    std::cout << " get motion corr params " << outputPrefix << std::endl;
    if( outputPrefix[0] == '0' && outputPrefix[1] == 'x' )
      {
      void* ptr;
      std::sscanf(outputPrefix.c_str(), "%p", (void **)&ptr);
      //      std::stringstream strstream;
      //      strstream << outputPrefix;
      //      void* ptr;
      //      strstream >> ptr;
      ( static_cast<std::pair<std::vector<std::string>, vnl_matrix<float> > *>( ptr ) )->first = ColumnHeaders;
      ( static_cast<std::pair<std::vector<std::string>, vnl_matrix<double> > *>( ptr ) )->second = param_values;
      std::cout << "motion-correction params written" << std::endl;
      }
    else
      {
      fnmp = outputPrefix + std::string("MOCOparams.csv");
      std::cout << " write " << fnmp << std::endl;
      writer->SetFileName( fnmp.c_str() );
      writer->SetColumnHeaders(ColumnHeaders);
      writer->SetInput( &param_values );
      writer->Write();
      }
    }

  return EXIT_SUCCESS;
}

void InitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
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
      + std::string( "are available for a single stage." );

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
        "use a fixed reference image instead of the neighor in the time series." );
    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "useFixedReferenceImage" );
    option->SetShortName( 'u' );
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
int antsMotionCorr( std::vector<std::string> args, std::ostream* out_stream = NULL )
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
  argv[argc] = 0;
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
    + std::string( "registration application meant to utilize ITKv4-only classes. The user can specify " )
    + std::string( "any number of \"stages\" where a stage consists of a transform; an image metric; " )
    + std::string( " and iterations, shrink factors, and smoothing sigmas for each level. " )
    + std::string(
      " Specialized for 4D time series data: fixed image is 3D, moving image should be the 4D time series. ")
    + std::string( " Fixed image is a reference space or time slice.");
  parser->SetCommandDescription( commandDescription );
  InitializeCommandLineOptions( parser );

  parser->Parse( argc, argv );

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
    std::cout << "Image dimensionality not specified.  See command line option --dimensionality" << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << std::endl << "Running " << argv[0] << "  for " << dimension << "-dimensional images." << std::endl
           << std::endl;

  switch( dimension )
    {
    case 2:
      {
      ants_motion<2>( parser );
      }
      break;
    case 3:
      {
      ants_motion<3>( parser );
      }
      break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return 0;
}
} // namespace ants
