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

#include "antsCommandLineParser.h"
#include "itkImageRegistrationMethodv4.h"
#include "itkSyNImageRegistrationMethod.h"
#include "itkTimeVaryingVelocityFieldImageRegistrationMethodv4.h"
#include "itkTimeVaryingBSplineVelocityFieldImageRegistrationMethod.h"

#include "itkANTSNeighborhoodCorrelationImageToImageMetricv4.h"
#include "itkDemonsImageToImageMetricv4.h"
#include "itkCorrelationImageToImageMetricv4.h"
#include "itkImageToImageMetricv4.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"

#include "itkAffineTransform.h"
#include "itkANTSAffine3DTransform.h"
#include "itkANTSCenteredAffine2DTransform.h"
#include "itkBSplineTransform.h"
#include "itkBSplineSmoothingOnUpdateDisplacementFieldTransform.h"
#include "itkCompositeTransform.h"
#include "itkDisplacementFieldTransform.h"
#include "itkGaussianSmoothingOnUpdateDisplacementFieldTransform.h"
#include "itkIdentityTransform.h"
#include "itkEuler2DTransform.h"
#include "itkEuler3DTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "itkQuaternionRigidTransform.h"
#include "itkSimilarity2DTransform.h"
#include "itkSimilarity3DTransform.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkTranslationTransform.h"
#include "itkTransform.h"
#include "itkTransformFactory.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"

#include "itkBSplineTransformParametersAdaptor.h"
#include "itkBSplineSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor.h"
#include "itkGaussianSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor.h"
#include "itkTimeVaryingVelocityFieldTransformParametersAdaptor.h"
#include "itkTimeVaryingBSplineVelocityFieldTransformParametersAdaptor.h"

#include "itkGradientDescentOptimizerv4.h"

#include "itkRegistrationParameterScalesFromShift.h"

#include "itkHistogramMatchingImageFilter.h"
#include "itkImageToHistogramFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkMacro.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkTimeProbe.h"
#include "itkVector.h"

#include <sstream>

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

    unsigned int currentLevel = 0;

    if( typeid( event ) == typeid( itk::IterationEvent ) )
      {
      currentLevel = filter->GetCurrentLevel() + 1;
      }
    if( currentLevel < this->m_NumberOfIterations.size() )
      {
      typename TFilter::ShrinkFactorsArrayType shrinkFactors = filter->GetShrinkFactorsPerLevel();
      typename TFilter::SmoothingSigmasArrayType smoothingSigmas = filter->GetSmoothingSigmasPerLevel();
      typename TFilter::TransformParametersAdaptorsContainerType adaptors =
        filter->GetTransformParametersAdaptorsPerLevel();

      std::cout << "  Current level = " << currentLevel << std::endl;
      std::cout << "    number of iterations = " << this->m_NumberOfIterations[currentLevel] << std::endl;
      std::cout << "    shrink factors = " << shrinkFactors[currentLevel] << std::endl;
      std::cout << "    smoothing sigmas = " << smoothingSigmas[currentLevel] << std::endl;
      std::cout << "    required fixed parameters = " << adaptors[currentLevel]->GetRequiredFixedParameters()
                << std::endl;

      typedef itk::GradientDescentOptimizerv4 GradientDescentOptimizerType;
      GradientDescentOptimizerType * optimizer = reinterpret_cast<GradientDescentOptimizerType *>(
          const_cast<typename TFilter::OptimizerType *>( filter->GetOptimizer() ) );
      optimizer->SetNumberOfIterations( this->m_NumberOfIterations[currentLevel] );
      }
  }

  void SetNumberOfIterations( std::vector<unsigned int> iterations )
  {
    this->m_NumberOfIterations = iterations;
  }

private:
  std::vector<unsigned int> m_NumberOfIterations;
};

//
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

void ConvertToLowerCase( std::string& str )
{
  std::transform( str.begin(), str.end(), str.begin(), tolower );
// You may need to cast the above line to (int(*)(int))
// tolower - this works as is on VC 7.1 but may not work on
// other compilers
}

template <class ImageType>
void PreprocessImages( ImageType *fixedImage, ImageType *movingImage,
                       itk::ants::CommandLineParser *parser )
{
  typedef typename itk::ants::CommandLineParser ParserType;
  typedef typename ParserType::OptionType       OptionType;

  std::string outputPreprocessingString = "";

  float lowerQuantile = 0.0;
  float upperQuantile = 1.0;

  bool doWinsorize = false;
  typename OptionType::Pointer winsorizeOption = parser->GetOption( "winsorize-image-intensities" );
  if( winsorizeOption && winsorizeOption->GetNumberOfParameters( 0 ) > 0 )
    {
    doWinsorize = true;

    lowerQuantile = 0.05;
    upperQuantile = 0.95;

    if( winsorizeOption->GetNumberOfParameters( 0 ) > 0 )
      {
      lowerQuantile = parser->Convert<float>( winsorizeOption->GetParameter( 0, 0 ) );
      }
    if( winsorizeOption->GetNumberOfParameters( 0 ) > 1 )
      {
      upperQuantile = parser->Convert<float>( winsorizeOption->GetParameter( 0, 1 ) );
      }
    outputPreprocessingString += std::string( "  preprocessing:  winsorizing the image intensities\n" );
    }

  bool doRescale = false;
  typename OptionType::Pointer rescaleOption = parser->GetOption( "rescale-images" );
  if( rescaleOption && rescaleOption->GetNumberOfValues() > 0 )
    {
    std::string histValue = rescaleOption->GetValue( 0 );
    ConvertToLowerCase( histValue );
    if( histValue.compare( "1" ) == 0 || histValue.compare( "true" ) == 0 )
      {
      doRescale = true;
      outputPreprocessingString += std::string( "  preprocessing:  rescaling the images to [0,1]\n" );
      }
    }

  typename ImageType::Pointer fixedWinsorizedRescaled = NULL;
  typename ImageType::Pointer movingWinsorizedRescaled = NULL;

  if( doWinsorize || doRescale )
    {
    typedef itk::Statistics::ImageToHistogramFilter<ImageType>   HistogramFilterType;
    typedef typename HistogramFilterType::InputBooleanObjectType InputBooleanObjectType;
    typedef typename HistogramFilterType::HistogramSizeType      HistogramSizeType;
    typedef typename HistogramFilterType::HistogramType          HistogramType;

    HistogramSizeType histogramSize( 1 );
    histogramSize[0] = 256;

    typename InputBooleanObjectType::Pointer autoMinMaxInputObject = InputBooleanObjectType::New();
    autoMinMaxInputObject->Set( true );

    typename HistogramFilterType::Pointer fixedHistogramFilter = HistogramFilterType::New();
    fixedHistogramFilter->SetInput( fixedImage );
    fixedHistogramFilter->SetAutoMinimumMaximumInput( autoMinMaxInputObject );
    fixedHistogramFilter->SetHistogramSize( histogramSize );
    fixedHistogramFilter->SetMarginalScale( 10.0 );
    fixedHistogramFilter->Update();

    float fixedLowerValue = fixedHistogramFilter->GetOutput()->Quantile( 0, lowerQuantile );
    float fixedUpperValue = fixedHistogramFilter->GetOutput()->Quantile( 0, upperQuantile );

    typename HistogramFilterType::Pointer movingHistogramFilter = HistogramFilterType::New();
    movingHistogramFilter->SetInput( movingImage );
    movingHistogramFilter->SetAutoMinimumMaximumInput( autoMinMaxInputObject );
    movingHistogramFilter->SetHistogramSize( histogramSize );
    movingHistogramFilter->SetMarginalScale( 10.0 );
    movingHistogramFilter->Update();

    float movingLowerValue = movingHistogramFilter->GetOutput()->Quantile( 0, lowerQuantile );
    float movingUpperValue = movingHistogramFilter->GetOutput()->Quantile( 0, upperQuantile );

    typedef itk::IntensityWindowingImageFilter<ImageType, ImageType> IntensityWindowingImageFilterType;

    typename IntensityWindowingImageFilterType::Pointer fixedWindowingFilter = IntensityWindowingImageFilterType::New();
    fixedWindowingFilter->SetInput( fixedImage );
    fixedWindowingFilter->SetWindowMinimum( fixedLowerValue );
    fixedWindowingFilter->SetWindowMaximum( fixedUpperValue );
    if( doRescale )
      {
      fixedWindowingFilter->SetOutputMinimum( 0.0 );
      fixedWindowingFilter->SetOutputMaximum( 1.0 );
      }
    else
      {
      fixedWindowingFilter->SetOutputMinimum( fixedLowerValue );
      fixedWindowingFilter->SetOutputMaximum( fixedUpperValue );
      }
    fixedWinsorizedRescaled = fixedWindowingFilter->GetOutput();
    fixedWinsorizedRescaled->Update();
    fixedWinsorizedRescaled->DisconnectPipeline();

    typename IntensityWindowingImageFilterType::Pointer movingWindowingFilter =
      IntensityWindowingImageFilterType::New();
    movingWindowingFilter->SetInput( movingImage );
    movingWindowingFilter->SetWindowMinimum( movingLowerValue );
    movingWindowingFilter->SetWindowMaximum( movingUpperValue );
    if( doRescale )
      {
      movingWindowingFilter->SetOutputMinimum( 0.0 );
      movingWindowingFilter->SetOutputMaximum( 1.0 );
      }
    else
      {
      movingWindowingFilter->SetOutputMinimum( movingLowerValue );
      movingWindowingFilter->SetOutputMaximum( movingUpperValue );
      }
    movingWinsorizedRescaled = movingWindowingFilter->GetOutput();
    movingWinsorizedRescaled->Update();
    movingWinsorizedRescaled->DisconnectPipeline();
    }

  typename OptionType::Pointer histOption = parser->GetOption( "use-histogram-matching" );
  if( histOption && histOption->GetNumberOfValues() > 0 )
    {
    std::string histValue = histOption->GetValue( 0 );
    ConvertToLowerCase( histValue );
    if( histValue.compare( "1" ) == 0 || histValue.compare( "true" ) == 0 )
      {
      outputPreprocessingString += std::string( "  preprocessing:  histogram matching the images\n" );

      typedef itk::HistogramMatchingImageFilter<ImageType, ImageType> HistogramMatchingFilterType;
      typename HistogramMatchingFilterType::Pointer matchingFilter = HistogramMatchingFilterType::New();
      if( doWinsorize || doRescale )
        {
        matchingFilter->SetSourceImage( movingWinsorizedRescaled );
        matchingFilter->SetReferenceImage( fixedWinsorizedRescaled );
        }
      else
        {
        matchingFilter->SetSourceImage( movingImage );
        matchingFilter->SetReferenceImage( fixedImage );
        }
      matchingFilter->SetNumberOfHistogramLevels( 256 );
      matchingFilter->SetNumberOfMatchPoints( 12 );
      matchingFilter->ThresholdAtMeanIntensityOn();
      matchingFilter->Update();

      movingImage = matchingFilter->GetOutput();
      }
    }

  std::cout << outputPreprocessingString << std::flush;
}

template <unsigned int ImageDimension>
int antsRegistration( itk::ants::CommandLineParser *parser )
{
  itk::TimeProbe totalTimer;

  totalTimer.Start();

  // We infer the number of stages by the number of transformations
  // specified by the user which should match the number of metrics.

  unsigned numberOfStages = 0;

  typedef typename itk::ants::CommandLineParser ParserType;
  typedef typename ParserType::OptionType       OptionType;

  typename OptionType::Pointer transformOption = parser->GetOption( "transform" );
  if( transformOption && transformOption->GetNumberOfValues() > 0 )
    {
    numberOfStages = transformOption->GetNumberOfValues();
    }
  else
    {
    std::cerr << "No transformations are specified." << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Registration using " << numberOfStages << " total stages." << std::endl;

  typename OptionType::Pointer metricOption = parser->GetOption( "metric" );
  if( !metricOption || metricOption->GetNumberOfValues() != numberOfStages  )
    {
    std::cerr << "The number of metrics specified does not match the number of stages." << std::endl;
    return EXIT_FAILURE;
    }

  typename OptionType::Pointer iterationsOption = parser->GetOption( "iterations" );
  if( !iterationsOption || iterationsOption->GetNumberOfValues() != numberOfStages  )
    {
    std::cerr << "The number of iteration sets specified does not match the number of stages." << std::endl;
    return EXIT_FAILURE;
    }

  typename OptionType::Pointer shrinkFactorsOption = parser->GetOption( "shrink-factors" );
  if( !shrinkFactorsOption || shrinkFactorsOption->GetNumberOfValues() != numberOfStages  )
    {
    std::cerr << "The number of shrinkFactor sets specified does not match the number of stages." << std::endl;
    return EXIT_FAILURE;
    }

  typename OptionType::Pointer smoothingSigmasOption = parser->GetOption( "smoothing-sigmas" );
  if( !smoothingSigmasOption || smoothingSigmasOption->GetNumberOfValues() != numberOfStages  )
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
  std::string outputPrefix = outputOption->GetValue( 0 );
  if( outputOption->GetNumberOfParameters( 0 ) > 0 )
    {
    outputPrefix = outputOption->GetParameter( 0, 0 );
    }

  typedef float                                 PixelType;
  typedef double                                RealType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::CompositeTransform<RealType, ImageDimension> CompositeTransformType;
  typename CompositeTransformType::Pointer compositeTransform = CompositeTransformType::New();

  // Register the matrix offset transform base class to the
  // transform factory for compatibility with the current ANTs.
  typedef itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension> MatrixOffsetTransformType;
  itk::TransformFactory<MatrixOffsetTransformType>::RegisterTransform();
  typedef itk::TranslationTransform<double, ImageDimension> TranslationTransformType;
  itk::TransformFactory<TranslationTransformType>::RegisterTransform();
  typedef itk::AffineTransform<double, ImageDimension> AffineTransformType;
  itk::TransformFactory<AffineTransformType>::RegisterTransform();
  typedef itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension> MatrixOffsetTransformType;
  itk::TransformFactory<MatrixOffsetTransformType>::RegisterTransform();
  typedef itk::TranslationTransform<double, ImageDimension> TranslationTransformType;
  itk::TransformFactory<TranslationTransformType>::RegisterTransform();
  typedef itk::AffineTransform<double, ImageDimension> AffineTransformType;
  itk::TransformFactory<AffineTransformType>::RegisterTransform();

  // Load an identity transform in case no transforms are loaded.
  typedef itk::IdentityTransform<RealType, ImageDimension> IdentityTransformType;
  typename IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();

  compositeTransform->AddTransform( identityTransform );

  // Load an initial initialTransform if requested
  unsigned int numberOfInitialTransforms = 0;

  typename itk::ants::CommandLineParser::OptionType::Pointer initialTransformOption =
    parser->GetOption( "initial-transform" );
  if( initialTransformOption && initialTransformOption->GetNumberOfValues() > 0 )
    {
    numberOfInitialTransforms = initialTransformOption->GetNumberOfValues();

    std::deque<std::string> initialTransformNames;
    std::deque<std::string> initialTransformTypes;
    for( unsigned int n = 0; n < initialTransformOption->GetNumberOfValues(); n++ )
      {
      std::string initialTransformName;
      std::string initialTransformType;

      typedef itk::Transform<double, ImageDimension, ImageDimension> TransformType;
      typename TransformType::Pointer initialTransform;

      bool hasTransformBeenRead = false;
      try
        {
        initialTransformName = initialTransformOption->GetValue( n );

        typedef itk::DisplacementFieldTransform<double, ImageDimension>
          DisplacementFieldTransformType;

        typedef typename DisplacementFieldTransformType::DisplacementFieldType
          DisplacementFieldType;

        typedef itk::ImageFileReader<DisplacementFieldType> DisplacementFieldReaderType;
        typename DisplacementFieldReaderType::Pointer fieldReader =
          DisplacementFieldReaderType::New();
        fieldReader->SetFileName( initialTransformName.c_str() );
        fieldReader->Update();

        typename DisplacementFieldTransformType::Pointer displacementFieldTransform =
          DisplacementFieldTransformType::New();
        displacementFieldTransform->SetDisplacementField( fieldReader->GetOutput() );
        initialTransform = dynamic_cast<TransformType *>( displacementFieldTransform.GetPointer() );

        hasTransformBeenRead = true;
        }
      catch( ... )
        {
        hasTransformBeenRead = false;
        }

      if( !hasTransformBeenRead )
        {
        try
          {
          typedef itk::TransformFileReader TransformReaderType;
          typename TransformReaderType::Pointer initialTransformReader
            = TransformReaderType::New();

          if( initialTransformOption->GetNumberOfParameters( n ) == 0 )
            {
            initialTransformName = initialTransformOption->GetValue( n );
            initialTransformReader->SetFileName( initialTransformName.c_str() );
            initialTransformReader->Update();
            initialTransform = dynamic_cast<TransformType *>(
                ( ( initialTransformReader->GetTransformList() )->front() ).GetPointer() );
            }
          else
            {
            initialTransformName = initialTransformOption->GetParameter( n, 0 );
            initialTransformReader->SetFileName( initialTransformName.c_str() );
            initialTransformReader->Update();
            initialTransform = dynamic_cast<TransformType *>(
                ( ( initialTransformReader->GetTransformList() )->front() ).GetPointer() );
            if( ( initialTransformOption->GetNumberOfParameters( n ) > 1 ) &&
                parser->Convert<bool>( initialTransformOption->GetParameter( n, 1 ) ) )
              {
              initialTransform = dynamic_cast<TransformType *>(
                  initialTransform->GetInverseTransform().GetPointer() );
              if( !initialTransform )
                {
                std::cerr << "Inverse does not exist for " << initialTransformName
                          << std::endl;
                return EXIT_FAILURE;
                }
              initialTransformName = std::string( "inverse of " ) + initialTransformName;
              }
            }
          }
        catch( const itk::ExceptionObject & e )
          {
          std::cerr << "Transform reader for "
                    << initialTransformName << " caught an ITK exception:\n";
          e.Print( std::cerr );
          return EXIT_FAILURE;
          }
        catch( const std::exception & e )
          {
          std::cerr << "Transform reader for "
                    << initialTransformName << " caught an exception:\n";
          std::cerr << e.what() << std::endl;
          return EXIT_FAILURE;
          }
        catch( ... )
          {
          std::cerr << "Transform reader for "
                    << initialTransformName << " caught an unknown exception!!!\n";
          return EXIT_FAILURE;
          }
        }
      compositeTransform->AddTransform( initialTransform );
      initialTransformNames.push_back( initialTransformName );
      initialTransformTypes.push_back( initialTransform->GetNameOfClass() );
      }
    std::cout << "Initializing with the following transforms " << "(in order): " << std::endl;
    for( unsigned int n = 0; n < initialTransformNames.size(); n++ )
      {
      std::cout << "  " << n + 1 << ". " << initialTransformNames[n] << " (type = "
                << initialTransformTypes[n] << ")" << std::endl;
      }
    }
  // We iterate backwards because the command line options are stored as a stack (first in last out)
  for( int currentStage = numberOfStages - 1; currentStage >= 0; currentStage-- )
    {
    itk::TimeProbe timer;
    timer.Start();

    typedef itk::ImageRegistrationMethodv4<ImageType, ImageType> AffineRegistrationType;

    std::cout << std::endl << "Stage "
              << ( numberOfInitialTransforms + numberOfStages - currentStage - 1 ) << std::endl;
    std::stringstream currentStageString;
    currentStageString << ( numberOfInitialTransforms + numberOfStages - currentStage - 1 );

    // Get the fixed and moving images

    std::string fixedImageFileName = metricOption->GetParameter( currentStage, 0 );
    std::string movingImageFileName = metricOption->GetParameter( currentStage, 1 );

    std::cout << "  fixed image: " << fixedImageFileName << std::endl;
    std::cout << "  moving image: " << movingImageFileName << std::endl;

    typedef itk::ImageFileReader<ImageType> ImageReaderType;
    typename ImageReaderType::Pointer fixedImageReader = ImageReaderType::New();
    fixedImageReader->SetFileName( fixedImageFileName.c_str() );
    fixedImageReader->Update();
    typename ImageType::Pointer fixedImage = fixedImageReader->GetOutput();
    try
      {
      fixedImage->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
      }
    fixedImage->DisconnectPipeline();

    typename ImageReaderType::Pointer movingImageReader = ImageReaderType::New();
    movingImageReader->SetFileName( movingImageFileName.c_str() );
    movingImageReader->Update();
    typename ImageType::Pointer movingImage = movingImageReader->GetOutput();
    try
      {
      movingImage->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
      }
    movingImage->DisconnectPipeline();

    // Preprocess images

    PreprocessImages<ImageType>( fixedImage, movingImage, parser );

    // Get the number of iterations and use that information to specify the number of levels

    std::vector<unsigned int> iterations = parser->ConvertVector<unsigned int>( iterationsOption->GetValue(
                                                                                  currentStage ) );
    std::cout << "  iterations = " << iterationsOption->GetValue( currentStage ) << std::endl;
    unsigned int numberOfLevels = iterations.size();
    std::cout << "  number of levels = " << numberOfLevels << std::endl;

    // Get shrink factors

    std::vector<unsigned int> factors = parser->ConvertVector<unsigned int>( shrinkFactorsOption->GetValue(
                                                                               currentStage ) );
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
      std::cout << "  shrink factors per level: " << shrinkFactorsPerLevel << std::endl;
      }

    // Get smoothing sigmas

    std::vector<float> sigmas = parser->ConvertVector<float>( smoothingSigmasOption->GetValue( currentStage ) );
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
      std::cout << "  smoothing sigmas per level: " << smoothingSigmasPerLevel << std::endl;
      }

    // Set up the image metric and scales estimator

    typedef itk::ImageToImageMetricv4<ImageType, ImageType> MetricType;
    typename MetricType::Pointer metric;

    std::string whichMetric = metricOption->GetValue( currentStage );
    ConvertToLowerCase( whichMetric );

    float samplingPercentage = 1.0;
    if( metricOption->GetNumberOfParameters( currentStage ) > 5 )
      {
      samplingPercentage = parser->Convert<float>( metricOption->GetParameter( currentStage, 5 ) );
      }

    std::string samplingStrategy = "";
    if( metricOption->GetNumberOfParameters( currentStage ) > 4 )
      {
      samplingStrategy = metricOption->GetParameter( currentStage, 4 );
      }
    ConvertToLowerCase( samplingStrategy );
    typename AffineRegistrationType::MetricSamplingStrategyType metricSamplingStrategy = AffineRegistrationType::NONE;
    if( std::strcmp( samplingStrategy.c_str(), "random" ) == 0 )
      {
      std::cout << "  random sampling (percentage = " << samplingPercentage << ")" << std::endl;
      metricSamplingStrategy = AffineRegistrationType::RANDOM;
      }
    if( std::strcmp( samplingStrategy.c_str(), "regular" ) == 0 )
      {
      std::cout << "  regular sampling (percentage = " << samplingPercentage << ")" << std::endl;
      metricSamplingStrategy = AffineRegistrationType::REGULAR;
      }

    if( std::strcmp( whichMetric.c_str(), "cc" ) == 0 )
      {
      unsigned int radiusOption = parser->Convert<unsigned int>( metricOption->GetParameter( currentStage, 3 ) );

      std::cout << "  using the CC metric (radius = " << radiusOption << ")" << std::endl;
      typedef itk::ANTSNeighborhoodCorrelationImageToImageMetricv4<ImageType, ImageType> CorrelationMetricType;
      typename CorrelationMetricType::Pointer correlationMetric = CorrelationMetricType::New();
      typename CorrelationMetricType::RadiusType radius;
      radius.Fill( radiusOption );
      correlationMetric->SetRadius( radius );
      correlationMetric->SetDoFixedImagePreWarp( false );
      correlationMetric->SetDoMovingImagePreWarp( false );
      correlationMetric->SetUseMovingImageGradientFilter( false );
      correlationMetric->SetUseFixedImageGradientFilter( false );

      metric = correlationMetric;
      }
    else if( std::strcmp( whichMetric.c_str(), "mi" ) == 0 )
      {
      unsigned int binOption = parser->Convert<unsigned int>( metricOption->GetParameter( currentStage, 3 ) );
      std::cout << "  using the Mattes MI metric (number of bins = " << binOption << ")" << std::endl;
      typedef itk::MattesMutualInformationImageToImageMetricv4<ImageType, ImageType> MutualInformationMetricType;
      typename MutualInformationMetricType::Pointer mutualInformationMetric = MutualInformationMetricType::New();
      mutualInformationMetric = mutualInformationMetric;
      mutualInformationMetric->SetNumberOfHistogramBins( binOption );
      mutualInformationMetric->SetDoFixedImagePreWarp( false );
      mutualInformationMetric->SetDoMovingImagePreWarp( false );
      mutualInformationMetric->SetUseMovingImageGradientFilter( false );
      mutualInformationMetric->SetUseFixedImageGradientFilter( false );
      mutualInformationMetric->SetUseFixedSampledPointSet( false );
      metric = mutualInformationMetric;
      }
    else if( std::strcmp( whichMetric.c_str(), "mi2" ) == 0 )
      {
      unsigned int binOption = parser->Convert<unsigned int>( metricOption->GetParameter( currentStage, 3 ) );

      std::cout << "  using the MI metric (number of bins = " << binOption << ")" << std::endl;
      typedef itk::JointHistogramMutualInformationImageToImageMetricv4<ImageType,
                                                                       ImageType> MutualInformationMetricType;
      typename MutualInformationMetricType::Pointer mutualInformationMetric = MutualInformationMetricType::New();
      mutualInformationMetric = mutualInformationMetric;
      mutualInformationMetric->SetNumberOfHistogramBins( binOption );
      mutualInformationMetric->SetDoFixedImagePreWarp( false );
      mutualInformationMetric->SetDoMovingImagePreWarp( false );
      mutualInformationMetric->SetUseMovingImageGradientFilter( false );
      mutualInformationMetric->SetUseFixedImageGradientFilter( false );
      mutualInformationMetric->SetUseFixedSampledPointSet( false );
      mutualInformationMetric->SetVarianceForJointPDFSmoothing( 1.0 );
      metric = mutualInformationMetric;
      }
    else if( std::strcmp( whichMetric.c_str(), "demons" ) == 0 )
      {
      std::cout << "  using the Demons metric." << std::endl;

      typedef itk::DemonsImageToImageMetricv4<ImageType, ImageType> DemonsMetricType;
      typename DemonsMetricType::Pointer demonsMetric = DemonsMetricType::New();
      demonsMetric = demonsMetric;
      demonsMetric->SetDoFixedImagePreWarp( false );
      demonsMetric->SetDoMovingImagePreWarp( false );

      metric = demonsMetric;
      }
    else if( std::strcmp( whichMetric.c_str(), "gc" ) == 0 )
      {
      std::cout << "  using the global correlation metric." << std::endl;
      typedef itk::CorrelationImageToImageMetricv4<ImageType, ImageType> corrMetricType;
      typename corrMetricType::Pointer corrMetric = corrMetricType::New();
      corrMetric->SetDoFixedImagePreWarp( false );
      corrMetric->SetDoMovingImagePreWarp( false );
      metric = corrMetric;
      }
    else
      {
      std::cerr << "ERROR: Unrecognized image metric: " << whichMetric << std::endl;
      }
    /** Can really impact performance */
    bool gaussian = false;
    metric->SetUseMovingImageGradientFilter( gaussian );
    metric->SetUseFixedImageGradientFilter( gaussian );

    // Set up the optimizer.  To change the iteration number for each level we rely
    // on the command observer.

    float learningRate = parser->Convert<float>( transformOption->GetParameter( currentStage, 0 ) );
    typedef itk::RegistrationParameterScalesFromShift<MetricType> ScalesEstimatorType;
    typename ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();
    scalesEstimator->SetMetric( metric );
    scalesEstimator->SetTransformForward( true );

    typedef itk::GradientDescentOptimizerv4 GradientDescentOptimizerType;
    typename GradientDescentOptimizerType::Pointer optimizer = GradientDescentOptimizerType::New();
    optimizer->SetLearningRate( learningRate );
    optimizer->SetMaximumStepSizeInPhysicalUnits( learningRate );
    optimizer->SetNumberOfIterations( iterations[0] );
    optimizer->SetScalesEstimator( scalesEstimator );
    //    optimizer->SetMinimumConvergenceValue( -1 );
    //    optimizer->SetConvergenceWindowSize( 10 );

    // Set up the image registration methods along with the transforms

    std::string whichTransform = transformOption->GetValue( currentStage );
    ConvertToLowerCase( whichTransform );
    if( std::strcmp( whichTransform.c_str(), "affine" ) == 0 )
      {
      typename AffineRegistrationType::Pointer affineRegistration = AffineRegistrationType::New();

      typedef itk::AffineTransform<double, ImageDimension> AffineTransformType;
      typename AffineTransformType::Pointer affineTransform = AffineTransformType::New();

      affineRegistration->SetFixedImage( fixedImage );
      affineRegistration->SetMovingImage( movingImage );
      affineRegistration->SetNumberOfLevels( numberOfLevels );
      affineRegistration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
      affineRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
      affineRegistration->SetMetric( metric );
      affineRegistration->SetMetricSamplingStrategy( metricSamplingStrategy );
      affineRegistration->SetMetricSamplingPercentage( samplingPercentage );
      affineRegistration->SetOptimizer( optimizer );
      affineRegistration->SetTransform( affineTransform );
      affineRegistration->SetCompositeTransform( compositeTransform );

      typedef CommandIterationUpdate<AffineRegistrationType> AffineCommandType;
      typename AffineCommandType::Pointer affineObserver = AffineCommandType::New();
      affineObserver->SetNumberOfIterations( iterations );

      affineRegistration->AddObserver( itk::IterationEvent(), affineObserver );

      try
        {
        std::cout << std::endl << "*** Running affine registration ***" << std::endl << std::endl;
        affineObserver->Execute( affineRegistration, itk::StartEvent() );
        affineRegistration->StartRegistration();
        }
      catch( itk::ExceptionObject & e )
        {
        std::cerr << "Exception caught: " << e << std::endl;
        return EXIT_FAILURE;
        }

      // Write out the affine transform

      std::string filename = outputPrefix + currentStageString.str() + std::string( "Affine.mat" );

      typedef itk::TransformFileWriter TransformWriterType;
      typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
      transformWriter->SetInput( affineRegistration->GetOutput()->Get() );
      transformWriter->SetFileName( filename.c_str() );
      transformWriter->Update();
      }
    else if( std::strcmp( whichTransform.c_str(), "rigid" ) == 0 )
      {
      typedef typename RigidTransformTraits<ImageDimension>::TransformType RigidTransformType;
      typename RigidTransformType::Pointer rigidTransform = RigidTransformType::New();

      typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, RigidTransformType> RigidRegistrationType;
      typename RigidRegistrationType::Pointer rigidRegistration = RigidRegistrationType::New();

      rigidRegistration->SetFixedImage( fixedImage );
      rigidRegistration->SetMovingImage( movingImage );
      rigidRegistration->SetNumberOfLevels( numberOfLevels );
      rigidRegistration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
      rigidRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
      rigidRegistration->SetMetric( metric );
      rigidRegistration->SetMetricSamplingStrategy(
        static_cast<typename RigidRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
      rigidRegistration->SetMetricSamplingPercentage( samplingPercentage );
      rigidRegistration->SetOptimizer( optimizer );
      rigidRegistration->SetTransform( rigidTransform );
      rigidRegistration->SetCompositeTransform( compositeTransform );

      typedef CommandIterationUpdate<RigidRegistrationType> RigidCommandType;
      typename RigidCommandType::Pointer rigidObserver = RigidCommandType::New();
      rigidObserver->SetNumberOfIterations( iterations );

      rigidRegistration->AddObserver( itk::IterationEvent(), rigidObserver );

      try
        {
        std::cout << std::endl << "*** Running rigid registration ***" << std::endl << std::endl;
        rigidObserver->Execute( rigidRegistration, itk::StartEvent() );
        rigidRegistration->StartRegistration();
        }
      catch( itk::ExceptionObject & e )
        {
        std::cerr << "Exception caught: " << e << std::endl;
        return EXIT_FAILURE;
        }
      // Write out the affine transform

      std::string filename = outputPrefix + currentStageString.str() + std::string( "Rigid.mat" );

      typedef itk::TransformFileWriter TransformWriterType;
      typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
      transformWriter->SetInput( rigidRegistration->GetOutput()->Get() );
      transformWriter->SetFileName( filename.c_str() );
      transformWriter->Update();
      }
    else if( std::strcmp( whichTransform.c_str(),
                          "compositeaffine" ) == 0 || std::strcmp( whichTransform.c_str(), "compaff" ) == 0 )
      {
      typedef typename CompositeAffineTransformTraits<ImageDimension>::TransformType CompositeAffineTransformType;
      typename CompositeAffineTransformType::Pointer compositeAffineTransform = CompositeAffineTransformType::New();
      typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, CompositeAffineTransformType> AffineRegistrationType;
      typename AffineRegistrationType::Pointer affineRegistration = AffineRegistrationType::New();

      affineRegistration->SetFixedImage( fixedImage );
      affineRegistration->SetMovingImage( movingImage );
      affineRegistration->SetNumberOfLevels( numberOfLevels );
      affineRegistration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
      affineRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
      affineRegistration->SetMetric( metric );
      affineRegistration->SetMetricSamplingStrategy(
        static_cast<typename AffineRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
      affineRegistration->SetMetricSamplingPercentage( samplingPercentage );
      affineRegistration->SetOptimizer( optimizer );
      affineRegistration->SetTransform( compositeAffineTransform );
      affineRegistration->SetCompositeTransform( compositeTransform );

      typedef CommandIterationUpdate<AffineRegistrationType> AffineCommandType;
      typename AffineCommandType::Pointer affineObserver = AffineCommandType::New();
      affineObserver->SetNumberOfIterations( iterations );

      affineRegistration->AddObserver( itk::IterationEvent(), affineObserver );

      try
        {
        std::cout << std::endl << "*** Running composite affine registration ***" << std::endl << std::endl;
        affineObserver->Execute( affineRegistration, itk::StartEvent() );
        affineRegistration->StartRegistration();
        }
      catch( itk::ExceptionObject & e )
        {
        std::cerr << "Exception caught: " << e << std::endl;
        return EXIT_FAILURE;
        }
      // Write out the affine transform

      std::string filename = outputPrefix + currentStageString.str() + std::string( "Affine.mat" );

      typedef itk::TransformFileWriter TransformWriterType;
      typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
      transformWriter->SetInput( affineRegistration->GetOutput()->Get() );
      transformWriter->SetFileName( filename.c_str() );
      transformWriter->Update();
      }
    else if( std::strcmp( whichTransform.c_str(), "similarity" ) == 0 )
      {
      typedef typename SimilarityTransformTraits<ImageDimension>::TransformType SimilarityTransformType;
      typename SimilarityTransformType::Pointer similarityTransform = SimilarityTransformType::New();

      typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, SimilarityTransformType> SimilarityRegistrationType;
      typename SimilarityRegistrationType::Pointer similarityRegistration = SimilarityRegistrationType::New();

      similarityRegistration->SetFixedImage( fixedImage );
      similarityRegistration->SetMovingImage( movingImage );
      similarityRegistration->SetNumberOfLevels( numberOfLevels );
      similarityRegistration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
      similarityRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
      similarityRegistration->SetMetric( metric );
      similarityRegistration->SetMetricSamplingStrategy(
        static_cast<typename SimilarityRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
      similarityRegistration->SetMetricSamplingPercentage( samplingPercentage );
      similarityRegistration->SetOptimizer( optimizer );
      similarityRegistration->SetTransform( similarityTransform );
      similarityRegistration->SetCompositeTransform( compositeTransform );

      typedef CommandIterationUpdate<SimilarityRegistrationType> SimilarityCommandType;
      typename SimilarityCommandType::Pointer similarityObserver = SimilarityCommandType::New();
      similarityObserver->SetNumberOfIterations( iterations );

      similarityRegistration->AddObserver( itk::IterationEvent(), similarityObserver );

      try
        {
        std::cout << std::endl << "*** Running similarity registration ***" << std::endl << std::endl;
        similarityObserver->Execute( similarityRegistration, itk::StartEvent() );
        similarityRegistration->StartRegistration();
        }
      catch( itk::ExceptionObject & e )
        {
        std::cerr << "Exception caught: " << e << std::endl;
        return EXIT_FAILURE;
        }
      // Write out the affine transform

      std::string filename = outputPrefix + currentStageString.str() + std::string( "Similarity.mat" );

      typedef itk::TransformFileWriter TransformWriterType;
      typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
      transformWriter->SetInput( similarityRegistration->GetOutput()->Get() );
      transformWriter->SetFileName( filename.c_str() );
      transformWriter->Update();
      }
    else if( std::strcmp( whichTransform.c_str(), "translation" ) == 0 )
      {
      typedef itk::TranslationTransform<RealType, ImageDimension> TranslationTransformType;
      typename TranslationTransformType::Pointer translationTransform = TranslationTransformType::New();

      typedef itk::ImageRegistrationMethodv4<ImageType, ImageType,
                                             TranslationTransformType> TranslationRegistrationType;
      typename TranslationRegistrationType::Pointer translationRegistration = TranslationRegistrationType::New();

      translationRegistration->SetFixedImage( fixedImage );
      translationRegistration->SetMovingImage( movingImage );
      translationRegistration->SetNumberOfLevels( numberOfLevels );
      translationRegistration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
      translationRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
      translationRegistration->SetMetric( metric );
      translationRegistration->SetMetricSamplingStrategy(
        static_cast<typename TranslationRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
      translationRegistration->SetMetricSamplingPercentage( samplingPercentage );
      translationRegistration->SetOptimizer( optimizer );
      translationRegistration->SetTransform( translationTransform );
      translationRegistration->SetCompositeTransform( compositeTransform );

      typedef CommandIterationUpdate<TranslationRegistrationType> TranslationCommandType;
      typename TranslationCommandType::Pointer translationObserver = TranslationCommandType::New();
      translationObserver->SetNumberOfIterations( iterations );

      translationRegistration->AddObserver( itk::IterationEvent(), translationObserver );

      try
        {
        std::cout << std::endl << "*** Running translation registration ***" << std::endl << std::endl;
        translationObserver->Execute( translationRegistration, itk::StartEvent() );
        translationRegistration->StartRegistration();
        }
      catch( itk::ExceptionObject & e )
        {
        std::cerr << "Exception caught: " << e << std::endl;
        return EXIT_FAILURE;
        }
      // Write out the affine transform

      std::string filename = outputPrefix + currentStageString.str() + std::string( "Translation.mat" );

      typedef itk::TransformFileWriter TransformWriterType;
      typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
      transformWriter->SetInput( translationRegistration->GetOutput()->Get() );
      transformWriter->SetFileName( filename.c_str() );
      transformWriter->Update();
      }
    else if( std::strcmp( whichTransform.c_str(),
                          "gaussiandisplacementfield" ) == 0 ||  std::strcmp( whichTransform.c_str(), "gdf" ) == 0 )
      {
      metric->SetDoFixedImagePreWarp( true );
      metric->SetDoMovingImagePreWarp( true );
      typedef itk::Vector<RealType, ImageDimension> VectorType;
      VectorType zeroVector( 0.0 );
      typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;
      typename DisplacementFieldType::Pointer displacementField = DisplacementFieldType::New();
      displacementField->CopyInformation( fixedImage );
      displacementField->SetRegions( fixedImage->GetBufferedRegion() );
      displacementField->Allocate();
      displacementField->FillBuffer( zeroVector );

      typedef itk::GaussianSmoothingOnUpdateDisplacementFieldTransform<RealType,
                                                                       ImageDimension> DisplacementFieldTransformType;

      typedef itk::ImageRegistrationMethodv4<ImageType, ImageType,
                                             DisplacementFieldTransformType> DisplacementFieldRegistrationType;

      // Create the transform adaptors

      typedef itk::DisplacementFieldTransformParametersAdaptor<DisplacementFieldTransformType>
        DisplacementFieldTransformAdaptorType;
      typename DisplacementFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;

      // Extract parameters

      RealType sigmaForUpdateField = parser->Convert<float>( transformOption->GetParameter( currentStage, 1 ) );
      RealType sigmaForTotalField = parser->Convert<float>( transformOption->GetParameter( currentStage, 2 ) );

      typedef itk::GaussianSmoothingOnUpdateDisplacementFieldTransform<RealType,
                                                                       ImageDimension>
        GaussianDisplacementFieldTransformType;
      typename GaussianDisplacementFieldTransformType::Pointer gaussianFieldTransform =
        GaussianDisplacementFieldTransformType::New();
      gaussianFieldTransform->SetGaussianSmoothingVarianceForTheUpdateField( sigmaForUpdateField );
      gaussianFieldTransform->SetGaussianSmoothingVarianceForTheTotalField( sigmaForTotalField );
      gaussianFieldTransform->SetDisplacementField( displacementField );
      // Create the transform adaptors
      // For the gaussian displacement field, the specified variances are in image spacing terms
      // and, in normal practice, we typically don't change these values at each level.  However,
      // if the user wishes to add that option, they can use the class
      // GaussianSmoothingOnUpdateDisplacementFieldTransformAdaptor
      for( unsigned int level = 0; level < numberOfLevels; level++ )
        {
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
        fieldTransformAdaptor->SetTransform( gaussianFieldTransform );

        adaptors.push_back( fieldTransformAdaptor.GetPointer() );
        }

      typename DisplacementFieldRegistrationType::Pointer displacementFieldRegistration =
        DisplacementFieldRegistrationType::New();
      displacementFieldRegistration->SetFixedImage( fixedImage );
      displacementFieldRegistration->SetMovingImage( movingImage );
      displacementFieldRegistration->SetNumberOfLevels( numberOfLevels );
      displacementFieldRegistration->SetCompositeTransform( compositeTransform );
      displacementFieldRegistration->SetTransform( gaussianFieldTransform );
      displacementFieldRegistration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
      displacementFieldRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
      displacementFieldRegistration->SetMetric( metric );
      displacementFieldRegistration->SetMetricSamplingStrategy(
        static_cast<typename DisplacementFieldRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
      displacementFieldRegistration->SetMetricSamplingPercentage( samplingPercentage );
      displacementFieldRegistration->SetOptimizer( optimizer );
      displacementFieldRegistration->SetTransformParametersAdaptorsPerLevel( adaptors );

      typedef CommandIterationUpdate<DisplacementFieldRegistrationType> DisplacementFieldCommandType;
      typename DisplacementFieldCommandType::Pointer displacementFieldRegistrationObserver =
        DisplacementFieldCommandType::New();
      displacementFieldRegistrationObserver->SetNumberOfIterations( iterations );

      displacementFieldRegistration->AddObserver( itk::IterationEvent(), displacementFieldRegistrationObserver );

      try
        {
        std::cout << std::endl << "*** Running gaussian displacement field registration (sigmaForUpdateField = "
                  << sigmaForUpdateField << ", sigmaForTotalField = " << sigmaForTotalField << ") ***" << std::endl
                  << std::endl;
        displacementFieldRegistrationObserver->Execute( displacementFieldRegistration, itk::StartEvent() );
        displacementFieldRegistration->StartRegistration();
        }
      catch( itk::ExceptionObject & e )
        {
        std::cerr << "Exception caught: " << e << std::endl;
        return EXIT_FAILURE;
        }

      // Write out the displacement field

      std::string filename = outputPrefix + currentStageString.str() + std::string( "Warp.nii.gz" );

      typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( const_cast<typename DisplacementFieldRegistrationType::TransformType *>(
                          displacementFieldRegistration->GetOutput()->Get() )->GetDisplacementField() );
      writer->SetFileName( filename.c_str() );
      writer->Update();
      }
    else if( std::strcmp( whichTransform.c_str(),
                          "bsplinedisplacementfield" ) == 0 || std::strcmp( whichTransform.c_str(), "dmffd" ) == 0 )
      {
      typedef itk::Vector<RealType, ImageDimension> VectorType;
      VectorType zeroVector( 0.0 );
      typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;
      typename DisplacementFieldType::Pointer displacementField = DisplacementFieldType::New();
      displacementField->CopyInformation( fixedImage );
      displacementField->SetRegions( fixedImage->GetBufferedRegion() );
      displacementField->Allocate();
      displacementField->FillBuffer( zeroVector );

      typedef itk::BSplineSmoothingOnUpdateDisplacementFieldTransform<RealType,
                                                                      ImageDimension> DisplacementFieldTransformType;
      typename DisplacementFieldTransformType::Pointer bsplineFieldTransform = DisplacementFieldTransformType::New();

      typedef itk::ImageRegistrationMethodv4<ImageType, ImageType,
                                             DisplacementFieldTransformType> DisplacementFieldRegistrationType;

      // Create the transform adaptors

      typedef itk::DisplacementFieldTransformParametersAdaptor<DisplacementFieldTransformType>
        DisplacementFieldTransformAdaptorType;
      typename DisplacementFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;

      // Extract parameters

      std::vector<unsigned int> meshSizeForTheUpdateField = parser->ConvertVector<unsigned int>( transformOption->GetParameter(
                                                                                                   currentStage,
                                                                                                   1 ) );
      std::vector<unsigned int> meshSizeForTheTotalField = parser->ConvertVector<unsigned int>( transformOption->GetParameter(
                                                                                                  currentStage,
                                                                                                  2 ) );

      if( meshSizeForTheUpdateField.size() != ImageDimension || meshSizeForTheTotalField.size() != ImageDimension )
        {
        std::cerr << "ERROR:  The mesh size(s) don't match the ImageDimension." << std::endl;
        return EXIT_FAILURE;
        }

      unsigned int splineOrder = 3;
      if( transformOption->GetNumberOfParameters( currentStage ) > 3 )
        {
        splineOrder = parser->Convert<unsigned int>( transformOption->GetParameter( currentStage, 3 ) );
        }

      typename DisplacementFieldTransformType::ArrayType updateMeshSize;
      typename DisplacementFieldTransformType::ArrayType totalMeshSize;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        updateMeshSize[d] = meshSizeForTheUpdateField[d];
        totalMeshSize[d] = meshSizeForTheTotalField[d];
        }

      bsplineFieldTransform->SetMeshSizeForTheUpdateField( updateMeshSize );
      bsplineFieldTransform->SetMeshSizeForTheTotalField( totalMeshSize );
      bsplineFieldTransform->SetSplineOrder( splineOrder );
      bsplineFieldTransform->SetDisplacementField( displacementField );
      // Create the transform adaptors specific to B-splines
      for( unsigned int level = 0; level < numberOfLevels; level++ )
        {
        // We use the shrink image filter to calculate the fixed parameters of the virtual
        // domain at each level.  To speed up calculation and avoid unnecessary memory
        // usage, we could calculate these fixed parameters directly.

        typedef itk::ShrinkImageFilter<DisplacementFieldType, DisplacementFieldType> ShrinkFilterType;
        typename ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
        shrinkFilter->SetShrinkFactors( shrinkFactorsPerLevel[level] );
        shrinkFilter->SetInput( displacementField );
        shrinkFilter->Update();

        typedef itk::BSplineSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor<DisplacementFieldTransformType>
          BSplineDisplacementFieldTransformAdaptorType;
        typename BSplineDisplacementFieldTransformAdaptorType::Pointer bsplineFieldTransformAdaptor =
          BSplineDisplacementFieldTransformAdaptorType::New();
        bsplineFieldTransformAdaptor->SetRequiredSpacing( shrinkFilter->GetOutput()->GetSpacing() );
        bsplineFieldTransformAdaptor->SetRequiredSize( shrinkFilter->GetOutput()->GetBufferedRegion().GetSize() );
        bsplineFieldTransformAdaptor->SetRequiredDirection( shrinkFilter->GetOutput()->GetDirection() );
        bsplineFieldTransformAdaptor->SetRequiredOrigin( shrinkFilter->GetOutput()->GetOrigin() );
        bsplineFieldTransformAdaptor->SetTransform( bsplineFieldTransform );

        // A good heuristic is to double the b-spline mesh resolution at each level
        typename DisplacementFieldTransformType::ArrayType newUpdateMeshSize = updateMeshSize;
        typename DisplacementFieldTransformType::ArrayType newTotalMeshSize = totalMeshSize;
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          newUpdateMeshSize[d] = newUpdateMeshSize[d] << ( level + 1 );
          newTotalMeshSize[d] = newTotalMeshSize[d] << ( level + 1 );
          }
        bsplineFieldTransformAdaptor->SetMeshSizeForTheUpdateField( newUpdateMeshSize );
        bsplineFieldTransformAdaptor->SetMeshSizeForTheTotalField( newTotalMeshSize );

        adaptors.push_back( bsplineFieldTransformAdaptor.GetPointer() );
        }

      typename DisplacementFieldRegistrationType::Pointer displacementFieldRegistration =
        DisplacementFieldRegistrationType::New();
      displacementFieldRegistration->SetFixedImage( fixedImage );
      displacementFieldRegistration->SetMovingImage( movingImage );
      displacementFieldRegistration->SetNumberOfLevels( numberOfLevels );
      displacementFieldRegistration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
      displacementFieldRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
      displacementFieldRegistration->SetCompositeTransform( compositeTransform );
      displacementFieldRegistration->SetTransform( bsplineFieldTransform );
      displacementFieldRegistration->SetMetric( metric );
      displacementFieldRegistration->SetMetricSamplingStrategy(
        static_cast<typename DisplacementFieldRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
      displacementFieldRegistration->SetMetricSamplingPercentage( samplingPercentage );
      displacementFieldRegistration->SetOptimizer( optimizer );
      displacementFieldRegistration->SetTransformParametersAdaptorsPerLevel( adaptors );

      typedef CommandIterationUpdate<DisplacementFieldRegistrationType> DisplacementFieldCommandType;
      typename DisplacementFieldCommandType::Pointer displacementFieldRegistrationObserver =
        DisplacementFieldCommandType::New();
      displacementFieldRegistrationObserver->SetNumberOfIterations( iterations );

      displacementFieldRegistration->AddObserver( itk::IterationEvent(), displacementFieldRegistrationObserver );

      try
        {
        std::cout << std::endl << "*** Running bspline displacement field registration (updateMeshSizeAtBaseLevel = "
                  << updateMeshSize << ", totalMeshSizeAtBaseLevel = " << totalMeshSize << ") ***" << std::endl
                  << std::endl;
        displacementFieldRegistrationObserver->Execute( displacementFieldRegistration, itk::StartEvent() );
        displacementFieldRegistration->StartRegistration();
        }
      catch( itk::ExceptionObject & e )
        {
        std::cerr << "Exception caught: " << e << std::endl;
        return EXIT_FAILURE;
        }

      // Write out the displacement field

      std::string filename = outputPrefix + currentStageString.str() + std::string( "Warp.nii.gz" );

      typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( const_cast<typename DisplacementFieldRegistrationType::TransformType *>(
                          displacementFieldRegistration->GetOutput()->Get() )->GetDisplacementField() );
      writer->SetFileName( filename.c_str() );
      writer->Update();
      }
    else if( std::strcmp( whichTransform.c_str(),
                          "bspline" ) == 0 || std::strcmp( whichTransform.c_str(), "ffd" ) == 0 )
      {
      const unsigned int SplineOrder = 3;
      typedef itk::BSplineTransform<RealType, ImageDimension, SplineOrder> BSplineTransformType;
      typename BSplineTransformType::Pointer bsplineTransform = BSplineTransformType::New();

      typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, BSplineTransformType> BSplineRegistrationType;

      std::vector<unsigned int> size =
        parser->ConvertVector<unsigned int>( transformOption->GetParameter( currentStage, 1 ) );

      typename BSplineTransformType::PhysicalDimensionsType physicalDimensions;
      typename BSplineTransformType::MeshSizeType meshSize;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        physicalDimensions[d] = fixedImage->GetSpacing()[d]
          * static_cast<RealType>( fixedImage->GetLargestPossibleRegion().GetSize()[d] - 1 );
        meshSize[d] = size[d];
        }

      bsplineTransform->SetTransformDomainOrigin( fixedImage->GetOrigin() );
      bsplineTransform->SetTransformDomainPhysicalDimensions( physicalDimensions );
      bsplineTransform->SetTransformDomainMeshSize( meshSize );
      bsplineTransform->SetTransformDomainDirection( fixedImage->GetDirection() );
      bsplineTransform->SetIdentity();

      // Create the transform adaptors

      typedef itk::BSplineTransformParametersAdaptor<BSplineTransformType> BSplineTransformAdaptorType;
      typename BSplineRegistrationType::TransformParametersAdaptorsContainerType adaptors;
      // Create the transform adaptors specific to B-splines
      for( unsigned int level = 0; level < numberOfLevels; level++ )
        {
        // A good heuristic is to double the b-spline mesh resolution at each level

        typename BSplineTransformType::MeshSizeType requiredMeshSize;
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          requiredMeshSize[d] = meshSize[d] << level;
          }

        typedef itk::BSplineTransformParametersAdaptor<BSplineTransformType> BSplineAdaptorType;
        typename BSplineAdaptorType::Pointer bsplineAdaptor = BSplineAdaptorType::New();
        bsplineAdaptor->SetTransform( bsplineTransform );
        bsplineAdaptor->SetRequiredTransformDomainMeshSize( requiredMeshSize );
        bsplineAdaptor->SetRequiredTransformDomainOrigin( bsplineTransform->GetTransformDomainOrigin() );
        bsplineAdaptor->SetRequiredTransformDomainDirection( bsplineTransform->GetTransformDomainDirection() );
        bsplineAdaptor->SetRequiredTransformDomainPhysicalDimensions(
          bsplineTransform->GetTransformDomainPhysicalDimensions() );

        adaptors.push_back( bsplineAdaptor.GetPointer() );
        }

      optimizer->SetScalesEstimator( NULL );

      typename BSplineRegistrationType::Pointer bsplineRegistration = BSplineRegistrationType::New();
      bsplineRegistration->SetFixedImage( fixedImage );
      bsplineRegistration->SetMovingImage( movingImage );
      bsplineRegistration->SetNumberOfLevels( numberOfLevels );
      bsplineRegistration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
      bsplineRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
      bsplineRegistration->SetCompositeTransform( compositeTransform );
      bsplineRegistration->SetTransform( bsplineTransform );
      bsplineRegistration->SetMetric( metric );
      bsplineRegistration->SetMetricSamplingStrategy(
        static_cast<typename BSplineRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
      bsplineRegistration->SetMetricSamplingPercentage( samplingPercentage );
      bsplineRegistration->SetOptimizer( optimizer );
      bsplineRegistration->SetTransformParametersAdaptorsPerLevel( adaptors );

      typedef CommandIterationUpdate<BSplineRegistrationType> BSplineCommandType;
      typename BSplineCommandType::Pointer bsplineObserver = BSplineCommandType::New();
      bsplineObserver->SetNumberOfIterations( iterations );

      bsplineRegistration->AddObserver( itk::IterationEvent(), bsplineObserver );

      try
        {
        std::cout << std::endl << "*** Running bspline registration (meshSizeAtBaseLevel = " << meshSize << ") ***"
                  << std::endl << std::endl;
        bsplineObserver->Execute( bsplineRegistration, itk::StartEvent() );
        bsplineRegistration->StartRegistration();
        }
      catch( itk::ExceptionObject & e )
        {
        std::cerr << "Exception caught: " << e << std::endl;
        return EXIT_FAILURE;
        }

      // Write out B-spline transform

      std::string filename = outputPrefix + currentStageString.str() + std::string( "BSpline.txt" );

      typedef itk::TransformFileWriter TransformWriterType;
      typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
      transformWriter->SetInput( bsplineRegistration->GetOutput()->Get() );
      transformWriter->SetFileName( filename.c_str() );
      transformWriter->Update();
      }
    else if( std::strcmp( whichTransform.c_str(),
                          "timevaryingvelocityfield" ) == 0 || std::strcmp( whichTransform.c_str(), "tvf" ) == 0 )
      {
//      metric->SetDoFixedImagePreWarp( true );
//      metric->SetDoMovingImagePreWarp( true );
      typedef itk::Vector<RealType, ImageDimension> VectorType;
      VectorType zeroVector( 0.0 );

      // Determine the parameters (size, spacing, etc) for the time-varying velocity field

      typedef itk::Image<VectorType, ImageDimension + 1> TimeVaryingVelocityFieldType;
      typename TimeVaryingVelocityFieldType::Pointer velocityField = TimeVaryingVelocityFieldType::New();

      typename TimeVaryingVelocityFieldType::IndexType velocityFieldIndex;
      typename TimeVaryingVelocityFieldType::SizeType velocityFieldSize;
      typename TimeVaryingVelocityFieldType::PointType velocityFieldOrigin;
      typename TimeVaryingVelocityFieldType::SpacingType velocityFieldSpacing;
      typename TimeVaryingVelocityFieldType::DirectionType velocityFieldDirection;
      typename TimeVaryingVelocityFieldType::RegionType velocityFieldRegion;

      typename ImageType::IndexType fixedImageIndex = fixedImage->GetBufferedRegion().GetIndex();
      typename ImageType::SizeType fixedImageSize = fixedImage->GetBufferedRegion().GetSize();
      typename ImageType::PointType fixedImageOrigin = fixedImage->GetOrigin();
      typename ImageType::SpacingType fixedImageSpacing = fixedImage->GetSpacing();
      typename ImageType::DirectionType fixedImageDirection = fixedImage->GetDirection();

      unsigned int numberOfTimeIndices = parser->Convert<unsigned int>( transformOption->GetParameter( 0, 1 ) );

      velocityFieldIndex.Fill( 0 );
      velocityFieldSize.Fill( numberOfTimeIndices );
      velocityFieldOrigin.Fill( 0.0 );
      velocityFieldSpacing.Fill( 1.0 );
      velocityFieldDirection.SetIdentity();
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        velocityFieldIndex[i] = fixedImageIndex[i];
        velocityFieldSize[i] = fixedImageSize[i];
        velocityFieldOrigin[i] = fixedImageOrigin[i];
        velocityFieldSpacing[i] = fixedImageSpacing[i];
        for( unsigned int j = 0; j < ImageDimension; j++ )
          {
          velocityFieldDirection[i][j] = fixedImageDirection[i][j];
          }
        }

      velocityFieldRegion.SetSize( velocityFieldSize );
      velocityFieldRegion.SetIndex( velocityFieldIndex );

      velocityField->SetOrigin( velocityFieldOrigin );
      velocityField->SetSpacing( velocityFieldSpacing );
      velocityField->SetDirection( velocityFieldDirection );
      velocityField->SetRegions( velocityFieldRegion );
      velocityField->Allocate();
      velocityField->FillBuffer( zeroVector );

      // Extract parameters

      RealType sigmaForUpdateField = parser->Convert<float>( transformOption->GetParameter( currentStage, 2 ) );
      RealType sigmaForUpdateFieldTime = parser->Convert<float>( transformOption->GetParameter( currentStage, 3 ) );
      RealType sigmaForTotalField = parser->Convert<float>( transformOption->GetParameter( currentStage, 4 ) );
      RealType sigmaForTotalFieldTime = parser->Convert<float>( transformOption->GetParameter( currentStage, 5 ) );

      typedef itk::TimeVaryingVelocityFieldImageRegistrationMethodv4<ImageType,
                                                                     ImageType> VelocityFieldRegistrationType;
      typename VelocityFieldRegistrationType::Pointer velocityFieldRegistration = VelocityFieldRegistrationType::New();
      velocityFieldRegistration->SetFixedImage( fixedImage );
      velocityFieldRegistration->SetMovingImage( movingImage );
      velocityFieldRegistration->SetNumberOfLevels( numberOfLevels );
      velocityFieldRegistration->SetCompositeTransform( compositeTransform );
      velocityFieldRegistration->SetMetric( metric );
      velocityFieldRegistration->SetMetricSamplingStrategy(
        static_cast<typename VelocityFieldRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
      velocityFieldRegistration->SetMetricSamplingPercentage( samplingPercentage );
      velocityFieldRegistration->SetLearningRate( learningRate );
      velocityFieldRegistration->GetTransform()->SetGaussianSpatialSmoothingVarianceForTheTotalField( sigmaForTotalField );
      velocityFieldRegistration->GetTransform()->SetGaussianSpatialSmoothingVarianceForTheUpdateField(
        sigmaForUpdateField );
      velocityFieldRegistration->GetTransform()->SetGaussianTemporalSmoothingVarianceForTheTotalField(
        sigmaForTotalFieldTime );
      velocityFieldRegistration->GetTransform()->SetGaussianTemporalSmoothingVarianceForTheUpdateField(
        sigmaForUpdateFieldTime );

      velocityFieldRegistration->GetTransform()->SetTimeVaryingVelocityField( velocityField );
      velocityFieldRegistration->GetTransform()->SetLowerTimeBound( 0.0 );
      velocityFieldRegistration->GetTransform()->SetUpperTimeBound( 1.0 );

      typename VelocityFieldRegistrationType::NumberOfIterationsArrayType numberOfIterationsPerLevel;
      numberOfIterationsPerLevel.SetSize( numberOfLevels );
      for( unsigned int d = 0; d < numberOfLevels; d++ )
        {
        numberOfIterationsPerLevel[d] = iterations[d];
        }
      velocityFieldRegistration->SetNumberOfIterationsPerLevel( numberOfIterationsPerLevel );
      velocityFieldRegistration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
      velocityFieldRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );

      typedef itk::TimeVaryingVelocityFieldTransformParametersAdaptor<typename VelocityFieldRegistrationType::
                                                                      TransformType> VelocityFieldTransformAdaptorType;

      typename VelocityFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;
      for( unsigned int level = 0; level < shrinkFactorsPerLevel.Size(); level++ )
        {
        typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkFilterType;
        typename ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
        shrinkFilter->SetShrinkFactors( shrinkFactorsPerLevel[level] );
        shrinkFilter->SetInput( fixedImage );
        shrinkFilter->Update();

        // Although we shrink the images for the given levels,
        // we keep the size in time the same

        velocityFieldSize.Fill( numberOfTimeIndices );
        velocityFieldOrigin.Fill( 0.0 );
        velocityFieldSpacing.Fill( 1.0 );
        velocityFieldDirection.SetIdentity();

        fixedImageSize = shrinkFilter->GetOutput()->GetBufferedRegion().GetSize();
        fixedImageOrigin = shrinkFilter->GetOutput()->GetOrigin();
        fixedImageSpacing = shrinkFilter->GetOutput()->GetSpacing();
        fixedImageDirection = shrinkFilter->GetOutput()->GetDirection();
        for( unsigned int i = 0; i < ImageDimension; i++ )
          {
          velocityFieldSize[i] = fixedImageSize[i];
          velocityFieldOrigin[i] = fixedImageOrigin[i];
          velocityFieldSpacing[i] = fixedImageSpacing[i];
          for( unsigned int j = 0; j < ImageDimension; j++ )
            {
            velocityFieldDirection[i][j] = fixedImageDirection[i][j];
            }
          }

        typename VelocityFieldTransformAdaptorType::Pointer fieldTransformAdaptor =
          VelocityFieldTransformAdaptorType::New();
        fieldTransformAdaptor->SetRequiredSpacing( velocityFieldSpacing );
        fieldTransformAdaptor->SetRequiredSize( velocityFieldSize );
        fieldTransformAdaptor->SetRequiredDirection( velocityFieldDirection );
        fieldTransformAdaptor->SetRequiredOrigin( velocityFieldOrigin );

        adaptors.push_back( fieldTransformAdaptor.GetPointer() );
        }

      velocityFieldRegistration->SetTransformParametersAdaptorsPerLevel( adaptors );

      typedef CommandIterationUpdate<VelocityFieldRegistrationType> VelocityFieldCommandType;
      typename VelocityFieldCommandType::Pointer velocityFieldRegistrationObserver = VelocityFieldCommandType::New();
      velocityFieldRegistrationObserver->SetNumberOfIterations( iterations );

      velocityFieldRegistration->AddObserver( itk::IterationEvent(), velocityFieldRegistrationObserver );

      try
        {
        std::cout << std::endl << "*** Running time-varying velocity field registration (sigmaForUpdateField = "
                  << sigmaForUpdateField << ", sigmaForTotalField = " << sigmaForTotalField
                  << ", sigmaForUpdateFieldTime = "
                  << sigmaForUpdateFieldTime << ", sigmaForTotalFieldTime = " << sigmaForTotalFieldTime << ") ***"
                  << std::endl << std::endl;
        velocityFieldRegistrationObserver->Execute( velocityFieldRegistration, itk::StartEvent() );
        velocityFieldRegistration->StartRegistration();
        }
      catch( itk::ExceptionObject & e )
        {
        std::cerr << "Exception caught: " << e << std::endl;
        return EXIT_FAILURE;
        }

      // Write out the displacement fields

      std::string filename = outputPrefix + currentStageString.str() + std::string( "Warp.nii.gz" );

      typedef typename VelocityFieldRegistrationType::TransformType::DisplacementFieldType DisplacementFieldType;

      typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( const_cast<typename VelocityFieldRegistrationType::TransformType *>( velocityFieldRegistration->
                                                                                             GetOutput()->Get() )->GetDisplacementField() );
      writer->SetFileName( filename.c_str() );
      writer->Update();

      std::string inverseFilename = outputPrefix + currentStageString.str() + std::string( "InverseWarp.nii.gz" );

      typedef itk::ImageFileWriter<DisplacementFieldType> InverseWriterType;
      typename InverseWriterType::Pointer inverseWriter = InverseWriterType::New();
      inverseWriter->SetInput( const_cast<typename VelocityFieldRegistrationType::TransformType *>(
                                 velocityFieldRegistration->GetOutput()->Get() )->GetInverseDisplacementField() );
      inverseWriter->SetFileName( inverseFilename.c_str() );
      inverseWriter->Update();
      }
    else if( std::strcmp( whichTransform.c_str(),
                          "timevaryingbsplinevelocityfield" ) == 0 ||
             std::strcmp( whichTransform.c_str(), "tvdmffd" ) == 0 )
      {
      typedef itk::Vector<RealType, ImageDimension> VectorType;
      VectorType zeroVector( 0.0 );

      // Determine the parameters (size, spacing, etc) for the time-varying velocity field control point lattice

      std::vector<unsigned int> meshSize = parser->ConvertVector<unsigned int>( transformOption->GetParameter( 0, 1 ) );
      if( meshSize.size() != ImageDimension + 1 )
        {
        std::cerr << "The transform domain mesh size does not have the correct number of elements."
                  << "For image dimension = " << ImageDimension << ", you need " << ImageDimension + 1
                  << "elements. " << std::endl;
        return EXIT_FAILURE;
        }

      unsigned int numberOfTimePointSamples = 4;
      if( transformOption->GetNumberOfParameters( currentStage ) > 2 )
        {
        numberOfTimePointSamples = parser->Convert<unsigned int>( transformOption->GetParameter( currentStage, 2 ) );
        }
      unsigned int splineOrder = 3;
      if( transformOption->GetNumberOfParameters( currentStage ) > 3 )
        {
        splineOrder = parser->Convert<unsigned int>( transformOption->GetParameter( currentStage, 3 ) );
        }

      typedef itk::Image<VectorType, ImageDimension + 1> TimeVaryingVelocityFieldControlPointLatticeType;
      typename TimeVaryingVelocityFieldControlPointLatticeType::Pointer velocityFieldLattice =
        TimeVaryingVelocityFieldControlPointLatticeType::New();

      typename ImageType::IndexType fixedImageIndex = fixedImage->GetBufferedRegion().GetIndex();
      typename ImageType::SizeType fixedImageSize = fixedImage->GetBufferedRegion().GetSize();
      typename ImageType::PointType fixedImageOrigin = fixedImage->GetOrigin();
      typename ImageType::SpacingType fixedImageSpacing = fixedImage->GetSpacing();
      typename ImageType::DirectionType fixedImageDirection = fixedImage->GetDirection();

      typename TimeVaryingVelocityFieldControlPointLatticeType::SizeType transformDomainMeshSize;
      typename TimeVaryingVelocityFieldControlPointLatticeType::PointType transformDomainOrigin;
      typename TimeVaryingVelocityFieldControlPointLatticeType::SpacingType transformDomainPhysicalDimensions;
      typename TimeVaryingVelocityFieldControlPointLatticeType::DirectionType transformDomainDirection;

      transformDomainDirection.SetIdentity();
      transformDomainOrigin.Fill( 0.0 );
      transformDomainPhysicalDimensions.Fill( 1.0 );
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        transformDomainOrigin[i] = fixedImageOrigin[i];
        transformDomainMeshSize[i] = 3;
        transformDomainPhysicalDimensions[i] = static_cast<double>( fixedImageSize[i] - 1 ) * fixedImageSpacing[i];
        for( unsigned int j = 0; j < ImageDimension; j++ )
          {
          transformDomainDirection[i][j] = fixedImageDirection[i][j];
          }
        }
      for( unsigned int i = 0; i < meshSize.size(); i++ )
        {
        transformDomainMeshSize[i] = meshSize[i];
        }
      typename TimeVaryingVelocityFieldControlPointLatticeType::SizeType initialTransformDomainMeshSize =
        transformDomainMeshSize;

      typedef itk::TimeVaryingBSplineVelocityFieldImageRegistrationMethod<ImageType,
                                                                          ImageType> VelocityFieldRegistrationType;
      typename VelocityFieldRegistrationType::Pointer velocityFieldRegistration = VelocityFieldRegistrationType::New();
      velocityFieldRegistration->SetFixedImage( fixedImage );
      velocityFieldRegistration->SetMovingImage( movingImage );
      velocityFieldRegistration->SetNumberOfLevels( numberOfLevels );
      velocityFieldRegistration->SetNumberOfTimePointSamples( numberOfTimePointSamples );
      velocityFieldRegistration->SetCompositeTransform( compositeTransform );
      velocityFieldRegistration->SetMetric( metric );
      velocityFieldRegistration->SetMetricSamplingStrategy(
        static_cast<typename VelocityFieldRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
      velocityFieldRegistration->SetMetricSamplingPercentage( samplingPercentage );
      velocityFieldRegistration->SetLearningRate( learningRate );
      velocityFieldRegistration->GetTransform()->SetSplineOrder( splineOrder );
      velocityFieldRegistration->GetTransform()->SetLowerTimeBound( 0.0 );
      velocityFieldRegistration->GetTransform()->SetUpperTimeBound( 1.0 );

      typedef typename VelocityFieldRegistrationType::TransformType TransformType;

      typedef itk::TimeVaryingBSplineVelocityFieldTransformParametersAdaptor<TransformType>
        VelocityFieldTransformAdaptorType;
      typename VelocityFieldTransformAdaptorType::Pointer initialFieldTransformAdaptor =
        VelocityFieldTransformAdaptorType::New();
      initialFieldTransformAdaptor->SetTransform( velocityFieldRegistration->GetTransform() );
      initialFieldTransformAdaptor->SetRequiredTransformDomainOrigin( transformDomainOrigin );
      initialFieldTransformAdaptor->SetRequiredTransformDomainPhysicalDimensions( transformDomainPhysicalDimensions );
      initialFieldTransformAdaptor->SetRequiredTransformDomainMeshSize( transformDomainMeshSize );
      initialFieldTransformAdaptor->SetRequiredTransformDomainDirection( transformDomainDirection );

      velocityFieldLattice->SetOrigin( initialFieldTransformAdaptor->GetRequiredControlPointLatticeOrigin() );
      velocityFieldLattice->SetSpacing( initialFieldTransformAdaptor->GetRequiredControlPointLatticeSpacing() );
      velocityFieldLattice->SetDirection( initialFieldTransformAdaptor->GetRequiredControlPointLatticeDirection() );
      velocityFieldLattice->SetRegions( initialFieldTransformAdaptor->GetRequiredControlPointLatticeSize() );
      velocityFieldLattice->Allocate();
      velocityFieldLattice->FillBuffer( zeroVector );

      typename TransformType::VelocityFieldPointType        sampledVelocityFieldOrigin;
      typename TransformType::VelocityFieldSpacingType      sampledVelocityFieldSpacing;
      typename TransformType::VelocityFieldSizeType         sampledVelocityFieldSize;
      typename TransformType::VelocityFieldDirectionType    sampledVelocityFieldDirection;

      sampledVelocityFieldOrigin.Fill( 0.0 );
      sampledVelocityFieldSpacing.Fill( 1.0 );
      sampledVelocityFieldSize.Fill( numberOfTimePointSamples );
      sampledVelocityFieldDirection.SetIdentity();
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        sampledVelocityFieldOrigin[i] = fixedImage->GetOrigin()[i];
        sampledVelocityFieldSpacing[i] = fixedImage->GetSpacing()[i];
        sampledVelocityFieldSize[i] = fixedImage->GetRequestedRegion().GetSize()[i];
        for( unsigned int j = 0; j < ImageDimension; j++ )
          {
          sampledVelocityFieldDirection[i][j] = fixedImage->GetDirection()[i][j];
          }
        }

      velocityFieldRegistration->GetTransform()->SetTimeVaryingVelocityFieldControlPointLattice( velocityFieldLattice );
      velocityFieldRegistration->GetTransform()->SetVelocityFieldOrigin( sampledVelocityFieldOrigin );
      velocityFieldRegistration->GetTransform()->SetVelocityFieldDirection( sampledVelocityFieldDirection );
      velocityFieldRegistration->GetTransform()->SetVelocityFieldSpacing( sampledVelocityFieldSpacing );
      velocityFieldRegistration->GetTransform()->SetVelocityFieldSize( sampledVelocityFieldSize );
//      velocityFieldRegistration->GetTransform()->IntegrateVelocityField();

      typename VelocityFieldRegistrationType::NumberOfIterationsArrayType numberOfIterationsPerLevel;
      numberOfIterationsPerLevel.SetSize( numberOfLevels );
      for( unsigned int d = 0; d < numberOfLevels; d++ )
        {
        numberOfIterationsPerLevel[d] = iterations[d];
        }
      velocityFieldRegistration->SetNumberOfIterationsPerLevel( numberOfIterationsPerLevel );
      velocityFieldRegistration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
      velocityFieldRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );

      typename VelocityFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;
      for( unsigned int level = 0; level < shrinkFactorsPerLevel.Size(); level++ )
        {
        typename VelocityFieldTransformAdaptorType::Pointer fieldTransformAdaptor =
          VelocityFieldTransformAdaptorType::New();
        fieldTransformAdaptor->SetTransform( velocityFieldRegistration->GetTransform() );
        fieldTransformAdaptor->SetRequiredTransformDomainOrigin( transformDomainOrigin );
        fieldTransformAdaptor->SetRequiredTransformDomainMeshSize( transformDomainMeshSize );
        fieldTransformAdaptor->SetRequiredTransformDomainDirection( transformDomainDirection );
        fieldTransformAdaptor->SetRequiredTransformDomainPhysicalDimensions( transformDomainPhysicalDimensions );

        adaptors.push_back( fieldTransformAdaptor.GetPointer() );
        for( unsigned int i = 0; i <= ImageDimension; i++ )
          {
          transformDomainMeshSize[i] <<= 1;
          }
        }
      velocityFieldRegistration->SetTransformParametersAdaptorsPerLevel( adaptors );

      typedef CommandIterationUpdate<VelocityFieldRegistrationType> VelocityFieldCommandType;
      typename VelocityFieldCommandType::Pointer velocityFieldRegistrationObserver = VelocityFieldCommandType::New();
      velocityFieldRegistrationObserver->SetNumberOfIterations( iterations );

      velocityFieldRegistration->AddObserver( itk::IterationEvent(), velocityFieldRegistrationObserver );

      try
        {
        std::cout << std::endl
                  << "*** Running time-varying b-spline velocity field registration (initial mesh size = "
                  << initialTransformDomainMeshSize << ") ***" << std::endl << std::endl;
        velocityFieldRegistrationObserver->Execute( velocityFieldRegistration, itk::StartEvent() );
        velocityFieldRegistration->StartRegistration();
        }
      catch( itk::ExceptionObject & e )
        {
        std::cerr << "Exception caught: " << e << std::endl;
        return EXIT_FAILURE;
        }

      // Write out the displacement fields

      std::string filename = outputPrefix + currentStageString.str() + std::string( "Warp.nii.gz" );

      typedef typename VelocityFieldRegistrationType::TransformType::DisplacementFieldType DisplacementFieldType;

      typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( const_cast<typename VelocityFieldRegistrationType::TransformType *>( velocityFieldRegistration->
                                                                                             GetOutput()->Get() )->GetDisplacementField() );
      writer->SetFileName( filename.c_str() );
      writer->Update();

      std::string inverseFilename = outputPrefix + currentStageString.str() + std::string( "InverseWarp.nii.gz" );

      typedef itk::ImageFileWriter<DisplacementFieldType> InverseWriterType;
      typename InverseWriterType::Pointer inverseWriter = InverseWriterType::New();
      inverseWriter->SetInput( const_cast<typename VelocityFieldRegistrationType::TransformType *>(
                                 velocityFieldRegistration->GetOutput()->Get() )->GetInverseDisplacementField() );
      inverseWriter->SetFileName( inverseFilename.c_str() );
      inverseWriter->Update();
      }
    else if( std::strcmp( whichTransform.c_str(),
                          "syn" ) == 0 ||  std::strcmp( whichTransform.c_str(), "symmetricnormalization" ) == 0 )
      {
//      metric->SetDoFixedImagePreWarp( true );
//      metric->SetDoMovingImagePreWarp( true );
      typedef itk::Vector<RealType, ImageDimension> VectorType;
      VectorType zeroVector( 0.0 );
      typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;
      typename DisplacementFieldType::Pointer displacementField = DisplacementFieldType::New();
      displacementField->CopyInformation( fixedImage );
      displacementField->SetRegions( fixedImage->GetBufferedRegion() );
      displacementField->Allocate();
      displacementField->FillBuffer( zeroVector );

      typename DisplacementFieldType::Pointer inverseDisplacementField = DisplacementFieldType::New();
      inverseDisplacementField->CopyInformation( fixedImage );
      inverseDisplacementField->SetRegions( fixedImage->GetBufferedRegion() );
      inverseDisplacementField->Allocate();
      inverseDisplacementField->FillBuffer( zeroVector );

      typedef itk::DisplacementFieldTransform<RealType, ImageDimension> DisplacementFieldTransformType;
      typename DisplacementFieldTransformType::Pointer displacementFieldTransform =
        DisplacementFieldTransformType::New();
      displacementFieldTransform->SetDisplacementField( displacementField );
      displacementFieldTransform->SetInverseDisplacementField( inverseDisplacementField );

      typedef itk::SyNImageRegistrationMethod<ImageType, ImageType,
                                              DisplacementFieldTransformType> DisplacementFieldRegistrationType;

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
        fieldTransformAdaptor->SetTransform( displacementFieldTransform );

        adaptors.push_back( fieldTransformAdaptor.GetPointer() );
        }

      // Extract parameters
      typename DisplacementFieldRegistrationType::NumberOfIterationsArrayType numberOfIterationsPerLevel;
      numberOfIterationsPerLevel.SetSize( numberOfLevels );
      for( unsigned int d = 0; d < numberOfLevels; d++ )
        {
        numberOfIterationsPerLevel[d] = iterations[d];
        }

      RealType sigmaForUpdateField = parser->Convert<float>( transformOption->GetParameter( currentStage, 1 ) );
      RealType sigmaForTotalField = parser->Convert<float>( transformOption->GetParameter( currentStage, 2 ) );

      typename DisplacementFieldRegistrationType::Pointer displacementFieldRegistration =
        DisplacementFieldRegistrationType::New();
      displacementFieldRegistration->SetDownsampleImagesForMetricDerivatives( true );
      displacementFieldRegistration->SetAverageMidPointGradients( false );
      displacementFieldRegistration->SetFixedImage( fixedImage );
      displacementFieldRegistration->SetMovingImage( movingImage );
      displacementFieldRegistration->SetNumberOfLevels( numberOfLevels );
      displacementFieldRegistration->SetCompositeTransform( compositeTransform );
      displacementFieldRegistration->SetTransform( displacementFieldTransform );
      displacementFieldRegistration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
      displacementFieldRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
      displacementFieldRegistration->SetMetric( metric );
      displacementFieldRegistration->SetLearningRate( learningRate );
      displacementFieldRegistration->SetNumberOfIterationsPerLevel( numberOfIterationsPerLevel );
      displacementFieldRegistration->SetTransformParametersAdaptorsPerLevel( adaptors );
      displacementFieldRegistration->SetGaussianSmoothingVarianceForTheUpdateField( sigmaForUpdateField );
      displacementFieldRegistration->SetGaussianSmoothingVarianceForTheTotalField( sigmaForTotalField );

      typedef CommandIterationUpdate<DisplacementFieldRegistrationType> DisplacementFieldCommandType;
      typename DisplacementFieldCommandType::Pointer displacementFieldRegistrationObserver =
        DisplacementFieldCommandType::New();
      displacementFieldRegistrationObserver->SetNumberOfIterations( iterations );

      displacementFieldRegistration->AddObserver( itk::IterationEvent(), displacementFieldRegistrationObserver );

      try
        {
        std::cout << std::endl << "*** Running SyN registration (sigmaForUpdateField = "
                  << sigmaForUpdateField << ", sigmaForTotalField = " << sigmaForTotalField << ") ***" << std::endl
                  << std::endl;
        displacementFieldRegistrationObserver->Execute( displacementFieldRegistration, itk::StartEvent() );
        displacementFieldRegistration->StartRegistration();
        }
      catch( itk::ExceptionObject & e )
        {
        std::cerr << "Exception caught: " << e << std::endl;
        return EXIT_FAILURE;
        }

      // Write out the displacement field and its inverse

      std::string filename = outputPrefix + currentStageString.str() + std::string( "Warp.nii.gz" );

      typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( const_cast<typename DisplacementFieldRegistrationType::TransformType *>(
                          displacementFieldRegistration->GetOutput()->Get() )->GetDisplacementField() );
      writer->SetFileName( filename.c_str() );
      writer->Update();

      filename = outputPrefix + currentStageString.str() + std::string( "InverseWarp.nii.gz" );

      typename WriterType::Pointer inverseWriter = WriterType::New();
      inverseWriter->SetInput( const_cast<typename DisplacementFieldRegistrationType::TransformType *>(
                                 displacementFieldRegistration->GetOutput()->Get() )->GetInverseDisplacementField() );
      inverseWriter->SetFileName( filename.c_str() );
      inverseWriter->Update();
      }
    else
      {
      std::cerr << "ERROR:  Unrecognized transform option - " << whichTransform << std::endl;
      return EXIT_FAILURE;
      }
    timer.Stop();
    std::cout << "  Elapsed time (stage "
              << ( numberOfStages - currentStage - 1 ) << "): " << timer.GetMeanTime() << std::endl << std::endl;
    }

  // Write out warped image(s), if requested.

  if( outputOption && outputOption->GetNumberOfParameters( 0 ) > 1 )
    {
    std::string fixedImageFileName = metricOption->GetParameter( 0, 0 );
    std::string movingImageFileName = metricOption->GetParameter( 0, 1 );

    std::cout << "Warping " << movingImageFileName << " to " << fixedImageFileName << std::endl;

    typedef itk::ImageFileReader<ImageType> ImageReaderType;
    typename ImageReaderType::Pointer fixedImageReader = ImageReaderType::New();
    fixedImageReader->SetFileName( fixedImageFileName.c_str() );
    fixedImageReader->Update();
    typename ImageType::Pointer fixedImage = fixedImageReader->GetOutput();
    fixedImage->Update();
    fixedImage->DisconnectPipeline();

    typename ImageReaderType::Pointer movingImageReader = ImageReaderType::New();
    movingImageReader->SetFileName( movingImageFileName.c_str() );
    movingImageReader->Update();
    typename ImageType::Pointer movingImage = movingImageReader->GetOutput();
    movingImage->Update();
    movingImage->DisconnectPipeline();

    typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
    typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetTransform( compositeTransform );
    resampler->SetInput( movingImage );
    resampler->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
    resampler->SetOutputOrigin(  fixedImage->GetOrigin() );
    resampler->SetOutputSpacing( fixedImage->GetSpacing() );
    resampler->SetOutputDirection( fixedImage->GetDirection() );
    resampler->SetDefaultPixelValue( 0 );
    resampler->Update();

    std::string fileName = outputOption->GetParameter( 0, 1 );

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( fileName.c_str() );
    writer->SetInput( resampler->GetOutput() );
    writer->Update();

    if( outputOption->GetNumberOfParameters( 0 ) > 2 && compositeTransform->GetInverseTransform() )
      {
      std::cout << "Warping " << fixedImageFileName << " to " << movingImageFileName << std::endl;

      typedef itk::ResampleImageFilter<ImageType, ImageType> InverseResampleFilterType;
      typename InverseResampleFilterType::Pointer inverseResampler = ResampleFilterType::New();
      inverseResampler->SetTransform( compositeTransform->GetInverseTransform() );
      inverseResampler->SetInput( fixedImage );
      inverseResampler->SetSize( movingImage->GetBufferedRegion().GetSize() );
      inverseResampler->SetOutputOrigin( movingImage->GetOrigin() );
      inverseResampler->SetOutputSpacing( movingImage->GetSpacing() );
      inverseResampler->SetOutputDirection( movingImage->GetDirection() );
      inverseResampler->SetDefaultPixelValue( 0 );
      inverseResampler->Update();

      std::string inverseFileName = outputOption->GetParameter( 0, 2 );

      typedef itk::ImageFileWriter<ImageType> InverseWriterType;
      typename InverseWriterType::Pointer inverseWriter = InverseWriterType::New();
      inverseWriter->SetFileName( inverseFileName.c_str() );
      inverseWriter->SetInput( inverseResampler->GetOutput() );
      inverseWriter->Update();
      }
    }

  totalTimer.Stop();
  std::cout << std::endl << "Total elapsed time: " << totalTimer.GetMeanTime() << std::endl;

  return EXIT_SUCCESS;
}

void InitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

    {
    std::string description =
      std::string( "This option forces the image to be treated as a specified-" )
      + std::string( "dimensional image.  If not specified, N4 tries to " )
      + std::string( "infer the dimensionality from the input image." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "dimensionality" );
    option->SetShortName( 'd' );
    option->SetUsageOption( 0, "2/3" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Specify the output transform prefix (output format is .nii.gz ). " )
      + std::string( "Optionally, one can choose to warp the moving image to the fixed space and, if the " )
      + std::string( "inverse transform exists, one can also output the warped fixed image." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "output" );
    option->SetShortName( 'o' );
    option->SetUsageOption( 0, "outputTransformPrefix" );
    option->SetUsageOption( 1, "[outputTransformPrefix,<outputWarpedImage>,<outputInverseWarpedImage>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Specify the initial transform(s) which get immediately " )
      + std::string( "incorporated into the composite transform.  The order of the " )
      + std::string( "transforms is stack-esque in that the last transform specified on " )
      + std::string( "the command line is the first to be applied.  See antsApplyTransforms " )
      + std::string( "for additional information." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "initial-transform" );
    option->SetShortName( 'r' );
    option->SetUsageOption( 0, "initialTransform" );
    option->SetUsageOption( 1, "[initialTransform,<useInverse>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "These image metrics are available--- " )
      + std::string( "CC:  ANTS neighborhood cross correlation, MI:  Mutual information, and " )
      + std::string( "Demons:  Thirion's Demons (modified mean-squares). " )
      + std::string( "GC, Global Correlation. " )
      + std::string( "Note that the metricWeight is currently not used.  " )
      + std::string( "Rather, it is a temporary place holder until multivariate metrics " )
      + std::string( "are available for a single stage. " )
      + std::string( "The metrics can also employ a sampling strategy defined by a " )
      + std::string( "sampling percentage. The sampling strategy defaults to dense, otherwise " )
      + std::string( "it defines a point set over which to optimize the metric. " )
      + std::string( "The point set can be on a regular lattice or a random lattice of points slightly " )
      + std::string( "perturbed to minimize aliasing artifacts. samplingPercentage defines the " )
      + std::string( "fraction of points to select from the domain. " );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "metric" );
    option->SetShortName( 'm' );
    option->SetUsageOption( 0,
                            "CC[fixedImage,movingImage,metricWeight,radius,<samplingStrategy={Regular,Random}>,<samplingPercentage=[0,1]>]" );
    option->SetUsageOption( 1,
                            "MI[fixedImage,movingImage,metricWeight,numberOfBins,<samplingStrategy={Regular,Random}>,<samplingPercentage=[0,1]>]" );
    option->SetUsageOption( 2,
                            "Mattes[fixedImage,movingImage,metricWeight,numberOfBins,<samplingStrategy={Regular,Random}>,<samplingPercentage=[0,1]>]" );
    option->SetUsageOption( 3,
                            "Demons[fixedImage,movingImage,metricWeight,radius,<samplingStrategy={Regular,Random}>,<samplingPercentage=[0,1]>]" );
    option->SetUsageOption( 4,
                            "GC[fixedImage,movingImage,metricWeight,radius,<samplingStrategy={Regular,Random}>,<samplingPercentage=[0,1]>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Several transform options are available.  The gradientStep or " )
      + std::string( "learningRate characterizes the gradient descent optimization and is scaled appropriately " )
      + std::string( "for each transform using the shift scales estimator.  Subsequent parameters are " )
      + std::string( "transform-specific and can be determined from the usage. " );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "transform" );
    option->SetShortName( 't' );
    option->SetUsageOption( 0, "Rigid[gradientStep]" );
    option->SetUsageOption( 1, "Affine[gradientStep]" );
    option->SetUsageOption( 2, "CompositeAffine[gradientStep]" );
    option->SetUsageOption( 3, "Similarity[gradientStep]" );
    option->SetUsageOption( 4, "Translation[gradientStep]" );
    option->SetUsageOption( 5, "BSpline[gradientStep,meshSizeAtBaseLevel]" );
    option->SetUsageOption( 6,
                            "GaussianDisplacementField[gradientStep,updateFieldSigmaInPhysicalSpace,totalFieldSigmaInPhysicalSpace]" );
    option->SetUsageOption( 7,
                            "BSplineDisplacementField[gradientStep,updateFieldMeshSizeAtBaseLevel,totalFieldMeshSizeAtBaseLevel,<splineOrder=3>]" );
    option->SetUsageOption( 8,
                            "TimeVaryingVelocityField[gradientStep,numberOfTimeIndices,updateFieldSigmaInPhysicalSpace,updateFieldTimeSigma,totalFieldSigmaInPhysicalSpace,totalFieldTimeSigma]" );
    option->SetUsageOption( 9,
                            "TimeVaryingBSplineVelocityField[gradientStep,velocityFieldMeshSize,<numberOfTimePointSamples=4>,<splineOrder=3>]" );
    option->SetUsageOption( 10, "SyN[gradientStep,updateFieldSigmaInPhysicalSpace,totalFieldSigmaInPhysicalSpace]" );
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
    option->SetLongName( "smoothing-sigmas" );
    option->SetShortName( 's' );
    option->SetUsageOption( 0, "MxNx0..." );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string(
        "Specify the shrink factor for the virtual domain (typically the fixed image) at each level." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "shrink-factors" );
    option->SetShortName( 'f' );
    option->SetUsageOption( 0, "MxNx0..." );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Histogram match the images before registration." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "use-histogram-matching" );
    option->SetShortName( 'u' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Rescale the images to [0,1] before registration." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "rescale-images" );
    option->SetShortName( 'c' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Winsorize data based on specified quantiles." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "winsorize-image-intensities" );
    option->SetShortName( 'w' );
    option->SetUsageOption( 0, "[lowerQuantile,upperQuantile]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Print the help menu (short version)." );

    OptionType::Pointer option = OptionType::New();
    option->SetShortName( 'h' );
    option->SetDescription( description );
    option->AddValue( std::string( "0" ) );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Print the help menu." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "help" );
    option->SetDescription( description );
    option->AddValue( std::string( "0" ) );
    parser->AddOption( option );
    }
}

int main( int argc, char *argv[] )
{
  itk::ants::CommandLineParser::Pointer parser = itk::ants::CommandLineParser::New();

  parser->SetCommand( argv[0] );

  std::string commandDescription = std::string( "This program is a user-level " )
    + std::string( "registration application meant to utilize ITKv4-only classes. The user can specify " )
    + std::string( "any number of \"stages\" where a stage consists of a transform; an image metric; " )
    + std::string( "and iterations, shrink factors, and smoothing sigmas for each level." );

  parser->SetCommandDescription( commandDescription );
  InitializeCommandLineOptions( parser );

  parser->Parse( argc, argv );

  if( argc < 2 || parser->Convert<bool>( parser->GetOption( "help" )->GetValue() ) )
    {
    parser->PrintMenu( std::cout, 5, false );
    exit( EXIT_FAILURE );
    }
  else if( parser->Convert<bool>( parser->GetOption( 'h' )->GetValue() ) )
    {
    parser->PrintMenu( std::cout, 5, true );
    exit( EXIT_FAILURE );
    }

  // Get dimensionality
  unsigned int dimension = 3;

  itk::ants::CommandLineParser::OptionType::Pointer dimOption = parser->GetOption( "dimensionality" );
  if( dimOption && dimOption->GetNumberOfValues() > 0 )
    {
    dimension = parser->Convert<unsigned int>( dimOption->GetValue() );
    }
  else
    {
    std::cerr << "Image dimensionality not specified.  See command line option --dimensionality" << std::endl;
    exit( EXIT_FAILURE );
    }

  std::cout << std::endl << "Running antsRegistration for " << dimension << "-dimensional images." << std::endl
            << std::endl;

  switch( dimension )
    {
    case 2:
      {
      antsRegistration<2>( parser );
      }
      break;
    case 3:
      {
      antsRegistration<3>( parser );
      }
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
}
