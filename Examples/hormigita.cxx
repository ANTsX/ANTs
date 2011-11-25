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

#include "itkSimpleImageRegistrationMethod.h"
#include "itkTimeVaryingVelocityFieldImageRegistrationMethod.h"

#include "itkANTSNeighborhoodCorrelationImageToImageObjectMetric.h"
#include "itkDemonsImageToImageObjectMetric.h"
#include "itkImageToImageObjectMetric.h"
#include "itkJointHistogramMutualInformationImageToImageObjectMetric.h"

#include "itkAffineTransform.h"
#include "itkBSplineTransform.h"
#include "itkBSplineSmoothingOnUpdateDisplacementFieldTransform.h"
#include "itkCompositeTransform.h"
#include "itkGaussianSmoothingOnUpdateDisplacementFieldTransform.h"
#include "itkIdentityTransform.h"
#include "itkEuler2DTransform.h"
#include "itkEuler3DTransform.h"
#include "itkTransform.h"

#include "itkBSplineTransformParametersAdaptor.h"
#include "itkBSplineSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor.h"
#include "itkGaussianSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor.h"
#include "itkTimeVaryingVelocityFieldTransformParametersAdaptor.h"

#include "itkGradientDescentObjectOptimizer.h"

#include "itkHistogramMatchingImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMacro.h"
#include "itkRegistrationParameterScalesFromShift.h"
#include "itkResampleImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkTimeProbe.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
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

    if( typeid( event ) != typeid( itk::IterationEvent ) )
      {
      return;
      }

    unsigned int currentLevel = filter->GetCurrentLevel();
    typename TFilter::ShrinkFactorsArrayType shrinkFactors = filter->GetShrinkFactorsPerLevel();
    typename TFilter::SmoothingSigmasArrayType smoothingSigmas = filter->GetSmoothingSigmasPerLevel();
    typename TFilter::TransformParametersAdaptorsContainerType adaptors =
      filter->GetTransformParametersAdaptorsPerLevel();

    std::cout << "  Current level = " << currentLevel << std::endl;
    std::cout << "    number of iterations = " << this->m_NumberOfIterations[currentLevel] << std::endl;
    std::cout << "    shrink factor = " << shrinkFactors[currentLevel] << std::endl;
    std::cout << "    smoothing sigma = " << smoothingSigmas[currentLevel] << std::endl;
    std::cout << "    required fixed parameters = " << adaptors[currentLevel]->GetRequiredFixedParameters()
              << std::endl;

    typedef itk::GradientDescentObjectOptimizer GradientDescentObjectOptimizerType;

    GradientDescentObjectOptimizerType * optimizer = reinterpret_cast<GradientDescentObjectOptimizerType *>(
        const_cast<typename TFilter::OptimizerType *>( filter->GetOptimizer() ) );
    optimizer->SetNumberOfIterations( this->m_NumberOfIterations[currentLevel] );
  }

  void SetNumberOfIterations( std::vector<unsigned int> iterations )
  {
    this->m_NumberOfIterations = iterations;
  }

private:

  std::vector<unsigned int> m_NumberOfIterations;
};

void ConvertToLowerCase( std::string& str )
{
  std::transform( str.begin(), str.end(), str.begin(), tolower );
// You may need to cast the above line to (int(*)(int))
// tolower - this works as is on VC 7.1 but may not work on
// other compilers
}

template <unsigned int ImageDimension>
int hormigita( itk::ants::CommandLineParser *parser )
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

  typename OptionType::Pointer shrinkFactorsOption = parser->GetOption( "shrinkFactors" );
  if( !shrinkFactorsOption || shrinkFactorsOption->GetNumberOfValues() != numberOfStages  )
    {
    std::cerr << "The number of shrinkFactor sets specified does not match the number of stages." << std::endl;
    return EXIT_FAILURE;
    }

  typename OptionType::Pointer smoothingSigmasOption = parser->GetOption( "smoothingSigmas" );
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
  std::string outputPrefix = outputOption->GetParameter( 0, 0 );

  typedef float                                 PixelType;
  typedef double                                RealType;
  typedef itk::Image<PixelType, ImageDimension> FixedImageType;
  typedef itk::Image<PixelType, ImageDimension> MovingImageType;

  typedef itk::CompositeTransform<RealType, ImageDimension> CompositeTransformType;
  typename CompositeTransformType::Pointer compositeTransform = CompositeTransformType::New();

  typedef itk::IdentityTransform<RealType, ImageDimension> IdentityTransformType;
  typename IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();

  compositeTransform->AddTransform( identityTransform );
  // We iterate backwards because the command line options are stored as a stack (first in last out)
  for( int currentStage = numberOfStages - 1; currentStage >= 0; currentStage-- )
    {
    itk::TimeProbe timer;
    timer.Start();

    typedef itk::SimpleImageRegistrationMethod<FixedImageType, MovingImageType> AffineRegistrationType;

    std::cout << std::endl << "Stage " << numberOfStages - currentStage << std::endl;
    std::stringstream currentStageString;
    currentStageString << currentStage;

    // Get the fixed and moving images

    std::string fixedImageFileName = metricOption->GetParameter( currentStage, 0 );
    std::string movingImageFileName = metricOption->GetParameter( currentStage, 1 );

    std::cout << "  fixed image: " << fixedImageFileName << std::endl;
    std::cout << "  moving image: " << movingImageFileName << std::endl;

    typedef itk::ImageFileReader<FixedImageType> ImageReaderType;
    typename ImageReaderType::Pointer fixedImageReader = ImageReaderType::New();
    fixedImageReader->SetFileName( fixedImageFileName.c_str() );
    fixedImageReader->Update();
    typename FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();
    fixedImage->Update();
    fixedImage->DisconnectPipeline();

    typename ImageReaderType::Pointer movingImageReader = ImageReaderType::New();
    movingImageReader->SetFileName( movingImageFileName.c_str() );
    movingImageReader->Update();
    typename MovingImageType::Pointer movingImage = movingImageReader->GetOutput();
    movingImage->Update();
    movingImage->DisconnectPipeline();

    // Histogram match images if requested by the user

    typename OptionType::Pointer histOption = parser->GetOption( "useHistogramMatching" );
    if( histOption && histOption->GetNumberOfValues() > 0 )
      {
      std::string histValue = histOption->GetValue( 0 );
      ConvertToLowerCase( histValue );
      if( histValue.compare( "1" ) == 0 || histValue.compare( "true" ) == 0 )
        {
        std::cout << "  (histogram matching the images)" << std::endl;

        typedef itk::HistogramMatchingImageFilter<MovingImageType, MovingImageType> HistogramMatchingFilterType;
        typename HistogramMatchingFilterType::Pointer matchingFilter = HistogramMatchingFilterType::New();
        matchingFilter->SetSourceImage( movingImage );
        matchingFilter->SetReferenceImage( fixedImage );
        matchingFilter->SetNumberOfHistogramLevels( 256 );
        matchingFilter->SetNumberOfMatchPoints( 12 );
        matchingFilter->ThresholdAtMeanIntensityOn();

        movingImage = matchingFilter->GetOutput();
        movingImage->Update();
        movingImage->DisconnectPipeline();
        }
      }

    // Get the number of iterations and use that information to specify the number of levels

    std::vector<unsigned int> iterations = parser->ConvertVector<unsigned int>( iterationsOption->GetValue(
                                                                                  currentStage ) );
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
        smoothingSigmasPerLevel[n] = factors[n];
        }
      std::cout << "  smoothing sigmas per level: " << smoothingSigmasPerLevel << std::endl;
      }

    // Set up the image metric and scales estimator

    typedef itk::ImageToImageObjectMetric<FixedImageType, MovingImageType> MetricType;
    typename MetricType::Pointer metric;

    std::string whichMetric = metricOption->GetValue( currentStage );
    ConvertToLowerCase( whichMetric );
    if( std::strcmp( whichMetric.c_str(), "cc" ) == 0 )
      {
      unsigned int radiusOption = parser->Convert<unsigned int>( metricOption->GetParameter( currentStage, 3 ) );

      std::cout << "  using the CC metric (radius = " << radiusOption << ")." << std::endl;
      typedef itk::ANTSNeighborhoodCorrelationImageToImageObjectMetric<FixedImageType,
                                                                       MovingImageType> CorrelationMetricType;
      typename CorrelationMetricType::Pointer correlationMetric = CorrelationMetricType::New();
      typename CorrelationMetricType::RadiusType radius;
      radius.Fill( radiusOption );
      correlationMetric->SetRadius( radius );
      correlationMetric->SetDoFixedImagePreWarp( true );
      correlationMetric->SetDoMovingImagePreWarp( true );
      correlationMetric->SetUseMovingImageGradientFilter( false );
      correlationMetric->SetUseFixedImageGradientFilter( false );

      metric = correlationMetric;
      }
    else if( std::strcmp( whichMetric.c_str(), "mi" ) == 0 )
      {
      unsigned int binOption = parser->Convert<unsigned int>( metricOption->GetParameter( currentStage, 3 ) );

      std::cout << "  using the MI metric (number of bins = " << binOption << ")" << std::endl;
      typedef itk::JointHistogramMutualInformationImageToImageObjectMetric<FixedImageType,
                                                                           MovingImageType> MutualInformationMetricType;
      typename MutualInformationMetricType::Pointer mutualInformationMetric = MutualInformationMetricType::New();
      mutualInformationMetric = mutualInformationMetric;
      mutualInformationMetric->SetNumberOfHistogramBins( binOption );
      mutualInformationMetric->SetDoFixedImagePreWarp( true );
      mutualInformationMetric->SetDoMovingImagePreWarp( true );
      mutualInformationMetric->SetUseMovingImageGradientFilter( false );
      mutualInformationMetric->SetUseFixedImageGradientFilter( false );
      mutualInformationMetric->SetUseFixedSampledPointSet( false );

      metric = mutualInformationMetric;
      }
    else if( std::strcmp( whichMetric.c_str(), "demons" ) == 0 )
      {
      std::cout << "  using the Demons metric." << std::endl;

      typedef itk::DemonsImageToImageObjectMetric<FixedImageType, MovingImageType> DemonsMetricType;
      typename DemonsMetricType::Pointer demonsMetric = DemonsMetricType::New();
      demonsMetric = demonsMetric;
      demonsMetric->SetDoFixedImagePreWarp( true );
      demonsMetric->SetDoMovingImagePreWarp( true );

      metric = demonsMetric;
      }
    else
      {
      std::cerr << "ERROR: Unrecognized image metric: " << whichMetric << std::endl;
      }

    // Set up the optimizer.  To change the iteration number for each level we rely
    // on the command observer.

    typedef itk::RegistrationParameterScalesFromShift<MetricType> ScalesEstimatorType;
    typename ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();
    scalesEstimator->SetMetric( metric );
    scalesEstimator->SetTransformForward( true );

    float learningRate = parser->Convert<float>( transformOption->GetParameter( currentStage, 0 ) );

    typedef itk::GradientDescentObjectOptimizer GradientDescentObjectOptimizerType;
    typename GradientDescentObjectOptimizerType::Pointer optimizer = GradientDescentObjectOptimizerType::New();
    optimizer->SetLearningRate( learningRate );
    optimizer->SetNumberOfIterations( iterations[0] );
    optimizer->SetScalesEstimator( scalesEstimator );

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
        affineRegistration->StartRegistration();
        }
      catch( itk::ExceptionObject & e )
        {
        std::cerr << "Exception caught: " << e << std::endl;
        return EXIT_FAILURE;
        }

      // Write out the affine transform

      std::string filename = outputPrefix + currentStageString.str() + std::string( "Affine.txt" );

      typedef itk::TransformFileWriter TransformWriterType;
      typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
      transformWriter->SetInput( affineRegistration->GetOutput()->Get() );
      transformWriter->SetFileName( filename.c_str() );
      transformWriter->Update();
      }
//    else if( std::strcmp( whichTransform.c_str(), "rigid" ) == 0 )
//      {
//      typedef itk::Euler2DTransform<double> RigidTransformType;
//      typename RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
//
//      typedef itk::SimpleImageRegistrationMethod<FixedImageType, MovingImageType, RigidTransformType>
// RigidRegistrationType;
//      typename RigidRegistrationType::Pointer rigidRegistration = RigidRegistrationType::New();
//
//      rigidRegistration->SetFixedImage( fixedImage );
//      rigidRegistration->SetMovingImage( movingImage );
//      rigidRegistration->SetNumberOfLevels( numberOfLevels );
//      rigidRegistration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
//      rigidRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
//      rigidRegistration->SetMetric( metric );
//      rigidRegistration->SetOptimizer( optimizer );
//      rigidRegistration->SetTransform( rigidTransform );
//      rigidRegistration->SetCompositeTransform( compositeTransform );
//
//      typedef CommandIterationUpdate<RigidRegistrationType> RigidCommandType;
//      typename RigidCommandType::Pointer rigidObserver = RigidCommandType::New();
//      rigidObserver->SetNumberOfIterations( iterations );
//
//      rigidRegistration->AddObserver( itk::IterationEvent(), rigidObserver );
//
//      try
//        {
//        std::cout << std::endl << "*** Running rigid registration ***" << std::endl << std::endl;
//        rigidRegistration->StartRegistration();
//        }
//      catch( itk::ExceptionObject &e )
//        {
//        std::cerr << "Exception caught: " << e << std::endl;
//        return EXIT_FAILURE;
//        }
//      }
    else if( std::strcmp( whichTransform.c_str(),
                          "gaussiandisplacementfield" ) == 0 ||  std::strcmp( whichTransform.c_str(), "gdf" ) == 0 )
      {
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

      typedef itk::SimpleImageRegistrationMethod<FixedImageType, MovingImageType,
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
      displacementFieldRegistration->SetOptimizer( optimizer );
      displacementFieldRegistration->SetTransformParametersAdaptorsPerLevel( adaptors );

      typedef CommandIterationUpdate<DisplacementFieldRegistrationType> DisplacementFieldCommandType;
      typename DisplacementFieldCommandType::Pointer dfObserver = DisplacementFieldCommandType::New();
      dfObserver->SetNumberOfIterations( iterations );

      displacementFieldRegistration->AddObserver( itk::IterationEvent(), dfObserver );

      try
        {
        std::cout << std::endl << "*** Running gaussian displacement field registration (sigmaForUpdateField = "
                  << sigmaForUpdateField << ", sigmaForTotalField = " << sigmaForTotalField << ") ***" << std::endl
                  << std::endl;
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

      typedef itk::SimpleImageRegistrationMethod<FixedImageType, MovingImageType,
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
      displacementFieldRegistration->SetOptimizer( optimizer );
      displacementFieldRegistration->SetTransformParametersAdaptorsPerLevel( adaptors );

      typedef CommandIterationUpdate<DisplacementFieldRegistrationType> DisplacementFieldCommandType;
      typename DisplacementFieldCommandType::Pointer dfObserver = DisplacementFieldCommandType::New();
      dfObserver->SetNumberOfIterations( iterations );

      displacementFieldRegistration->AddObserver( itk::IterationEvent(), dfObserver );

      try
        {
        std::cout << std::endl << "*** Running bspline displacement field registration (updateMeshSizeAtBaseLevel = "
                  << updateMeshSize << ", totalMeshSizeAtBaseLevel = " << totalMeshSize << ") ***" << std::endl
                  << std::endl;
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

      typedef itk::SimpleImageRegistrationMethod<FixedImageType, MovingImageType,
                                                 BSplineTransformType> BSplineRegistrationType;

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
                          "TimeVaryingVelocityField" ) == 0 || std::strcmp( whichTransform.c_str(), "tvf" ) == 0 )
      {
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

      typename FixedImageType::IndexType fixedImageIndex = fixedImage->GetBufferedRegion().GetIndex();
      typename FixedImageType::SizeType fixedImageSize = fixedImage->GetBufferedRegion().GetSize();
      typename FixedImageType::PointType fixedImageOrigin = fixedImage->GetOrigin();
      typename FixedImageType::SpacingType fixedImageSpacing = fixedImage->GetSpacing();
      typename FixedImageType::DirectionType fixedImageDirection = fixedImage->GetDirection();

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

      typedef itk::TimeVaryingVelocityFieldImageRegistrationMethod<FixedImageType,
                                                                   MovingImageType> VelocityFieldRegistrationType;
      typename VelocityFieldRegistrationType::Pointer velocityFieldRegistration = VelocityFieldRegistrationType::New();
      velocityFieldRegistration->SetFixedImage( fixedImage );
      velocityFieldRegistration->SetMovingImage( movingImage );
      velocityFieldRegistration->SetNumberOfLevels( numberOfLevels );
      velocityFieldRegistration->SetCompositeTransform( compositeTransform );
      velocityFieldRegistration->SetMetric( metric );
      velocityFieldRegistration->SetLearningRate( learningRate );
      velocityFieldRegistration->SetNumberOfIntegrationStepsPerTimeIndex( 5 );
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

      typename VelocityFieldRegistrationType::ShrinkFactorsArrayType numberOfIterationsPerLevel;
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
        typedef itk::ShrinkImageFilter<FixedImageType, FixedImageType> ShrinkFilterType;
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

      try
        {
        std::cout << std::endl << "*** Running time-varying velocity field registration (sigmaForUpdateField = "
                  << sigmaForUpdateField << ", sigmaForTotalField = " << sigmaForTotalField
                  << ", sigmaForUpdateFieldTime = "
                  << sigmaForUpdateFieldTime << ", sigmaForTotalFieldTime = " << sigmaForTotalFieldTime << ") ***"
                  << std::endl << std::endl;
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
      inverseWriter->SetFileName( filename.c_str() );
      inverseWriter->Update();
      }
    else
      {
      std::cerr << "ERROR:  Unrecognized transform option - " << whichTransform << std::endl;
      return EXIT_FAILURE;
      }
    timer.Stop();
    std::cout << "  Elapsed time (stage " << currentStage + 1 << "): " << timer.GetMeanTime() << std::endl << std::endl;
    }

  // Write out warped image(s), if requested.

  if( outputOption && outputOption->GetNumberOfParameters( 0 ) > 1 )
    {
    std::string fixedImageFileName = metricOption->GetParameter( 0, 0 );
    std::string movingImageFileName = metricOption->GetParameter( 0, 1 );

    std::cout << "Warping " << movingImageFileName << " to " << fixedImageFileName << std::endl;

    typedef itk::ImageFileReader<FixedImageType> ImageReaderType;
    typename ImageReaderType::Pointer fixedImageReader = ImageReaderType::New();
    fixedImageReader->SetFileName( fixedImageFileName.c_str() );
    fixedImageReader->Update();
    typename FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();
    fixedImage->Update();
    fixedImage->DisconnectPipeline();

    typename ImageReaderType::Pointer movingImageReader = ImageReaderType::New();
    movingImageReader->SetFileName( movingImageFileName.c_str() );
    movingImageReader->Update();
    typename MovingImageType::Pointer movingImage = movingImageReader->GetOutput();
    movingImage->Update();
    movingImage->DisconnectPipeline();

    typedef itk::ResampleImageFilter<MovingImageType, FixedImageType> ResampleFilterType;
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

    typedef itk::ImageFileWriter<FixedImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( fileName.c_str() );
    writer->SetInput( resampler->GetOutput() );
    writer->Update();

    if( outputOption->GetNumberOfParameters( 0 ) > 2 && compositeTransform->GetInverseTransform() )
      {
      std::cout << "Warping " << fixedImageFileName << " to " << movingImageFileName << std::endl;

      typedef itk::ResampleImageFilter<FixedImageType, MovingImageType> InverseResampleFilterType;
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

      typedef itk::ImageFileWriter<MovingImageType> InverseWriterType;
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
    std::string description = std::string( "Three image metrics are available--- " )
      + std::string( "CC:  ANTS neighborhood cross correlation, MI:  Mutual information, and " )
      + std::string( "Demons:  Thirion's Demons (modified mean-squares). " )
      + std::string( "Note that the metricWeight is currently not used.  " )
      + std::string( "Rather, it is a temporary place holder until multivariate metrics " )
      + std::string( "are available for a single stage." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "metric" );
    option->SetShortName( 'm' );
    option->SetUsageOption( 0, "CC[fixedImage,movingImage,metricWeight,radius]" );
    option->SetUsageOption( 1, "MI[fixedImage,movingImage,metricWeight,numberOfBins]" );
    option->SetUsageOption( 2, "Demons[fixedImage,movingImage,metricWeight]" );
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
    option->SetUsageOption( 0, "Rigid[gradientStep]" );
    option->SetUsageOption( 1, "Affine[gradientStep]" );
    option->SetUsageOption( 2, "BSpline[gradientStep,meshSizeAtBaseLevel]" );
    option->SetUsageOption( 3,
                            "GaussianDisplacementField[gradientStep,updateFieldSigmaInPhysicalSpace,totalFieldSigmaInPhysicalSpace]" );
    option->SetUsageOption( 4,
                            "BSplineDisplacementField[gradientStep,updateFieldMeshSizeAtBaseLevel,totalFieldMeshSizeAtBaseLevel,<splineOrder=3>]" );
    option->SetUsageOption( 5,
                            "TimeVaryingVelocityField[gradientStep,numberOfTimeIndices,updateFieldSigmaInPhysicalSpace,updateFieldTimeSigma,totalFieldSigmaInPhysicalSpace,totalFieldTimeSigma]" );
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
    std::string description = std::string( "Histogram match the images before registration." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "useHistogramMatching" );
    option->SetShortName( 'u' );
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
    option->SetUsageOption( 0, "[outputTransformPrefix,<outputWarpedImage>,<outputInverseWarpedImage>]" );
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

  std::string commandDescription = std::string( "hormigita---little ant.  This program is a user-level " )
    + std::string( "registration application meant to utilize ITKv4-only classes. The user can specify " )
    + std::string( "any number of \"stages\" where a stage consists of a transform; an image metric; " )
    + std::string( " and iterations, shrink factors, and smoothing sigmas for each level." );

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

  std::cout << std::endl << "Running hormigita for " << dimension << "-dimensional images." << std::endl << std::endl;

  switch( dimension )
    {
    case 2:
      {
      hormigita<2>( parser );
      }
      break;
    case 3:
      {
      hormigita<3>( parser );
      }
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
}
