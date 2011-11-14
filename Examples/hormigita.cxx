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

#include "itkANTSNeighborhoodCorrelationImageToImageObjectMetric.h"
#include "itkDemonsImageToImageObjectMetric.h"
#include "itkImageToImageObjectMetric.h"
#include "itkJointHistogramMutualInformationImageToImageObjectMetric.h"

#include "itkAffineTransform.h"
#include "itkBSplineSmoothingOnUpdateDisplacementFieldTransform.h"
#include "itkCompositeTransform.h"
#include "itkGaussianSmoothingOnUpdateDisplacementFieldTransform.h"
#include "itkIdentityTransform.h"
#include "itkRigid2DTransform.h"
#include "itkRigid3DTransform.h"
#include "itkTransform.h"

#include "itkBSplineSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor.h"
#include "itkGaussianSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor.h"

#include "itkGradientDescentObjectOptimizer.h"

#include "itkHistogramMatchingImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMacro.h"
#include "itkRegistrationParameterScalesFromShift.h"
#include "itkResampleImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkVector.h"

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

  typename OptionType::Pointer regularizationOption = parser->GetOption( "regularization" );
  if( !regularizationOption || regularizationOption->GetNumberOfValues() != numberOfStages  )
    {
    std::cerr << "The number of regularizations specified does not match the number of stages." << std::endl;
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

  typedef float                                 PixelType;
  typedef double                                RealType;
  typedef itk::Image<PixelType, ImageDimension> FixedImageType;
  typedef itk::Image<PixelType, ImageDimension> MovingImageType;

  typedef itk::CompositeTransform<RealType, ImageDimension> CompositeTransformType;
  typename CompositeTransformType::Pointer compositeTransform = CompositeTransformType::New();
  // We iterate backwards because the command line options are stored as a stack (first in last out)
  for( int currentStage = numberOfStages - 1; currentStage >= 0; currentStage-- )
    {
    typedef itk::SimpleImageRegistrationMethod<FixedImageType, MovingImageType> AffineRegistrationType;

    std::cout << std::endl << "Stage " << numberOfStages - currentStage << std::endl;

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
      }
    else if( std::strcmp( whichTransform.c_str(),
                          "displacementfield" ) == 0 ||  std::strcmp( whichTransform.c_str(), "df" ) == 0 )
      {
      typedef itk::Vector<RealType, ImageDimension> VectorType;
      VectorType zeroVector( 0.0 );
      typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;
      typename DisplacementFieldType::Pointer displacementField = DisplacementFieldType::New();
      displacementField->CopyInformation( fixedImage );
      displacementField->SetRegions( fixedImage->GetBufferedRegion() );
      displacementField->Allocate();
      displacementField->FillBuffer( zeroVector );

      std::string whichRegularization = regularizationOption->GetValue( currentStage );
      ConvertToLowerCase( whichRegularization );
      if( std::strcmp( whichRegularization.c_str(),
                       "gauss" ) == 0 || std::strcmp( whichRegularization.c_str(), "gaussian" ) == 0 )
        {
        typedef itk::GaussianSmoothingOnUpdateDisplacementFieldTransform<RealType,
                                                                         ImageDimension> DisplacementFieldTransformType;

        typedef itk::SimpleImageRegistrationMethod<FixedImageType, MovingImageType,
                                                   DisplacementFieldTransformType> DisplacementFieldRegistrationType;

        // Create the transform adaptors

        typedef itk::DisplacementFieldTransformParametersAdaptor<DisplacementFieldTransformType>
          DisplacementFieldTransformAdaptorType;
        typename DisplacementFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;

        // Extract parameters

        RealType sigmaForUpdateField = parser->Convert<float>( regularizationOption->GetParameter( currentStage, 0 ) );
        RealType sigmaForTotalField = parser->Convert<float>( regularizationOption->GetParameter( currentStage, 1 ) );

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
        }
      else if( std::strcmp( whichRegularization.c_str(),
                            "bspline" ) == 0 || std::strcmp( whichRegularization.c_str(), "dmffd" ) == 0 )
        {
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

        std::vector<unsigned int> meshSizeForTheUpdateField = parser->ConvertVector<unsigned int>( regularizationOption->GetParameter(
                                                                                                     currentStage,
                                                                                                     0 ) );
        std::vector<unsigned int> meshSizeForTheTotalField = parser->ConvertVector<unsigned int>( regularizationOption->GetParameter(
                                                                                                    currentStage,
                                                                                                    1 ) );

        if( meshSizeForTheUpdateField.size() != ImageDimension || meshSizeForTheTotalField.size() != ImageDimension )
          {
          std::cerr << "ERROR:  The mesh size(s) don't match the ImageDimension." << std::endl;
          return EXIT_FAILURE;
          }

        unsigned int splineOrder = 3;
        if( regularizationOption->GetNumberOfParameters( currentStage ) > 2 )
          {
          splineOrder = parser->Convert<unsigned int>( regularizationOption->GetParameter( currentStage, 2 ) );
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

          typedef itk::BSplineSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor<
              DisplacementFieldTransformType> BSplineDisplacementFieldTransformAdaptorType;
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
        }
      else
        {
        std::cerr << "ERROR:  Unrecognized regularization option - " << whichRegularization << std::endl;
        return EXIT_FAILURE;
        }
      }
    }

  // For debugging purposes, we warp the fixed and moving images from the first stage

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

  std::string outputFileName = "hormigita.nii.gz";
  typename OptionType::Pointer outputOption = parser->GetOption( "output" );
  if( outputOption && outputOption->GetNumberOfValues() > 0 )
    {
    outputFileName = outputOption->GetValue();
    }

  typedef itk::ImageFileWriter<FixedImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputFileName.c_str() );
  writer->SetInput( resampler->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}

void InitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

    {
    std::string description = std::string( "" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "dimensionality" );
    option->SetShortName( 'd' );
    option->SetUsageOption( 0, "2/3" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "" );

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
    std::string description = std::string( "" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "transform" );
    option->SetShortName( 't' );
    option->SetUsageOption( 0, "Affine[gradientStep]" );
    option->SetUsageOption( 1, "Field[gradientStep]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "regularization" );
    option->SetShortName( 'r' );
    option->SetUsageOption( 0, "none" );
    option->SetUsageOption( 1, "Gaussian[sigmaInPhysicalSpace]" );
    option->SetUsageOption( 2, "DMFFD[updateFieldMeshSize,totalFieldMeshSize,splineOrder]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "iterations" );
    option->SetShortName( 'i' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "smoothingSigmas" );
    option->SetShortName( 's' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "shrinkFactors" );
    option->SetShortName( 'f' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "useHistogramMatching" );
    option->SetShortName( 'u' );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "" );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "output" );
    option->SetShortName( 'o' );
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

  std::string commandDescription =
    std::string( "hormigita---little ant." );

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
