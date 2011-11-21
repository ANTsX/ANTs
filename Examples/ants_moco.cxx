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
#include "itkCSVNumericObjectFileWriter.h"
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
#include "itkExtractImageFilter.h"

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

template <class T>
inline std::string ants_moco_to_string(const T& t)
{
  std::stringstream ss;

  ss << t;
  return ss.str();
}

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
int ants_moco( itk::ants::CommandLineParser *parser )
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

  typedef float                                     PixelType;
  typedef double                                    RealType;
  typedef itk::Image<PixelType, ImageDimension>     FixedImageType;
  typedef itk::Image<PixelType, ImageDimension + 1> MovingImageType;
  typedef vnl_matrix<double>                        vMatrix;
  vMatrix metric_values;
  vMatrix param_values;

  typedef itk::CompositeTransform<RealType, ImageDimension> CompositeTransformType;
  unsigned int   nparams = 0;
  itk::TimeProbe totalTimer;
  totalTimer.Start();
  // We iterate backwards because the command line options are stored as a stack (first in last out)
  for( int currentStage = numberOfStages - 1; currentStage >= 0; currentStage-- )
    {
    typedef itk::SimpleImageRegistrationMethod<FixedImageType, FixedImageType> AffineRegistrationType;

    std::cout << std::endl << "Stage " << numberOfStages - currentStage << std::endl;
    std::stringstream currentStageString;
    currentStageString << currentStage;

    // Get the fixed and moving images

    std::string fixedImageFileName = metricOption->GetParameter( currentStage, 0 );
    std::string movingImageFileName = metricOption->GetParameter( currentStage, 1 );

    std::cout << "  fixed image: " << fixedImageFileName << std::endl;
    std::cout << "  moving image: " << movingImageFileName << std::endl;

    typedef itk::ImageFileReader<FixedImageType> FixedImageReaderType;
    typename FixedImageReaderType::Pointer fixedImageReader = FixedImageReaderType::New();
    fixedImageReader->SetFileName( fixedImageFileName.c_str() );
    fixedImageReader->Update();
    typename FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();
    fixedImage->Update();
    fixedImage->DisconnectPipeline();

    typedef itk::ImageFileReader<MovingImageType> MovingImageReaderType;
    typename MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
    movingImageReader->SetFileName( movingImageFileName.c_str() );
    movingImageReader->Update();
    typename MovingImageType::Pointer movingImage = movingImageReader->GetOutput();
    movingImage->Update();
    movingImage->DisconnectPipeline();

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

    // the fixed image is a reference image in 3D while the moving is a 4D image
    // loop over every time point and register image_i+1 to image_i
    //
    // Set up the image metric and scales estimator
    unsigned int timedims = movingImage->GetLargestPossibleRegion().GetSize()[ImageDimension - 1];
    metric_values.set_size(timedims, 2);  metric_values.fill(0);
    for( unsigned int timedim = 0;  timedim < timedims - 1;  timedim++ )
      {
      typename CompositeTransformType::Pointer compositeTransform = CompositeTransformType::New();
      typedef itk::IdentityTransform<RealType, ImageDimension> IdentityTransformType;
      typename IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();
      //
      typedef itk::ExtractImageFilter<MovingImageType, FixedImageType> ExtractFilterType;
      typename MovingImageType::RegionType extractRegion = movingImage->GetLargestPossibleRegion();
      extractRegion.SetSize(ImageDimension, 0);
      extractRegion.SetIndex(ImageDimension, timedim );
      typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
      extractFilter->SetInput( movingImage );
      extractFilter->SetDirectionCollapseToSubmatrix();
      extractFilter->SetExtractionRegion( extractRegion );
      extractFilter->Update();
      typename FixedImageType::Pointer fixed_time_slice = extractFilter->GetOutput();

      extractRegion.SetIndex(ImageDimension, timedim + 1 );
      typename ExtractFilterType::Pointer extractFilter2 = ExtractFilterType::New();
      extractFilter2->SetInput( movingImage );
      extractFilter2->SetDirectionCollapseToSubmatrix();
      extractFilter2->SetExtractionRegion( extractRegion );
      extractFilter2->Update();
      typename FixedImageType::Pointer moving_time_slice = extractFilter->GetOutput();

      typedef itk::ImageToImageObjectMetric<FixedImageType, FixedImageType> MetricType;
      typename MetricType::Pointer metric;

      std::string whichMetric = metricOption->GetValue( currentStage );
      ConvertToLowerCase( whichMetric );
      if( std::strcmp( whichMetric.c_str(), "cc" ) == 0 )
        {
        unsigned int radiusOption = parser->Convert<unsigned int>( metricOption->GetParameter( currentStage, 3 ) );

        std::cout << "  using the CC metric (radius = " << radiusOption << ")." << std::endl;
        typedef itk::ANTSNeighborhoodCorrelationImageToImageObjectMetric<FixedImageType,
                                                                         FixedImageType> CorrelationMetricType;
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
        unsigned int npoints_to_skip = parser->Convert<unsigned int>( metricOption->GetParameter( currentStage, 4 ) );
        typedef itk::JointHistogramMutualInformationImageToImageObjectMetric<FixedImageType,
                                                                             FixedImageType> MutualInformationMetricType;
        typename MutualInformationMetricType::Pointer mutualInformationMetric = MutualInformationMetricType::New();
        mutualInformationMetric = mutualInformationMetric;
        mutualInformationMetric->SetNumberOfHistogramBins( binOption );
        mutualInformationMetric->SetDoFixedImagePreWarp( true );
        mutualInformationMetric->SetDoMovingImagePreWarp( true );
        mutualInformationMetric->SetUseMovingImageGradientFilter( false );
        mutualInformationMetric->SetUseFixedImageGradientFilter( false );
        typedef typename MutualInformationMetricType::FixedSampledPointSetType PointSetType;
        typedef typename PointSetType::PointType                               PointType;
        typename PointSetType::Pointer                    pset(PointSetType::New() );
        unsigned long                                     ind = 0, ct = 0;
        itk::ImageRegionIteratorWithIndex<FixedImageType> It(fixed_time_slice,
                                                             fixed_time_slice->GetLargestPossibleRegion() );
        for( It.GoToBegin(); !It.IsAtEnd(); ++It )
          {
          // take every N^th point
          if( ct % npoints_to_skip == 0  ) // about a factor of 5 speed-up over dense
            {
            PointType pt;
            fixed_time_slice->TransformIndexToPhysicalPoint( It.GetIndex(), pt);
            pset->SetPoint(ind, pt);
            ind++;
            }
          ct++;
          }
        mutualInformationMetric->SetFixedSampledPointSet( pset );
        mutualInformationMetric->SetUseFixedSampledPointSet( true );
        metric = mutualInformationMetric;
        }
      else if( std::strcmp( whichMetric.c_str(), "demons" ) == 0 )
        {
        std::cout << "  using the Demons metric." << std::endl;

        typedef itk::DemonsImageToImageObjectMetric<FixedImageType, FixedImageType> DemonsMetricType;
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
      //    typedef itk::JointHistogramMutualInformationImageToImageObjectMetric<FixedImageType, FixedImageType>
      // MutualInformationMetricType;
      typedef itk::RegistrationParameterScalesFromShift<MetricType> ScalesEstimatorType;
      typename ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();
      scalesEstimator->SetMetric( metric );
      scalesEstimator->SetTransformForward( true );
      scalesEstimator->SetSamplingStrategy(ScalesEstimatorType::CornerSampling);

      float learningRate = parser->Convert<float>( transformOption->GetParameter( currentStage, 0 ) );

      typedef itk::GradientDescentObjectOptimizer GradientDescentObjectOptimizerType;
      typename GradientDescentObjectOptimizerType::Pointer optimizer = GradientDescentObjectOptimizerType::New();
      optimizer->SetLearningRate( learningRate );
      optimizer->SetNumberOfIterations( iterations[0] );
      optimizer->SetScalesEstimator( scalesEstimator );
      double small_step = 0;
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        small_step += fixed_time_slice->GetSpacing()[i] * fixed_time_slice->GetSpacing()[i];
        }
      optimizer->SetMaximumStepSizeInPhysicalSpaceUnits(sqrt(small_step) * 0.1);

      // Set up the image registration methods along with the transforms
      std::string whichTransform = transformOption->GetValue( currentStage );
      ConvertToLowerCase( whichTransform );
      if( std::strcmp( whichTransform.c_str(), "affine" ) == 0 )
        {
        typename AffineRegistrationType::Pointer affineRegistration = AffineRegistrationType::New();
        typedef itk::AffineTransform<double, ImageDimension> AffineTransformType;
        typename AffineTransformType::Pointer affineTransform = AffineTransformType::New();
        nparams = affineTransform->GetNumberOfParameters();
        affineTransform->SetIdentity();
        affineRegistration->SetFixedImage( fixed_time_slice );
        affineRegistration->SetMovingImage( moving_time_slice );
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
          std::cout << std::endl << "*** Running affine registration ***" << timedim << std::endl << std::endl;
          affineRegistration->StartRegistration();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cerr << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }

        // Write out the affine transform
        std::string filename = outputPrefix + std::string("TimeSlice") + ants_moco_to_string<unsigned int>(timedim)
          + std::string( "Affine.txt" );
        typedef itk::TransformFileWriter TransformWriterType;
        typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
        transformWriter->SetInput( affineRegistration->GetOutput()->Get() );
        transformWriter->SetFileName( filename.c_str() );
        transformWriter->Update();
        if( timedim == 0 )
          {
          param_values.set_size(timedims, nparams);
          param_values.fill(0);
          }
        for( unsigned int i = 0; i < nparams; i++ )
          {
          param_values(timedim, i) = affineRegistration->GetOutput()->Get()->GetParameters()[i];
          }
        }
      else if( std::strcmp( whichTransform.c_str(), "rigid" ) == 0 )
        {
        typedef itk::Euler3DTransform<double>
                                                                       RigidTransformType;
        typedef itk::SimpleImageRegistrationMethod<FixedImageType, FixedImageType,
                                                   RigidTransformType> RigidRegistrationType;
        typename RigidRegistrationType::Pointer rigidRegistration = RigidRegistrationType::New();
        typename RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
        nparams = rigidTransform->GetNumberOfParameters();
        rigidTransform->SetIdentity();
        rigidRegistration->SetFixedImage( fixed_time_slice );
        rigidRegistration->SetMovingImage( moving_time_slice );
        rigidRegistration->SetNumberOfLevels( numberOfLevels );
        rigidRegistration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );
        rigidRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
        rigidRegistration->SetMetric( metric );
        rigidRegistration->SetOptimizer( optimizer );
        rigidRegistration->SetTransform( rigidTransform );
        rigidRegistration->SetCompositeTransform( compositeTransform );
        typedef CommandIterationUpdate<RigidRegistrationType> RigidCommandType;
        typename RigidCommandType::Pointer rigidObserver = RigidCommandType::New();
        rigidObserver->SetNumberOfIterations( iterations );
        rigidRegistration->AddObserver( itk::IterationEvent(), rigidObserver );
        try
          {
          std::cout << std::endl << "*** Running rigid registration ***" << timedim  << std::endl << std::endl;
          rigidRegistration->StartRegistration();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cerr << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }
        // Write out the rigid transform
        std::string filename = outputPrefix + std::string("TimeSlice") + ants_moco_to_string<unsigned int>(timedim)
          + std::string( "Rigid.txt" );
        typedef itk::TransformFileWriter TransformWriterType;
        typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
        transformWriter->SetInput( rigidRegistration->GetOutput()->Get() );
        transformWriter->SetFileName( filename.c_str() );
        transformWriter->Update();
        if( timedim == 0 )
          {
          param_values.set_size(timedims, nparams);
          param_values.fill(0);
          }
        for( unsigned int i = 0; i < nparams; i++ )
          {
          param_values(timedim, i) = rigidRegistration->GetOutput()->Get()->GetParameters()[i];
          }
        }
      else
        {
        std::cerr << "ERROR:  Unrecognized transform option - " << whichTransform << std::endl;
        return EXIT_FAILURE;
        }
      metric_values(timedim, 0) = metric->GetValue();
      }
    }
  // Write out warped image(s), if requested.
  totalTimer.Stop();
  metric_values(metric_values.rows() - 1, 0) = metric_values(metric_values.rows() - 2, 0);
  for( unsigned int i = 0; i < nparams; i++ )
    {
    param_values(metric_values.rows() - 1, i) = param_values(metric_values.rows() - 2, i);
    }
  std::cout << std::endl << "Total elapsed time: " << totalTimer.GetMeanTime() << std::endl;
    {
    std::vector<std::string> ColumnHeaders;
    for( unsigned int nv = 0; nv < nparams; nv++ )
      {
      std::string colname = std::string("MOCOparam") + ants_moco_to_string<unsigned int>(nv);
      ColumnHeaders.push_back( colname );
      }
    typedef itk::CSVNumericObjectFileWriter<double> WriterType;
    WriterType::Pointer writer = WriterType::New();
    std::string         fnmp = outputPrefix + std::string("MOCOparams.csv");
    std::cout << " write " << fnmp << std::endl;
    writer->SetFileName( fnmp.c_str() );
    writer->SetColumnHeaders(ColumnHeaders);
    writer->SetInput( &param_values );
    writer->Write();
    }
    {
    std::vector<std::string> ColumnHeaders;
    std::string              colname = std::string("MetricPost");
    ColumnHeaders.push_back( colname );
    colname = std::string("MetricPre");
    ColumnHeaders.push_back( colname );
    typedef itk::CSVNumericObjectFileWriter<double> WriterType;
    WriterType::Pointer writer = WriterType::New();
    std::string         fnmp = outputPrefix + std::string("Metric.csv");
    std::cout << " write " << fnmp << std::endl;
    writer->SetFileName( fnmp.c_str() );
    writer->SetColumnHeaders(ColumnHeaders);
    writer->SetInput( &metric_values );
    writer->Write();
    }
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
    option->SetUsageOption( 1, "MI[fixedImage,movingImage,metricWeight,numberOfBins,n_points_to_skip]" );
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
    option->SetUsageOption( 0, "Affine[gradientStep]" );
    option->SetUsageOption( 1, "Rigid[gradientStep]" );
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

  std::string commandDescription = std::string( "ants_moco = motion correction.  This program is a user-level " )
    + std::string( "registration application meant to utilize ITKv4-only classes. The user can specify " )
    + std::string( "any number of \"stages\" where a stage consists of a transform; an image metric; " )
    + std::string( " and iterations, shrink factors, and smoothing sigmas for each level. " )
    + std::string(
      " Specialized for 4D time series data: fixed image is 3D, moving image should be the 4D time series. ")
    + std::string( " Fixed image is a reference space or time slice.");
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
    //   case 2:
    //  hormigita<2>( parser );
    //  break;
    case 3:
      {
      ants_moco<3>( parser );
      }
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
}
