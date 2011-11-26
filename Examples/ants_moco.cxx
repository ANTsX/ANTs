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
#include "itkQuasiNewtonObjectOptimizer.h"

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

    typedef itk::GradientDescentObjectOptimizer OptimizerType;
    typedef itk::QuasiNewtonObjectOptimizer     OptimizerType2;

    OptimizerType * optimizer = reinterpret_cast<OptimizerType *>(
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

template <class TImageIn, class TImageOut>
void
AverageTimeImages( typename TImageIn::Pointer image_in,  typename TImageOut::Pointer image_avg, unsigned int time_dims )
{
  std::cout << " averaging images " << std::endl;

  typedef TImageIn  ImageType;
  typedef TImageOut OutImageType;
  enum { ImageDimension = ImageType::ImageDimension };
  typedef float                                           PixelType;
  typedef itk::ImageRegionIteratorWithIndex<OutImageType> Iterator;
  image_avg->FillBuffer(0);
  Iterator vfIter2(  image_avg, image_avg->GetLargestPossibleRegion() );
  for(  vfIter2.GoToBegin(); !vfIter2.IsAtEnd(); ++vfIter2 )
    {
    typename OutImageType::PixelType  fval = 0;
    typename ImageType::IndexType ind;
    typename OutImageType::IndexType spind = vfIter2.GetIndex();
    for( unsigned int xx = 0; xx < time_dims; xx++ )
      {
      for( unsigned int yy = 0; yy < ImageDimension - 1; yy++ )
        {
        ind[yy] = spind[yy];
        }
      ind[ImageDimension - 1] = xx;
      fval += image_in->GetPixel(ind);
      }
    fval /= (double)time_dims;
    image_avg->SetPixel(spind, fval);
    }
  std::cout << " averaging images done " << std::endl;
  return;
}

template <unsigned int ImageDimension, class TRigid>
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
  vMatrix param_values;
  typedef itk::CompositeTransform<RealType, ImageDimension> CompositeTransformType;
  std::vector<typename CompositeTransformType::Pointer> CompositeTransformVector;

  unsigned int   nparams = 2;
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
    typename FixedImageType::Pointer fixed_time_slice = NULL;
    typename FixedImageType::Pointer moving_time_slice = NULL;

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

    typename MovingImageReaderType::Pointer outputImageReader = MovingImageReaderType::New();
    outputImageReader->SetFileName( movingImageFileName.c_str() );
    outputImageReader->Update();
    typename MovingImageType::Pointer outputImage = outputImageReader->GetOutput();
    outputImage->Update();
    outputImage->DisconnectPipeline();

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
        smoothingSigmasPerLevel[n] = sigmas[n];
        }
      std::cout << "  smoothing sigmas per level: " << smoothingSigmasPerLevel << std::endl;
      }

    // the fixed image is a reference image in 3D while the moving is a 4D image
    // loop over every time point and register image_i+1 to image_i
    //
    // Set up the image metric and scales estimator
    unsigned int timedims = movingImage->GetLargestPossibleRegion().GetSize()[ImageDimension];
    for( unsigned int timedim = 0;  timedim < timedims;  timedim++ )
      {
      typename CompositeTransformType::Pointer compositeTransform = NULL;
      if(  currentStage == (numberOfStages - 1) )
        {
        compositeTransform = CompositeTransformType::New();
        CompositeTransformVector.push_back(compositeTransform);
        }
      else if( CompositeTransformVector.size() == timedims && !CompositeTransformVector[timedim].IsNull() )
        {
        compositeTransform = CompositeTransformVector[timedim];
        std::cout << " use existing transform " << compositeTransform->GetParameters() << std::endl;
        }
      typedef itk::IdentityTransform<RealType, ImageDimension> IdentityTransformType;
      typename IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();
      //
      typedef itk::ExtractImageFilter<MovingImageType, FixedImageType> ExtractFilterType;
      typename MovingImageType::RegionType extractRegion = movingImage->GetLargestPossibleRegion();
      extractRegion.SetSize(ImageDimension, 0);
      bool maptoneighbor = true;
      typename OptionType::Pointer fixedOption = parser->GetOption( "useFixedReferenceImage" );
      if( fixedOption && fixedOption->GetNumberOfValues() > 0 )
        {
        std::string fixedValue = fixedOption->GetValue( 0 );
        ConvertToLowerCase( fixedValue );
        if( fixedValue.compare( "1" ) == 0 || fixedValue.compare( "true" ) == 0 )
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
        extractFilter2->SetExtractionRegion( extractRegion );
        extractFilter2->Update();
        moving_time_slice = extractFilter2->GetOutput();
        }

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
      // scalesEstimator->SetSamplingStrategy(ScalesEstimatorType::CornerSampling);

      float learningRate = parser->Convert<float>( transformOption->GetParameter( currentStage, 0 ) );

      typedef itk::GradientDescentObjectOptimizer OptimizerType;
      typedef itk::QuasiNewtonObjectOptimizer     OptimizerType2;
      typename OptimizerType::Pointer optimizer = OptimizerType::New();
      optimizer->SetLearningRate( learningRate );
      optimizer->SetNumberOfIterations( iterations[0] );
      typename OptionType::Pointer scalesOption = parser->GetOption( "useScalesEstimator" );
      if( scalesOption && scalesOption->GetNumberOfValues() > 0 )
        {
        std::string scalesValue = scalesOption->GetValue( 0 );
        ConvertToLowerCase( scalesValue );
        if( scalesValue.compare( "1" ) == 0 || scalesValue.compare( "true" ) == 0 )
          {
          std::cout << " employing scales estimator " << std::endl;
          optimizer->SetScalesEstimator( scalesEstimator );
          }
        else
          {
          std::cout << " not employing scales estimator " << scalesValue << std::endl;
          }
        }
      double small_step = 0;
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        small_step += fixed_time_slice->GetSpacing()[i] * fixed_time_slice->GetSpacing()[i];
        }
      optimizer->SetMaximumStepSizeInPhysicalUnits(sqrt(small_step) * learningRate);
      //    optimizer->SetMaximumNewtonStepSizeInPhysicalUnits(sqrt(small_step)*learningR);

      // Set up the image registration methods along with the transforms
      std::string whichTransform = transformOption->GetValue( currentStage );
      ConvertToLowerCase( whichTransform );
      if( std::strcmp( whichTransform.c_str(), "affine" ) == 0 )
        {
        typename AffineRegistrationType::Pointer affineRegistration = AffineRegistrationType::New();
        typedef itk::AffineTransform<double, ImageDimension> AffineTransformType;
        typename AffineTransformType::Pointer affineTransform = AffineTransformType::New();
        affineTransform->SetIdentity();
        nparams = affineTransform->GetNumberOfParameters() + 2;
        typename ScalesEstimatorType::ScalesType scales(affineTransform->GetNumberOfParameters() );
        metric->SetFixedImage( fixed_time_slice );
        metric->SetVirtualDomainImage( fixed_time_slice );
        metric->SetMovingImage( moving_time_slice );
        metric->SetTransform( affineTransform );
        scalesEstimator->SetMetric(metric);
        scalesEstimator->EstimateScales(scales);
        optimizer->SetScales(scales);
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
        typedef TRigid
                                                                       RigidTransformType;
        typedef itk::SimpleImageRegistrationMethod<FixedImageType, FixedImageType,
                                                   RigidTransformType> RigidRegistrationType;
        typename RigidRegistrationType::Pointer rigidRegistration = RigidRegistrationType::New();
        typename RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
        rigidTransform->SetIdentity();
        nparams = rigidTransform->GetNumberOfParameters() + 2;
        typename ScalesEstimatorType::ScalesType scales(rigidTransform->GetNumberOfParameters() );
        metric->SetFixedImage( fixed_time_slice );
        metric->SetVirtualDomainImage( fixed_time_slice );
        metric->SetMovingImage( moving_time_slice );
        metric->SetTransform( rigidTransform );
        scalesEstimator->SetMetric(metric);
        scalesEstimator->EstimateScales(scales);
        optimizer->SetScales(scales);

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
        typedef itk::Vector<RealType, ImageDimension> VectorType;
        VectorType zeroVector( 0.0 );
        typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;
        typename DisplacementFieldType::Pointer displacementField = DisplacementFieldType::New();
        displacementField->CopyInformation( fixed_time_slice );
        displacementField->SetRegions(  fixed_time_slice->GetBufferedRegion() );
        displacementField->Allocate();
        displacementField->FillBuffer( zeroVector );
        typedef itk::GaussianSmoothingOnUpdateDisplacementFieldTransform<RealType,
                                                                         ImageDimension> DisplacementFieldTransformType;

        typedef itk::SimpleImageRegistrationMethod<FixedImageType, FixedImageType,
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
          fieldTransformAdaptor->SetTransform( gaussianFieldTransform );

          adaptors.push_back( fieldTransformAdaptor.GetPointer() );
          }

        typename DisplacementFieldRegistrationType::Pointer displacementFieldRegistration =
          DisplacementFieldRegistrationType::New();
        displacementFieldRegistration->SetFixedImage( fixed_time_slice );
        displacementFieldRegistration->SetMovingImage( moving_time_slice );
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
                    << sigmaForUpdateField << ", sigmaForTotalField = " << sigmaForTotalField << ") ***"
                    << " timedim " << timedim << std::endl << std::endl;
          displacementFieldRegistration->StartRegistration();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cerr << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }
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
      if( std::strcmp( whichTransform.c_str(),
                       "gaussiandisplacementfield" ) != 0 && std::strcmp( whichTransform.c_str(), "gdf" ) != 0 )
        {
        param_values(timedim, 1) = metric->GetValue();
        }

      // resample the moving image and then put it in its place
      typedef itk::ResampleImageFilter<FixedImageType, FixedImageType> ResampleFilterType;
      typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
      resampler->SetTransform( compositeTransform );
      resampler->SetInput( moving_time_slice );
      resampler->SetSize(  moving_time_slice->GetLargestPossibleRegion().GetSize() );
      resampler->SetOutputOrigin(   moving_time_slice->GetOrigin() );
      resampler->SetOutputSpacing(  moving_time_slice->GetSpacing() );
      resampler->SetOutputDirection(  moving_time_slice->GetDirection() );
      resampler->SetDefaultPixelValue( 0 );
      std::cout << " resampling " << std::endl;
      resampler->Update();
      std::cout << " done resampling " << std::endl;
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
        ind[ImageDimension] = timedim;
        outputImage->SetPixel(ind, fval);
        }
      }
    if( outputOption && outputOption->GetNumberOfParameters( 0 ) > 1 )
      {
      std::string fileName = outputOption->GetParameter( 0, 1 );
      typedef itk::ImageFileWriter<MovingImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( fileName.c_str() );
      writer->SetInput( outputImage );
      writer->Update();
      }
    if( outputOption && outputOption->GetNumberOfParameters( 0 ) > 2 && outputImage )
      {
      std::string fileName = outputOption->GetParameter( 0, 2 );
      typename FixedImageType::Pointer avgImage;
      typedef itk::ExtractImageFilter<MovingImageType, FixedImageType> ExtractFilterType;
      typename MovingImageType::RegionType extractRegion = movingImage->GetLargestPossibleRegion();
      extractRegion.SetSize(ImageDimension, 0);
      typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
      extractFilter->SetInput( movingImage );
      extractFilter->SetDirectionCollapseToSubmatrix();
      unsigned int td = 0;
      extractRegion.SetIndex(ImageDimension, td );
      extractFilter->SetExtractionRegion( extractRegion );
      extractFilter->Update();
      avgImage = extractFilter->GetOutput();
      AverageTimeImages<MovingImageType, FixedImageType>( outputImage, avgImage, timedims );
      typedef itk::ImageFileWriter<FixedImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( fileName.c_str() );
      writer->SetInput( avgImage );
      writer->Update();
      std::cout << " done writing avg image " << std::endl;
      }
    }
  totalTimer.Stop();
  std::cout << std::endl << "Total elapsed time: " << totalTimer.GetMeanTime() << std::endl;
    {
    std::vector<std::string> ColumnHeaders;
    std::string              colname;
    colname = std::string("MetricPre");
    ColumnHeaders.push_back( colname );
    colname = std::string("MetricPost");
    ColumnHeaders.push_back( colname );
    for( unsigned int nv = 2; nv < nparams; nv++ )
      {
      std::string colname = std::string("MOCOparam") + ants_moco_to_string<unsigned int>(nv - 2);
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
    option->SetUsageOption( 2,
                            "GaussianDisplacementField[gradientStep,updateFieldSigmaInPhysicalSpace,totalFieldSigmaInPhysicalSpace]" );
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
  typedef itk::Euler2DTransform<double> Euler2D;
  typedef itk::Euler3DTransform<double> Euler3D;

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

  std::cout << std::endl << "Running " << argv[0] << "  for " << dimension << "-dimensional images." << std::endl
            << std::endl;

  switch( dimension )
    {
    case 2:
      {
      ants_moco<2, Euler2D>( parser );
      }
      break;
    case 3:
      {
      ants_moco<3, Euler3D>( parser );
      }
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
}
