/*
   Instantiating all 4 combinations of 2D,3D and float and double in one file
   was causing object files that were too big to be linked with gcc44, and
   is anticipated to cause problems on Windows machines that require relatively
   small object files as well.

   This file will use explicit template instantiation to make the overall size of
   each object smaller.
 */
#ifndef  __ANTSREGISTRATIONTEMPLATEHEADER_H__
#define  __ANTSREGISTRATIONTEMPLATEHEADER_H__

#include "antsUtilities.h"
#include "itkantsRegistrationHelper.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkGaussianInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkLabelImageGaussianInterpolateImageFunction.h"
#include "include/antsRegistration.h"

namespace ants {

extern const char * RegTypeToFileName(const std::string & type, bool & writeInverse, bool & writeVelocityField);

template <class TComputeType, unsigned VImageDimension>
int
DoRegistration(typename ParserType::Pointer & parser)
{
  typedef TComputeType RealType;
  typedef typename ants::RegistrationHelper<TComputeType, VImageDimension>   RegistrationHelperType;
  typedef typename RegistrationHelperType::ImageType                         ImageType;
  typedef typename RegistrationHelperType::CompositeTransformType            CompositeTransformType;

  typename RegistrationHelperType::Pointer regHelper = RegistrationHelperType::New();

  OptionType::Pointer transformOption = parser->GetOption( "transform" );

  OptionType::Pointer metricOption = parser->GetOption( "metric" );

  OptionType::Pointer convergenceOption = parser->GetOption( "convergence" );

  OptionType::Pointer shrinkFactorsOption = parser->GetOption( "shrink-factors" );

  OptionType::Pointer smoothingSigmasOption = parser->GetOption( "smoothing-sigmas" );

  OptionType::Pointer restrictDeformationOption = parser->GetOption( "restrict-deformation" );

  OptionType::Pointer outputOption = parser->GetOption( "output" );

  OptionType::Pointer maskOption = parser->GetOption( "masks" );

  OptionType::Pointer compositeOutputOption = parser->GetOption( "write-composite-transform" );

  OptionType::Pointer collapseOutputTransformsOption = parser->GetOption( "collapse-output-transforms" );

  if( !outputOption || outputOption->GetNumberOfFunctions() == 0 )
    {
    std::cout << "Output option not specified." << std::endl;
    return EXIT_FAILURE;
    }

  OptionType::Pointer collapseLinearTransforms =
    parser->GetOption( "collapse-linear-transforms-to-fixed-image-header" );
  if( collapseLinearTransforms && parser->Convert<bool>( collapseLinearTransforms->GetFunction( 0 )->GetName() ) )
    {
    regHelper->SetApplyLinearTransformsToFixedImageHeader( true );
    }
  else
    {
    regHelper->SetApplyLinearTransformsToFixedImageHeader( false );
    }

  OptionType::Pointer printSimilarityMeasureInterval = parser->GetOption( "print-similarity-measure-interval" );
  if( printSimilarityMeasureInterval && printSimilarityMeasureInterval->GetNumberOfFunctions() )
    {
    unsigned int intervalLength = parser->Convert<unsigned int>( printSimilarityMeasureInterval->GetFunction(
                                                                   0 )->GetName() );
    regHelper->SetPrintSimilarityMeasureInterval(intervalLength);
    }
  else
    {
    regHelper->SetPrintSimilarityMeasureInterval( 0 );
    }

  OptionType::Pointer writeIntervalVolumes = parser->GetOption( "write-interval-volumes" );
  if( writeIntervalVolumes && writeIntervalVolumes->GetNumberOfFunctions() )
    {
    unsigned int LengthOfIntervals = parser->Convert<unsigned int>( writeIntervalVolumes->GetFunction( 0 )->GetName() );
    regHelper->SetWriteIntervalVolumes(LengthOfIntervals);
    }
  else
    {
    regHelper->SetPrintSimilarityMeasureInterval( 0 );
    }

  std::string outputPrefix = outputOption->GetFunction( 0 )->GetName();
  if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
    {
    outputPrefix = outputOption->GetFunction( 0 )->GetParameter( 0 );
    }
  std::string outputWarpedImageName;
  if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
    {
    outputWarpedImageName = outputOption->GetFunction( 0 )->GetParameter( 1 );
    }

  std::string outputInverseWarpedImageName;
  if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
    {
    outputInverseWarpedImageName = outputOption->GetFunction( 0 )->GetParameter( 2 );
    }

  ParserType::OptionType::Pointer initialMovingTransformOption = parser->GetOption( "initial-moving-transform" );

  if( initialMovingTransformOption && initialMovingTransformOption->GetNumberOfFunctions() )
    {
    std::vector<bool> isDerivedInitialMovingTransform;
    typename CompositeTransformType::Pointer compositeTransform =
      GetCompositeTransformFromParserOption<TComputeType, VImageDimension>( parser, initialMovingTransformOption,
                                                              isDerivedInitialMovingTransform );

    if( compositeTransform.IsNull() )
      {
      return EXIT_FAILURE;
      }
    regHelper->SetMovingInitialTransform( compositeTransform );

    // Write out initial derived transforms only if we're not collapsing them in the output
    if( !parser->Convert<bool>( collapseOutputTransformsOption->GetFunction( 0 )->GetName() ) )
      {
      for( unsigned int n = 0; n < isDerivedInitialMovingTransform.size(); n++ )
        {
        std::stringstream curFileName;
        curFileName << outputPrefix << n << "DerivedInitialMovingTranslation.mat";

        typename RegistrationHelperType::CompositeTransformType::TransformTypePointer curTransform =
          compositeTransform->GetNthTransform( n );
        itk::ants::WriteTransform<TComputeType, VImageDimension>( curTransform, curFileName.str() );
        }
      }
    }

  ParserType::OptionType::Pointer initialFixedTransformOption = parser->GetOption( "initial-fixed-transform" );

  if( initialFixedTransformOption && initialFixedTransformOption->GetNumberOfFunctions() )
    {
    std::vector<bool> isDerivedInitialFixedTransform;
    typename CompositeTransformType::Pointer compositeTransform =
      GetCompositeTransformFromParserOption<TComputeType, VImageDimension>( parser, initialFixedTransformOption,
                                                              isDerivedInitialFixedTransform );
    if( compositeTransform.IsNull() )
      {
      return EXIT_FAILURE;
      }
    regHelper->SetFixedInitialTransform( compositeTransform );

    // Write out initial derived transforms only if we're not collapsing them in the output
    if( !parser->Convert<bool>( collapseOutputTransformsOption->GetFunction( 0 )->GetName() ) )
      {
      for( unsigned int n = 0; n < isDerivedInitialFixedTransform.size(); n++ )
        {
        std::stringstream curFileName;
        curFileName << outputPrefix << n << "DerivedInitialFixedTranslation.mat";

        typename RegistrationHelperType::CompositeTransformType::TransformTypePointer curTransform =
          compositeTransform->GetNthTransform( n );
        itk::ants::WriteTransform<TComputeType, VImageDimension>( curTransform, curFileName.str() );
        }
      }
    }

  if( maskOption.IsNotNull() && maskOption->GetNumberOfFunctions() )
    {
    typedef typename RegistrationHelperType::MaskImageType MaskImageType;
    typedef itk::ImageFileReader<MaskImageType>            ImageReaderType;
    for( unsigned m = 0; m < maskOption->GetFunction( 0 )->GetNumberOfParameters(); m++ )
      {
      std::string fname = maskOption->GetFunction( 0 )->GetParameter( m );

      typename MaskImageType::Pointer maskImage;
      typename ImageReaderType::Pointer reader = ImageReaderType::New();

      reader->SetFileName( fname.c_str() );
      try
        {
        reader->Update();
        maskImage = reader->GetOutput();
        }
      catch( itk::ExceptionObject & err )
        {
        std::cout << "Can't read specified mask image " << fname.c_str() << std::endl;
        std::cout << "Exception Object caught: " << std::endl;
        std::cout << err << std::endl;
        return EXIT_FAILURE;
        }
      if( m == 0 )
        {
        regHelper->SetFixedImageMask( maskImage );
        }
      else if( m == 1 )
        {
        regHelper->SetMovingImageMask( maskImage );
        }
      }
    }

  // The misc. options
  //  * winsorize image intensities
  //  * use histogram matching
  //  * estimate learning rate
  // are currently specified once on the command line and then apply to all
  // stages.  Advanced parameter specification might require us to rewrite
  // this in the future.

  float lowerQuantile = 0.0;
  float upperQuantile = 1.0;

  bool doWinsorize = false;

  OptionType::Pointer winsorizeOption = parser->GetOption( "winsorize-image-intensities" );
  if( winsorizeOption && winsorizeOption->GetNumberOfFunctions() )
    {
    doWinsorize = true;
    if( winsorizeOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      lowerQuantile = parser->Convert<float>( winsorizeOption->GetFunction( 0 )->GetParameter( 0 ) );
      }
    if( winsorizeOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
      {
      upperQuantile = parser->Convert<float>( winsorizeOption->GetFunction( 0 )->GetParameter( 1 ) );
      }
    }
  regHelper->SetWinsorizeImageIntensities( doWinsorize, lowerQuantile, upperQuantile );

  bool doHistogramMatch = false;

  OptionType::Pointer histOption = parser->GetOption( "use-histogram-matching" );
  if( histOption && histOption->GetNumberOfFunctions() )
    {
    std::string histFunction = histOption->GetFunction( 0 )->GetName();
    ConvertToLowerCase( histFunction );
    if( histFunction.compare( "1" ) == 0 || histFunction.compare( "true" ) == 0 )
      {
      doHistogramMatch = true;
      }
    }
  regHelper->SetUseHistogramMatching( doHistogramMatch );

  bool doEstimateLearningRateAtEachIteration = true;

  OptionType::Pointer rateOption = parser->GetOption( "use-estimate-learning-rate-once" );
  if( rateOption && rateOption->GetNumberOfFunctions() )
    {
    std::string rateFunction = rateOption->GetFunction( 0 )->GetName();
    ConvertToLowerCase( rateFunction );
    if( rateFunction.compare( "1" ) == 0 || rateFunction.compare( "true" ) == 0 )
      {
      doEstimateLearningRateAtEachIteration = false;
      }
    }
  regHelper->SetDoEstimateLearningRateAtEachIteration( doEstimateLearningRateAtEachIteration );

  if( restrictDeformationOption && restrictDeformationOption->GetNumberOfFunctions() )
    {
    std::vector<RealType> restrictDeformationWeights =
      parser->ConvertVector<RealType>( restrictDeformationOption->GetFunction( 0 )->GetName() );
    if( restrictDeformationWeights.size() != VImageDimension )
      {
      std::cout << "The restrict deformation weights vector should be the same as the "
        << "number of local parameters (=ImageDimension)." << std::endl;
      }

    regHelper->SetRestrictDeformationOptimizerWeights( restrictDeformationWeights );
    }

  // We find both the number of transforms and the number of metrics

  unsigned int numberOfTransforms = transformOption->GetNumberOfFunctions();
  if( transformOption.IsNull() || numberOfTransforms == 0 )
    {
    std::cout << "No transformations are specified." << std::endl;
    return EXIT_FAILURE;
    }

  std::vector<std::vector<unsigned int> > iterationList;
  std::vector<RealType>                   convergenceThresholdList;
  std::vector<unsigned int>               convergenceWindowSizeList;
  std::vector<std::vector<unsigned int> > shrinkFactorsList;
  std::vector<std::vector<float> >        smoothingSigmasList;
  std::vector<bool>                       smoothingSigmasAreInPhysicalUnitsList;
  std::deque<std::string>                 TransformTypeNames;
  // Each image registration "stage" is characterized by
  //   * a transform
  //   * a set of convergence criteria (number of iterations, convergence threshold,
  //     and/or convergence window)
  //   * a set of shrink factors
  //   * a set of smoothing factors (specified in physical space or voxel space)
  // Note that the set of the number of iterations, the set of shrink factors, and the
  // set of smoothing factors imply the number of levels that will be used for that stage
  // and thus they all need to be the same vector length, e.g. "-c 100x50x10 -f 4x2x1 -s 2x1x0".
  // We use the number of transforms to implicitly guess how many stages are being used in
  // current registration call.  We don't add the metrics in this loop as there could be
  // more than one metric per stage.
  //
  // Also, we iterate backwards because the command line options are stored as a stack (first
  // in last out).
  for( int currentStage = numberOfTransforms - 1; currentStage >= 0; currentStage-- )
    {

    // Get the number of iterations and use that information to specify the number of levels
    std::vector<unsigned int> iterations;
    RealType                  convergenceThreshold = 1e-6;
    unsigned int              convergenceWindowSize = 10;

    if( convergenceOption.IsNotNull() && convergenceOption->GetNumberOfFunctions() )
      {
      if( convergenceOption->GetFunction( currentStage )->GetNumberOfParameters() == 0 )
        {
        iterations = parser->ConvertVector<unsigned int>( convergenceOption->GetFunction( currentStage )->GetName() );
        }
      else if( convergenceOption->GetFunction( currentStage )->GetNumberOfParameters() > 0 )
        {
        iterations =
          parser->ConvertVector<unsigned int>( convergenceOption->GetFunction( currentStage )->GetParameter( 0 ) );
        }
      if( convergenceOption->GetFunction( currentStage )->GetNumberOfParameters() > 1 )
        {
        convergenceThreshold = parser->Convert<RealType>( convergenceOption->GetFunction( currentStage )->GetParameter(
                                                          1 ) );
        }
      if( convergenceOption->GetFunction( currentStage )->GetNumberOfParameters() > 2 )
        {
        convergenceWindowSize = parser->Convert<unsigned int>( convergenceOption->GetFunction(
                                                                 currentStage )->GetParameter( 2 ) );
        const unsigned int minAllowedconvergenceWindowSize = 2; // The BSplineScatteredDataPoints requires at least 2
                                                                // points for interpolation.
        if( convergenceWindowSize < minAllowedconvergenceWindowSize )
          {
          std::cout << "Convergence Window Size must be greater than or equal to " << minAllowedconvergenceWindowSize
                   << std::endl;
          }
        }
      }
    else
      {
      std::cout << "No convergence criteria are specified." << std::endl;
      return EXIT_FAILURE;
      }

    iterationList.push_back( iterations );
    convergenceThresholdList.push_back( convergenceThreshold );
    convergenceWindowSizeList.push_back( convergenceWindowSize );

    unsigned int numberOfLevels = iterations.size();
    std::cout << "  number of levels = " << numberOfLevels << std::endl;

    // Get the first metricOption for the currentStage (for use with the B-spline transforms)
				unsigned int numberOfMetrics = metricOption->GetNumberOfFunctions();
				std::string fixedImageFileName;
				for( int currentMetricNumber = numberOfMetrics - 1; currentMetricNumber >= 0; currentMetricNumber-- )
						{
						// Get the fixed filename to read later in the case of a B-spline transform
						unsigned int stageID = metricOption->GetFunction( currentMetricNumber )->GetStageID();
      if( stageID == static_cast<unsigned int>( currentStage ) )
        {
    				fixedImageFileName = metricOption->GetFunction( currentMetricNumber )->GetParameter( 0 );
        break;
        }
      }

    // Get shrink factors
    std::vector<unsigned int> factors =
      parser->ConvertVector<unsigned int>( shrinkFactorsOption->GetFunction( currentStage )->GetName() );
    shrinkFactorsList.push_back( factors );

    // Get smoothing sigmas
    std::string  smoothingSigmasString = smoothingSigmasOption->GetFunction( currentStage )->GetName();
    const size_t mmPosition = smoothingSigmasString.find( "mm" );
    const size_t voxPosition = smoothingSigmasString.find( "vox" );
    if( mmPosition != std::string::npos )
      {
      smoothingSigmasString.replace( mmPosition, 2, "" );
      smoothingSigmasAreInPhysicalUnitsList.push_back( true );
      }
    else if( voxPosition != std::string::npos )
      {
      smoothingSigmasString.replace( voxPosition, 3, "" );
      smoothingSigmasAreInPhysicalUnitsList.push_back( false );
      }
    else
      {
      smoothingSigmasAreInPhysicalUnitsList.push_back( false );
      }

    std::vector<float> sigmas = parser->ConvertVector<float>( smoothingSigmasString );
    if( sigmas.size() == 1 )
      {
      sigmas.resize( numberOfLevels, sigmas[0] );
      }
    smoothingSigmasList.push_back( sigmas );

    // Set up the optimizer.  To change the iteration number for each level we rely
    // on the command observer.

    float learningRate = parser->Convert<float>( transformOption->GetFunction( currentStage )->GetParameter( 0 ) );

    std::string whichTransform = transformOption->GetFunction( currentStage )->GetName();
    ConvertToLowerCase( whichTransform );

    TransformTypeNames.push_back( whichTransform );

    typename RegistrationHelperType::XfrmMethod xfrmMethod = regHelper->StringToXfrmMethod( whichTransform );

    switch( xfrmMethod )
      {
      case RegistrationHelperType::Affine:
        {
        regHelper->AddAffineTransform( learningRate );
        }
        break;
      case RegistrationHelperType::Rigid:
        {
        regHelper->AddRigidTransform( learningRate );
        }
        break;
      case RegistrationHelperType::CompositeAffine:
        {
        regHelper->AddCompositeAffineTransform( learningRate );
        }
        break;
      case RegistrationHelperType::Similarity:
        {
        regHelper->AddSimilarityTransform( learningRate );
        }
        break;
      case RegistrationHelperType::Translation:
        {
        regHelper->AddTranslationTransform( learningRate );
        }
        break;
      case RegistrationHelperType::GaussianDisplacementField:
        {
        const float varianceForUpdateField = parser->Convert<float>( transformOption->GetFunction(
                                                                       currentStage )->GetParameter( 1 ) );
        const float varianceForTotalField = parser->Convert<float>( transformOption->GetFunction(
                                                                      currentStage )->GetParameter( 2 ) );
        regHelper->AddGaussianDisplacementFieldTransform(learningRate, varianceForUpdateField, varianceForTotalField);
        }
        break;
      case RegistrationHelperType::BSplineDisplacementField:
        {
        unsigned int splineOrder = 3;
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 3 )
          {
          splineOrder =
            parser->Convert<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 3 ) );
          }

        std::vector<unsigned int> meshSizeForTheUpdateField =
          parser->ConvertVector<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 1 ) );

        if( meshSizeForTheUpdateField.size() == 1 )
          {
										typename ImageType::Pointer fixedImage;
										ReadImage<ImageType>( fixedImage, fixedImageFileName.c_str() );
										fixedImage->DisconnectPipeline();

          meshSizeForTheUpdateField = regHelper->CalculateMeshSizeForSpecifiedKnotSpacing(
            fixedImage, meshSizeForTheUpdateField[0], splineOrder );
          }

        std::vector<unsigned int> meshSizeForTheTotalField;
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 2 )
          {
          meshSizeForTheTotalField =
            parser->ConvertVector<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 2 ) );
										if( meshSizeForTheTotalField.size() == 1 )
												{
												typename ImageType::Pointer fixedImage;
												ReadImage<ImageType>( fixedImage, fixedImageFileName.c_str() );
												fixedImage->DisconnectPipeline();

												meshSizeForTheTotalField = regHelper->CalculateMeshSizeForSpecifiedKnotSpacing(
														fixedImage, meshSizeForTheTotalField[0], splineOrder );
												}
          }
        else
          {
          for( unsigned int d = 0; d < VImageDimension; d++ )
            {
            meshSizeForTheTotalField.push_back( 0 );
            }
          }


        regHelper->AddBSplineDisplacementFieldTransform( learningRate, meshSizeForTheUpdateField,
                                                         meshSizeForTheTotalField,
                                                         splineOrder );
        }
        break;
      case RegistrationHelperType::BSpline:
        {
        std::vector<unsigned int> meshSizeAtBaseLevel =
          parser->ConvertVector<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 1 ) );
								if( meshSizeAtBaseLevel.size() == 1 )
										{
										typename ImageType::Pointer fixedImage;
										ReadImage<ImageType>( fixedImage, fixedImageFileName.c_str() );
										fixedImage->DisconnectPipeline();

										meshSizeAtBaseLevel = regHelper->CalculateMeshSizeForSpecifiedKnotSpacing(
												fixedImage, meshSizeAtBaseLevel[0], 3 );
										}
        regHelper->AddBSplineTransform( learningRate, meshSizeAtBaseLevel );
        }
        break;
      case RegistrationHelperType::TimeVaryingVelocityField:
        {
        unsigned int numberOfTimeIndices = parser->Convert<unsigned int>( transformOption->GetFunction(
                                                                            0 )->GetParameter( 1 ) );

        const float varianceForUpdateField = parser->Convert<float>( transformOption->GetFunction(
                                                                       currentStage )->GetParameter( 2 ) );
        const float varianceForUpdateFieldTime = parser->Convert<float>( transformOption->GetFunction(
                                                                           currentStage )->GetParameter( 3 ) );
        const float varianceForTotalField = parser->Convert<float>( transformOption->GetFunction(
                                                                      currentStage )->GetParameter( 4 ) );
        const float varianceForTotalFieldTime = parser->Convert<float>( transformOption->GetFunction(
                                                                          currentStage )->GetParameter( 5 ) );
        regHelper->AddTimeVaryingVelocityFieldTransform( learningRate,
                                                         numberOfTimeIndices,
                                                         varianceForUpdateField,
                                                         varianceForUpdateFieldTime,
                                                         varianceForTotalField,
                                                         varianceForTotalFieldTime );
        }
        break;
      case RegistrationHelperType::TimeVaryingBSplineVelocityField:
        {
        std::vector<unsigned int> meshSize = parser->ConvertVector<unsigned int>( transformOption->GetFunction(
                                                                                    0 )->GetParameter( 1 ) );

        unsigned int numberOfTimePointSamples = 4;

        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 2 )
          {
          numberOfTimePointSamples = parser->Convert<unsigned int>( transformOption->GetFunction(
                                                                      currentStage )->GetParameter( 2 ) );
          }
        unsigned int splineOrder = 3;
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 3 )
          {
          splineOrder =
            parser->Convert<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 3 ) );
          }
        regHelper->AddTimeVaryingBSplineVelocityFieldTransform( learningRate,
                                                                meshSize,
                                                                numberOfTimePointSamples,
                                                                splineOrder );
        }
        break;
      case RegistrationHelperType::SyN:
        {
        float varianceForUpdateField = 3.0;
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 1 )
          {
          varianceForUpdateField =
            parser->Convert<float>( transformOption->GetFunction( currentStage )->GetParameter( 1 ) );
          }
        float varianceForTotalField = 0.0;
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 2 )
          {
          varianceForTotalField = parser->Convert<float>( transformOption->GetFunction( currentStage )->GetParameter( 2 ) );
          }
        regHelper->AddSyNTransform( learningRate, varianceForUpdateField, varianceForTotalField );
        }
        break;
      case RegistrationHelperType::BSplineSyN:
        {
        unsigned int splineOrder = 3;
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 3 )
          {
          splineOrder =
            parser->Convert<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 3 ) );
          }

        std::vector<unsigned int> meshSizeForTheUpdateField =
          parser->ConvertVector<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 1 ) );
        if( meshSizeForTheUpdateField.size() == 1 )
          {
										typename ImageType::Pointer fixedImage;
										ReadImage<ImageType>( fixedImage, fixedImageFileName.c_str() );
										fixedImage->DisconnectPipeline();

          meshSizeForTheUpdateField = regHelper->CalculateMeshSizeForSpecifiedKnotSpacing(
            fixedImage, meshSizeForTheUpdateField[0], splineOrder );
          }

        std::vector<unsigned int> meshSizeForTheTotalField;
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 2 )
          {
          meshSizeForTheTotalField =
            parser->ConvertVector<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 2 ) );
										if( meshSizeForTheTotalField.size() == 1 )
												{
												typename ImageType::Pointer fixedImage;
												ReadImage<ImageType>( fixedImage, fixedImageFileName.c_str() );
												fixedImage->DisconnectPipeline();

												meshSizeForTheTotalField = regHelper->CalculateMeshSizeForSpecifiedKnotSpacing(
														fixedImage, meshSizeForTheTotalField[0], splineOrder );
												}
          }
        else
          {
          for( unsigned int d = 0; d < VImageDimension; d++ )
            {
            meshSizeForTheTotalField.push_back( 0 );
            }
          }

        regHelper->AddBSplineSyNTransform( learningRate, meshSizeForTheUpdateField,
                                           meshSizeForTheTotalField,
                                           splineOrder );
        }
        break;
      case RegistrationHelperType::Exponential:
        {
        const float varianceForUpdateField = parser->Convert<float>( transformOption->GetFunction(
                                                                       currentStage )->GetParameter( 1 ) );
        const float varianceForVelocityField = parser->Convert<float>( transformOption->GetFunction(
                                                                         currentStage )->GetParameter( 2 ) );
        unsigned int numberOfIntegrationSteps = 0;  // If the number of integration steps = 0, compute steps
                                                    // automatically
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 3 )
          {
          numberOfIntegrationSteps = parser->Convert<unsigned int>( transformOption->GetFunction(
                                                                      currentStage )->GetParameter( 3 ) );
          }
        regHelper->AddExponentialTransform( learningRate, varianceForUpdateField, varianceForVelocityField,
                                            numberOfIntegrationSteps );
        }
        break;
      case RegistrationHelperType::BSplineExponential:
        {
        unsigned int splineOrder = 3;
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 4 )
          {
          splineOrder =
            parser->Convert<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 4 ) );
          }

        std::vector<unsigned int> meshSizeForTheUpdateField =
          parser->ConvertVector<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 1 ) );
								if( meshSizeForTheUpdateField.size() == 1 )
										{
										typename ImageType::Pointer fixedImage;
										ReadImage<ImageType>( fixedImage, fixedImageFileName.c_str() );
										fixedImage->DisconnectPipeline();

										meshSizeForTheUpdateField = regHelper->CalculateMeshSizeForSpecifiedKnotSpacing(
												fixedImage, meshSizeForTheUpdateField[0], splineOrder );
										}

        std::vector<unsigned int> meshSizeForTheVelocityField;
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 2 )
          {
          meshSizeForTheVelocityField =
            parser->ConvertVector<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 2 ) );
										if( meshSizeForTheVelocityField.size() == 1 )
												{
												typename ImageType::Pointer fixedImage;
												ReadImage<ImageType>( fixedImage, fixedImageFileName.c_str() );
												fixedImage->DisconnectPipeline();

												meshSizeForTheVelocityField = regHelper->CalculateMeshSizeForSpecifiedKnotSpacing(
														fixedImage, meshSizeForTheVelocityField[0], splineOrder );
												}
          }
        else
          {
          for( unsigned int d = 0; d < VImageDimension; d++ )
            {
            meshSizeForTheVelocityField.push_back( 0 );
            }
          }

        unsigned int numberOfIntegrationSteps = 0;  // If the number of integration steps = 0, compute steps
                                                    // automatically
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 3 )
          {
          numberOfIntegrationSteps = parser->Convert<unsigned int>( transformOption->GetFunction(
                                                                      currentStage )->GetParameter( 3 ) );
          }

        regHelper->AddBSplineExponentialTransform( learningRate, meshSizeForTheUpdateField,
                                                   meshSizeForTheVelocityField,
                                                   numberOfIntegrationSteps, splineOrder );
        }
        break;
      default:
        {
        std::cout << "Unknown registration method " << "\"" << whichTransform << "\"" << std::endl;
        }
        break;
      }
    }

  // set the vector-vector parameters accumulated
  regHelper->SetIterations( iterationList );
  regHelper->SetConvergenceWindowSizes( convergenceWindowSizeList );
  regHelper->SetConvergenceThresholds( convergenceThresholdList );
  regHelper->SetSmoothingSigmas( smoothingSigmasList );
  regHelper->SetSmoothingSigmasAreInPhysicalUnits( smoothingSigmasAreInPhysicalUnitsList );
  regHelper->SetShrinkFactors( shrinkFactorsList );

  // We iterate through each of the metric "functions" specified on the command
  // line and add it the registration helper.  We also need to assign the stage
  // ID to the added metric.  Multiple metrics for a single stage are specified
  // on the command line by being specified adjacently.

  unsigned int numberOfMetrics = metricOption->GetNumberOfFunctions();
  for( int currentMetricNumber = numberOfMetrics - 1; currentMetricNumber >= 0; currentMetricNumber-- )
    {
    // Get the fixed and moving images

    std::string fixedImageFileName = metricOption->GetFunction( currentMetricNumber )->GetParameter( 0 );
    std::string movingImageFileName = metricOption->GetFunction( currentMetricNumber )->GetParameter( 1 );
    std::cout << "  fixed image: " << fixedImageFileName << std::endl;
    std::cout << "  moving image: " << movingImageFileName << std::endl;

    typename ImageType::Pointer fixedImage;
    typename ImageType::Pointer movingImage;
    ReadImage<ImageType>( fixedImage,  fixedImageFileName.c_str() );
    ReadImage<ImageType>( movingImage, movingImageFileName.c_str() );
    fixedImage->DisconnectPipeline();
    movingImage->DisconnectPipeline();

    // Get the stage ID
    unsigned int stageID = metricOption->GetFunction( currentMetricNumber )->GetStageID();

    // We check the last stage ID (first iteration) to ensure that the number of stages
    // (as determined by the number of transforms) is equal to the number of stages (as
    // determined by the metrics command line specification).
    if( currentMetricNumber == static_cast<int>( numberOfMetrics - 1 ) )
      {
      if( stageID != numberOfTransforms - 1 )
        {
        std::cout << "\n\n\n"
                         << "Error:  The number of stages does not match up with the metrics." << std::endl
                         << "The number of transforms is " << numberOfTransforms << " and the last stage ID "
                         << " as determined by the metrics is " << stageID << "." << std::endl;
        return EXIT_FAILURE;
        }
      }

    std::string whichMetric = metricOption->GetFunction( currentMetricNumber )->GetName();
    ConvertToLowerCase( whichMetric );

    float metricWeighting = 1.0;
    if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 2 )
      {
      metricWeighting = parser->Convert<float>( metricOption->GetFunction( currentMetricNumber )->GetParameter( 2 ) );
      }

    float samplingPercentage = 1.0;
    if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 5 )
      {
      samplingPercentage =
        parser->Convert<float>( metricOption->GetFunction( currentMetricNumber )->GetParameter( 5 ) );
      }

    std::string strategy = "none";
    if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 4 )
      {
      strategy = metricOption->GetFunction( currentMetricNumber )->GetParameter( 4 );
      }
    ConvertToLowerCase( strategy );

    typename RegistrationHelperType::SamplingStrategy samplingStrategy = RegistrationHelperType::invalid;
    if( strategy == "random" )
      {
      samplingStrategy = RegistrationHelperType::random;
      }
    else if( strategy == "regular" )
      {
      samplingStrategy = RegistrationHelperType::regular;
      }
    else if( ( strategy == "none" ) || ( strategy == "" ) )
      {
      samplingStrategy = RegistrationHelperType::none;
      }
    else
      {
      samplingStrategy = RegistrationHelperType::invalid;
      std::cout << "ERROR: invalid sampling strategy specified: " << strategy << std::endl;
      return EXIT_FAILURE;
      }

    typename RegistrationHelperType::MetricEnumeration curMetric = regHelper->StringToMetricType( whichMetric );

    switch( curMetric )
      {
      case RegistrationHelperType::CC:
        {
        unsigned int radiusOption = parser->Convert<unsigned int>( metricOption->GetFunction(
                                                                     currentMetricNumber )->GetParameter( 3 ) );
        regHelper->AddMetric( curMetric,
                              fixedImage,
                              movingImage,
                              stageID,
                              metricWeighting,
                              samplingStrategy,
                              1,
                              radiusOption,
                              samplingPercentage );
        }
        break;
      case RegistrationHelperType::GC:
      case RegistrationHelperType::MeanSquares:
      case RegistrationHelperType::Demons:
        {
        regHelper->AddMetric( curMetric,
                              fixedImage,
                              movingImage,
                              stageID,
                              metricWeighting,
                              samplingStrategy,
                              1,
                              1,
                              samplingPercentage );
        }
        break;
      case RegistrationHelperType::Mattes:
      case RegistrationHelperType::MI:
        {
        unsigned int binOption = parser->Convert<unsigned int>( metricOption->GetFunction(
                                                                  currentMetricNumber )->GetParameter( 3 ) );
        regHelper->AddMetric( curMetric,
                              fixedImage,
                              movingImage,
                              stageID,
                              metricWeighting,
                              samplingStrategy,
                              binOption,
                              1,
                              samplingPercentage );
        }
        break;
      default:
        std::cout << "ERROR: Unrecognized image metric: " << whichMetric << std::endl;
        return EXIT_FAILURE;
      }
    }

  // Perform the registration

  if( regHelper->DoRegistration() == EXIT_FAILURE )
    {
    return EXIT_FAILURE;
    }
  //
  // write out transforms stored in the composite
  typename CompositeTransformType::Pointer resultTransform = regHelper->GetModifiableCompositeTransform();

  if( parser->Convert<bool>( compositeOutputOption->GetFunction( 0 )->GetName() ) )
    {
    std::string compositeTransformFileName = outputPrefix + std::string( "Composite.h5" );
    std::string inverseCompositeTransformFileName = outputPrefix + std::string( "InverseComposite.h5" );

    typename RegistrationHelperType::CompositeTransformType::TransformTypePointer compositeTransform =
      resultTransform.GetPointer();
    itk::ants::WriteTransform<TComputeType, VImageDimension>( compositeTransform, compositeTransformFileName.c_str() );

    typename RegistrationHelperType::CompositeTransformType::TransformTypePointer inverseCompositeTransform =
      compositeTransform->GetInverseTransform();
    if( inverseCompositeTransform.IsNotNull() )
      {
      itk::ants::WriteTransform<TComputeType, VImageDimension>( inverseCompositeTransform,
                                                  inverseCompositeTransformFileName.c_str() );
      }
    }
  unsigned int numTransforms = resultTransform->GetNumberOfTransforms();

  // write out transforms actually computed, so skip any initial transforms unless
  // we're collapsing the output transforms.

  typedef typename RegistrationHelperType::CompositeTransformType         CompositeTransformType;
  typedef typename CompositeTransformType::Pointer                        CompositeTransformPointer;
  typedef typename RegistrationHelperType::DisplacementFieldTransformType DisplacementFieldTransformType;
  typedef typename RegistrationHelperType::TransformType                  TransformType;

  unsigned int              startIndex = initialMovingTransformOption->GetNumberOfFunctions();
  CompositeTransformPointer collapsedResultTransform;
  if( parser->Convert<bool>( collapseOutputTransformsOption->GetFunction( 0 )->GetName() ) )
    {
    collapsedResultTransform = regHelper->CollapseCompositeTransform( resultTransform );
    numTransforms = collapsedResultTransform->GetNumberOfTransforms();
    startIndex = 0;
    TransformTypeNames.clear();
    for( unsigned int i = 0; i < numTransforms; i++ )
      {
      if( collapsedResultTransform->GetNthTransform( i )->GetTransformCategory() == TransformType::Linear )
        {
        // The output type must be Affine, not matrixoffset!  TransformTypeNames.push_back( "matrixoffset" );
        TransformTypeNames.push_back( "genericaffine" );
        }
      else if( collapsedResultTransform->GetNthTransform( i )->GetTransformCategory() ==
               TransformType::DisplacementField )
        {
        typename DisplacementFieldTransformType::Pointer nthTransform =
          dynamic_cast<DisplacementFieldTransformType *>( collapsedResultTransform->GetNthTransform( i ).GetPointer() );

        // We don't know what set of displacement field transforms were optimized.
        // All we know is whether or not an inverse displacement field exists.  If so,
        // we simply pass a transform name which either does have an inverse or does
        // not.
        if( nthTransform && nthTransform->GetInverseDisplacementField() )
          {
          TransformTypeNames.push_back( "syn" );
          }
        else
          {
          TransformTypeNames.push_back( "gdf" );
          }
        }
      else if( collapsedResultTransform->GetNthTransform( i )->GetTransformCategory() == TransformType::BSpline )
        {
        TransformTypeNames.push_back( "bspline" );
        }
      }
    }
  for( unsigned int i = startIndex; i < numTransforms; ++i )
    {
    typename CompositeTransformType::TransformTypePointer curTransform;
    if( parser->Convert<bool>( collapseOutputTransformsOption->GetFunction( 0 )->GetName() ) )
      {
      curTransform = collapsedResultTransform->GetNthTransform( i );
      }
    else
      {
      curTransform = resultTransform->GetNthTransform( i );
      }

    //
    // only registrations not part of the initial transforms in the
    // TransformTypeNames list.
    const std::string curTransformType = TransformTypeNames.front();
    TransformTypeNames.pop_front();

    bool writeInverse;
    bool writeVelocityField;

    std::string transformTemplateName = RegTypeToFileName( curTransformType, writeInverse, writeVelocityField );

    std::stringstream curFileName;
    curFileName << outputPrefix << i << transformTemplateName;

    // WriteTransform will spit all sorts of error messages if it
    // fails, and we want to keep going even if it does so ignore its
    // return value.
    itk::ants::WriteTransform<TComputeType, VImageDimension>( curTransform, curFileName.str() );

    typedef typename RegistrationHelperType::DisplacementFieldTransformType DisplacementFieldTransformType;
    typedef typename DisplacementFieldTransformType::DisplacementFieldType  DisplacementFieldType;
    typename DisplacementFieldTransformType::Pointer dispTransform =
      dynamic_cast<DisplacementFieldTransformType *>(curTransform.GetPointer() );
    // write inverse transform file
    if( writeInverse && dispTransform.IsNotNull() )
      {
      typename DisplacementFieldType::ConstPointer inverseDispField = dispTransform->GetInverseDisplacementField();
      if( inverseDispField.IsNotNull() )
        {
        std::stringstream curInverseFileName;
        curInverseFileName << outputPrefix << i << "InverseWarp.nii.gz";
        typedef itk::ImageFileWriter<DisplacementFieldType> InverseWriterType;
        typename InverseWriterType::Pointer inverseWriter = InverseWriterType::New();
        inverseWriter->SetInput( dispTransform->GetInverseDisplacementField() );
        inverseWriter->SetFileName( curInverseFileName.str().c_str() );
        try
          {
          inverseWriter->Update();
          }
        catch( itk::ExceptionObject & err )
          {
          std::cout << "Can't write transform file " << curInverseFileName.str().c_str() << std::endl;
          std::cout << "Exception Object caught: " << std::endl;
          std::cout << err << std::endl;
          }
        }
      }
    if( writeVelocityField )
      {
      // write velocity field (if applicable)
      typedef typename RegistrationHelperType::TimeVaryingVelocityFieldTransformType
        VelocityFieldTransformType;

      typedef itk::Image<itk::Vector<TComputeType, VImageDimension>, VImageDimension + 1> VelocityFieldType;
      typename VelocityFieldTransformType::Pointer velocityFieldTransform =
        dynamic_cast<VelocityFieldTransformType *>(curTransform.GetPointer() );
      if( !velocityFieldTransform.IsNull() )
        {
        std::stringstream curVelocityFieldFileName;
        curVelocityFieldFileName << outputPrefix << i << "VelocityField.nii.gz";

        typedef itk::ImageFileWriter<VelocityFieldType> VelocityFieldWriterType;
        typename VelocityFieldWriterType::Pointer velocityFieldWriter = VelocityFieldWriterType::New();
        velocityFieldWriter->SetInput( velocityFieldTransform->GetTimeVaryingVelocityField() );
        velocityFieldWriter->SetFileName( curVelocityFieldFileName.str().c_str() );
        try
          {
          velocityFieldWriter->Update();
          }
        catch( itk::ExceptionObject & err )
          {
          std::cout << "Can't write velocity field transform file " << curVelocityFieldFileName.str().c_str()
                   << std::endl;
          std::cout << "Exception Object caught: " << std::endl;
          std::cout << err << std::endl;
          }
        }
      }
    }

  std::string whichInterpolator( "linear" );
  typename itk::ants::CommandLineParser::OptionType::Pointer interpolationOption = parser->GetOption( "interpolation" );
  if( interpolationOption && interpolationOption->GetNumberOfFunctions() )
    {
    whichInterpolator = interpolationOption->GetFunction( 0 )->GetName();
    ConvertToLowerCase( whichInterpolator );
    }

  typename ImageType::SpacingType cache_spacing_for_smoothing_sigmas
    (itk::NumericTraits<typename ImageType::SpacingType::ValueType>::Zero);
  if( !std::strcmp( whichInterpolator.c_str(), "gaussian" )
      ||   !std::strcmp( whichInterpolator.c_str(), "multilabel" )
      )
    {
#if 1
    // HACK:: This can just be cached when reading the fixedImage from above!!
    //
    std::string fixedImageFileName = metricOption->GetFunction( numberOfTransforms - 1 )->GetParameter( 0 );

    typedef itk::ImageFileReader<ImageType> ImageReaderType;
    typename ImageReaderType::Pointer fixedImageReader = ImageReaderType::New();

    fixedImageReader->SetFileName( fixedImageFileName.c_str() );
    fixedImageReader->Update();
    typename ImageType::Pointer fixedImage = fixedImageReader->GetOutput();
#endif
    cache_spacing_for_smoothing_sigmas = fixedImage->GetSpacing();
    }

#include "make_interpolator_snip.tmpl"
  regHelper->SetInterpolator( interpolator );

  typename ImageType::Pointer warpedImage = regHelper->GetWarpedImage();
  typedef itk::ImageFileWriter<ImageType> WarpedImageWriterType;
  if( !outputWarpedImageName.empty() )
    {
    WriteImage<ImageType>( warpedImage, outputWarpedImageName.c_str()  );
    }

  if( !outputInverseWarpedImageName.empty() )
    {
    typename ImageType::Pointer inverseWarpedImage = regHelper->GetInverseWarpedImage();
    if( inverseWarpedImage.IsNotNull() )
      {
      WriteImage<ImageType>( inverseWarpedImage, outputInverseWarpedImageName.c_str()  );
      }
    }

  return EXIT_SUCCESS;
}

extern int antsRegistration2DDouble(ParserType::Pointer & parser);
extern int antsRegistration3DDouble(ParserType::Pointer & parser);
extern int antsRegistration2DFloat(ParserType::Pointer & parser);
extern int antsRegistration3DFloat(ParserType::Pointer & parser);

} //End namespace

#endif // __ANTSREGISTRATIONTEMPLATEHEADER_H__
