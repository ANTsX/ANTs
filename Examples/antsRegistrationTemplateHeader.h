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
#include "itkLabelImageGenericInterpolateImageFunction.h"
#include "include/antsRegistration.h"
#include "ReadWriteData.h"

namespace ants
{

extern const char * RegTypeToFileName( const std::string & type, bool & writeInverse, bool & writeVelocityField, bool minc );

template <typename TComputeType, unsigned VImageDimension>
int
DoRegistration(typename ParserType::Pointer & parser)
{
  typedef TComputeType                                                     RealType;
  typedef typename ants::RegistrationHelper<TComputeType, VImageDimension> RegistrationHelperType;
  typedef typename RegistrationHelperType::ImageType                       ImageType;
  typedef typename RegistrationHelperType::MaskImageType                   MaskImageType;
  typedef typename RegistrationHelperType::LabeledPointSetType             LabeledPointSetType;
  typedef typename RegistrationHelperType::IntensityPointSetType           IntensityPointSetType;
  typedef typename RegistrationHelperType::CompositeTransformType          CompositeTransformType;

  typename RegistrationHelperType::Pointer regHelper = RegistrationHelperType::New();

  OptionType::Pointer useMincFormatOption = parser->GetOption( "minc" );
  const bool useMincFormat = parser->Convert<bool>( useMincFormatOption->GetFunction( 0 )->GetName() );

  bool verbose = false;
  typename itk::ants::CommandLineParser::OptionType::Pointer verboseOption =
    parser->GetOption( "verbose" );
  if( verboseOption && verboseOption->GetNumberOfFunctions() )
    {
    verbose = parser->Convert<bool>( verboseOption->GetFunction( 0 )->GetName() );
    }

  nullStream cnul;
  if( ! verbose )
    {
    regHelper->SetLogStream( cnul );
    }

  OptionType::Pointer fixRandomSeed = parser->GetOption( "random-seed" );
  if( fixRandomSeed && fixRandomSeed->GetNumberOfFunctions() )
    {
    int randomSeed = parser->Convert<int>( fixRandomSeed->GetFunction(0)->GetName() );
    regHelper->SetRegistrationRandomSeed(randomSeed);
    }
  else
    {
    char* randomSeedEnv = getenv( "ANTS_RANDOM_SEED" );
    if ( randomSeedEnv != nullptr )
      {
      regHelper->SetRegistrationRandomSeed( std::stoi( randomSeedEnv ) );
      }
    else
      {
      regHelper->SetRegistrationRandomSeed(0);
      }
    }

  OptionType::Pointer transformOption = parser->GetOption( "transform" );
  if( !transformOption || transformOption->GetNumberOfFunctions() == 0 )
    {
    if( verbose )
      {
      std::cerr << "ERROR: the transform option ('-t') must be specified.  See help menu." << std::endl;
      }
    return EXIT_FAILURE;
    }

  OptionType::Pointer metricOption = parser->GetOption( "metric" );
  if( !metricOption || metricOption->GetNumberOfFunctions() == 0 )
    {
    if( verbose )
      {
      std::cerr << "ERROR: the metric option ('-m') must be specified.  See help menu." << std::endl;
      }
    return EXIT_FAILURE;
    }

  OptionType::Pointer convergenceOption = parser->GetOption( "convergence" );
  if( !convergenceOption || convergenceOption->GetNumberOfFunctions() == 0 )
    {
    if( verbose )
      {
      std::cerr << "ERROR: the convergence option ('-c') must be specified.  See help menu." << std::endl;
      }
    return EXIT_FAILURE;
    }

  OptionType::Pointer shrinkFactorsOption = parser->GetOption( "shrink-factors" );
  if( !shrinkFactorsOption || shrinkFactorsOption->GetNumberOfFunctions() == 0 )
    {
    if( verbose )
      {
      std::cerr << "ERROR: the shrink factors option ('-f') must be specified.  See help menu." << std::endl;
      }
    return EXIT_FAILURE;
    }

  OptionType::Pointer smoothingSigmasOption = parser->GetOption( "smoothing-sigmas" );
  if( !smoothingSigmasOption || smoothingSigmasOption->GetNumberOfFunctions() == 0 )
    {
    if( verbose )
      {
      std::cerr << "ERROR: the smoothing sigmas option ('-s') must be specified.  See help menu." << std::endl;
      }
    return EXIT_FAILURE;
    }

  OptionType::Pointer restrictDeformationOption = parser->GetOption( "restrict-deformation" );

  OptionType::Pointer outputOption = parser->GetOption( "output" );
  if( !outputOption || outputOption->GetNumberOfFunctions() == 0 )
    {
    if( verbose )
      {
      std::cerr << "ERROR: the output option ('-o') must be specified.  See help menu." << std::endl;
      }
    return EXIT_FAILURE;
    }

  OptionType::Pointer maskOption = parser->GetOption( "masks" );

  OptionType::Pointer compositeOutputOption = parser->GetOption( "write-composite-transform" );
  const bool writeCompositeTransform = parser->Convert<bool>( compositeOutputOption->GetFunction( 0 )->GetName() );

  OptionType::Pointer saveStateOption = parser->GetOption( "save-state" );

  OptionType::Pointer collapseOutputTransformsOption = parser->GetOption( "collapse-output-transforms" );
  const bool shouldCollapseBeDone = parser->Convert<bool>( collapseOutputTransformsOption->GetFunction( 0 )->GetName() );

  OptionType::Pointer initializeTransformsPerStageOption = parser->GetOption( "initialize-transforms-per-stage" );
  if( initializeTransformsPerStageOption && parser->Convert<bool>( initializeTransformsPerStageOption->GetFunction( 0 )->GetName() ) )
    {
    if( shouldCollapseBeDone )
      {
      if( verbose )
        {
        std::cerr << "ERROR: initialize-transforms-per-stage & collapse-output-transforms options are mutually exclusive."
                  << std::endl;
        }
      return EXIT_FAILURE;
      }
    regHelper->SetInitializeTransformsPerStage( true );
    }
  else
    {
    regHelper->SetInitializeTransformsPerStage( false );
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
    if( !shouldCollapseBeDone )
      {
      for( unsigned int n = 0; n < isDerivedInitialMovingTransform.size(); n++ )
        {
        std::stringstream currentFileName;
        if( useMincFormat )
          {
          currentFileName << outputPrefix << n << "DerivedInitialMovingTranslation.xfm";
          }
        else
          {
          currentFileName << outputPrefix << n << "DerivedInitialMovingTranslation.mat";
          }

        typename RegistrationHelperType::CompositeTransformType::TransformTypePointer currentTransform =
          compositeTransform->GetNthTransform( n );
        if( currentTransform->IsLinear() && isDerivedInitialMovingTransform[n] )
          {
          itk::ants::WriteTransform<TComputeType, VImageDimension>( currentTransform, currentFileName.str() );
          }
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
    if( !shouldCollapseBeDone )
      {
      for( unsigned int n = 0; n < isDerivedInitialFixedTransform.size(); n++ )
        {
        std::stringstream currentFileName;
        if( useMincFormat )
          {
          currentFileName << outputPrefix << n << "DerivedInitialFixedTranslation.xfm";
          }
        else
          {
          currentFileName << outputPrefix << n << "DerivedInitialFixedTranslation.mat";
          }

        typename RegistrationHelperType::CompositeTransformType::TransformTypePointer currentTransform =
          compositeTransform->GetNthTransform( n );
        if( currentTransform->IsLinear() && isDerivedInitialFixedTransform[n] )
          {
          itk::ants::WriteTransform<TComputeType, VImageDimension>( currentTransform, currentFileName.str() );
          }
        }
      }
    }

  ParserType::OptionType::Pointer restoreStateOption = parser->GetOption( "restore-state" );

  if( restoreStateOption && restoreStateOption->GetNumberOfFunctions() )
    {
    if( initialMovingTransformOption->GetNumberOfFunctions() || initialFixedTransformOption->GetNumberOfFunctions() )
      {
      if( verbose )
        {
        std::cerr << "restore-state option is mutually exclusive with "
                  << "initial-moving-transform & initial-fixed-transform options." << std::endl;
        }
      return EXIT_FAILURE;
      }

    if( verbose )
      {
      std::cout << "Restoring previous registration state" << std::endl;
      }
    std::vector<bool> isDerivedInitialMovingTransform;
    typename CompositeTransformType::Pointer compositeTransform =
      GetCompositeTransformFromParserOption<TComputeType, VImageDimension>( parser, restoreStateOption,
                                                                           isDerivedInitialMovingTransform );
    if( verbose )
      {
      std::cout << "+" << std::endl;
      }
    if( compositeTransform.IsNull() )
      {
      return EXIT_FAILURE;
      }
    regHelper->SetRestoreStateTransform( compositeTransform );
    if( verbose )
      {
      std::cout << "+" << std::endl;
      }
    }

  if( maskOption && maskOption->GetNumberOfFunctions() )
    {
    if( verbose )
      {
      std::cout << "  Reading mask(s)." << std::endl;
      }
    for( int l = maskOption->GetNumberOfFunctions() - 1; l >= 0; l-- )
      {
      if( verbose )
        {
        std::cout << "    Registration stage " << ( maskOption->GetNumberOfFunctions() - l - 1 ) << std::endl;
        }
      if( maskOption->GetFunction( l )->GetNumberOfParameters() > 0 )
        {
        for( unsigned m = 0; m < maskOption->GetFunction( l )->GetNumberOfParameters(); m++ )
          {
          std::string fname = maskOption->GetFunction( l )->GetParameter( m );
          typename MaskImageType::Pointer maskImage;
          ReadImage<MaskImageType>( maskImage, fname.c_str() );
          if( m == 0 )
            {
            regHelper->AddFixedImageMask( maskImage );
            if( verbose )
              {
              if( maskImage.IsNotNull() )
                {
                std::cout << "      Fixed mask = " << fname.c_str() << std::endl;
                }
              else
                {
                std::cout << "      No fixed mask" << std::endl;
                }
              }
            }
          else if( m == 1 )
            {
            regHelper->AddMovingImageMask( maskImage );
            if( verbose )
              {
              if( maskImage.IsNotNull() )
                {
                std::cout << "      Moving mask = " << fname << std::endl;
                }
              else
                {
                std::cout << "      No moving mask" << std::endl;
                }
              }
            }
          }
        }
      else
        {
        std::string fname = maskOption->GetFunction( l )->GetName();
        typename MaskImageType::Pointer maskImage;
        ReadImage<MaskImageType>( maskImage, fname.c_str() );
        regHelper->AddFixedImageMask( maskImage );
        if( verbose )
          {
          if( maskImage.IsNotNull() )
            {
            std::cout << "      Fixed mask = " << fname << std::endl;
            }
          else
            {
            std::cout << "      No fixed mask" << std::endl;
            }
          }
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

  // We find both the number of transforms and the number of metrics

  unsigned int numberOfTransforms = transformOption->GetNumberOfFunctions();
  if( transformOption.IsNull() || numberOfTransforms == 0 )
    {
    if( verbose )
      {
      std::cerr << "No transformations are specified." << std::endl;
      }
    return EXIT_FAILURE;
    }

  std::vector<std::vector<unsigned int> > iterationList;
  std::vector<std::vector<RealType> >     restrictDeformationWeightsList;
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
        constexpr unsigned int minAllowedconvergenceWindowSize = 2; // The BSplineScatteredDataPoints requires at least 2
                                                                // points for interpolation.
        if( convergenceWindowSize < minAllowedconvergenceWindowSize )
          {
          if( verbose )
            {
            std::cerr << "Convergence Window Size must be greater than or equal to " << minAllowedconvergenceWindowSize << std::endl;
            }
          return EXIT_FAILURE;
          }
        }
      }
    else
      {
      if( verbose )
        {
        std::cerr << "No convergence criteria are specified." << std::endl;
        }
      return EXIT_FAILURE;
      }

    iterationList.push_back( iterations );
    convergenceThresholdList.push_back( convergenceThreshold );
    convergenceWindowSizeList.push_back( convergenceWindowSize );

    if( restrictDeformationOption.IsNotNull() && restrictDeformationOption->GetNumberOfFunctions() > static_cast<unsigned int>( currentStage ) )
      {
      std::vector<RealType> restrictDeformationWeights =
        parser->ConvertVector<RealType>( restrictDeformationOption->GetFunction( currentStage )->GetName() );
      restrictDeformationWeightsList.push_back( restrictDeformationWeights );
      }

    unsigned int numberOfLevels = iterations.size();
    if( verbose )
      {
      std::cout << "  number of levels = " << numberOfLevels << std::endl;
      }

    // Get the first metricOption for the currentStage (for use with the B-spline transforms)
    unsigned int numberOfMetrics = metricOption->GetNumberOfFunctions();
    std::string  fixedImageFileName;
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
          varianceForTotalField = parser->Convert<float>( transformOption->GetFunction( currentStage )->GetParameter(
                                                            2 ) );
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

        std::vector<unsigned int> meshSizeForTheUpdateField;

        std::vector<float> meshSizeForTheUpdateFieldFloat =
          parser->ConvertVector<float>( transformOption->GetFunction( currentStage )->GetParameter( 1 ) );
        if( meshSizeForTheUpdateFieldFloat.size() == 1 )
          {
          typename ImageType::Pointer fixedImage;
          ReadImage<ImageType>( fixedImage, fixedImageFileName.c_str() );
          fixedImage->DisconnectPipeline();

          meshSizeForTheUpdateField = regHelper->CalculateMeshSizeForSpecifiedKnotSpacing(
              fixedImage, meshSizeForTheUpdateFieldFloat[0], splineOrder );
          }
        else
          {
          meshSizeForTheUpdateField =
            parser->ConvertVector<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 1 ) );
          }

        std::vector<unsigned int> meshSizeForTheTotalField;
        if( transformOption->GetFunction( currentStage )->GetNumberOfParameters() > 2 )
          {
          std::vector<float> meshSizeForTheTotalFieldFloat =
            parser->ConvertVector<float>( transformOption->GetFunction( currentStage )->GetParameter( 2 ) );
          if( meshSizeForTheTotalFieldFloat.size() == 1 )
            {
            typename ImageType::Pointer fixedImage;
            ReadImage<ImageType>( fixedImage, fixedImageFileName.c_str() );
            fixedImage->DisconnectPipeline();

            meshSizeForTheTotalField = regHelper->CalculateMeshSizeForSpecifiedKnotSpacing(
                fixedImage, meshSizeForTheTotalFieldFloat[0], splineOrder );
            }
          else
            {
            meshSizeForTheTotalField =
              parser->ConvertVector<unsigned int>( transformOption->GetFunction( currentStage )->GetParameter( 2 ) );
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
        if( verbose )
          {
          std::cerr << "Unknown registration method " << "\"" << whichTransform << "\"" << std::endl;
          }
        return EXIT_FAILURE;
        }
        break;
      }
    }

  // set the vector-vector parameters accumulated
  regHelper->SetIterations( iterationList );
  regHelper->SetRestrictDeformationOptimizerWeights( restrictDeformationWeightsList );
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
    // Get the stage ID
    unsigned int stageID = metricOption->GetFunction( currentMetricNumber )->GetStageID();

    // We check the last stage ID (first iteration) to ensure that the number of stages
    // (as determined by the number of transforms) is equal to the number of stages (as
    // determined by the metrics command line specification).
    if( currentMetricNumber == static_cast<int>( numberOfMetrics - 1 ) )
      {
      if( stageID != numberOfTransforms - 1 )
        {
        if( verbose )
          {
          std::cerr << "\n\n\n"
                    << "Error:  The number of stages does not match up with the metrics." << std::endl
                    << "The number of transforms is " << numberOfTransforms << " and the last stage ID "
                    << " as determined by the metrics is " << stageID << "." << std::endl;
          }
        return EXIT_FAILURE;
        }
      }

    std::string whichMetric = metricOption->GetFunction( currentMetricNumber )->GetName();
    ConvertToLowerCase( whichMetric );
    typename RegistrationHelperType::MetricEnumeration currentMetric = regHelper->StringToMetricType( whichMetric );

    // Get the fixed and moving images or point sets

    typename ImageType::Pointer fixedImage = nullptr;
    typename ImageType::Pointer movingImage = nullptr;
    typename LabeledPointSetType::Pointer fixedLabeledPointSet = nullptr;
    typename LabeledPointSetType::Pointer movingLabeledPointSet = nullptr;
    typename IntensityPointSetType::Pointer fixedIntensityPointSet = nullptr;
    typename IntensityPointSetType::Pointer movingIntensityPointSet = nullptr;

    float metricWeighting = 1.0;
    if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 2 )
      {
      metricWeighting = parser->Convert<float>( metricOption->GetFunction( currentMetricNumber )->GetParameter( 2 ) );
      }

    // assign default image metric variables
    typename RegistrationHelperType::SamplingStrategy samplingStrategy = RegistrationHelperType::none;
    unsigned int numberOfBins = 32;
    unsigned int radius = 4;

    // assign default point-set variables

    //   labeled point sets
    bool useBoundaryPointsOnly = false;
    RealType pointSetSigma = itk::NumericTraits<RealType>::OneValue();
    unsigned int evaluationKNeighborhood = 50;
    RealType alpha = static_cast<RealType>( 1.1 );
    bool useAnisotropicCovariances = false;
    RealType samplingPercentage = itk::NumericTraits<RealType>::OneValue();
    //   intensity point sets
    RealType intensityDistanceSigma = itk::NumericTraits<RealType>::ZeroValue();
    RealType euclideanDistanceSigma = itk::NumericTraits<RealType>::ZeroValue();

    if( !regHelper->IsPointSetMetric( currentMetric ) )
      {
      if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 5 )
        {
        samplingPercentage =
          parser->Convert<float>( metricOption->GetFunction( currentMetricNumber )->GetParameter( 5 ) );
        }
      std::string fixedFileName = metricOption->GetFunction( currentMetricNumber )->GetParameter( 0 );
      std::string movingFileName = metricOption->GetFunction( currentMetricNumber )->GetParameter( 1 );

      if( verbose )
        {
        std::cout << "  fixed image: " << fixedFileName << std::endl;
        std::cout << "  moving image: " << movingFileName << std::endl;
        }

      ReadImage<ImageType>( fixedImage, fixedFileName.c_str() );
      ReadImage<ImageType>( movingImage, movingFileName.c_str() );
      fixedImage->DisconnectPipeline();
      movingImage->DisconnectPipeline();

      std::string strategy = "none";
      if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 4 )
        {
        strategy = metricOption->GetFunction( currentMetricNumber )->GetParameter( 4 );
        }
      ConvertToLowerCase( strategy );

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
        if( verbose )
          {
          std::cerr << "ERROR: invalid sampling strategy specified: " << strategy << std::endl;
          }
        return EXIT_FAILURE;
        }
      if( currentMetric == RegistrationHelperType::CC )
        {
        if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 3 )
          {
          radius = parser->Convert<unsigned int>( metricOption->GetFunction( currentMetricNumber )->GetParameter( 3 ) );
          }
        }
      else if( currentMetric == RegistrationHelperType::Mattes || currentMetric == RegistrationHelperType::MI )
        {
        if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 3 )
          {
          numberOfBins = parser->Convert<unsigned int>( metricOption->GetFunction(
                                                                  currentMetricNumber )->GetParameter( 3 ) );
          }
        }
      }
    else
      {
      if( whichMetric == "igdm" )
        {
        if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() < 5 )
          {
          if( verbose )
            {
            std::cerr << "The expected number of parameters aren't specified.  Please see help menu." << std::endl;
            }
          return EXIT_FAILURE;
          }

        std::string fixedFileName = metricOption->GetFunction( currentMetricNumber )->GetParameter( 0 );
        std::string movingFileName = metricOption->GetFunction( currentMetricNumber )->GetParameter( 1 );

        if( verbose )
          {
          std::cout << "  fixed intensity point set: " << fixedFileName << std::endl;
          std::cout << "  moving intensity point set: " << movingFileName << std::endl;
          }

        std::string fixedPointSetMaskFile = metricOption->GetFunction( currentMetricNumber )->GetParameter( 3 );
        std::string movingPointSetMaskFile = metricOption->GetFunction( currentMetricNumber )->GetParameter( 4 );

        if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 6 )
          {
          intensityDistanceSigma =
            parser->Convert<RealType>( metricOption->GetFunction( currentMetricNumber )->GetParameter( 6 ) );
          }
        if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 7 )
          {
          euclideanDistanceSigma =
            parser->Convert<RealType>( metricOption->GetFunction( currentMetricNumber )->GetParameter( 7 ) );
          }
        evaluationKNeighborhood = 1;
        if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 8 )
          {
          evaluationKNeighborhood =
            parser->Convert<unsigned int>( metricOption->GetFunction( currentMetricNumber )->GetParameter( 8 ) );
          }

        RealType gradientPointSetSigma = itk::NumericTraits<RealType>::OneValue();
        if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 9 )
          {
          gradientPointSetSigma = parser->Convert<RealType>(
            metricOption->GetFunction( currentMetricNumber )->GetParameter( 9 ) );
          }
        std::vector<unsigned int> neighborhoodRadius;
        if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 5 )
          {
          neighborhoodRadius = parser->ConvertVector<unsigned int>(
            metricOption->GetFunction( currentMetricNumber )->GetParameter( 5 ) );
          }
        else
          {
          for( unsigned d = 0; d < VImageDimension; d++ )
            {
            neighborhoodRadius.push_back( 0 );
            }
          }

        if( neighborhoodRadius.size() != VImageDimension )
          {
          if( verbose )
            {
            std::cerr << "The neighborhood size must equal the dimension." << std::endl;
            }
          return EXIT_FAILURE;
          }

        ReadImageIntensityPointSet<ImageType, MaskImageType, IntensityPointSetType>(
          fixedIntensityPointSet, fixedFileName.c_str(), fixedPointSetMaskFile.c_str(),
          neighborhoodRadius, gradientPointSetSigma );

        ReadImageIntensityPointSet<ImageType, MaskImageType, IntensityPointSetType>(
          movingIntensityPointSet, movingFileName.c_str(), movingPointSetMaskFile.c_str(),
          neighborhoodRadius, gradientPointSetSigma );

        fixedIntensityPointSet->DisconnectPipeline();
        movingIntensityPointSet->DisconnectPipeline();
        }
      else
        {
        if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 3 )
          {
          samplingPercentage =
            parser->Convert<RealType>( metricOption->GetFunction( currentMetricNumber )->GetParameter( 3 ) );
          }
        if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 4 )
          {
          useBoundaryPointsOnly =
            parser->Convert<bool>( metricOption->GetFunction( currentMetricNumber )->GetParameter( 4 ) );
          }
        if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 5 )
          {
          pointSetSigma =
            parser->Convert<RealType>( metricOption->GetFunction( currentMetricNumber )->GetParameter( 5 ) );
          }
        if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 6 )
          {
          evaluationKNeighborhood =
            parser->Convert<unsigned int>( metricOption->GetFunction( currentMetricNumber )->GetParameter( 6 ) );
          }
        if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 7 )
          {
          alpha =
            parser->Convert<RealType>( metricOption->GetFunction( currentMetricNumber )->GetParameter( 7 ) );
          }
        if( metricOption->GetFunction( currentMetricNumber )->GetNumberOfParameters() > 8 )
          {
          useAnisotropicCovariances =
            parser->Convert<bool>( metricOption->GetFunction( currentMetricNumber )->GetParameter( 8 ) );
          }
        std::string fixedFileName = metricOption->GetFunction( currentMetricNumber )->GetParameter( 0 );
        std::string movingFileName = metricOption->GetFunction( currentMetricNumber )->GetParameter( 1 );
        if( verbose )
          {
          std::cout << "  fixed labeled point set: " << fixedFileName << std::endl;
          std::cout << "  moving labeled point set: " << movingFileName << std::endl;
          }

        ReadLabeledPointSet<LabeledPointSetType>( fixedLabeledPointSet, fixedFileName.c_str(), useBoundaryPointsOnly, samplingPercentage );
        ReadLabeledPointSet<LabeledPointSetType>( movingLabeledPointSet, movingFileName.c_str(), useBoundaryPointsOnly, samplingPercentage );

        fixedLabeledPointSet->DisconnectPipeline();
        movingLabeledPointSet->DisconnectPipeline();
        }
      }

    regHelper->AddMetric( currentMetric,
                          fixedImage,
                          movingImage,
                          fixedLabeledPointSet,
                          movingLabeledPointSet,
                          fixedIntensityPointSet,
                          movingIntensityPointSet,
                          stageID,
                          metricWeighting,
                          samplingStrategy,
                          numberOfBins,
                          radius,
                          useBoundaryPointsOnly,
                          pointSetSigma,
                          evaluationKNeighborhood,
                          alpha,
                          useAnisotropicCovariances,
                          samplingPercentage,
                          intensityDistanceSigma,
                          euclideanDistanceSigma
                          );
    }

  // Perform the registration

  if( regHelper->DoRegistration() == EXIT_FAILURE )
    {
    return EXIT_FAILURE;
    }

  // write out transforms stored in the composite
  typename CompositeTransformType::Pointer resultTransform = regHelper->GetModifiableCompositeTransform();
  unsigned int numTransforms = resultTransform->GetNumberOfTransforms();

  ////
  typedef typename RegistrationHelperType::CompositeTransformType         CompositeTransformType;
  typedef typename CompositeTransformType::Pointer                        CompositeTransformPointer;
  typedef typename RegistrationHelperType::DisplacementFieldTransformType DisplacementFieldTransformType;
  typedef typename RegistrationHelperType::TransformType                  TransformType;

  if( saveStateOption && saveStateOption->GetNumberOfFunctions() )
    {
    CompositeTransformPointer savedStateTx =
      dynamic_cast<CompositeTransformType *>( regHelper->GetModifiableRegistrationState() );
    if( savedStateTx.IsNull() )
      {
      return EXIT_FAILURE;
      }
    unsigned int numStateComponents = savedStateTx->GetNumberOfTransforms();
    // If the last two transforms are displacement field transforms, we add their inverse displacement field to the saved state composite.
    if( savedStateTx->GetNthTransform( numStateComponents-1 )->GetTransformCategory() == TransformType::TransformCategoryEnum::DisplacementField
       && savedStateTx->GetNthTransform( numStateComponents-2 )->GetTransformCategory() == TransformType::TransformCategoryEnum::DisplacementField )
      {
      typename DisplacementFieldTransformType::Pointer oneToEndTransform =
        dynamic_cast<DisplacementFieldTransformType *>( savedStateTx->GetNthTransform( numStateComponents-2 ).GetPointer() );
      typename DisplacementFieldTransformType::Pointer endTransform =
        dynamic_cast<DisplacementFieldTransformType *>( savedStateTx->GetNthTransform( numStateComponents-1 ).GetPointer() );
      if( oneToEndTransform && oneToEndTransform->GetInverseDisplacementField()
         && endTransform && endTransform->GetInverseDisplacementField() )
        {
        savedStateTx->RemoveTransform();
        savedStateTx->AddTransform( oneToEndTransform->GetInverseTransform() );
        savedStateTx->AddTransform( endTransform );
        savedStateTx->AddTransform( endTransform->GetInverseTransform() );
        }
      }
    typename RegistrationHelperType::CompositeTransformType::TransformTypePointer savedStateCompositeTransform =
      savedStateTx.GetPointer();
    const std::string saveStateFileName = saveStateOption->GetFunction( 0 )->GetName();

    // The savedState includes:
    // output linear transforms
    //  + SyN FixedToMiddle displacement field + SyN FixedToMiddle inverse displacement field
    //  + SyN MovingToMiddle displacement field + SyN MovingToMiddle inverse displacement field
    //
    itk::ants::WriteTransform<TComputeType, VImageDimension>( savedStateCompositeTransform,
                                                             saveStateFileName.c_str() );
    }

  // write out transforms actually computed, so skip any initial transforms unless
  // we're collapsing the output transforms.

    CompositeTransformPointer transformToWrite;
    if( shouldCollapseBeDone )
      {
      transformToWrite = regHelper->CollapseCompositeTransform( resultTransform );

      numTransforms = transformToWrite->GetNumberOfTransforms();
      TransformTypeNames.clear();
      for( unsigned int i = 0; i < numTransforms; i++ )
        {
        if( transformToWrite->GetNthTransform( i )->GetTransformCategory() == TransformType::TransformCategoryEnum::Linear )
          {
          // The output type must be Affine, not matrixoffset!  TransformTypeNames.push_back( "matrixoffset" );
          TransformTypeNames.emplace_back("genericaffine" );
          }
        else if( transformToWrite->GetNthTransform( i )->GetTransformCategory() ==
          TransformType::TransformCategoryEnum::DisplacementField )
          {
          typename DisplacementFieldTransformType::Pointer nthTransform =
            dynamic_cast<DisplacementFieldTransformType *>( transformToWrite->GetNthTransform( i ).GetPointer() );

          // We don't know what set of displacement field transforms were optimized.
          // All we know is whether or not an inverse displacement field exists.  If so,
          // we simply pass a transform name which either does have an inverse or does
          // not.
          if( nthTransform && nthTransform->GetInverseDisplacementField() )
            {
            TransformTypeNames.emplace_back("syn" );
            }
          else
            {
            TransformTypeNames.emplace_back("gdf" );
            }
          }
        else if( transformToWrite->GetNthTransform( i )->GetTransformCategory() == TransformType::TransformCategoryEnum::BSpline )
          {
          TransformTypeNames.emplace_back("bspline" );
          }
        }
      }
    else
      {
      transformToWrite = resultTransform.GetPointer();
      }
    if( writeCompositeTransform )
      {
      std::string compositeTransformFileName = outputPrefix;
      if( useMincFormat )
        {
        compositeTransformFileName += std::string( ".xfm" );
        }
      else
        {
        compositeTransformFileName += std::string( "Composite.h5" );
        }

      std::string inverseCompositeTransformFileName = outputPrefix;
      if( useMincFormat )
        {
        inverseCompositeTransformFileName += std::string( "_inverse.xfm" );
        }
      else
        {
        inverseCompositeTransformFileName += std::string( "InverseComposite.h5" );
        }

      typename RegistrationHelperType::CompositeTransformType::TransformTypePointer compositeTransform =
        transformToWrite.GetPointer();
      itk::ants::WriteTransform<TComputeType, VImageDimension>( compositeTransform, compositeTransformFileName.c_str() );

      typename RegistrationHelperType::CompositeTransformType::TransformTypePointer inverseCompositeTransform =
        compositeTransform->GetInverseTransform();
      if( inverseCompositeTransform.IsNotNull() )
        {
        itk::ants::WriteTransform<TComputeType, VImageDimension>( inverseCompositeTransform,
          inverseCompositeTransformFileName.c_str() );
        }
      }
    else //Write out each individual transform
      {
      const unsigned int startIndex = ( shouldCollapseBeDone ) ? 0 : initialMovingTransformOption->GetNumberOfFunctions();
      for( unsigned int i = startIndex; i < numTransforms; ++i )
        {
        typename CompositeTransformType::TransformTypePointer currentTransform = transformToWrite->GetNthTransform( i );

        // only registrations not part of the initial transforms in the
        // TransformTypeNames list.
        const std::string currentTransformType = TransformTypeNames.front();
        TransformTypeNames.pop_front();

        bool writeInverse;
        bool writeVelocityField;

        std::string transformTemplateName = RegTypeToFileName( currentTransformType, writeInverse, writeVelocityField, useMincFormat );

        std::stringstream currentFileName;
        currentFileName << outputPrefix << i << transformTemplateName;

        // WriteTransform will spit all sorts of error messages if it
        // fails, and we want to keep going even if it does so ignore its
        // return value.
        itk::ants::WriteTransform<TComputeType, VImageDimension>( currentTransform, currentFileName.str() );

        typename DisplacementFieldTransformType::Pointer dispTransform =
          dynamic_cast<DisplacementFieldTransformType *>( currentTransform.GetPointer() );
        if( writeInverse && dispTransform.IsNotNull() )
          {
          std::stringstream currentInverseFileName;
          if( useMincFormat )
            {
            currentInverseFileName << outputPrefix << i << "_inverse" << transformTemplateName;
            }
          else
            {
            currentInverseFileName << outputPrefix << i << "Inverse" << transformTemplateName;
            }

          // write inverse transform file
          itk::ants::WriteInverseTransform<TComputeType, VImageDimension>( dispTransform, currentInverseFileName.str() );
          }
        if( writeVelocityField )
          {
          // write velocity field (if applicable)
          typedef typename RegistrationHelperType::TimeVaryingVelocityFieldTransformType
            TimeVaryingVelocityFieldTransformType;

          typedef itk::GaussianExponentialDiffeomorphicTransform<TComputeType, VImageDimension>
            GaussianDisplacementFieldTransformType;

          typename TimeVaryingVelocityFieldTransformType::Pointer tvVelocityFieldTransform =
            dynamic_cast<TimeVaryingVelocityFieldTransformType *>( currentTransform.GetPointer() );
          typename GaussianDisplacementFieldTransformType::Pointer constVelocityFieldTransform =
            dynamic_cast<GaussianDisplacementFieldTransformType *>( currentTransform.GetPointer() );

          std::stringstream currentVelocityFieldFileName;
          if( useMincFormat )
            {
            currentVelocityFieldFileName << outputPrefix << i << "_VelocityField.mnc";
            }
          else
            {
            currentVelocityFieldFileName << outputPrefix << i << "VelocityField.nii.gz";
            }

          try
            {
            if( !tvVelocityFieldTransform.IsNull() )
              {

              typedef itk::Image<itk::Vector<TComputeType, VImageDimension>, VImageDimension + 1> VelocityFieldType;
              typedef itk::ImageFileWriter<VelocityFieldType> VelocityFieldWriterType;
              typename VelocityFieldWriterType::Pointer velocityFieldWriter = VelocityFieldWriterType::New();

              velocityFieldWriter->SetInput( tvVelocityFieldTransform->GetTimeVaryingVelocityField() );
              velocityFieldWriter->SetFileName( currentVelocityFieldFileName.str().c_str() );
                velocityFieldWriter->Update();
              }
            else if( !constVelocityFieldTransform.IsNull() )
              {
              typedef itk::Image<itk::Vector<TComputeType, VImageDimension>, VImageDimension> VelocityFieldType;
              typedef itk::ImageFileWriter<VelocityFieldType> VelocityFieldWriterType;
              typename VelocityFieldWriterType::Pointer velocityFieldWriter = VelocityFieldWriterType::New();

              velocityFieldWriter->SetInput( constVelocityFieldTransform->GetModifiableConstantVelocityField() );
              velocityFieldWriter->SetFileName( currentVelocityFieldFileName.str().c_str() );
                velocityFieldWriter->Update();
              }
            }
          catch( itk::ExceptionObject & err )
            {
            if( verbose )
              {
              std::cerr << "Can't write velocity field transform file " << currentVelocityFieldFileName.str().c_str()
                << std::endl;
              std::cerr << "Exception Object caught: " << std::endl;
              std::cerr << err << std::endl;
              }
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
    (itk::NumericTraits<typename ImageType::SpacingType::ValueType>::ZeroValue());
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

  if( !outputWarpedImageName.empty() )
    {
    typename ImageType::Pointer warpedImage = regHelper->GetWarpedImage();
    if( warpedImage.IsNotNull() )
      {
      WriteImage<ImageType>( warpedImage, outputWarpedImageName.c_str()  );
     }
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

extern int antsRegistration4DDouble(ParserType::Pointer & parser);

extern int antsRegistration2DFloat(ParserType::Pointer & parser);

extern int antsRegistration3DFloat(ParserType::Pointer & parser);

extern int antsRegistration4DFloat(ParserType::Pointer & parser);

} // End namespace

#endif // __ANTSREGISTRATIONTEMPLATEHEADER_H__
