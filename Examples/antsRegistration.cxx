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
#include "itkantsRegistrationHelper.h"

void ConvertToLowerCase( std::string& str )
{
  std::transform( str.begin(), str.end(), str.begin(), tolower );
// You may need to cast the above line to (int(*)(int))
// tolower - this works as is on VC 7.1 but may not work on
// other compilers
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
      + std::string( "MeanSquares:  Thirion's MeanSquares (modified mean-squares). " )
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
    option->SetUsageOption(
      0,
      "CC[fixedImage,movingImage,metricWeight,radius,<samplingStrategy={Regular,Random}>,<samplingPercentage=[0,1]>]" );
    option->SetUsageOption(
      1,
      "MI[fixedImage,movingImage,metricWeight,numberOfBins,<samplingStrategy={Regular,Random}>,<samplingPercentage=[0,1]>]" );
    option->SetUsageOption(
      2,
      "Mattes[fixedImage,movingImage,metricWeight,numberOfBins,<samplingStrategy={Regular,Random}>,<samplingPercentage=[0,1]>]" );
    option->SetUsageOption(
      3,
      "MeanSquares[fixedImage,movingImage,metricWeight,radius,<samplingStrategy={Regular,Random}>,<samplingPercentage=[0,1]>]" );
    option->SetUsageOption(
      4,
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
    option->SetUsageOption(
      6, "GaussianDisplacementField[gradientStep,updateFieldSigmaInPhysicalSpace,totalFieldSigmaInPhysicalSpace]" );
    option->SetUsageOption(
      7,
      "BSplineDisplacementField[gradientStep,updateFieldMeshSizeAtBaseLevel,totalFieldMeshSizeAtBaseLevel,<splineOrder=3>]" );
    option->SetUsageOption(
      8,
      "TimeVaryingVelocityField[gradientStep,numberOfTimeIndices,updateFieldSigmaInPhysicalSpace,updateFieldTimeSigma,totalFieldSigmaInPhysicalSpace,totalFieldTimeSigma]" );
    option->SetUsageOption(
      9,
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

typedef itk::ants::CommandLineParser ParserType;
typedef ParserType::OptionType       OptionType;

template <unsigned VDimension>
int
DoRegistration(typename ParserType::Pointer & parser)
{
  typedef typename itk::ants::RegistrationHelper<VDimension> RegistrationHelperType;
  typename RegistrationHelperType::Pointer regHelper =
    RegistrationHelperType::New();

  regHelper->SetWriteOutputs(true);

  OptionType::Pointer transformOption = parser->GetOption( "transform" );

  OptionType::Pointer metricOption = parser->GetOption( "metric" );

  OptionType::Pointer iterationsOption = parser->GetOption( "iterations" );

  OptionType::Pointer shrinkFactorsOption = parser->GetOption( "shrink-factors" );

  OptionType::Pointer smoothingSigmasOption = parser->GetOption( "smoothing-sigmas" );

  OptionType::Pointer outputOption = parser->GetOption( "output" );

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
  regHelper->SetOutputTransformPrefix(outputPrefix);
  if( outputOption->GetNumberOfParameters(0) > 1 )
    {
    std::string outputWarpedImageName = outputOption->GetParameter( 0, 1 );
    regHelper->SetOutputWarpedImageName(outputWarpedImageName);
    }
  if( outputOption->GetNumberOfParameters(0) > 2 )
    {
    std::string outputInverseWarpedImageName = outputOption->GetParameter( 0, 2 );
    regHelper->SetOutputInverseWarpedImageName(outputInverseWarpedImageName);
    }

  ParserType::OptionType::Pointer initialTransformOption = parser->GetOption( "initial-transform" );

  if( initialTransformOption && initialTransformOption->GetNumberOfValues() > 0 )
    {
    for( unsigned int n = 0; n < initialTransformOption->GetNumberOfValues(); n++ )
      {
      std::string initialTransformName;

      bool useInverse(false);

      if( initialTransformOption->GetNumberOfParameters(n) == 0 )
        {
        initialTransformName = initialTransformOption->GetValue( n );
        }
      else
        {
        initialTransformName = initialTransformOption->GetParameter( n, 0 );
        if( initialTransformOption->GetNumberOfParameters( n ) > 1 )
          {
          useInverse = parser->Convert<bool>( initialTransformOption->GetParameter( n, 1  ) );
          }
        }
      regHelper->AddInitialTransform(initialTransformName, useInverse);
      }
    }

  unsigned int numberOfStages;

  if( transformOption.IsNull() || ( numberOfStages = transformOption->GetNumberOfValues() ) == 0 )
    {
    std::cerr << "No transformations are specified." << std::endl;
    return EXIT_FAILURE;
    }
  std::vector<std::vector<unsigned int> > iterationList;
  std::vector<std::vector<unsigned int> > shrinkFactorsList;
  std::vector<std::vector<float> >        smoothingSigmasList;
  for( unsigned int currentStage = 0; currentStage < numberOfStages; currentStage++ )
    {
    // Get the fixed and moving images

    std::string fixedImageFileName = metricOption->GetParameter( currentStage, 0 );
    std::string movingImageFileName = metricOption->GetParameter( currentStage, 1 );
    std::string whichMetric = metricOption->GetValue( currentStage );
    ConvertToLowerCase( whichMetric );

    typename RegistrationHelperType::MetricEnumeration curMetric =
      regHelper->StringToMetricType(whichMetric);

    float lowerQuantile = 0.0;
    float upperQuantile = 1.0;

    OptionType::Pointer winsorizeOption = parser->GetOption( "winsorize-image-intensities" );
    bool                doWinsorize(false);

    if( winsorizeOption && winsorizeOption->GetNumberOfParameters( 0 ) > 0 )
      {
      doWinsorize = true;
      if( winsorizeOption->GetNumberOfParameters( 0 ) > 0 )
        {
        lowerQuantile = parser->Convert<float>( winsorizeOption->GetParameter( 0, 0 ) );
        }
      if( winsorizeOption->GetNumberOfParameters( 0 ) > 1 )
        {
        upperQuantile = parser->Convert<float>( winsorizeOption->GetParameter( 0, 1 ) );
        }
      }
    regHelper->SetWinsorizeImageIntensities(doWinsorize, lowerQuantile, upperQuantile);

    bool                doHistogramMatch = false;
    OptionType::Pointer histOption = parser->GetOption( "use-histogram-matching" );
    if( histOption && histOption->GetNumberOfValues() > 0 )
      {
      std::string histValue = histOption->GetValue( 0 );
      ConvertToLowerCase( histValue );
      if( histValue.compare( "1" ) == 0 || histValue.compare( "true" ) == 0 )
        {
        doHistogramMatch = true;
        }
      }

    // Get the number of iterations and use that information to specify the number of levels

    std::vector<unsigned int> iterations =
      parser->ConvertVector<unsigned int>( iterationsOption->GetValue( currentStage ) );

    iterationList.push_back(iterations);

    unsigned int numberOfLevels = iterations.size();
    std::cout << "  number of levels = " << numberOfLevels << std::endl;

    // Get shrink factors

    std::vector<unsigned int> factors =
      parser->ConvertVector<unsigned int>( shrinkFactorsOption->GetValue( currentStage ) );
    shrinkFactorsList.push_back(factors);

    // Get smoothing sigmas
    std::vector<float> sigmas = parser->ConvertVector<float>( smoothingSigmasOption->GetValue( currentStage ) );
    smoothingSigmasList.push_back(sigmas);

    float samplingPercentage = 1.0;
    if( metricOption->GetNumberOfParameters( currentStage ) > 5 )
      {
      samplingPercentage = parser->Convert<float>( metricOption->GetParameter( currentStage, 5 ) );
      }

    std::string Strategy = "";
    if( metricOption->GetNumberOfParameters( currentStage ) > 4 )
      {
      Strategy = metricOption->GetParameter( currentStage, 4 );
      }
    ConvertToLowerCase( Strategy );

    typename RegistrationHelperType::SamplingStrategy samplingStrategy =
      RegistrationHelperType::regular;
    if( Strategy == "random" )
      {
      samplingStrategy = RegistrationHelperType::random;
      }

    switch( curMetric )
      {
      case RegistrationHelperType::CC:
        {
        unsigned int radiusOption = parser->Convert<unsigned int>( metricOption->GetParameter( currentStage, 3 ) );
        regHelper->AddMetric(curMetric,
                             fixedImageFileName,
                             movingImageFileName,
                             1.0,
                             samplingStrategy,
                             1,
                             radiusOption,
                             samplingPercentage);
        }
        break;
      case RegistrationHelperType::GC:
      case RegistrationHelperType::MeanSquares:
        {
        regHelper->AddMetric(curMetric,
                             fixedImageFileName,
                             movingImageFileName,
                             1.0,
                             samplingStrategy,
                             1,
                             1,
                             samplingPercentage);
        }
        break;
      case RegistrationHelperType::Mattes:
      case RegistrationHelperType::MI:
        {
        unsigned int binOption = parser->Convert<unsigned int>( metricOption->GetParameter( currentStage, 3 ) );
        regHelper->AddMetric(curMetric,
                             fixedImageFileName,
                             movingImageFileName,
                             1.0,
                             samplingStrategy,
                             binOption,
                             1,
                             samplingPercentage);
        }
        break;
      default:
        std::cerr << "ERROR: Unrecognized image metric: " << whichMetric << std::endl;
        return EXIT_FAILURE;
      }

    // Set up the optimizer.  To change the iteration number for each level we rely
    // on the command observer.

    float learningRate = parser->Convert<float>( transformOption->GetParameter( currentStage, 0 ) );

    std::string whichTransform = transformOption->GetValue( currentStage );
    ConvertToLowerCase( whichTransform );

    typename RegistrationHelperType::XfrmMethod xfrmMethod = regHelper->StringToXfrmMethod(whichTransform);

    switch( xfrmMethod )
      {
      case RegistrationHelperType::Affine:
        {
        regHelper->AddAffineTransform(learningRate);
        }
        break;
      case RegistrationHelperType::Rigid:
        {
        regHelper->AddRigidTransform(learningRate);
        }
        break;
      case RegistrationHelperType::CompositeAffine:
        {
        regHelper->AddCompositeAffineTransform(learningRate);
        }
        break;
      case RegistrationHelperType::Similarity:
        {
        regHelper->AddSimilarityTransform(learningRate);
        }
        break;
      case RegistrationHelperType::Translation:
        {
        regHelper->AddTranslationTransform(learningRate);
        }
      case RegistrationHelperType::GaussianDisplacementField:
        {
        const float varianceForUpdateField = parser->Convert<float>( transformOption->GetParameter( currentStage, 1 ) );
        const float varianceForTotalField = parser->Convert<float>( transformOption->GetParameter( currentStage, 2 ) );
        regHelper->AddGaussianDisplacementFieldTransform(learningRate, varianceForUpdateField, varianceForTotalField);
        }
      case RegistrationHelperType::BSplineDisplacementField:
        {
        std::vector<unsigned int> meshSizeForTheUpdateField = parser->ConvertVector<unsigned int>(
            transformOption->GetParameter( currentStage, 1 ) );
        std::vector<unsigned int> meshSizeForTheTotalField = parser->ConvertVector<unsigned int>(
            transformOption->GetParameter( currentStage, 2 ) );

        unsigned int splineOrder = 3;
        if( transformOption->GetNumberOfParameters( currentStage ) > 3 )
          {
          splineOrder = parser->Convert<unsigned int>( transformOption->GetParameter( currentStage, 3 ) );
          }

        regHelper->AddBSplineDisplacementFieldTransform(learningRate, meshSizeForTheUpdateField,
                                                        meshSizeForTheTotalField,
                                                        splineOrder);
        }
        break;
      case RegistrationHelperType::TimeVaryingVelocityField:
        {
        unsigned int numberOfTimeIndices = parser->Convert<unsigned int>( transformOption->GetParameter( 0, 1 ) );

        const float varianceForUpdateField = parser->Convert<float>( transformOption->GetParameter( currentStage, 2 ) );
        const float varianceForUpdateFieldTime = parser->Convert<float>( transformOption->GetParameter( currentStage,
                                                                                                        3 ) );
        const float varianceForTotalField = parser->Convert<float>( transformOption->GetParameter( currentStage, 4 ) );
        const float varianceForTotalFieldTime =
          parser->Convert<float>( transformOption->GetParameter( currentStage, 5 ) );
        regHelper->AddTimeVaryingVelocityFieldTransform(learningRate,
                                                        numberOfTimeIndices,
                                                        varianceForUpdateField,
                                                        varianceForUpdateFieldTime,
                                                        varianceForTotalField,
                                                        varianceForTotalFieldTime);
        }
        break;
      case RegistrationHelperType::TimeVaryingBSplineVelocityField:
        {
        std::vector<unsigned int> meshSize =
          parser->ConvertVector<unsigned int>( transformOption->GetParameter( 0, 1 ) );
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
        regHelper->AddTimeVaryingBSplineVelocityFieldTransform(learningRate,
                                                               meshSize,
                                                               numberOfTimePointSamples,
                                                               splineOrder);
        }
      case RegistrationHelperType::SyN:
        {
        const float varianceForUpdateField = parser->Convert<float>( transformOption->GetParameter( currentStage, 1 ) );
        const float varianceForTotalField = parser->Convert<float>( transformOption->GetParameter( currentStage, 2 ) );
        regHelper->AddSyNTransform(learningRate, varianceForUpdateField, varianceForTotalField);
        }
        break;
      default:
        {
        std::cerr << "Unknown registration method " << whichTransform << std::endl;
        }
        break;
      }
    }

  // set the vector-vector parameters accumulated
  regHelper->SetIterations(iterationList);
  regHelper->SetSmoothingSigmas(smoothingSigmasList);
  regHelper->SetShrinkFactors(shrinkFactorsList);

  if( regHelper->DoRegistration() == EXIT_FAILURE )
    {
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;

  // Write out warped image(s), if requested.
}

int main( int argc, char *argv[] )
{
  ParserType::Pointer parser = ParserType::New();

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
    return EXIT_FAILURE;
    }
  else if( parser->Convert<bool>( parser->GetOption( 'h' )->GetValue() ) )
    {
    parser->PrintMenu( std::cout, 5, true );
    return EXIT_FAILURE;
    }
  unsigned int dimension = 3;

  ParserType::OptionType::Pointer dimOption = parser->GetOption( "dimensionality" );
  if( dimOption && dimOption->GetNumberOfValues() > 0 )
    {
    dimension = parser->Convert<unsigned int>( dimOption->GetValue() );
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
      return DoRegistration<2>(parser);
      }
    case 3:
      {
      return DoRegistration<3>(parser);
      }
    default:
      std::cerr << "bad image dimension " << dimension << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
