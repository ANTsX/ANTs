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
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

typedef itk::ants::CommandLineParser ParserType;
typedef ParserType::OptionType       OptionType;

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
    option->SetUsageOption( 11,
                            "BSplineSyN[gradientStep,updateFieldMeshSizeAtBaseLevel,totalFieldMeshSizeAtBaseLevel,<splineOrder=3>]" );
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
    std::string         description = "Image masks to limit voxels considered by the metric.";
    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "masks" );
    option->SetShortName( 'x' );
    option->SetUsageOption( 0, "[fixedImageMask,movingImageMask]" );
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

static
const char *
RegTypeToFileName(const std::string & type, bool & writeInverse, bool & writeVelocityField)
{
  std::string str(type);

  ConvertToLowerCase(str);
  if( str == "syn" ||
      str == "symmetricnormalization" ||
      str == "bsplinesyn" ||
      str == "timevaryingbsplinevelocityfield" ||
      str == "tvdmffd" ||
      str == "timevaryingvelocityfield" ||
      str == "tvf" )
    {
    writeInverse = true;
    }
  else
    {
    writeInverse = false;
    }

  if( str == "timevaryingvelocityfield" ||
      str == "tvf" )
    {
    writeVelocityField = true;
    }
  else
    {
    writeVelocityField = false;
    }

  if( str == "rigid" )
    {
    return "Rigid.mat";
    }
  else if( str == "affine" ||
           str == "compositeaffine" || str == "compaff" )
    {
    return "Affine.mat";
    }
  else if( str == "similarity" )
    {
    return "Similarity.mat";
    }
  else if( str == "translation" )
    {
    return "Translation.mat";
    }
  else if( str == "bspline" ||
           str == "ffd" )
    {
    return "BSpline.txt";
    }
  else if( str == "gaussiandisplacementfield" ||
           str == "gdf" ||
           str == "bsplinedisplacementfield" ||
           str == "dmffd" ||
           str == "syn" ||
           str == "symmetricnormalization" ||
           str == "bsplinesyn" ||
           str == "timevaryingvelocityfield" ||
           str == "tvf" ||
           str == "timevaryingbsplinevelocityfield" ||
           str == "tvdmffd" )
    {
    return "Warp.nii.gz";
    }
  return "BOGUS.XXXX";
}

namespace
{
bool MatOffRegistered(false);
}

template <unsigned VImageDimension>
int
AddInitialTransform(
  typename itk::ants::RegistrationHelper<VImageDimension>::CompositeTransformType::Pointer & compositeTransform,
  const std::string & filename,
  bool useInverse)
{
  typedef itk::ants::RegistrationHelper<VImageDimension>                  RegistrationHelperType;
  typedef typename RegistrationHelperType::DisplacementFieldTransformType DisplacementFieldTransformType;

  if( !MatOffRegistered )
    {
    MatOffRegistered = true;
    // Register the matrix offset transform base class to the
    // transform factory for compatibility with the current ANTs.
    typedef itk::MatrixOffsetTransformBase<double, VImageDimension,
                                           VImageDimension> MatrixOffsetTransformType;
    itk::TransformFactory<MatrixOffsetTransformType>::RegisterTransform();
    }

  typedef typename RegistrationHelperType::TransformType TransformType;
  typename TransformType::Pointer initialTransform;

  bool hasTransformBeenRead = false;

  typedef typename DisplacementFieldTransformType::DisplacementFieldType DisplacementFieldType;

  typedef itk::ImageFileReader<DisplacementFieldType> DisplacementFieldReaderType;
  typename DisplacementFieldReaderType::Pointer fieldReader =
    DisplacementFieldReaderType::New();
  try
    {
    fieldReader->SetFileName( filename.c_str() );
    fieldReader->Update();
    hasTransformBeenRead = true;
    }
  catch( ... )
    {
    hasTransformBeenRead = false;
    }

  if( hasTransformBeenRead )
    {
    typename DisplacementFieldTransformType::Pointer displacementFieldTransform =
      DisplacementFieldTransformType::New();
    displacementFieldTransform->SetDisplacementField( fieldReader->GetOutput() );
    initialTransform = dynamic_cast<TransformType *>( displacementFieldTransform.GetPointer() );
    }
  else
    {
    typename itk::TransformFileReader::Pointer initialTransformReader
      = itk::TransformFileReader::New();

    initialTransformReader->SetFileName( filename.c_str() );
    try
      {
      initialTransformReader->Update();
      }
    catch( const itk::ExceptionObject & e )
      {
      std::cerr << "Transform reader for "
                << filename << " caught an ITK exception:\n";
      e.Print( std::cerr );
      return EXIT_FAILURE;
      }
    catch( const std::exception & e )
      {
      std::cerr << "Transform reader for "
                << filename << " caught an exception:\n";
      std::cerr << e.what() << std::endl;
      return EXIT_FAILURE;
      }
    catch( ... )
      {
      std::cerr << "Transform reader for "
                << filename << " caught an unknown exception!!!\n";
      return EXIT_FAILURE;
      }

    initialTransform =
      dynamic_cast<TransformType *>( ( ( initialTransformReader->GetTransformList() )->front() ).GetPointer() );

    if( useInverse )
      {
      initialTransform =
        dynamic_cast<TransformType *>(initialTransform->GetInverseTransform().GetPointer() );
      if( initialTransform.IsNull() )
        {
        std::cerr << "Inverse does not exist for " << filename
                  << std::endl;
        return EXIT_FAILURE;
        }
      }
    }
  compositeTransform->AddTransform( initialTransform );
  return EXIT_SUCCESS;
}

template <unsigned VDimension>
int
DoRegistration(typename ParserType::Pointer & parser)
{
  typedef typename itk::ants::RegistrationHelper<VDimension>      RegistrationHelperType;
  typedef typename RegistrationHelperType::ImageType              ImageType;
  typedef typename RegistrationHelperType::CompositeTransformType CompositeTransformType;

  typename RegistrationHelperType::Pointer regHelper =
    RegistrationHelperType::New();

  OptionType::Pointer transformOption = parser->GetOption( "transform" );

  OptionType::Pointer metricOption = parser->GetOption( "metric" );

  OptionType::Pointer iterationsOption = parser->GetOption( "iterations" );

  OptionType::Pointer shrinkFactorsOption = parser->GetOption( "shrink-factors" );

  OptionType::Pointer smoothingSigmasOption = parser->GetOption( "smoothing-sigmas" );

  OptionType::Pointer outputOption = parser->GetOption( "output" );

  OptionType::Pointer maskOption = parser->GetOption( "masks" );

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
  std::string outputWarpedImageName;
  if( outputOption->GetNumberOfParameters( 0 ) > 1 )
    {
    outputWarpedImageName = outputOption->GetParameter( 0, 1 );
    }

  std::string outputInverseWarpedImageName;
  if( outputOption->GetNumberOfParameters(0) > 2 )
    {
    outputInverseWarpedImageName = outputOption->GetParameter( 0, 2 );
    }

  ParserType::OptionType::Pointer initialTransformOption = parser->GetOption( "initial-transform" );

  if( initialTransformOption && initialTransformOption->GetNumberOfValues() > 0 )
    {
    typename CompositeTransformType::Pointer compositeTransform = CompositeTransformType::New();
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

      if( AddInitialTransform<VDimension>(compositeTransform, initialTransformName, useInverse) != EXIT_SUCCESS )
        {
        std::cerr << "Can't read initialTransform " << initialTransformName << std::endl;
        return EXIT_FAILURE;
        }
      }
    regHelper->SetInitialTransform(compositeTransform);
    }

  if( maskOption.IsNotNull() )
    {
    typedef typename RegistrationHelperType::MaskImageType MaskImageType;
    typedef itk::ImageFileReader<MaskImageType>            ImageReaderType;
    for( unsigned m = 0; m < maskOption->GetNumberOfParameters(); m++ )
      {
      std::string fname = maskOption->GetParameter(0, m);

      typename MaskImageType::Pointer maskImage;
      typename ImageReaderType::Pointer reader = ImageReaderType::New();

      reader->SetFileName(fname.c_str() );
      try
        {
        reader->Update();
        maskImage = reader->GetOutput();
        }
      catch( itk::ExceptionObject & err )
        {
        std::cerr << "Can't read specified mask image " << fname.c_str() << std::endl;
        std::cerr << "Exception Object caught: " << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
        }
      if( m == 0 )
        {
        regHelper->SetFixedImageMask(maskImage);
        }
      else if( m == 1 )
        {
        regHelper->SetMovingImageMask(maskImage);
        }
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
  std::deque<std::string>                 TransformTypeNames;
  // We iterate backwards because the command line options are stored as a stack (first in last out)
  for( int currentStage = numberOfStages - 1; currentStage >= 0; currentStage-- )
    {
    // Get the fixed and moving images

    std::string fixedImageFileName = metricOption->GetParameter( currentStage, 0 );
    std::string movingImageFileName = metricOption->GetParameter( currentStage, 1 );
    std::cout << "  fixed image: " << fixedImageFileName << std::endl;
    std::cout << "  moving image: " << movingImageFileName << std::endl;

    typename ImageType::Pointer fixedImage;
    typename ImageType::Pointer movingImage;

    typedef itk::ImageFileReader<ImageType> ImageReaderType;
    typename ImageReaderType::Pointer fixedImageReader = ImageReaderType::New();

    fixedImageReader->SetFileName( fixedImageFileName.c_str() );
    fixedImageReader->Update();
    fixedImage = fixedImageReader->GetOutput();
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
    movingImage = movingImageReader->GetOutput();
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

    bool                doHistogramMatch(false);
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
    regHelper->SetUseHistogramMatching(doHistogramMatch);

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

    typename RegistrationHelperType::SamplingStrategy samplingStrategy = RegistrationHelperType::none;
    if( Strategy == "random" )
      {
      samplingStrategy = RegistrationHelperType::random;
      }
    else if( Strategy == "regular" )
      {
      samplingStrategy = RegistrationHelperType::regular;
      }

    switch( curMetric )
      {
      case RegistrationHelperType::CC:
        {
        unsigned int radiusOption = parser->Convert<unsigned int>( metricOption->GetParameter( currentStage, 3 ) );
        regHelper->AddMetric(curMetric,
                             fixedImage,
                             movingImage,
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
                             fixedImage,
                             movingImage,
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
                             fixedImage,
                             movingImage,
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

    TransformTypeNames.push_back(whichTransform);

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
        break;
      case RegistrationHelperType::GaussianDisplacementField:
        {
        const float varianceForUpdateField = parser->Convert<float>( transformOption->GetParameter( currentStage, 1 ) );
        const float varianceForTotalField = parser->Convert<float>( transformOption->GetParameter( currentStage, 2 ) );
        regHelper->AddGaussianDisplacementFieldTransform(learningRate, varianceForUpdateField, varianceForTotalField);
        }
        break;
      case RegistrationHelperType::BSpline:
        {
        std::vector<unsigned int> MeshSizeAtBaseLevel =
          parser->ConvertVector<unsigned int>( transformOption->GetParameter( currentStage, 1 ) );
        regHelper->AddBSplineTransform(learningRate, MeshSizeAtBaseLevel);
        }
        break;
      case RegistrationHelperType::BSplineDisplacementField:
        {
        std::vector<unsigned int> meshSizeForTheUpdateField =
          parser->ConvertVector<unsigned int>( transformOption->GetParameter( currentStage, 1 ) );
        std::vector<unsigned int> meshSizeForTheTotalField =
          parser->ConvertVector<unsigned int>( transformOption->GetParameter( currentStage, 2 ) );

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
        break;
      case RegistrationHelperType::SyN:
        {
        const float varianceForUpdateField = parser->Convert<float>( transformOption->GetParameter( currentStage, 1 ) );
        const float varianceForTotalField = parser->Convert<float>( transformOption->GetParameter( currentStage, 2 ) );
        regHelper->AddSyNTransform(learningRate, varianceForUpdateField, varianceForTotalField);
        }
        break;
      case RegistrationHelperType::BSplineSyN:
        {
        std::vector<unsigned int> meshSizeForTheUpdateField =
          parser->ConvertVector<unsigned int>( transformOption->GetParameter( currentStage, 1 ) );
        std::vector<unsigned int> meshSizeForTheTotalField =
          parser->ConvertVector<unsigned int>( transformOption->GetParameter( currentStage, 2 ) );

        unsigned int splineOrder = 3;
        if( transformOption->GetNumberOfParameters( currentStage ) > 3 )
          {
          splineOrder = parser->Convert<unsigned int>( transformOption->GetParameter( currentStage, 3 ) );
          }

        regHelper->AddBSplineSyNTransform(learningRate, meshSizeForTheUpdateField,
                                          meshSizeForTheTotalField,
                                          splineOrder);
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
  //
  // write out transforms stored in the composite
  typename CompositeTransformType::Pointer resultTransform =
    regHelper->GetCompositeTransform();
  unsigned int numTransforms = resultTransform->GetNumberOfTransforms();
  // write out transforms actually computed, so skip any initial transforms.
  for( unsigned int i = initialTransformOption->GetNumberOfValues(); i < numTransforms; ++i )
    {
    typename RegistrationHelperType::CompositeTransformType::TransformTypePointer curTransform =
      resultTransform->GetNthTransform(i);

    //
    // only registrations not part of the initial transforms in the
    // TransformTypeNames list.
    std::string curTransformType = TransformTypeNames.front();
    TransformTypeNames.pop_front();

    bool writeInverse;
    bool writeVelocityField;

    std::string transformTemplateName =
      RegTypeToFileName(curTransformType, writeInverse, writeVelocityField);

    std::stringstream curFileName;
    curFileName << outputPrefix << i << transformTemplateName;

    // write transform file
    if( transformTemplateName != "Warp.nii.gz" )
      {
      typedef itk::TransformFileWriter TransformWriterType;
      typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
      transformWriter->SetInput(curTransform);
      transformWriter->SetFileName(curFileName.str().c_str() );
      try
        {
        transformWriter->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        std::cerr << "Can't write transform file " << curFileName.str().c_str() << std::endl;
        std::cerr << "Exception Object caught: " << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
        }
      }
    else
    // write displacement field
      {
      typedef typename RegistrationHelperType::DisplacementFieldTransformType
        DisplacementFieldTransformType;

      typedef itk::Image<itk::Vector<double, VDimension>, VDimension> DisplacementFieldType;
      typename DisplacementFieldTransformType::Pointer dispTransform =
        dynamic_cast<DisplacementFieldTransformType *>(curTransform.GetPointer() );
      if( dispTransform.IsNull() )
        {
        std::cerr << "Cast failed on transform of type "
                  << curTransform->GetNameOfClass()
                  << " to write " << curFileName.str().c_str()
                  << std::endl;
        continue;
        }

      typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( dispTransform->GetDisplacementField() );
      writer->SetFileName( curFileName.str().c_str() );
      try
        {
        writer->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        std::cerr << "Can't write transform file " << curFileName.str().c_str() << std::endl;
        std::cerr << "Exception Object caught: " << std::endl;
        std::cerr << err << std::endl;
        }

      if( writeInverse )
        {
        typename DisplacementFieldType::Pointer inverseDispField =
          dispTransform->GetInverseDisplacementField();
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
            std::cerr << "Can't write transform file " << curInverseFileName.str().c_str() << std::endl;
            std::cerr << "Exception Object caught: " << std::endl;
            std::cerr << err << std::endl;
            }
          }
        }
      if( writeVelocityField )
        {
        // write velocity field (if applicable)
        typedef typename RegistrationHelperType::TimeVaryingVelocityFieldTransformType
          VelocityFieldTransformType;

        typedef itk::Image<itk::Vector<double, VDimension>, VDimension + 1> VelocityFieldType;
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
            std::cerr << "Can't write velocity field transform file " << curVelocityFieldFileName.str().c_str()
                      << std::endl;
            std::cerr << "Exception Object caught: " << std::endl;
            std::cerr << err << std::endl;
            }
          }
        }
      }
    }
  typename ImageType::Pointer warpedImage = regHelper->GetWarpedImage();

  typedef itk::ImageFileWriter<ImageType> WarpedImageWriterType;

  if( !outputWarpedImageName.empty() )
    {
    typename WarpedImageWriterType::Pointer writer = WarpedImageWriterType::New();
    writer->SetFileName( outputWarpedImageName.c_str() );
    writer->SetInput( warpedImage );
    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "Can't write warped image " << outputWarpedImageName << std::endl;
      std::cerr << "Exception Object caught: " << std::endl;
      std::cerr << err << std::endl;
      }
    }

  if( !outputInverseWarpedImageName.empty() )
    {
    typename ImageType::Pointer inverseWarpedImage = regHelper->GetInverseWarpedImage();
    if( inverseWarpedImage.IsNotNull() )
      {
      typename WarpedImageWriterType::Pointer inverseWriter = WarpedImageWriterType::New();
      inverseWriter->SetFileName( outputInverseWarpedImageName.c_str() );
      inverseWriter->SetInput( inverseWarpedImage );
      try
        {
        inverseWriter->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        std::cerr << "Can't write inverse warped image " << outputInverseWarpedImageName << std::endl;
        std::cerr << "Exception Object caught: " << std::endl;
        std::cerr << err << std::endl;
        }
      }
    }

  return EXIT_SUCCESS;
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
