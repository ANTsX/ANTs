
#include "antsUtilities.h"
#include <algorithm>

#include "antsCommandLineParser.h"

#include "itkDiReCTImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkImage.h"
#include "ReadWriteData.h"
#include "itkTimeProbe.h"

#include <string>
#include <algorithm>
#include <vector>

namespace ants
{
template <typename TFilter>
class CommandIterationUpdate final : public itk::Command
{
public:
  typedef CommandIterationUpdate  Self;
  typedef itk::Command            Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() = default;
public:

  void Execute(itk::Object *caller, const itk::EventObject & event) override
  {
    Execute( (const itk::Object *) caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event) override
  {
    const auto * filter =
      dynamic_cast<const TFilter *>( object );

    if( typeid( event ) != typeid( itk::IterationEvent ) )
      {
      return;
      }
    std::cout << "  Iteration " << filter->GetElapsedIterations()
             << " (of "
             << filter->GetMaximumNumberOfIterations()
             << ").  ";
    std::cout << "Current energy = " << filter->GetCurrentEnergy() << ".  ";
    if( filter->GetElapsedIterations() >= filter->GetConvergenceWindowSize() )
      {
      std::cout << "(convergence value = "
               << filter->GetCurrentConvergenceMeasurement()
               << ", threshold = " << filter->GetConvergenceThreshold()
               << ")";
      }
    std::cout << std::endl;
  }
};

template <unsigned int ImageDimension>
int DiReCT( itk::ants::CommandLineParser *parser )
{
  typedef float        RealType;
  typedef unsigned int LabelType;

  typedef itk::Image<LabelType, ImageDimension> LabelImageType;
  typename LabelImageType::Pointer segmentationImage;

  typedef itk::Image<RealType, ImageDimension> ImageType;
  typename ImageType::Pointer grayMatterProbabilityImage;
  typename ImageType::Pointer whiteMatterProbabilityImage;
  typename ImageType::Pointer thicknessPriorImage;

  typedef itk::DiReCTImageFilter<LabelImageType, ImageType> DiReCTFilterType;
  typename DiReCTFilterType::Pointer direct = DiReCTFilterType::New();
  typedef typename DiReCTFilterType::LabelType DirectLabelType;

  bool verbose = false;
  typename itk::ants::CommandLineParser::OptionType::Pointer verboseOption =
    parser->GetOption( "verbose" );
  if( verboseOption && verboseOption->GetNumberOfFunctions() )
    {
    verbose = parser->Convert<bool>( verboseOption->GetFunction( 0 )->GetName() );
    }

  if( verbose )
    {
    std::cout << "Running DiReCT for " << ImageDimension <<
      "-dimensional images." << std::endl << std::endl;
    }

  //
  // debugging information
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer
  debugOption = parser->GetOption( "print-debug-information" );
  if( debugOption && debugOption->GetNumberOfFunctions() )
    {
    std::string value = debugOption->GetFunction()->GetName();
    ConvertToLowerCase( value );
    if( std::strcmp( value.c_str(), "true" ) || parser->Convert<int>( value ) != 0 )
      {
      direct->DebugOn();
      }
    }
  //
  // segmentation image
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer
  segmentationImageOption = parser->GetOption( "segmentation-image" );
  if( segmentationImageOption && segmentationImageOption->GetNumberOfFunctions() )
    {
    if( segmentationImageOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
      {
      std::string inputFile = segmentationImageOption->GetFunction( 0 )->GetName();
      ReadImage<LabelImageType>( segmentationImage, inputFile.c_str()   );
      }
    else if( segmentationImageOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      std::string inputFile = segmentationImageOption->GetFunction( 0 )->GetParameter( 0 );
      ReadImage<LabelImageType>( segmentationImage, inputFile.c_str()   );
      if( segmentationImageOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        auto grayMatterValue = parser->Convert<DirectLabelType>(
          segmentationImageOption->GetFunction( 0 )->GetParameter( 1 ) );
        direct->SetGrayMatterLabel( grayMatterValue );
        }
      if( segmentationImageOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
        {
        auto whiteMatterValue = parser->Convert<DirectLabelType>(
          segmentationImageOption->GetFunction( 0 )->GetParameter( 2 ) );
        direct->SetWhiteMatterLabel( whiteMatterValue );
        }
      }
    }
  else
    {
    if( verbose )
      {
      std::cerr << "Segmentation image not specified." << std::endl;
      }
    return EXIT_FAILURE;
    }
  direct->SetSegmentationImage( segmentationImage );
  //
  // gray matter probability image
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer
  grayMatterOption = parser->GetOption( "gray-matter-probability-image" );
  if( grayMatterOption && grayMatterOption->GetNumberOfFunctions() )
    {
    std::string gmFile = grayMatterOption->GetFunction()->GetName();
    ReadImage<ImageType>( grayMatterProbabilityImage, gmFile.c_str()   );
    }
  else
    {
    if( verbose )
      {
      std::cout << "  Grey matter probability image not specified. "
               << "Creating one from the segmentation image." << std::endl;
      }

    typedef itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType> ThresholderType;
    typename ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetInput( segmentationImage );
    thresholder->SetLowerThreshold( direct->GetGrayMatterLabel() );
    thresholder->SetUpperThreshold( direct->GetGrayMatterLabel() );
    thresholder->SetInsideValue( 1 );
    thresholder->SetOutsideValue( 0 );

    typedef itk::DiscreteGaussianImageFilter<LabelImageType, ImageType> SmootherType;
    typename SmootherType::Pointer smoother = SmootherType::New();
    smoother->SetVariance( 1.0 );
    smoother->SetUseImageSpacingOn();
    smoother->SetMaximumError( 0.01 );
    smoother->SetInput( thresholder->GetOutput() );
    smoother->Update();

    grayMatterProbabilityImage = smoother->GetOutput();
    }
  direct->SetGrayMatterProbabilityImage( grayMatterProbabilityImage );
  //
  // white matter probability image
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer
  whiteMatterOption = parser->GetOption( "white-matter-probability-image" );
  if( whiteMatterOption && whiteMatterOption->GetNumberOfFunctions() )
    {
    std::string wmFile = whiteMatterOption->GetFunction( 0 )->GetName();
    ReadImage<ImageType>( whiteMatterProbabilityImage, wmFile.c_str()   );
    }
  else
    {
    if( verbose )
      {
      std::cout << "  White matter probability image not specified. "
               << "Creating one from the segmentation image." << std::endl << std::endl;
      }

    typedef itk::BinaryThresholdImageFilter<LabelImageType, ImageType>
      ThresholderType;
    typename ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetInput( segmentationImage );
    thresholder->SetLowerThreshold( direct->GetWhiteMatterLabel() );
    thresholder->SetUpperThreshold( direct->GetWhiteMatterLabel() );
    thresholder->SetInsideValue( 1 );
    thresholder->SetOutsideValue( 0 );

    typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> SmootherType;
    typename SmootherType::Pointer smoother = SmootherType::New();
    smoother->SetVariance( 1.0 );
    smoother->SetUseImageSpacingOn();
    smoother->SetMaximumError( 0.01 );
    smoother->SetInput( thresholder->GetOutput() );
    smoother->Update();

    whiteMatterProbabilityImage = smoother->GetOutput();
    }
  direct->SetWhiteMatterProbabilityImage( whiteMatterProbabilityImage );

  //
  // label priors
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer
  tpOption = parser->GetOption( "thickness-prior-image" );
  if( tpOption && tpOption->GetNumberOfFunctions() )
    {
    std::string labFile = tpOption->GetFunction( 0 )->GetName();
    ReadImage<ImageType>( thicknessPriorImage, labFile.c_str()   );
    direct->SetThicknessPriorImage( thicknessPriorImage );
    }

  //
  // convergence options
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer convergenceOption =
    parser->GetOption( "convergence" );
  if( convergenceOption && convergenceOption->GetNumberOfFunctions() )
    {
    if( convergenceOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      direct->SetMaximumNumberOfIterations( parser->Convert<unsigned int>(
        convergenceOption->GetFunction( 0 )->GetParameter( 0 ) ) );
      }
    if( convergenceOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
      {
      direct->SetConvergenceThreshold( parser->Convert<float>(
        convergenceOption->GetFunction( 0 )->GetParameter( 1 ) ) );
      }
    if( convergenceOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
      {
      direct->SetConvergenceWindowSize( parser->Convert<unsigned int>(
        convergenceOption->GetFunction( 0 )->GetParameter( 2 ) ) );
      }
    }

  //
  // thickness prior estimate
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer
    thicknessPriorOption = parser->GetOption( "thickness-prior-estimate" );
  if( thicknessPriorOption && thicknessPriorOption->GetNumberOfFunctions() )
    {
    direct->SetThicknessPriorEstimate( parser->Convert<RealType>(
      thicknessPriorOption->GetFunction( 0 )->GetName() ) );
    }

  //
  // gradient step
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer
    gradientStepOption = parser->GetOption( "gradient-step" );
  if( gradientStepOption && gradientStepOption->GetNumberOfFunctions() )
    {
    direct->SetInitialGradientStep( parser->Convert<RealType>(
      gradientStepOption->GetFunction( 0 )->GetName() ) );
    }

  //
  // do B-spline smoothing?
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer
    bsplineSmoothingOption = parser->GetOption( "use-bspline-smoothing" );
  if( bsplineSmoothingOption && bsplineSmoothingOption->GetNumberOfFunctions() )
    {
    direct->SetUseBSplineSmoothing( parser->Convert<bool>(
      bsplineSmoothingOption->GetFunction( 0 )->GetName() ) );
    }


    //
    // do matrix-based smoothing?
    //
    typename itk::ants::CommandLineParser::OptionType::Pointer
      maskedSmoothingOption = parser->GetOption( "use-masked-smoothing" );
    if( maskedSmoothingOption && maskedSmoothingOption->GetNumberOfFunctions() )
      {
      direct->SetUseMaskedSmoothing( parser->Convert<bool>(
        maskedSmoothingOption->GetFunction( 0 )->GetName() ) );
      }

    //
    // time points
    //
    typename itk::ants::CommandLineParser::OptionType::Pointer
      timePointsOption = parser->GetOption( "time-points" );
    if( timePointsOption && timePointsOption->GetNumberOfFunctions() )
      {
      direct->SetTimePoints( parser->ConvertVector<RealType>(
        timePointsOption->GetFunction( 0 )->GetParameter( 0 ) ) );
      if( timePointsOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        direct->SetTimeSmoothingVariance( parser->Convert<float>(
          timePointsOption->GetFunction( 0 )->GetParameter( 1 ) ) );
        }
      }

  typename itk::ants::CommandLineParser::OptionType::Pointer
    restrictOption = parser->GetOption( "restrict-deformation" );
  if( restrictOption && restrictOption->GetNumberOfFunctions() )
    {
    direct->SetRestrictDeformation( parser->Convert<bool>(
      restrictOption->GetFunction( 0 )->GetName() ) );
    }

  //
  // smoothing parameter for the velocity field
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer
    smoothingVelocityFieldParameterOption =
    parser->GetOption( "smoothing-velocity-field-parameter" );
  if( smoothingVelocityFieldParameterOption &&
    smoothingVelocityFieldParameterOption->GetNumberOfFunctions() )
    {
    if( direct->GetUseBSplineSmoothing() )
      {
      direct->SetBSplineSmoothingIsotropicMeshSpacing( parser->Convert<RealType>(
        smoothingVelocityFieldParameterOption->GetFunction( 0 )->GetName() ) );
      }
    else
      {
      direct->SetSmoothingVelocityFieldVariance( parser->Convert<RealType>(
        smoothingVelocityFieldParameterOption->GetFunction( 0 )->GetName() ) );
      }
    }

  //
  // smoothing variance for the hit and total images
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer
    smoothingVarianceOption = parser->GetOption( "smoothing-variance" );
  if( smoothingVarianceOption && smoothingVarianceOption->GetNumberOfFunctions() )
    {
    direct->SetSmoothingVariance( parser->Convert<RealType>(
      smoothingVarianceOption->GetFunction( 0 )->GetName() ) );
    }

  //
  // number of integration points
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer
    numberOfIntegrationPointsOption = parser->GetOption( "number-of-integration-points" );
  if( numberOfIntegrationPointsOption &&
    numberOfIntegrationPointsOption->GetNumberOfFunctions() )
    {
    direct->SetNumberOfIntegrationPoints( parser->Convert<unsigned int>(
      numberOfIntegrationPointsOption->GetFunction( 0 )->GetName() ) );
    }

  //
  // number of invert displacement field iterations
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer
    numberOfInvertDisplacementFieldIterationsOption = parser->GetOption(
    "maximum-number-of-invert-displacement-field-iterations" );
  if( numberOfInvertDisplacementFieldIterationsOption &&
    numberOfInvertDisplacementFieldIterationsOption->GetNumberOfFunctions() )
    {
    direct->SetMaximumNumberOfInvertDisplacementFieldIterations(
      parser->Convert<unsigned int>(
      numberOfInvertDisplacementFieldIterationsOption->GetFunction( 0 )->GetName() ) );
    }

  if( verbose )
    {
    typedef CommandIterationUpdate<DiReCTFilterType> CommandType;
    typename CommandType::Pointer observer = CommandType::New();
    direct->AddObserver( itk::IterationEvent(), observer );
    }

  itk::TimeProbe timer;
  timer.Start();
  try
    {
    direct->Update(); // causes problems with ANTsR , unknown reason
    }
  catch( itk::ExceptionObject & e )
    {
    if( verbose )
      {
      std::cerr << "Exception caught: " << e << std::endl;
      }
    return EXIT_FAILURE;
    }
  timer.Stop();

  if( verbose )
    {
    direct->Print( std::cout, 3 );
    std::cout << "DiReCT elapsed time: " << timer.GetMean() << std::endl;
    }

  /**
   * output
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer outputOption =
    parser->GetOption( "output" );

  if( outputOption && outputOption->GetNumberOfFunctions() > 0 )
    {
    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
      {
      WriteImage<ImageType>( direct->GetOutput(), ( outputOption->GetFunction( 0 )->GetName() ).c_str() );
      }
    else if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      WriteImage<ImageType>( direct->GetOutput(), ( outputOption->GetFunction( 0 )->GetParameter() ).c_str() );
      if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        WriteImage<ImageType>( direct->GetOutput( 1 ), ( outputOption->GetFunction( 0 )->GetParameter( 1 ) ).c_str() );
        }
      }
    }

  if( segmentationImageOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
    {
    std::string inputFile = segmentationImageOption->GetFunction( 0 )->GetName();
    ReadImage<LabelImageType>( segmentationImage, inputFile.c_str()   );
    }
  else if( segmentationImageOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )

    {
    }
  if( outputOption && outputOption->GetNumberOfFunctions() > 1 )
    {
    WriteImage<ImageType>(  direct->GetOutput( 1 ), ( outputOption->GetFunction( 1 )->GetName() ).c_str() );
    }

  return EXIT_SUCCESS;
}

void KellyKapowskiInitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

  {
  std::string description =
    std::string( "This option forces the image to be treated as a specified-" )
    + std::string( "dimensional image.  If not specified, DiReCT tries to " )
    + std::string( "infer the dimensionality from the input image." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "image-dimensionality" );
  option->SetShortName( 'd' );
  option->SetUsageOption( 0, "2/3" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "A segmentation image must be supplied labeling the gray" )
    + std::string( "and white matters.  Default values = 2 and 3, respectively." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "segmentation-image" );
  option->SetShortName( 's' );
  option->SetUsageOption( 0, "imageFilename" );
  option->SetUsageOption( 1, "[imageFilename,<grayMatterLabel=2>,<whiteMatterLabel=3>]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "In addition to the segmentation image, a gray matter " )
    + std::string( "probability image can be used. If no such image is " )
    + std::string( "supplied, one is created using the segmentation image " )
    + std::string( "and a variance of 1.0 mm." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "gray-matter-probability-image" );
  option->SetShortName( 'g' );
  option->SetUsageOption( 0, "imageFilename" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "In addition to the segmentation image, a white matter " )
    + std::string( "probability image can be used. If no such image is " )
    + std::string( "supplied, one is created using the segmentation image " )
    + std::string( "and a variance of 1.0 mm." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "white-matter-probability-image" );
  option->SetShortName( 'w' );
  option->SetUsageOption( 0, "imageFilename" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Convergence is determined by fitting a line to the normalized energy " )
    + std::string( "profile of the last N iterations (where N is specified by " )
    + std::string( "the window size) and determining the slope which is then " )
    + std::string( "compared with the convergence threshold." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "convergence" );
  option->SetShortName( 'c' );
  option->SetUsageOption( 0, "[<numberOfIterations=50>,<convergenceThreshold=0.001>,<convergenceWindowSize=10>]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Provides a prior constraint on the final thickness measurement. Default = 10 mm." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "thickness-prior-estimate" );
  option->SetShortName( 't' );
  option->SetUsageOption( 0, "thicknessPriorEstimate" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "An image containing spatially varying prior thickness values." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "thickness-prior-image" );
  option->SetShortName( 'a' );
  option->SetUsageOption( 0, "thicknessPriorFileName" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Gradient step size for the optimization.  Default = 0.025." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "gradient-step" );
  option->SetShortName( 'r' );
  option->SetUsageOption( 0, "stepSize" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Defines the Gaussian smoothing of the hit and total images.  Default = 1.0." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "smoothing-variance" );
  option->SetShortName( 'l' );
  option->SetUsageOption( 0, "variance" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Defines the Gaussian smoothing of the velocity field (default = 1.5)." )
    + std::string( "If the b-spline smoothing option is chosen, then this " )
    + std::string( "defines the isotropic mesh spacing for the smoothing spline (default = 15)." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "smoothing-velocity-field-parameter" );
  option->SetShortName( 'm' );
  option->SetUsageOption( 0, "variance" );
  option->SetUsageOption( 1, "isotropicMeshSpacing" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Sets the option for B-spline smoothing of the velocity field." )
    + std::string( "Default = false." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "use-bspline-smoothing" );
  option->SetShortName( 'b' );
  option->SetUsageOption( 0, "1/(0)" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Sets the option for masked-based smoothing of the velocity field." )
    + std::string( "Default = false." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "use-masked-smoothing" );
  option->SetShortName( 'x' );
  option->SetUsageOption( 0, "1/(0)" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Time points for irregularly spaced time samples and " ) +
    std::string( "time-variance with which to compute distance metric. " ) +
    std::string( "The user specifies [0.0x1.2x4.5,3] for input with 3 time " ) +
    std::string( "slices where the vector of numeric value defines the time " ) +
    std::string( "of sampling e.g. in years and the scalar value (here \'3\')" ) +
    std::string( "defines the variance." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "time-points" );
  option->SetShortName( 'p' );
  option->SetUsageOption( 0, "1" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Restrict the last dimension's deformation.  Meant for use " ) +
    std::string( "with multiple time points.  Default = false." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "restrict-deformation" );
  option->SetShortName( 'e' );
  option->SetUsageOption( 0, "1/(0)" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Number of compositions of the diffeomorphism per iteration.  Default = 10." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "number-of-integration-points" );
  option->SetShortName( 'n' );
  option->SetUsageOption( 0, "numberOfPoints" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Maximum number of iterations for estimating the invert displacement field.  Default = 20." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "maximum-number-of-invert-displacement-field-iterations" );
  option->SetShortName( 'p' );
  option->SetUsageOption( 0, "numberOfIterations" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "The output consists of a thickness map defined in the " )
    + std::string( "segmented gray matter. " );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "output" );
  option->SetShortName( 'o' );
  option->SetUsageOption( 0, "imageFileName" );
  option->SetUsageOption( 1, "[imageFileName,warpedWhiteMatterImageFileName]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Verbose output." );

  OptionType::Pointer option = OptionType::New();
  option->SetShortName( 'v' );
  option->SetLongName( "verbose" );
  option->SetUsageOption( 0, "(0)/1" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Print the help menu (short version)." );

  OptionType::Pointer option = OptionType::New();
  option->SetShortName( 'h' );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Print the help menu." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "help" );
  option->SetDescription( description );
  parser->AddOption( option );
  }
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int KellyKapowski( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "KellyKapowski" );

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
  argv[argc] = nullptr;
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

  itk::ants::CommandLineParser::Pointer parser =
    itk::ants::CommandLineParser::New();

  parser->SetCommand( argv[0] );

  std::string commandDescription =
    std::string( "DiReCT is a registration based estimate of cortical " )
    + std::string( "thickness.  It was published in S. R. Das, B. B. " )
    + std::string( "Avants, M. Grossman, and J. C. Gee, Registration based " )
    + std::string( "cortical thickness measurement, Neuroimage 2009, " )
    + std::string( "45:867--879.  See also N. J. Tustison, P. A. Cook, " )
    + std::string( "A. Klein, G. Song, S. R. Das, J. T. Duda, B M. Kandel, " )
    + std::string( "N. van Strien, J. R. Stone, J. C. Gee, and B. B. Avants. " )
    + std::string( "Large-Scale Evaluation of ANTs and FreeSurfer Cortical " )
    + std::string( "Thickness Measurements. NeuroImage, 99:166-179, Oct 2014." );

  parser->SetCommandDescription( commandDescription );
  KellyKapowskiInitializeCommandLineOptions( parser );

  if( parser->Parse( argc, argv ) == EXIT_FAILURE )
    {
    return EXIT_FAILURE;
    }

  if( argc == 1 )
    {
    parser->PrintMenu( std::cout, 5, false );
    return EXIT_FAILURE;
    }
  else if( parser->GetOption( "help" )->GetFunction() && parser->Convert<bool>( parser->GetOption( "help" )->GetFunction()->GetName() ) )
    {
    parser->PrintMenu( std::cout, 5, false );
    return EXIT_SUCCESS;
    }
  else if( parser->GetOption( 'h' )->GetFunction() && parser->Convert<bool>( parser->GetOption( 'h' )->GetFunction()->GetName() ) )
    {
    parser->PrintMenu( std::cout, 5, true );
    return EXIT_SUCCESS;
    }

  // Get dimensionality
  unsigned int dimension = 3;

  itk::ants::CommandLineParser::OptionType::Pointer dimOption =
    parser->GetOption( "image-dimensionality" );
  if( dimOption && dimOption->GetNumberOfFunctions() > 0 )
    {
    dimension = parser->Convert<unsigned int>( dimOption->GetFunction( 0 )->GetName() );
    }
  else
    {
    // Read in the first intensity image to get the image dimension.
    std::string filename;

    itk::ants::CommandLineParser::OptionType::Pointer imageOption =
      parser->GetOption( "input-image" );
    if( imageOption && imageOption->GetNumberOfFunctions() > 0 )
      {
      if( imageOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
        {
        filename = imageOption->GetFunction( 0 )->GetParameter( 0 );
        }
      else
        {
        filename = imageOption->GetFunction( 0 )->GetName();
        }
      }
    else
      {
      std::cout << "No input images were specified.  Specify an input "
               << " segmentation image with the -s option" << std::endl;
      return EXIT_FAILURE;
      }
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
        filename.c_str(), itk::ImageIOFactory::FileModeEnum::ReadMode );
    dimension = imageIO->GetNumberOfDimensions();
    }

  switch( dimension )
    {
    case 2:
      {
      return DiReCT<2>( parser );
      }
      break;
    case 3:
      {
      return DiReCT<3>( parser );
      }
      break;
    case 4:
      {
      return DiReCT<4>( parser );
      }
      break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
