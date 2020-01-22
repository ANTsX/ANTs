#include "antsAllocImage.h"
#include "antsCommandLineParser.h"
#include "antsUtilities.h"

#include "ReadWriteData.h"

#include "itkNonLocalSuperresolutionImageFilter.h"

#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkGaussianInterpolateImageFunction.h"
#include "itkInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkLabelImageGaussianInterpolateImageFunction.h"
#include "itkLabelImageGenericInterpolateImageFunction.h"

#include "itkTimeProbe.h"

#include "ANTsVersion.h"

namespace ants
{


template <typename TFilter>
class CommandProgressUpdate : public itk::Command
{
public:
  typedef  CommandProgressUpdate                      Self;
  typedef  itk::Command                               Superclass;
  typedef  itk::SmartPointer<CommandProgressUpdate>  Pointer;
  itkNewMacro( CommandProgressUpdate );
protected:

  CommandProgressUpdate() : m_CurrentProgress( 0 ) {};

  typedef TFilter FilterType;

  unsigned int m_CurrentProgress;

public:

  void Execute(itk::Object *caller, const itk::EventObject & event) override
    {
    const auto * filter = dynamic_cast<const TFilter *>( caller );

    if( typeid( event ) == typeid( itk::IterationEvent ) )
      {
      if( filter->GetCurrentIteration() > 0 )
        {
        std::cout << "(epsilon value = " << filter->GetCurrentEpsilon() << ")." << std::endl;
        }
      std::cout << "Level " << filter->GetCurrentIteration() << ": " << std::flush;
      this->m_CurrentProgress = 0;
      }

    auto *po = dynamic_cast<itk::ProcessObject *>( caller );
    if (! po) return;
//    std::cout << po->GetProgress() << std::endl;
    if( typeid( event ) == typeid ( itk::ProgressEvent )  )
      {
      if( this->m_CurrentProgress < 99 )
        {
        this->m_CurrentProgress++;
        if( this->m_CurrentProgress % 10 == 0 )
          {
          std::cout << this->m_CurrentProgress << std::flush;
          }
        else
          {
          std::cout << "*" << std::flush;
          }
        }
      }
    }

  void Execute(const itk::Object * object, const itk::EventObject & event) override
    {
    const auto * filter = dynamic_cast<const TFilter *>( object );

    if( typeid( event ) == typeid( itk::IterationEvent ) )
      {
      if( filter->GetCurrentIteration() > 0 )
        {
        std::cout << "(epsilon value = " << filter->GetCurrentEpsilon() << ")." << std::endl;
        }
      std::cout << "Level " << filter->GetCurrentIteration() << ": " << std::flush;
      this->m_CurrentProgress = 0;
      }

    auto *po = dynamic_cast<itk::ProcessObject *>(
      const_cast<itk::Object *>( object ) );
    if (! po) return;

    if( typeid( event ) == typeid ( itk::ProgressEvent )  )
      {
      if( this->m_CurrentProgress < 99 )
        {
        this->m_CurrentProgress++;
        if( this->m_CurrentProgress % 10 == 0 )
          {
          std::cout << this->m_CurrentProgress << std::flush;
          }
        else
          {
          std::cout << "*" << std::flush;
          }
        }
      }
    }
};

template <unsigned int ImageDimension>
int NonLocalSuperResolution( itk::ants::CommandLineParser *parser )
{
  typedef float RealType;

  typedef typename itk::ants::CommandLineParser::OptionType   OptionType;

  bool verbose = false;
  typename itk::ants::CommandLineParser::OptionType::Pointer verboseOption =
    parser->GetOption( "verbose" );
  if( verboseOption && verboseOption->GetNumberOfFunctions() )
    {
    verbose = parser->Convert<bool>( verboseOption->GetFunction( 0 )->GetName() );
    }

  if( verbose )
    {
    std::cout << std::endl << "Running for "
             << ImageDimension << "-dimensional images." << std::endl << std::endl;
    }

  typedef itk::Image<RealType, ImageDimension> ImageType;
  typename ImageType::Pointer inputImage = nullptr;

  typename OptionType::Pointer inputImageOption = parser->GetOption( "input-image" );
  if( inputImageOption && inputImageOption->GetNumberOfFunctions() )
    {
    std::string inputFile = inputImageOption->GetFunction( 0 )->GetName();
    ReadImage<ImageType>( inputImage, inputFile.c_str() );
    }
  else
    {
    if( verbose )
      {
      std::cerr << "Input image not specified." << std::endl;
      }
    return EXIT_FAILURE;
    }

  typename ImageType::Pointer referenceImage = nullptr;
  typename OptionType::Pointer referenceImageOption = parser->GetOption( "reference-image" );
  typename ImageType::Pointer interpolatedImage = nullptr;
  typename OptionType::Pointer interpolatedImageOption = parser->GetOption( "interpolated-image" );

  if( referenceImageOption && referenceImageOption->GetNumberOfFunctions() )
    {
    std::string inputFile = referenceImageOption->GetFunction( 0 )->GetName();

    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( inputFile.c_str() );

    referenceImage = reader->GetOutput();
    referenceImage->Update();
    referenceImage->DisconnectPipeline();
    }
  else if( interpolatedImageOption && interpolatedImageOption->GetNumberOfFunctions() )
      {
      std::string inputFile = interpolatedImageOption->GetFunction( 0 )->GetName();

      typedef itk::ImageFileReader<ImageType> ReaderType;
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( inputFile.c_str() );

      interpolatedImage = reader->GetOutput();
      interpolatedImage->Update();
      interpolatedImage->DisconnectPipeline();
      }
  else
    {
    if( verbose )
      {
      std::cerr << "Reference image or interpolated image not specified." << std::endl;
      }
    return EXIT_FAILURE;
    }

  typedef itk::NonLocalSuperresolutionImageFilter<ImageType, ImageType> SuperresoluterType;
  typename SuperresoluterType::Pointer superresoluter = SuperresoluterType::New();

  superresoluter->SetLowResolutionInputImage( inputImage );
  if( referenceImage )
    {
    superresoluter->SetHighResolutionReferenceImage( referenceImage );
    superresoluter->SetPerformInitialMeanCorrection( false );
    }
  else if( interpolatedImage )
    {
    superresoluter->SetHighResolutionReferenceImage( interpolatedImage );
    superresoluter->SetPerformInitialMeanCorrection( true );
    }

  typename SuperresoluterType::NeighborhoodRadiusType neighborhoodPatchRadius;
  typename SuperresoluterType::NeighborhoodRadiusType neighborhoodSearchRadius;

  neighborhoodPatchRadius.Fill( 1 );
  neighborhoodSearchRadius.Fill( 3 );

  // Get the search and patch radii
  typename OptionType::Pointer searchRadiusOption = parser->GetOption( "search-radius" );
  if( searchRadiusOption && searchRadiusOption->GetNumberOfFunctions() )
    {
    std::string searchRadiusString = searchRadiusOption->GetFunction( 0 )->GetName();

    std::vector<unsigned int> searchRadius;
    searchRadius.push_back( 3 );
    if( searchRadiusOption && searchRadiusOption->GetNumberOfFunctions() )
      {
      searchRadius = parser->ConvertVector<unsigned int>( searchRadiusString );
      }
    if( searchRadius.size() == 1 )
      {
      for( unsigned int d = 1; d < ImageDimension; d++ )
        {
        searchRadius.push_back( searchRadius[0] );
        }
      }
    if( searchRadius.size() != ImageDimension )
      {
      if( verbose )
        {
        std::cerr << "Search radius specified incorrectly.  Please see usage options." << std::endl;
        }
      return EXIT_FAILURE;
      }
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      neighborhoodSearchRadius[d] = searchRadius[d];
      }
    }
  superresoluter->SetNeighborhoodSearchRadius( neighborhoodSearchRadius );

  typename OptionType::Pointer patchRadiusOption = parser->GetOption( "patch-radius" );
  if( patchRadiusOption && patchRadiusOption->GetNumberOfFunctions() )
    {
    std::vector<unsigned int> patchRadius;
    patchRadius.push_back( 1 );
    patchRadius = parser->ConvertVector<unsigned int>( patchRadiusOption->GetFunction( 0 )->GetName() );

    if( patchRadius.size() == 1 )
      {
      for( unsigned int d = 1; d < ImageDimension; d++ )
        {
        patchRadius.push_back( patchRadius[0] );
        }
      }
    if( patchRadius.size() != ImageDimension )
      {
      if( verbose )
        {
        std::cerr << "Patch radius specified incorrectly.  Please see usage options." << std::endl;
        }
      return EXIT_FAILURE;
      }
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      neighborhoodPatchRadius[d] = patchRadius[d];
      }
    }
  superresoluter->SetNeighborhoodPatchRadius( neighborhoodPatchRadius );

  RealType intensitySigma = 1.0;

  typename OptionType::Pointer intensitySigmaOption = parser->GetOption( "intensity-difference-sigma" );
  if( intensitySigmaOption && intensitySigmaOption->GetNumberOfFunctions() )
    {
    intensitySigma = parser->Convert<RealType>( intensitySigmaOption->GetFunction( 0 )->GetName() );
    }
  superresoluter->SetIntensityDifferenceSigma( intensitySigma );

  RealType patchSimilaritySigma = 1.0;

  typename OptionType::Pointer patchSimilaritySigmaOption = parser->GetOption( "patch-similarity-sigma" );
  if( patchSimilaritySigmaOption && patchSimilaritySigmaOption->GetNumberOfFunctions() )
    {
    patchSimilaritySigma = parser->Convert<RealType>( patchSimilaritySigmaOption->GetFunction( 0 )->GetName() );
    }
  superresoluter->SetPatchSimilaritySigma( patchSimilaritySigma );


  std::vector<RealType> scaleLevels;
  scaleLevels.push_back( 32.0 );
  scaleLevels.push_back( 16.0 );
  scaleLevels.push_back(  8.0 );
  scaleLevels.push_back(  4.0 );
  scaleLevels.push_back(  2.0 );
  scaleLevels.push_back(  1.0 );

  typename OptionType::Pointer scaleLevelsOption = parser->GetOption( "scale-levels" );
  if( scaleLevelsOption && scaleLevelsOption->GetNumberOfFunctions() )
    {
    scaleLevels = parser->ConvertVector<RealType>( scaleLevelsOption->GetFunction( 0 )->GetName() );
    }
  superresoluter->SetScaleLevels( scaleLevels );

  // Get the interpolator and possible parameters
  std::string whichInterpolator( "linear" );
  typename itk::ants::CommandLineParser::OptionType::Pointer interpolationOption = parser->GetOption( "interpolation" );
  if( interpolationOption && interpolationOption->GetNumberOfFunctions() )
    {
    whichInterpolator = interpolationOption->GetFunction( 0 )->GetName();
    ConvertToLowerCase( whichInterpolator );
    }
  if( !std::strcmp( whichInterpolator.c_str(), "multilabel" ) || !std::strcmp( whichInterpolator.c_str(), "genericlabel" ) )
    {
    if( verbose )
      {
      std::cerr << "A label-based interpolator is not appropriate for this application." << std::endl;
      }
    return EXIT_FAILURE;
    }

  const size_t VImageDimension = ImageDimension;
  typename ImageType::SpacingType
    cache_spacing_for_smoothing_sigmas(itk::NumericTraits<typename ImageType::SpacingType::ValueType>::ZeroValue());
  if( !std::strcmp( whichInterpolator.c_str(), "gaussian" ) )
    {
    cache_spacing_for_smoothing_sigmas = referenceImage->GetSpacing();
    }

#include "make_interpolator_snip.tmpl"

  superresoluter->SetInterpolator( interpolator );

  itk::TimeProbe timer;
  timer.Start();

  if( verbose )
    {
    typedef CommandProgressUpdate<SuperresoluterType> CommandType;
    typename CommandType::Pointer observer = CommandType::New();
    superresoluter->AddObserver( itk::ProgressEvent(), observer );
    superresoluter->AddObserver( itk::IterationEvent(), observer );
    }

  try
    {
    // superresoluter->DebugOn();
    superresoluter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    if( verbose )
      {
      std::cerr << "Exception caught: " << e << std::endl;
      }
    return EXIT_FAILURE;
    }

  if( verbose )
    {
    std::cout << std::endl << std::endl;
    superresoluter->Print( std::cout, 3 );
    }

  timer.Stop();
  if( verbose )
    {
    std::cout << "Elapsed time: " << timer.GetMean() << std::endl;
    }

  /**
   * output
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer outputOption =
    parser->GetOption( "output" );
  if( outputOption && outputOption->GetNumberOfFunctions() )
    {
    WriteImage<ImageType>( superresoluter->GetOutput(),  ( outputOption->GetFunction( 0 )->GetName() ).c_str() );
    }

  return EXIT_SUCCESS;
}

void InitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

  {
  std::string description =
      std::string( "This option forces the image to be treated as a specified-" )
      + std::string( "dimensional image.  If not specified, the program tries to " )
      + std::string( "infer the dimensionality from the input image." );
  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "image-dimensionality" );
  option->SetShortName( 'd' );
  option->SetUsageOption( 0, "2/3/4" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "A low-resolution image input image to be superresoluted.  " );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "input-image" );
  option->SetShortName( 'i' );
  option->SetUsageOption( 0, "inputImageFilename" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "An interpolated version of the low-resolution image (such as B-spline). " ) +
    std::string( "One should specify either this option as a secondary input or a high-resolution " ) +
    std::string( "multi-modal counterpart (cf the -k option)." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "interpolated-image" );
  option->SetShortName( 'j' );
  option->SetUsageOption( 0, "inputImageFilename" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "A high resolution reference multi-modal image.  Assumed to be in the same " ) +
    std::string( "space as the low-resolution input image (i.e., registered)." ) +
    std::string( "One should specify either this option as a secondary input or an interpolated " ) +
    std::string( "version (cf the -j option)." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "reference-image" );
  option->SetShortName( 'k' );
  option->SetUsageOption( 0, "inputImageFilename" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Patch radius.  Default = 1x1x1" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "patch-radius" );
  option->SetShortName( 'p' );
  option->SetUsageOption( 0, "1" );
  option->SetUsageOption( 1, "1x1x1" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Search radius.  Default = 3x3x3." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "search-radius" );
  option->SetShortName( 'r' );
  option->SetUsageOption( 0, "3" );
  option->SetUsageOption( 1, "3x3x3" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Intensity difference sigma.  Default = 1.0" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "intensity-difference-sigma" );
  option->SetShortName( 'g' );
  option->SetUsageOption( 0, "1.0" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Patch similarity sigma.  Default = 1.0" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "patch-similarity-sigma" );
  option->SetShortName( 't' );
  option->SetUsageOption( 0, "1.0" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Scale levels.  Default = 32x16x8x2x1" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "scale-levels" );
  option->SetShortName( 's' );
  option->SetUsageOption( 0, "32x16x8x2x1" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Several interpolation options are available in ITK. " )
    + std::string( "These have all been made available." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "interpolation" );
  option->SetShortName( 'n' );
  option->SetUsageOption( 0, "Linear" );
  option->SetUsageOption( 1, "NearestNeighbor" );
  option->SetUsageOption( 2, "Gaussian[<sigma=imageSpacing>,<alpha=1.0>]" );
  option->SetUsageOption( 3, "BSpline[<order=3>]" );
  option->SetUsageOption( 4, "CosineWindowedSinc" );
  option->SetUsageOption( 5, "WelchWindowedSinc" );
  option->SetUsageOption( 6, "HammingWindowedSinc" );
  option->SetUsageOption( 7, "LanczosWindowedSinc" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "The output consists of the noise corrected version of the " )
    + std::string( "input image.  Optionally, one can also output the estimated " )
    + std::string( "noise image." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "output" );
  option->SetShortName( 'o' );
  option->SetUsageOption( 0, "outputImage" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Get Version Information." );
  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "version" );
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
int NonLocalSuperResolution( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "NonLocalSuperResolution" );

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
    std::string( "Non-local super resolution described in the following papers:  " )
    + std::string( "1) JV Manjon, P Coupe, A Buades, V Fonov, DL Collins, and Montserrat Robles. " )
    + std::string( "Non-local MRI Upsampling." )
    + std::string( "Medical Image Analysis, 14:784-792, 2010 and" )
    + std::string( "2) JV Manjon, P Coupe, A Buades, DL Collins, and Montserrat Robles. " )
    + std::string( "MRI Superresolution Using Self-Similarity and Image Priors." )
    + std::string( "International Journal of Biomedical Imaging, 2010." );

  parser->SetCommandDescription( commandDescription );
  InitializeCommandLineOptions( parser );

  if( parser->Parse( argc, argv ) == EXIT_FAILURE )
    {
    return EXIT_FAILURE;
    }

  if( argc == 1 )
    {
    parser->PrintMenu( std::cerr, 5, false );
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
  // Show automatic version
  itk::ants::CommandLineParser::OptionType::Pointer versionOption = parser->GetOption( "version" );
  if( versionOption && versionOption->GetNumberOfFunctions() )
    {
    std::string versionFunction = versionOption->GetFunction( 0 )->GetName();
    ConvertToLowerCase( versionFunction );
    if( versionFunction.compare( "1" ) == 0 || versionFunction.compare( "true" ) == 0 )
      {
      //Print Version Information
      std::cout << ANTs::Version::ExtendedVersionString() << std::endl;
      return EXIT_SUCCESS;
      }
    }
  // Get dimensionality
  unsigned int dimension = 3;

  itk::ants::CommandLineParser::OptionType::Pointer dimOption =
    parser->GetOption( "image-dimensionality" );
  if( dimOption && dimOption->GetNumberOfFunctions() )
    {
    dimension = parser->Convert<unsigned int>( dimOption->GetFunction( 0 )->GetName() );
    }
  else
    {
    // Read in the first intensity image to get the image dimension.
    std::string filename;

    itk::ants::CommandLineParser::OptionType::Pointer imageOption =
      parser->GetOption( "input-image" );
    itk::ants::CommandLineParser::OptionType::Pointer interpolatedImageOption =
      parser->GetOption( "interpolated-image" );
    itk::ants::CommandLineParser::OptionType::Pointer referenceImageOption =
      parser->GetOption( "reference-image" );
    if( imageOption && imageOption->GetNumberOfFunctions() > 0 &&
        (
          ( interpolatedImageOption && interpolatedImageOption->GetNumberOfFunctions() > 0 ) ||
          ( referenceImageOption && referenceImageOption->GetNumberOfFunctions() > 0 )
        )
      )
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
      std::cerr << "Not enough input images were specified.  Specify an input image"
               << " with the -i option and a corresponding high-resoution image.  Either"
               << " an interpolated version (-j) or multi-modal counterpart (-k)." << std::endl;
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
      return NonLocalSuperResolution<2>( parser );
      }
      break;
    case 3:
      {
      return NonLocalSuperResolution<3>( parser );
      }
      break;
    case 4:
      {
      return NonLocalSuperResolution<4>( parser );
      }
      break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
