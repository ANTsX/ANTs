#include "antsCommandLineParser.h"
#include "antsUtilities.h"
#include "ReadWriteImage.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

namespace ants
{
int antsSurfaceFunction( itk::ants::CommandLineParser *parser )
{
  const unsigned int ImageDimension = 3;

  typedef float                                      RealType;

  typedef RealType                                   PixelType;
  typedef unsigned char                              RgbComponentType;
  typedef itk::RGBPixel<RgbComponentType>            RgbPixelType;

  typedef itk::Image<PixelType, ImageDimension>      ImageType;

  typedef itk::Image<RgbPixelType, ImageDimension>   RgbImageType;

  return EXIT_SUCCESS;
}

void InitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

    {
    std::string description =
      std::string( "Main input binary image for 3-D rendering.  " );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "surface-image" );
    option->SetShortName( 's' );
    option->SetUsageOption( 0, "surfaceImageFilename" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "A functional overlay can be specified using both " )
      + std::string( "and rgb image and a mask specifying where that " )
      + std::string( "rgb image should be applied.  Both images must " )
      + std::string( "have the same image geometry as the input image. " )
      + std::string( "Optionally, an alpha parameter can be specified." )
      + std::string( "Note that more than one functional overlays can " )
      + std::string( "be rendered." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "functional-overlay" );
    option->SetShortName( 'f' );
    option->SetUsageOption( 0, "[rgbImageFileName,maskImageFileName,<alpha=1>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "The output consists of a." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "output" );
    option->SetShortName( 'o' );
    option->SetUsageOption( 0, "surface.vtk" );
    option->SetUsageOption( 1, "[surface.vtk,<output2dImageFileFormat>]" );
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
int antsSurfaceFunction( std::vector<std::string> args, std::ostream* /*out_stream = NULL */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "antsSurfaceFunction" );

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
  argv[argc] = 0;
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
    std::string( "Produce a 3-D surface rendering with optional RGB overlay." );

  parser->SetCommandDescription( commandDescription );
  InitializeCommandLineOptions( parser );

  parser->Parse( argc, argv );

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

  std::string filename;

  itk::ants::CommandLineParser::OptionType::Pointer imageOption =
    parser->GetOption( "surface-image" );
  if( imageOption && imageOption->GetNumberOfFunctions() > 0 )
    {
    filename = imageOption->GetFunction( 0 )->GetName();

    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
        filename.c_str(), itk::ImageIOFactory::ReadMode );
    unsigned int dimension = imageIO->GetNumberOfDimensions();

    if( dimension == 3 )
      {
      antsSurfaceFunction( parser );
      }
    else
      {
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
      }
    }
  else
    {
    std::cout << "Input image not specified." << std::endl;
    return EXIT_FAILURE;
    }


  return EXIT_SUCCESS;
}
} // namespace ants
