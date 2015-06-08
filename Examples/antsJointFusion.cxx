#include "antsCommandLineParser.h"
#include "antsUtilities.h"
#include "antsAllocImage.h"
#include "ReadWriteData.h"

#include "itkWeightedVotingFusionImageFilter.h"
#include "itkTimeProbe.h"


#include <string>
#include <algorithm>
#include <vector>

#include "ANTsVersion.h"

namespace ants
{

// template <class TFilter>
// class CommandIterationUpdate : public itk::Command
// {
// public:
//   typedef CommandIterationUpdate  Self;
//   typedef itk::Command            Superclass;
//   typedef itk::SmartPointer<Self> Pointer;
//   itkNewMacro( Self );
// protected:
//   CommandIterationUpdate()
//   {
//   };
// public:
//
//   void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
//   {
//     Execute( (const itk::Object *) caller, event);
//   }
//
//   void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
//   {
//     const TFilter * filter =
//       dynamic_cast<const TFilter *>( object );
//
//     if( typeid( event ) != typeid( itk::IterationEvent ) )
//       {
//       return;
//       }
//     if( filter->GetElapsedIterations() == 1 )
//       {
//       std::cout << "Current level = " << filter->GetCurrentLevel() + 1
//                << std::endl;
//       }
//     std::cout << "  Iteration " << filter->GetElapsedIterations()
//              << " (of "
//              << filter->GetMaximumNumberOfIterations()[filter->GetCurrentLevel()]
//              << ").  ";
//     std::cout << " Current convergence value = "
//              << filter->GetCurrentConvergenceMeasurement()
//              << " (threshold = " << filter->GetConvergenceThreshold()
//              << ")" << std::endl;
//   }
// };

template <unsigned int ImageDimension>
int antsJointFusion( itk::ants::CommandLineParser *parser )
{
  typedef float RealType;

  typedef itk::Image<RealType, ImageDimension> ImageType;
  typename ImageType::Pointer inputImage = ITK_NULLPTR;

  typedef itk::Image<unsigned, ImageDimension> MaskImageType;
  typename MaskImageType::Pointer maskImage = ITK_NULLPTR;

  bool verbose = false;
  typename itk::ants::CommandLineParser::OptionType::Pointer verboseOption =
    parser->GetOption( "verbose" );
  if( verboseOption && verboseOption->GetNumberOfFunctions() )
    {
    verbose = parser->Convert<bool>( verboseOption->GetFunction( 0 )->GetName() );
    }

  if( verbose )
    {
    std::cout << std::endl << "Running antsJointFusion for "
             << ImageDimension << "-dimensional images." << std::endl << std::endl;
    }

  typedef itk::WeightedVotingFusionImageFilter<ImageType, ImageType> FusionFilterType;
  typename FusionFilterType::Pointer fusionFilter = FusionFilterType::New();


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
    std::string( "The target image (or multimodal target images) assumed to be " )
    + std::string( "aligned to a common image domain." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "target-image" );
  option->SetShortName( 't' );
  option->SetUsageOption( 0, "targetImage" );
  option->SetUsageOption( 1, "[targetImageModality0,targetImageModality1,...,targetImageModalityN]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "The atlas image (or multimodal atlas images) assumed to be " )
    + std::string( "aligned to a common image domain." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "atlas-image" );
  option->SetShortName( 'g' );
  option->SetUsageOption( 0, "atlasImage" );
  option->SetUsageOption( 1, "[atlasImageModality0,atlasImageModality1,...,atlasImageModalityN]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "The atlas segmentation images.  For performing label fusion the number of " )
    + std::string( "specified segmentations should be identical to the number of atlas image sets." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "atlas-segmentation" );
  option->SetShortName( 'l' );
  option->SetUsageOption( 0, "atlasSegmentation" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Regularization term added to matrix Mx for calculating the inverse.  Default = 0.1" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "alpha" );
  option->SetShortName( 'a' );
  option->SetUsageOption( 0, "0.1" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Exponent for mapping intensity difference to the joint error.  Default = 2.0" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "beta" );
  option->SetShortName( 'b' );
  option->SetUsageOption( 0, "2.0" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Exponent for mapping intensity difference to the joint error.  Default = 2.0" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "beta" );
  option->SetShortName( 'b' );
  option->SetUsageOption( 0, "2.0" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Patch radius for similarity measures.  Default = 2x2x2" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "patch-radius" );
  option->SetShortName( 'p' );
  option->SetUsageOption( 0, "2" );
  option->SetUsageOption( 1, "2x2x2" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Search radius for similarity measures.  Default = 3x3x3" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "search-radius" );
  option->SetShortName( 's' );
  option->SetUsageOption( 0, "3" );
  option->SetUsageOption( 1, "3x3x3" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Specify an exclusion region for the given label." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "exclusion-image" );
  option->SetShortName( 'x' );
  option->SetUsageOption( 0, "label[exclusion-image]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "The output is the intensity and/or label fusion image.  Additional " )
    + std::string( "optional outputs include the label posterior probability images " )
    + std::string( "and the atlas voting weight images." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "output" );
  option->SetShortName( 'o' );
  option->SetUsageOption( 0, "intensityFusionImage" );
  option->SetUsageOption( 1, "labelFusionImage" );
  option->SetUsageOption( 2, "[intensityFusionImage,labelFusionImage,<labelPosteriorProbabilityImageFileNameFormat>,<atlasVotingWeightImageFileNameFormat>]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }


  {
  std::string description = std::string( "Get version information." );
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
int antsJointFusion( std::vector<std::string> args, std::ostream* /*out_stream = NULL */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "antsJointFusion" );

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
  argv[argc] = ITK_NULLPTR;
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
    std::string( "antsJointFusion is an image fusion algorithm developed by Hongzhi Wang and " )
    + std::string( "Paul Yushkevich which won segmentation challenges at MICCAI 2012 and MICCAI 2013. " )
    + std::string( "The original label fusion framework was extended to accommodate intensities by " )
    + std::string( "Brain Avants.  This implementation is based on Paul's original ITK-style " )
    + std::string( "implementation and Brian's ANTsR implementation.  References include  1) H. Wang, " )
    + std::string( "J. W. Suh, S. Das, J. Pluta, C. Craige, P. Yushkevich, Multi-atlas " )
    + std::string( "segmentation with joint label fusion IEEE Trans. on Pattern " )
    + std::string( "Analysis and Machine Intelligence, 35(3), 611-623, 2013. and 2) " )
    + std::string( "H. Wang and P. A. Yushkevich, Multi-atlas segmentation with joint " )
    + std::string( "label fusion and corrective learning--an open source implementation, " )
    + std::string( "Front. Neuroinform., 2013. " );

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
      parser->GetOption( "target-image" );
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
      std::cerr << "No input images were specified.  Specify an input image"
               << " with the -t option" << std::endl;
      return EXIT_FAILURE;
      }
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
        filename.c_str(), itk::ImageIOFactory::ReadMode );
    dimension = imageIO->GetNumberOfDimensions();
    }

  switch( dimension )
    {
    case 2:
      {
      antsJointFusion<2>( parser );
      }
      break;
    case 3:
      {
      antsJointFusion<3>( parser );
      }
      break;
    case 4:
      {
      antsJointFusion<4>( parser );
      }
      break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
