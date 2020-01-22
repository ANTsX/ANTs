#include "antsUtilities.h"
#include "antsAllocImage.h"
#include "itkantsRegistrationHelper.h"
#include "ReadWriteData.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

#include "itkAffineTransform.h"
#include "itkCompositeTransform.h"
#include "itkDisplacementFieldTransform.h"
#include "itkIdentityTransform.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkTransformFactory.h"
#include "itkTransformFileWriter.h"
#include "itkTransformToDisplacementFieldFilter.h"

#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkGaussianInterpolateImageFunction.h"
#include "itkInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkLabelImageGaussianInterpolateImageFunction.h"

namespace ants
{
template <unsigned int Dimension>
int antsAlignOriginImplementation( itk::ants::CommandLineParser::Pointer & parser, unsigned int inputImageType )
{
  if( inputImageType != 0 )
    {
    std::cerr << "inputImageType is not used, therefore only mode 0 is supported at the momemnt." << std::endl;
    return EXIT_FAILURE;
    }
  typedef double                           PixelType;

  typedef itk::Image<PixelType, Dimension> ImageType;

  typename ImageType::Pointer inputImage;
  typename ImageType::Pointer outputImage;

  /**
   * Input object option - for now, we're limiting this to images.
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer inputOption = parser->GetOption( "input" );
  typename itk::ants::CommandLineParser::OptionType::Pointer outputOption = parser->GetOption( "output" );

  if( inputOption && inputOption->GetNumberOfFunctions() > 0 )
    {
    if( inputOption->GetFunction()->GetNumberOfParameters() > 1 &&
        parser->Convert<unsigned int>( outputOption->GetFunction( 0 )->GetParameter( 1 ) ) == 0 )
      {
      std::cerr << "An input image is required." << std::endl;
      return EXIT_FAILURE;
      }

    std::cout << "Input image: " << inputOption->GetFunction()->GetName() << std::endl;

    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( ( inputOption->GetFunction()->GetName() ).c_str() );
    reader->Update();
    inputImage = reader->GetOutput();
    }

  std::string outputTransform;
  std::string outputWarpedImageName;

  if( outputOption && outputOption->GetNumberOfFunctions() > 0 )
    {
    outputTransform = outputOption->GetFunction( 0 )->GetName();
    if( outputOption->GetFunction()->GetNumberOfParameters() > 0 )
      {
      outputTransform = outputOption->GetFunction( 0 )->GetParameter( 0 );
      }

    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
      {
      outputWarpedImageName = outputOption->GetFunction( 0 )->GetParameter( 1 );
      }

    std::cout << "Output transform: " << outputTransform << std::endl;
    std::cout << "Output image: " << outputWarpedImageName << std::endl;
    }

  /**
    * Reference image option
    */

  // read in the image as char since we only need the header information.
  typedef itk::Image<char, Dimension> ReferenceImageType;
  typename ReferenceImageType::Pointer referenceImage;

  typename itk::ants::CommandLineParser::OptionType::Pointer referenceOption =
    parser->GetOption( "reference-image" );
  if( referenceOption && referenceOption->GetNumberOfFunctions() > 0 )
    {
    std::cout << "Reference image: " << referenceOption->GetFunction()->GetName() << std::endl;

    // read in the image as char since we only need the header information.
    typedef itk::ImageFileReader<ReferenceImageType> ReferenceReaderType;
    typename ReferenceReaderType::Pointer referenceReader =
      ReferenceReaderType::New();
    referenceReader->SetFileName( ( referenceOption->GetFunction()->GetName() ).c_str() );

    referenceImage = referenceReader->GetOutput();
    referenceImage->Update();
    referenceImage->DisconnectPipeline();
    }
  else
    {
    std::cerr << "Error:  No reference image specified." << std::endl;
    return EXIT_FAILURE;
    }

  typename ImageType::PointType::VectorType translation = inputImage->GetOrigin() - referenceImage->GetOrigin();
  translation = referenceImage->GetDirection() * inputImage->GetDirection() * translation;
  std::cout << "offset = " << translation << std::endl;

  typedef itk::MatrixOffsetTransformBase<double, Dimension, Dimension> TransformType;
  typename TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();
  transform->SetTranslation( translation );

  typename itk::TransformFileWriter::Pointer transform_writer = itk::TransformFileWriter::New();
  transform_writer->SetFileName( outputTransform );
  transform_writer->SetInput( transform );
#if ITK_VERSION_MAJOR >= 5
  transform_writer->SetUseCompression(true);
#endif
  transform_writer->Update();

  return EXIT_SUCCESS;
}

static void InitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{

    {
    std::string description =
      std::string( "This option forces the image to be treated as a specified-" )
      + std::string( "dimensional image.  If not specified, antsWarp tries to " )
      + std::string( "infer the dimensionality from the input image." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "dimensionality" );
    option->SetShortName( 'd' );
    option->SetUsageOption( 0, "2/3" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Currently, the only input objects supported are image " )
      + std::string( "objects.  However, the current framework allows for " )
      + std::string( "warping of other objects such as meshes and point sets. ");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "input" );
    option->SetShortName( 'i' );
    option->SetUsageOption( 0, "inputFileName" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "For warping input images, the reference image defines the " )
      + std::string( "spacing, origin, size, and direction of the output warped " )
      + std::string( "image. ");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "reference-image" );
    option->SetShortName( 'r' );
    option->SetUsageOption( 0, "imageFileName" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "One can either output the warped image or, if the boolean " )
      + std::string( "is set, one can print out the displacement field based on the" )
      + std::string( "composite transform and the reference image." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "output" );
    option->SetShortName( 'o' );
    option->SetUsageOption( 0, "warpedOutputFileName" );
    option->SetUsageOption( 1, "[transform,alignedImage]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Print the help menu (short version)." );

    OptionType::Pointer option = OptionType::New();
    option->SetShortName( 'h' );
    option->SetDescription( description );
    option->AddFunction( std::string( "0" ) );
    parser->AddOption( option );
    }

    {
    std::string description = std::string( "Print the help menu." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "help" );
    option->SetDescription( description );
    option->AddFunction( std::string( "0" ) );
    parser->AddOption( option );
    }
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int antsAlignOrigin( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "antsAlignOrigin" );
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
    std::string( "antsAlignOrigin, applied to an input image, transforms it " )
    + std::string( "according to a reference image and a transform " )
    + std::string( "(or a set of transforms)." );

  parser->SetCommandDescription( commandDescription );
  InitializeCommandLineOptions( parser );

  if( parser->Parse( argc, argv ) == EXIT_FAILURE )
    {
    return EXIT_FAILURE;
    }

  if( argc < 2 || ( parser->GetOption( "help" ) &&
                    ( parser->Convert<bool>( parser->GetOption( "help" )->GetFunction()->GetName() ) ) ) )
    {
    parser->PrintMenu( std::cout, 5, false );
    if( argc < 2 )
      {
      return EXIT_FAILURE;
      }
    return EXIT_SUCCESS;
    }
  else if( parser->GetOption( 'h' ) &&
           ( parser->Convert<bool>( parser->GetOption( 'h' )->GetFunction()->GetName() ) ) )
    {
    parser->PrintMenu( std::cout, 5, true );
    return EXIT_SUCCESS;
    }

  // Read in the first intensity image to get the image dimension.
  std::string filename;

  itk::ants::CommandLineParser::OptionType::Pointer inputOption =
    parser->GetOption( "reference-image" );
  if( inputOption && inputOption->GetNumberOfFunctions() > 0 )
    {
    if( inputOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      filename = inputOption->GetFunction( 0 )->GetParameter( 0 );
      }
    else
      {
      filename = inputOption->GetFunction( 0 )->GetName();
      }
    }
  else
    {
    std::cerr << "No reference image was specified." << std::endl;
    return EXIT_FAILURE;
    }

  itk::ants::CommandLineParser::OptionType::Pointer inputImageTypeOption =
    parser->GetOption( "input-image-type" );

  itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
      filename.c_str(), itk::ImageIOFactory::FileModeEnum::ReadMode );
  unsigned int dimension = imageIO->GetNumberOfDimensions();

  itk::ants::CommandLineParser::OptionType::Pointer dimOption =
    parser->GetOption( "dimensionality" );
  if( dimOption && dimOption->GetNumberOfFunctions() > 0 )
    {
    dimension = parser->Convert<unsigned int>( dimOption->GetFunction()->GetName() );
    }

  switch( dimension )
    {
    case 2:
      {
      if( inputImageTypeOption )
        {
        std::string inputImageType = inputImageTypeOption->GetFunction()->GetName();

        if( !std::strcmp( inputImageType.c_str(), "scalar" ) || !std::strcmp( inputImageType.c_str(), "0" ) )
          {
          return antsAlignOriginImplementation<2>( parser, 0 );
          }
        else if( !std::strcmp( inputImageType.c_str(), "vector" ) || !std::strcmp( inputImageType.c_str(), "1" ) )
          {
          return antsAlignOriginImplementation<2>( parser, 1 );
          }
        else if( !std::strcmp( inputImageType.c_str(), "tensor" ) || !std::strcmp( inputImageType.c_str(), "2" ) )
          {
          std::cerr << "antsApplyTransforms is not implemented for 2-D tensor images." << std::endl;
          }
        else
          {
          std::cerr << "Unrecognized input image type (cf --input-image-type option)." << std::endl;
          return EXIT_FAILURE;
          }
        }
      else
        {
        return antsAlignOriginImplementation<2>( parser, 0 );
        }
      }
      break;
    case 3:
      {
      if( inputImageTypeOption )
        {
        std::string inputImageType = inputImageTypeOption->GetFunction()->GetName();

        if( !std::strcmp( inputImageType.c_str(), "scalar" ) || !std::strcmp( inputImageType.c_str(), "0" ) )
          {
          return antsAlignOriginImplementation<3>( parser, 0 );
          }
        else if( !std::strcmp( inputImageType.c_str(), "vector" ) || !std::strcmp( inputImageType.c_str(), "1" ) )
          {
          return antsAlignOriginImplementation<3>( parser, 1 );
          }
        else if( !std::strcmp( inputImageType.c_str(), "tensor" ) || !std::strcmp( inputImageType.c_str(), "2" ) )
          {
          return antsAlignOriginImplementation<3>( parser, 2 );
          }
        else
          {
          std::cerr << "Unrecognized input image type (cf --input-image-type option)." << std::endl;
          return EXIT_FAILURE;
          }
        }
      else
        {
        return antsAlignOriginImplementation<3>( parser, 0 );
        }
      }
      break;
    case 4:
      {
      if( inputImageTypeOption )
        {
        std::string inputImageType = inputImageTypeOption->GetFunction()->GetName();

        if( !std::strcmp( inputImageType.c_str(), "scalar" ) || !std::strcmp( inputImageType.c_str(), "0" ) )
          {
          return antsAlignOriginImplementation<4>( parser, 0 );
          }
        else if( !std::strcmp( inputImageType.c_str(), "vector" ) || !std::strcmp( inputImageType.c_str(), "1" ) )
          {
          return antsAlignOriginImplementation<4>( parser, 1 );
          }
        else if( !std::strcmp( inputImageType.c_str(), "tensor" ) || !std::strcmp( inputImageType.c_str(), "2" ) )
          {
          std::cerr << "antsApplyTransforms is not implemented for 4-D tensor images." << std::endl;
          }
        else
          {
          std::cerr << "Unrecognized input image type (cf --input-image-type option)." << std::endl;
          return EXIT_FAILURE;
          }
        }
      else
        {
        return antsAlignOriginImplementation<3>( parser, 0 );
        }
      }
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
