
#include "antsUtilities.h"
#include <algorithm>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "antsAtroposSegmentationImageFilter.h"
#include "antsCommandLineOption.h"
#include "antsCommandLineParser.h"
#include "antsJointHistogramParzenShapeAndOrientationListSampleFunction.h"

#include <string>
#include <algorithm>
#include <vector>

namespace ants
{
template <unsigned int ImageDimension>
int AtroposSegmentation( itk::ants::CommandLineParser *parser )
{
  typedef float                                 PixelType;
  typedef float                                 RealType;
  typedef itk::Image<PixelType, ImageDimension> InputImageType;

  typedef unsigned int                          LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;

  typedef  itk::ants::AtroposSegmentationImageFilter
    <InputImageType, LabelImageType> SegmentationFilterType;
  typename SegmentationFilterType::Pointer segmenter
    = SegmentationFilterType::New();

  /**
   * mask image
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer maskOption =
    parser->GetOption( "mask-image" );
  if( maskOption && maskOption->GetNumberOfFunctions() > 0 )
    {
    try
      {
      typedef  itk::ImageFileReader<LabelImageType> ReaderType;
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( ( maskOption->GetFunction( 0 )->GetName() ).c_str() );
      reader->Update();

      segmenter->SetMaskImage( reader->GetOutput() );
      }
    catch( ... )
      {
      }
    }
  else
    {
    std::cout << "An image mask is required.  Specify a mask image"
             << " with the -x option." << std::endl;
    return EXIT_FAILURE;
    }

  /**
   * intensity images
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer imageOption =
    parser->GetOption( "intensity-image" );
  if( imageOption && imageOption->GetNumberOfFunctions() > 0 )
    {
    unsigned int count = 0;
    for( int n = imageOption->GetNumberOfFunctions() - 1; n >= 0; n-- )
      {
      typedef itk::ImageFileReader<InputImageType> ReaderType;
      typename ReaderType::Pointer reader = ReaderType::New();
      if( imageOption->GetFunction( n )->GetNumberOfParameters() > 0 )
        {
        reader->SetFileName( imageOption->GetFunction( n )->GetParameter( 0 ) );
        }
      else
        {
        reader->SetFileName( imageOption->GetFunction( n )->GetName() );
        }
      reader->Update();

      segmenter->SetIntensityImage( count, reader->GetOutput() );
      if( imageOption->GetFunction( count )->GetNumberOfParameters() > 1 )
        {
        segmenter->SetAdaptiveSmoothingWeight( count, parser->Convert<float>(
                                                 imageOption->GetFunction( count )->GetParameter( 1 ) ) );
        }
      else
        {
        segmenter->SetAdaptiveSmoothingWeight( count, 0.0 );
        }
      count++;
      }
    }
  else
    {
    std::cout << "No input images were specified.  Specify an input image"
             << " with the -a option." << std::endl;
    return EXIT_FAILURE;
    }

  /**
   * likelihood
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer likelihoodOption =
    parser->GetOption( "likelihood-model" );
  if( likelihoodOption && likelihoodOption->GetNumberOfFunctions() > 0 )
    {
    std::string likelihoodModel = likelihoodOption->GetFunction( 0 )->GetName();
    ConvertToLowerCase( likelihoodModel );
    if( !likelihoodModel.compare( std::string( "jointshapeandorientationprobability" ) ) )
      {
      if( segmenter->GetNumberOfIntensityImages() !=
          static_cast<unsigned int>( ImageDimension * ( ImageDimension + 1 ) / 2 ) )
        {
        std::cout << " Expect images in upper triangular order " << std::endl;
        std::cout << " xx xy xz yy yz zz " << std::endl;
        std::cout << "Incorrect number of intensity images specified." << std::endl;
        return EXIT_FAILURE;
        }
      typedef typename SegmentationFilterType::SampleType SampleType;
      typedef itk::ants::Statistics::
        JointHistogramParzenShapeAndOrientationListSampleFunction
        <SampleType, float, float> LikelihoodType;

      float shapeSigma = 1.0;
      if( likelihoodOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
        {
        shapeSigma = parser->Convert<float>(
            likelihoodOption->GetFunction( 0 )->GetParameter( 0 ) );
        }
      unsigned int numberOfShapeBins = 32;
      if( likelihoodOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        numberOfShapeBins = parser->Convert<unsigned int>(
            likelihoodOption->GetFunction( 0 )->GetParameter( 1 ) );
        }
      float orientationSigma = 2.0;
      if( likelihoodOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
        {
        orientationSigma = parser->Convert<float>(
            likelihoodOption->GetFunction( 0 )->GetParameter( 2 ) );
        }
      unsigned int numberOfOrientationBins = 64;
      if( likelihoodOption->GetFunction( 0 )->GetNumberOfParameters() > 3 )
        {
        numberOfOrientationBins = parser->Convert<unsigned int>(
            likelihoodOption->GetFunction( 0 )->GetParameter(3) );
        }
      for( unsigned int n = 0; n < segmenter->GetNumberOfTissueClasses(); n++ )
        {
        typename LikelihoodType::Pointer hpwLikelihood =
          LikelihoodType::New();
        hpwLikelihood->SetShapeSigma( shapeSigma );
        hpwLikelihood->SetOrientationSigma( orientationSigma);
        hpwLikelihood->SetNumberOfShapeJointHistogramBins( numberOfShapeBins );
        hpwLikelihood->SetNumberOfOrientationJointHistogramBins( numberOfOrientationBins);
        segmenter->SetLikelihoodFunction( n, hpwLikelihood );
        }
      }
    else
      {
      std::cout << "Unrecognized likelihood model request." << std::endl;
      return EXIT_FAILURE;
      }
    }

//  std::cout << std::endl << "Writing output:" << std::endl;
//  typename itk::ants::CommandLineParser::OptionType::Pointer outputOption =
//    parser->GetOption( "output" );
//  if( outputOption && outputOption->GetNumberOfFunctions() > 0 )
//    {
//    typedef  itk::ImageFileWriter<ImageType> WriterType;
//    typename WriterType::Pointer writer = WriterType::New();
//    writer->SetInput( segmenter->GetOutput() );
//    writer->SetFileName( ( outputOption->GetFunction( 0 )->GetName() ).c_str() );
//    writer->Update();
//    }

  std::cout << std::endl;
  segmenter->Print( std::cout, 2 );

  return EXIT_SUCCESS;
}

void InitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

    {
    std::string description =
      std::string( "This option forces the image to be treated as a specified-" )
      + std::string( "dimensional image.  If not specified, Atropos tries to " )
      + std::string( "infer the dimensionality from the first input image." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "image-dimensionality" );
    option->SetShortName( 'd' );
    option->SetUsageOption( 0, "2/3/4" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "One or more scalar images is specified for segmentation " )
      + std::string( "using the -a/--intensity-image option.  For segmentation " )
      + std::string( "scenarios with no prior information, the first scalar " )
      + std::string( "image encountered on the command line is used to order " )
      + std::string( "labelings such that the class with the smallest intensity " )
      + std::string( "signature is class \'1\' through class \'N\' which represents " )
      + std::string( "the voxels with the largest intensity values.  The " )
      + std::string( "optional adaptive smoothing weight parameter is applicable " )
      + std::string( "only when using prior label or probability images.  This " )
      + std::string( "scalar parameter is to be specified between [0,1] which " )
      + std::string( "smooths each labeled region separately and modulates the " )
      + std::string( "intensity measurement at each voxel in each intensity image " )
      + std::string( "between the original intensity and its smoothed " )
      + std::string( "counterpart.  The smoothness parameters are governed by the " )
      + std::string( "-b/--bspline option." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "intensity-image" );
    option->SetShortName( 'a' );
    option->SetUsageOption( 0, "[intensityImage,<adaptiveSmoothingWeight>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "The image mask (which is required) defines the region which " )
      + std::string( "is to be labeled by the Atropos algorithm." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "mask-image" );
    option->SetShortName( 'x' );
    option->SetUsageOption( 0, "maskImageFilename" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Both parametric and non-parametric options exist in Atropos. " )
      + std::string( "The Gaussian parametric option is commonly used " )
      + std::string( "(e.g. SPM & FAST) where the mean and standard deviation " )
      + std::string( "for the Gaussian of each class is calculated at each " )
      + std::string( "iteration.  Other groups use non-parametric approaches " )
      + std::string( "exemplified by option 2.  We recommend using options 1 " )
      + std::string( "or 2 as they are fairly standard and the " )
      + std::string( "default parameters work adequately." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "likelihood-model" );
    option->SetShortName( 'k' );
    option->SetUsageOption( 0, "JointShapeAndOrientationProbability[<sigma=1.0>,<numberOfBins=32>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "The output consists of a labeled image where each voxel " )
      + std::string( "in the masked region is assigned a label from 1, 2, " )
      + std::string( "..., N.  Optionally, one can also output the posterior " )
      + std::string( "probability images specified in the same format as the " )
      + std::string( "prior probability images, e.g. posterior%02d.nii.gz " )
      + std::string( "(C-style file name formatting)." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "output" );
    option->SetShortName( 'o' );
    option->SetUsageOption( 0, "[classifiedImage,<posteriorProbabilityImageFileNameFormat>]" );
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
int AtroposMin( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "AtroposMin" );
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
    std::string( "A finite mixture modeling (FMM) segmentation approach " )
    + std::string( "with possibilities for specifying prior constraints. " )
    + std::string( "These prior constraints include the specification " )
    + std::string( "of a prior label image, prior probability images " )
    + std::string( "(one for each class), and/or an MRF prior to " )
    + std::string( "enforce spatial smoothing of the labels.  Similar algorithms " )
    + std::string( "include FAST and SPM.  " );

  parser->SetCommandDescription( commandDescription );
  InitializeCommandLineOptions( parser );

  parser->Parse( argc, argv );

  if( argc < 2 || parser->Convert<bool>(
        parser->GetOption( "help" )->GetFunction()->GetName() ) )
    {
    parser->PrintMenu( std::cout, 5, false );
    if( argc < 2 )
      {
      return EXIT_FAILURE;
      }
    return EXIT_SUCCESS;
    }
  else if( parser->GetOption( 'h' ) &&
           parser->Convert<bool>( parser->GetOption( 'h' )->GetFunction()->GetName() ) )
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
      parser->GetOption( "intensity-image" );
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
      std::cout << "No input images were specified.  Specify an input image"
               << " with the -a option" << std::endl;
      return EXIT_FAILURE;
      }
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
        filename.c_str(), itk::ImageIOFactory::ReadMode );
    dimension = imageIO->GetNumberOfDimensions();
    }

  std::cout << std::endl << "Running Atropos for "
           << dimension << "-dimensional images." << std::endl;

  switch( dimension )
    {
    case 2:
      {
      AtroposSegmentation<2>( parser );
      }
      break;
    case 3:
      {
      AtroposSegmentation<3>( parser );
      }
      break;
    case 4:
      {
      AtroposSegmentation<4>( parser );
      }
      break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
