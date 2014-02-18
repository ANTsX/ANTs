#include "antsCommandLineParser.h"
#include "antsUtilities.h"
#include "ReadWriteImage.h"

#include <algorithm>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkConstantPadImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkTileImageFilter.h"

namespace ants
{
int CreateMosaic( itk::ants::CommandLineParser *parser )
{
  const unsigned int ImageDimension = 3;

  typedef float                                   PixelType;
  typedef float                                   RealType;
  typedef itk::Image<PixelType, ImageDimension>   ImageType;
  typedef itk::Image<PixelType, ImageDimension-1> SliceType;

  typedef itk::RGBPixel<unsigned char> RGBPixelType;
//  typedef itk::RGBAPixel<unsigned char> RGBPixelType;
  typedef itk::Image<RGBPixelType, ImageDimension> RGBImageType;

  // Read in input image

  ImageType::Pointer inputImage = NULL;

  itk::ants::CommandLineParser::OptionType::Pointer inputImageOption =
    parser->GetOption( "input-image" );
  if( inputImageOption && inputImageOption->GetNumberOfFunctions() )
    {
    std::string inputFile = inputImageOption->GetFunction( 0 )->GetName();
    ReadImage<ImageType>( inputImage, inputFile.c_str() );
    inputImage->Update();
    inputImage->DisconnectPipeline();
    }
  else
    {
    std::cout << "Input image not specified." << std::endl;
    return EXIT_FAILURE;
    }
  ImageType::SpacingType spacing = inputImage->GetSpacing();
  ImageType::SizeType size = inputImage->GetRequestedRegion().GetSize();

  // Get direction.  If not specified, pick direction with coarsest spacing.

  unsigned int direction = 0;

  itk::ants::CommandLineParser::OptionType::Pointer directionOption =
    parser->GetOption( "direction" );
  if( directionOption && directionOption->GetNumberOfFunctions() )
    {
    direction = parser->Convert<unsigned int>( directionOption->GetFunction( 0 )->GetName() );
    if( direction > 2 )
      {
      std::cerr << "The direction value must be 0, 1, or 2." << std::endl;
      return EXIT_FAILURE;
      }
    }
  else
    {
    float maxSpacing = spacing[0];
    unsigned int maxIndex = 0;
    for( unsigned int d = 1; d < ImageDimension; d++ )
      {
      if( spacing[d] > maxSpacing )
        {
        maxSpacing = spacing[d];
        maxIndex = d;
        }
      }
    direction = maxIndex;
    }

  // Get padding/cropping options.

  int paddingType = 0;
  RealType padValue = 0;

  unsigned long lowerBound[ImageDimension];
  unsigned long upperBound[ImageDimension];

  SliceType::RegionType              croppedSliceRegion;
  SliceType::RegionType::SizeType    croppedSliceSize;
  SliceType::RegionType::IndexType   croppedSliceIndex;

  itk::ants::CommandLineParser::OptionType::Pointer paddingOption =
    parser->GetOption( "pad-or-crop" );
  if( paddingOption && paddingOption->GetNumberOfFunctions() )
    {
    if( paddingOption->GetFunction( 0 )->GetNumberOfParameters() == 3 )
      {
      std::vector<int> lowerBoundVector = parser->ConvertVector<int>(
        paddingOption->GetFunction( 0 )->GetParameter( 0 ) );
      std::vector<int> upperBoundVector = parser->ConvertVector<int>(
        paddingOption->GetFunction( 0 )->GetParameter( 1 ) );

      if( lowerBoundVector.size() != 2 || upperBoundVector.size() != 2 )
        {
        std::cerr << "Incorrect padding specification." << std::endl;
        return EXIT_FAILURE;
        }
      int lowerBoundProduct = lowerBoundVector[0] * lowerBoundVector[1];
      int upperBoundProduct = upperBoundVector[0] * upperBoundVector[1];

      if( lowerBoundProduct < 0 || upperBoundProduct < 0 || upperBoundProduct * lowerBoundProduct < 0 )
        {
        std::cerr << "Current capabilities do not include mixing of cropping and padding,"
                  << " i.e. negative and positive pad values, respectively"  << std::endl;
        return EXIT_FAILURE;
        }

      if( lowerBoundVector[0] < 0 )
        {
        paddingType = -1;
        unsigned int count = 0;
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          if( d != direction )
            {
            croppedSliceSize[count] = size[d] - ( vnl_math_abs( lowerBoundVector[count] )
              + vnl_math_abs( upperBoundVector[count] ) );
            croppedSliceIndex[count] = vnl_math_abs( lowerBoundVector[count] );
            count++;
            }
          }
        croppedSliceRegion.SetSize( croppedSliceSize );
        croppedSliceRegion.SetIndex( croppedSliceIndex );
        }
      else
        {
        paddingType = 1;
        for( unsigned int d = 0; d < ImageDimension - 1; d++ )
          {
          lowerBound[d] = lowerBoundVector[d];
          upperBound[d] = upperBoundVector[d];
          }
        }
      }
    else
      {
      int padWidth = 0;
      if( paddingOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
        {
        padWidth = parser->Convert<int>( paddingOption->GetFunction( 0 )->GetName() );
        }
      else if( paddingOption->GetFunction( 0 )->GetNumberOfParameters() <= 2 )
        {
        padWidth = parser->Convert<int>( paddingOption->GetFunction( 0 )->GetParameter( 0 ) );
        padValue = parser->Convert<int>( paddingOption->GetFunction( 0 )->GetParameter( 1 ) );
        }
      if( padWidth < 0 )
        {
        paddingType = -1;
        unsigned int count = 0;
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          if( d != direction )
            {
            croppedSliceSize[count] = size[d] - 2 * vnl_math_abs( padWidth );
            croppedSliceIndex[count] = vnl_math_abs( padWidth );
            count++;
            }
          }
        croppedSliceRegion.SetSize( croppedSliceSize );
        croppedSliceRegion.SetIndex( croppedSliceIndex );
        }
      else
        {
        paddingType = 1;
        for( unsigned int d = 0; d < ImageDimension - 1; d++ )
          {
          lowerBound[d] = padWidth;
          upperBound[d] = padWidth;
          }
        }
      }
    }

  // Get the slices

  std::vector<unsigned int> whichSlices;
  for( unsigned int n = 0; n < size[direction]; n++ )
    {
    whichSlices.push_back( n );
    }

  itk::ants::CommandLineParser::OptionType::Pointer slicesOption =
    parser->GetOption( "slices" );
  if( slicesOption && slicesOption->GetNumberOfFunctions() )
    {
    unsigned int numberOfSlicesToIncrement = 1;
    unsigned int startingSlice = 0;
    unsigned int endSlice = size[direction] - 1;

    bool readSlices = false;

    if( slicesOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
      {
      whichSlices = parser->ConvertVector<unsigned int>( slicesOption->GetFunction( 0 )->GetName() );
      if( whichSlices.size() == 1 )
        {
        numberOfSlicesToIncrement = whichSlices[0];
        }
      else
        {
        readSlices = true;
        }
      }

    if( !readSlices )
      {
      if( slicesOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
        {
        numberOfSlicesToIncrement = parser->Convert<unsigned int>( slicesOption->GetFunction( 0 )->GetParameter( 0 ) );
        }
      if( numberOfSlicesToIncrement == 0 )
        {
        std::cerr << "Need greater than 0 slices for incrementing." << std::endl;
        return EXIT_FAILURE;
        }

      if( slicesOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        startingSlice = parser->Convert<unsigned int>( slicesOption->GetFunction( 0 )->GetParameter( 1 ) );
        }
      if( slicesOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
        {
        endSlice = parser->Convert<unsigned int>( slicesOption->GetFunction( 0 )->GetParameter( 2 ) );
        }

      whichSlices.clear();
      for( unsigned int n = startingSlice; n < endSlice; n += numberOfSlicesToIncrement )
        {
        whichSlices.push_back( n );
        }
      }
    }

  // Get tile geometry.

  unsigned long numberOfSlices = whichSlices.size();

  int numberOfRows = 0;
  int numberOfColumns = 0;

  itk::ants::CommandLineParser::OptionType::Pointer tileGeometryOption =
    parser->GetOption( "tile-geometry" );
  if( tileGeometryOption && tileGeometryOption->GetNumberOfFunctions() )
    {
    std::vector<int> layout = parser->ConvertVector<int>( tileGeometryOption->GetFunction( 0 )->GetName() );
    if( layout.size() > 2 )
      {
      std::cerr << "Tile geometry is specified as numberOfRowsxnumberOfColumns" << std::endl;
      return EXIT_FAILURE;
      }
    numberOfRows = vnl_math_min( static_cast<int>( layout[0] ), static_cast<int>( numberOfSlices ) );
    numberOfColumns = vnl_math_min( static_cast<int>( layout[1] ), static_cast<int>( numberOfSlices ) );
    }

  if( numberOfRows <= 0 && numberOfColumns > 0 )
    {
    numberOfRows = vcl_ceil( static_cast<float>( numberOfSlices ) / static_cast<float>( numberOfColumns ) );
    }
  else if( numberOfColumns <= 0 && numberOfRows > 0 )
    {
    numberOfColumns = vcl_ceil( static_cast<float>( numberOfSlices ) / static_cast<float>( numberOfRows ) );
    }
  else if( numberOfColumns <= 0 && numberOfRows <= 0 )
    {
    numberOfRows = static_cast<int>( vcl_sqrt( static_cast<float>( numberOfSlices ) ) );
    numberOfColumns = vcl_ceil( static_cast<float>( numberOfSlices ) / static_cast<float>( numberOfRows ) );
    }

  // Read in optional RGB image

  RGBImageType::Pointer rgbImage = NULL;

  itk::ants::CommandLineParser::OptionType::Pointer rgbImageOption =
    parser->GetOption( "rgb-image" );
  if( rgbImageOption && rgbImageOption->GetNumberOfFunctions() )
    {
    std::cerr << "Whoa there, Cowboy.  I haven't gotten to this part yet." << std::endl;
    return EXIT_FAILURE;

    std::string rgbFile = rgbImageOption->GetFunction( 0 )->GetName();
    ReadImage<RGBImageType>( rgbImage, rgbFile.c_str() );
    rgbImage->Update();
    rgbImage->DisconnectPipeline();
    }

  RealType alpha = 1.0;

  itk::ants::CommandLineParser::OptionType::Pointer alphaOption =
    parser->GetOption( "alpha" );
  if( alphaOption && alphaOption->GetNumberOfFunctions() )
    {
    alpha = parser->Convert<RealType>( alphaOption->GetFunction( 0 )->GetName() );
    }

  // Now do the tiling

  std::cout << "Slices[" << direction << "]: " << whichSlices.size() << std::endl;
  std::cout << "Rows:  " << numberOfRows << std::endl;
  std::cout << "Columns:  " << numberOfColumns << std::endl;

  typedef itk::TileImageFilter<SliceType, SliceType> FilterType;
  FilterType::LayoutArrayType array;

  array[0] = numberOfColumns;
  array[1] = numberOfRows;

  ImageType::RegionType region;
  size[direction] = 0;

  FilterType::Pointer filter = FilterType::New();
  filter->SetLayout( array );

  for( unsigned int i = 0; i < whichSlices.size(); i++ )
    {
    unsigned int n = whichSlices[i];

    ImageType::IndexType index;
    index.Fill( 0 );
    index[direction] = static_cast<int>( n );
    region.SetIndex( index );
    region.SetSize( size );

    typedef itk::ExtractImageFilter<ImageType, SliceType> ExtracterType;
    ExtracterType::Pointer extracter = ExtracterType::New();
    extracter->SetInput( inputImage );
    extracter->SetExtractionRegion( region );
    extracter->SetDirectionCollapseToIdentity();

    SliceType::Pointer outputSlice = NULL;

    if( paddingType == -1 )
      {
      typedef itk::ExtractImageFilter<SliceType, SliceType> ExtracterType2;
      ExtracterType2::Pointer extracter2 = ExtracterType2::New();
      extracter2->SetInput( extracter->GetOutput() );
      extracter2->SetExtractionRegion( croppedSliceRegion );
      extracter2->SetDirectionCollapseToIdentity();

      outputSlice = extracter2->GetOutput();
      outputSlice->Update();
      outputSlice->DisconnectPipeline();
      }
    else if( paddingType == 1 )
      {
      typedef itk::ConstantPadImageFilter<SliceType, SliceType> PadderType;
      PadderType::Pointer padder = PadderType::New();
      padder->SetInput( extracter->GetOutput() );
      padder->SetPadLowerBound( lowerBound );
      padder->SetPadUpperBound( upperBound );
      padder->SetConstant( static_cast<PixelType>( padValue ) );

      outputSlice = padder->GetOutput();
      outputSlice->Update();
      outputSlice->DisconnectPipeline();
      }
    else // paddingType == 0
      {
      outputSlice = extracter->GetOutput();
      outputSlice->Update();
      outputSlice->DisconnectPipeline();
      }

    filter->SetInput( i, outputSlice );
    }
  filter->Update();

  itk::ants::CommandLineParser::OptionType::Pointer outputOption =
    parser->GetOption( "output" );
  if( outputOption && outputOption->GetNumberOfFunctions() )
    {
    std::string outputFile = outputOption->GetFunction( 0 )->GetName();
    WriteImage<SliceType>( filter->GetOutput(), outputFile.c_str() );
    }
  else
    {
    std::cerr << "No output filename specified." << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

void InitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

    {
    std::string description =
      std::string( "Main input is a 3-D grayscale image.  " );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "input-image" );
    option->SetShortName( 'i' );
    option->SetUsageOption( 0, "inputImageFilename" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "An optional RGB image can be added as an overlay.  ")
      + std::string( "It must have the same image geometry as the input " )
      + std::string( "grayscale image." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "rgb-image" );
    option->SetShortName( 'r' );
    option->SetUsageOption( 0, "rgbImageFilename" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "If an RGB image is provided, render the overlay using the specified " )
      + std::string( "alpha parameter." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "alpha" );
    option->SetShortName( 'a' );
    option->SetUsageOption( 0, "value" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "The output consists of the tiled mosaic image." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "output" );
    option->SetShortName( 'o' );
    option->SetUsageOption( 0, "tiledMosaicImage" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "The tile geometry specifies the number of rows and columns " )
      + std::string( "in the output image.  For example, if the user specifies " )
      + std::string( "\'5x10\', then 5 rows by 10 columns of slices are rendered. " )
      + std::string( "If R < 0 and C > 0 (or vice versa), the negative value is " )
      + std::string( "selected based on direction." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "tile-geometry" );
    option->SetShortName( 't' );
    option->SetUsageOption( 0, "RxC" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Specifies the direction of the slices.  If no direction is specified, " )
      + std::string( "the direction with the coarsest spacing is chosen." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "direction" );
    option->SetShortName( 'd' );
    option->SetUsageOption( 0, "0/1/2" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "The user can specify whether to pad or crop a specified voxel-width " )
      + std::string( "boundary of each individual slice.  For this program, cropping is " )
      + std::string( "simply padding with negative voxel-widths.  If one pads (+), the " )
      + std::string( "user can also specify a constant pad value (default = 0). " );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "pad-or-crop" );
    option->SetShortName( 'p' );
    option->SetUsageOption( 0, "padVoxelWidth" );
    option->SetUsageOption( 1, "[padVoxelWidth,<constantValue=0>]" );
    option->SetUsageOption( 2, "[lowerPadding[0]xlowerPadding[1],upperPadding[0]xupperPadding[1],constantValue]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "This option gives the user more control over what slices " )
      + std::string( "to use for rendering.  The user can specify specific slices " )
      + std::string( "for a particular order.  Alternatively the user can specify " )
      + std::string( "the number slices to skip with the optional specification of " )
      + std::string( "which slices to start and end the sequence." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "slices" );
    option->SetShortName( 's' );
    option->SetUsageOption( 0, "Slice1xSlice2xSlice3..." );
    option->SetUsageOption( 1, "numberOfSlicesToIncrement" );
    option->SetUsageOption( 2, "[numberOfSlicesToIncrement,<startingSlice=0>,<endingSlice=lastSlice>]" );
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
int CreateTiledMosaic( std::vector<std::string> args, std::ostream* /*out_stream = NULL */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "CreateTiledMosaic" );

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
    std::string( "Render a 3-D image volume with optional RGB overlay." );

  parser->SetCommandDescription( commandDescription );
  InitializeCommandLineOptions( parser );

  parser->Parse( argc, argv );

  if( argc < 2 || parser->Convert<bool>( parser->GetOption( "help" )->GetFunction( 0 )->GetName() ) )
    {
    parser->PrintMenu( std::cout, 5, false );
    if( argc < 2 )
      {
      return EXIT_FAILURE;
      }
    return EXIT_SUCCESS;
    }
  else if( parser->Convert<bool>( parser->GetOption( 'h' )->GetFunction( 0 )->GetName() ) )
    {
    parser->PrintMenu( std::cout, 5, true );
    return EXIT_SUCCESS;
    }

  // Get dimensionality

  std::string filename;

  itk::ants::CommandLineParser::OptionType::Pointer imageOption =
    parser->GetOption( "input-image" );
  if( imageOption && imageOption->GetNumberOfFunctions() > 0 )
    {
    filename = imageOption->GetFunction( 0 )->GetName();

    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
        filename.c_str(), itk::ImageIOFactory::ReadMode );
    unsigned int dimension = imageIO->GetNumberOfDimensions();

    if( dimension == 3 )
      {
      CreateMosaic( parser );
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
