#include "antsCommandLineParser.h"
#include "antsUtilities.h"
#include "ReadWriteData.h"

#include <algorithm>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkCastImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkTileImageFilter.h"

namespace ants
{
int CreateMosaic( itk::ants::CommandLineParser *parser )
{
  const unsigned int ImageDimension = 3;

  typedef float                                      RealType;

  typedef RealType                                   PixelType;
  typedef unsigned char                              RgbComponentType;
  typedef itk::RGBPixel<RgbComponentType>            RgbPixelType;

  typedef itk::Image<PixelType, ImageDimension>      ImageType;
  typedef itk::Image<PixelType, ImageDimension-1>    SliceType;
  typedef itk::Image<RgbPixelType, ImageDimension-1> RgbSliceType;

  typedef itk::Image<RgbPixelType, ImageDimension> RgbImageType;

  // Read in input image

  ImageType::Pointer inputImage = NULL;

  itk::ants::CommandLineParser::OptionType::Pointer inputImageOption =
    parser->GetOption( "input-image" );
  if( inputImageOption && inputImageOption->GetNumberOfFunctions() )
    {
    std::string inputFile = inputImageOption->GetFunction( 0 )->GetName();
    ReadImage<ImageType>( inputImage, inputFile.c_str() );
    }
  else
    {
    std::cout << "Input image not specified." << std::endl;
    return EXIT_FAILURE;
    }
  ImageType::SpacingType spacing = inputImage->GetSpacing();
  ImageType::SizeType size = inputImage->GetRequestedRegion().GetSize();

  // Read in optional mask image

  ImageType::Pointer maskImage = NULL;
  ImageType::RegionType maskRegion;

  itk::ants::CommandLineParser::OptionType::Pointer maskImageOption =
    parser->GetOption( "mask-image" );
  if( maskImageOption && maskImageOption->GetNumberOfFunctions() )
    {
    std::string maskFile = maskImageOption->GetFunction( 0 )->GetName();
    ReadImage<ImageType>( maskImage, maskFile.c_str() );

    typedef itk::Image<unsigned short, ImageDimension>      ShortImageType;
    typedef itk::CastImageFilter<ImageType, ShortImageType> CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput( maskImage );
    caster->Update();

    typedef itk::LabelStatisticsImageFilter<ShortImageType, ShortImageType> StatsFilterType;
    StatsFilterType::Pointer stats = StatsFilterType::New();
    stats->SetLabelInput( caster->GetOutput() );
    stats->SetInput( caster->GetOutput() );
    stats->Update();

    maskRegion = stats->GetRegion( 1 );
    }

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
      std::string padWidthString;
      if( paddingOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
        {
        padWidthString = paddingOption->GetFunction( 0 )->GetName();
        }
      else if( paddingOption->GetFunction( 0 )->GetNumberOfParameters() <= 2 )
        {
        padWidthString = paddingOption->GetFunction( 0 )->GetParameter( 0 );
        padValue = parser->Convert<int>( paddingOption->GetFunction( 0 )->GetParameter( 1 ) );
        }

      if( padWidthString.find( std::string( "mask" ) ) != std::string::npos )
        {
        if( !maskImage )
          {
          std::cerr << "Mask image is not specified." << std::endl;
          return EXIT_FAILURE;
          }

        int offset = 0;
        if( padWidthString.find( std::string( "+" ) ) != std::string::npos )
          {
          std::string offsetString = padWidthString.substr( padWidthString.find( std::string( "+" ) ) + 1 );
          offset = parser->Convert<int>( offsetString );
          }
        if( padWidthString.find( std::string( "-" ) ) != std::string::npos )
          {
          std::string offsetString = padWidthString.substr( padWidthString.find( std::string( "-" ) ) + 1 );
          offset = parser->Convert<int>( offsetString );
          offset *= -1;
          }

        paddingType = -1;
        unsigned int count = 0;
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          if( d != direction )
            {
            croppedSliceSize[count] = maskRegion.GetSize()[d] + 2 * offset;
            croppedSliceIndex[count] = maskRegion.GetIndex()[d] - offset;
            count++;
            }
          }
        croppedSliceRegion.SetSize( croppedSliceSize );
        croppedSliceRegion.SetIndex( croppedSliceIndex );
        }
      else
        {
        padWidth = parser->Convert<int>( padWidthString );

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
    int numberOfSlicesToIncrement = 1;
    int startingSlice = 0;
    int endSlice = size[direction] - 1;

    bool readSlices = false;
    bool reverseOrder = false;

    if( slicesOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
      {
      std::vector<int> slicesVector = parser->ConvertVector<int>( slicesOption->GetFunction( 0 )->GetName() );
      if( slicesVector.size() == 1 )
        {
        numberOfSlicesToIncrement = slicesVector[0];
        if( numberOfSlicesToIncrement < 0 )
          {
          reverseOrder = true;
          }
        }
      else
        {
        whichSlices = parser->ConvertVector<unsigned int>( slicesOption->GetFunction( 0 )->GetName() );
        readSlices = true;
        }
      }

    if( !readSlices )
      {
      if( slicesOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
        {
        numberOfSlicesToIncrement = parser->Convert<int>( slicesOption->GetFunction( 0 )->GetParameter( 0 ) );
        if( numberOfSlicesToIncrement < 0 )
          {
          reverseOrder = true;
          }
        }
      if( numberOfSlicesToIncrement == 0 )
        {
        std::cerr << "Need greater than 0 slices for incrementing." << std::endl;
        return EXIT_FAILURE;
        }

      std::ostringstream stream;
      stream << startingSlice;
      std::string startingSliceString = stream.str();
      stream << endSlice;
      std::string endSliceString = stream.str();

      if( slicesOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        startingSliceString = slicesOption->GetFunction( 0 )->GetParameter( 1 );
        }
      if( slicesOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
        {
        endSliceString = slicesOption->GetFunction( 0 )->GetParameter( 2 );
        }

      bool isStartingSliceMaskDependent =
        startingSliceString.find( std::string( "mask" ) ) != std::string::npos;
      bool isEndSliceMaskDependent =
        endSliceString.find( std::string( "mask" ) ) != std::string::npos;

      if( isStartingSliceMaskDependent || isEndSliceMaskDependent )
        {
        if( !maskImage )
          {
          std::cerr << "Mask image is not specified." << std::endl;
          return EXIT_FAILURE;
          }

        if( isStartingSliceMaskDependent )
          {
          startingSlice = maskRegion.GetIndex()[direction];
          }
        else
          {
          startingSlice = parser->Convert<unsigned int>( startingSliceString );
          }

        if( startingSliceString.find( std::string( "+" ) ) != std::string::npos )
          {
          std::string offsetString = startingSliceString.substr( startingSliceString.find( std::string( "+" ) ) + 1 );
          int offset = parser->Convert<int>( offsetString );
          startingSlice += offset;
          }
        else if( startingSliceString.find( std::string( "-" ) ) != std::string::npos )
          {
          std::string offsetString = startingSliceString.substr( startingSliceString.find( std::string( "-" ) ) + 1 );
          int offset = parser->Convert<int>( offsetString );
          startingSlice -= offset;
          }

        if( isEndSliceMaskDependent )
          {
          endSlice = maskRegion.GetIndex()[direction] + maskRegion.GetSize()[direction];
          }
        else
          {
          endSlice = parser->Convert<unsigned int>( endSliceString );
          }

        if( endSliceString.find( std::string( "+" ) ) != std::string::npos )
          {
          std::string offsetString = endSliceString.substr( endSliceString.find( std::string( "+" ) ) + 1 );
          int offset = parser->Convert<int>( offsetString );
          endSlice += offset;
          }
        else if( endSliceString.find( std::string( "-" ) ) != std::string::npos )
          {
          std::string offsetString = endSliceString.substr( endSliceString.find( std::string( "-" ) ) + 1 );
          int offset = parser->Convert<int>( offsetString );
          endSlice -= offset;
          }

        }
      else
        {
        startingSlice = parser->Convert<unsigned int>( startingSliceString );
        endSlice = parser->Convert<unsigned int>( endSliceString );
        }

      startingSlice = vnl_math_max( itk::NumericTraits<int>::Zero, startingSlice );
      startingSlice = vnl_math_min( startingSlice, static_cast<int>( size[direction] - 1 ) );

      endSlice = vnl_math_max( itk::NumericTraits<int>::Zero, endSlice );
      endSlice = vnl_math_min( endSlice, static_cast<int>( size[direction] - 1 ) );

      whichSlices.clear();
      if( reverseOrder )
        {
        for( int n = endSlice; n >= startingSlice; n -= vnl_math_abs( numberOfSlicesToIncrement ) )
          {
          whichSlices.push_back( n );
          }
        }
      else
        {
        for( int n = startingSlice; n <= endSlice; n += numberOfSlicesToIncrement )
          {
          whichSlices.push_back( n );
          }
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

  itk::ants::CommandLineParser::OptionType::Pointer flipOption =
    parser->GetOption( "flip-slice" );

  bool doFlipHorizontally = false;
  bool doFlipVertically = false;

  if( flipOption && flipOption->GetNumberOfFunctions() )
    {
    std::vector<bool> layout = parser->ConvertVector<bool>( flipOption->GetFunction( 0 )->GetName() );
    if( layout.size() > 2 )
      {
      std::cerr << "Flip layout is specified as doFlipXxdoFlipY" << std::endl;
      return EXIT_FAILURE;
      }
    doFlipHorizontally = layout[0];
    doFlipVertically = layout[1];
    }

  itk::ants::CommandLineParser::OptionType::Pointer permuteOption =
    parser->GetOption( "permute-axes" );

  bool doPermute = false;

  if( permuteOption && permuteOption->GetNumberOfFunctions() )
    {
    doPermute = parser->Convert<bool>( permuteOption->GetFunction( 0 )->GetName() );
    }

  // Read in optional Rgb image

  RgbImageType::Pointer rgbImage = NULL;

  itk::ants::CommandLineParser::OptionType::Pointer rgbImageOption =
    parser->GetOption( "rgb-image" );
  if( rgbImageOption && rgbImageOption->GetNumberOfFunctions() )
    {
    std::string rgbFile = rgbImageOption->GetFunction( 0 )->GetName();
    ReadImage<RgbImageType>( rgbImage, rgbFile.c_str() );
    }

  RealType minIntensityValue = 0.0;
  RealType maxIntensityValue = 1.0;
  if( rgbImage )
    {
    typedef itk::StatisticsImageFilter<ImageType> StatisticsImageFilterType;
    StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();
    statisticsImageFilter->SetInput( inputImage );
    statisticsImageFilter->Update();

    minIntensityValue = statisticsImageFilter->GetMinimum();
    maxIntensityValue = statisticsImageFilter->GetMaximum();
    }

  RealType alpha = 1.0;

  itk::ants::CommandLineParser::OptionType::Pointer alphaOption =
    parser->GetOption( "alpha" );
  if( alphaOption && alphaOption->GetNumberOfFunctions() )
    {
    alpha = parser->Convert<RealType>( alphaOption->GetFunction( 0 )->GetName() );

    if( alpha < 0 || alpha > 1.0 )
      {
      std::cerr << "The alpha parameter must be between 0 and 1." << std::endl;
      return EXIT_FAILURE;
      }
    }

  // Now do the tiling

  std::cout << "Slices[" << direction << "]: " << whichSlices.size() << std::endl;
  std::cout << "Rows:  " << numberOfRows << std::endl;
  std::cout << "Columns:  " << numberOfColumns << std::endl;


  typedef itk::TileImageFilter<SliceType, SliceType> TileFilterType;
  TileFilterType::LayoutArrayType array;

  array[0] = numberOfColumns;
  array[1] = numberOfRows;

  ImageType::RegionType region;
  size[direction] = 0;

  TileFilterType::Pointer tileFilter = TileFilterType::New();
  tileFilter->SetLayout( array );

  typedef itk::TileImageFilter<RgbSliceType, RgbSliceType> RgbTileFilterType;
  RgbTileFilterType::Pointer rgbTileFilter = RgbTileFilterType::New();
  rgbTileFilter->SetLayout( array );

  for( unsigned int i = 0; i < whichSlices.size(); i++ )
    {
    unsigned int n = whichSlices[i];

    std::cout << "Processing slice " << n << std::endl;

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
    SliceType::Pointer outputSlice2 = NULL;

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

    typedef itk::FlipImageFilter<SliceType> FlipFilterType;
    FlipFilterType::Pointer flipper = FlipFilterType::New();
    FlipFilterType::FlipAxesArrayType flipArray;
    flipArray[0] = doFlipHorizontally;
    flipArray[1] = doFlipVertically;

    flipper->SetInput( outputSlice );
    flipper->SetFlipAxes( flipArray );

    typedef itk::PermuteAxesImageFilter<SliceType> PermuteAxesImageFilterType;
    itk::FixedArray<unsigned int, 2> order;
    order[0] = 0;
    order[1] = 1;
    if( doPermute )
      {
      order[0] = 1;
      order[1] = 0;
      }
    PermuteAxesImageFilterType::Pointer permuteAxesFilter = PermuteAxesImageFilterType::New();
    permuteAxesFilter->SetInput( flipper->GetOutput() );
    permuteAxesFilter->SetOrder( order );

    outputSlice2 = permuteAxesFilter->GetOutput();
    outputSlice2->Update();
    outputSlice2->DisconnectPipeline();

    if( rgbImage )
      {

      SliceType::Pointer outputMaskSlice = NULL;
      SliceType::Pointer outputMaskSlice2 = NULL;

      if( maskImage )
        {
        ExtracterType::Pointer maskExtracter = ExtracterType::New();
        maskExtracter->SetInput( maskImage );
        maskExtracter->SetExtractionRegion( region );
        maskExtracter->SetDirectionCollapseToIdentity();

        if( paddingType == -1 )
          {
          typedef itk::ExtractImageFilter<SliceType, SliceType> ExtracterType2;
          ExtracterType2::Pointer maskExtracter2 = ExtracterType2::New();
          maskExtracter2->SetInput( maskExtracter->GetOutput() );
          maskExtracter2->SetExtractionRegion( croppedSliceRegion );
          maskExtracter2->SetDirectionCollapseToIdentity();

          outputMaskSlice = maskExtracter2->GetOutput();
          outputMaskSlice->Update();
          outputMaskSlice->DisconnectPipeline();
          }
        else if( paddingType == 1 )
          {
          typedef itk::ConstantPadImageFilter<SliceType, SliceType> PadderType;
          PadderType::Pointer maskPadder = PadderType::New();
          maskPadder->SetInput( maskExtracter->GetOutput() );
          maskPadder->SetPadLowerBound( lowerBound );
          maskPadder->SetPadUpperBound( upperBound );
          maskPadder->SetConstant( 0 );

          outputMaskSlice = maskPadder->GetOutput();
          outputMaskSlice->Update();
          outputMaskSlice->DisconnectPipeline();
          }
        else // paddingType == 0
          {
          outputMaskSlice = maskExtracter->GetOutput();
          outputMaskSlice->Update();
          outputMaskSlice->DisconnectPipeline();
          }

        FlipFilterType::Pointer maskFlipper = FlipFilterType::New();
        maskFlipper->SetInput( outputMaskSlice );
        maskFlipper->SetFlipAxes( flipArray );

        PermuteAxesImageFilterType::Pointer maskPermuteAxesFilter = PermuteAxesImageFilterType::New();
        maskPermuteAxesFilter->SetInput( maskFlipper->GetOutput() );
        maskPermuteAxesFilter->SetOrder( order );

        outputMaskSlice2 = maskPermuteAxesFilter->GetOutput();
        outputMaskSlice2->Update();
        outputMaskSlice2->DisconnectPipeline();
        }

      RgbSliceType::Pointer outputRgbSlice = NULL;
      RgbSliceType::Pointer outputRgbSlice2 = NULL;

      typedef itk::ExtractImageFilter<RgbImageType, RgbSliceType> RgbExtracterType;
      RgbExtracterType::Pointer rgbExtracter = RgbExtracterType::New();
      rgbExtracter->SetInput( rgbImage );
      rgbExtracter->SetExtractionRegion( region );
      rgbExtracter->SetDirectionCollapseToIdentity();

      if( paddingType == -1 )
        {
        typedef itk::ExtractImageFilter<RgbSliceType, RgbSliceType> RgbExtracterType2;
        RgbExtracterType2::Pointer rgbExtracter2 = RgbExtracterType2::New();
        rgbExtracter2->SetInput( rgbExtracter->GetOutput() );
        rgbExtracter2->SetExtractionRegion( croppedSliceRegion );
        rgbExtracter2->SetDirectionCollapseToIdentity();

        outputRgbSlice = rgbExtracter2->GetOutput();
        outputRgbSlice->Update();
        outputRgbSlice->DisconnectPipeline();
        }
      else if( paddingType == 1 )
        {
        typedef itk::ConstantPadImageFilter<RgbSliceType, RgbSliceType> RgbPadderType;
        RgbPadderType::Pointer rgbPadder = RgbPadderType::New();
        rgbPadder->SetInput( rgbExtracter->GetOutput() );
        rgbPadder->SetPadLowerBound( lowerBound );
        rgbPadder->SetPadUpperBound( upperBound );
        rgbPadder->SetConstant( static_cast<PixelType>( padValue ) );

        outputRgbSlice = rgbPadder->GetOutput();
        outputRgbSlice->Update();
        outputRgbSlice->DisconnectPipeline();
        }
      else // paddingType == 0
        {
        outputRgbSlice = rgbExtracter->GetOutput();
        outputRgbSlice->Update();
        outputRgbSlice->DisconnectPipeline();
        }

      typedef itk::FlipImageFilter<RgbSliceType> RgbFlipFilterType;
      RgbFlipFilterType::Pointer rgbFlipper = RgbFlipFilterType::New();
      RgbFlipFilterType::FlipAxesArrayType rgbFlipArray;
      rgbFlipArray[0] = doFlipHorizontally;
      rgbFlipArray[1] = doFlipVertically;

      rgbFlipper->SetInput( outputRgbSlice );
      rgbFlipper->SetFlipAxes( rgbFlipArray );

      typedef itk::PermuteAxesImageFilter<RgbSliceType> RgbPermuteAxesImageFilterType;
      itk::FixedArray<unsigned int, 2> rgbOrder;
      rgbOrder[0] = 0;
      rgbOrder[1] = 1;
      if( doPermute )
        {
        rgbOrder[0] = 1;
        rgbOrder[1] = 0;
        }
      RgbPermuteAxesImageFilterType::Pointer rgbPermuteAxesFilter = RgbPermuteAxesImageFilterType::New();
      rgbPermuteAxesFilter->SetInput( rgbFlipper->GetOutput() );
      rgbPermuteAxesFilter->SetOrder( rgbOrder );

      outputRgbSlice2 = rgbPermuteAxesFilter->GetOutput();
      outputRgbSlice2->Update();
      outputRgbSlice2->DisconnectPipeline();

      // combine grayscale slice and rgb slice
      itk::ImageRegionConstIteratorWithIndex<SliceType> It( outputSlice2,
        outputSlice2->GetRequestedRegion() );
      itk::ImageRegionIterator<RgbSliceType> ItRgb( outputRgbSlice2,
        outputRgbSlice2->GetRequestedRegion() );
      for( It.GoToBegin(), ItRgb.GoToBegin(); !It.IsAtEnd(); ++It, ++ItRgb )
        {
        RgbPixelType rgbPixel = ItRgb.Get();
        PixelType pixel = ( It.Get() - minIntensityValue ) / ( maxIntensityValue - minIntensityValue )
          * itk::NumericTraits<RgbComponentType>::max();

        if( outputMaskSlice2 && outputMaskSlice2->GetPixel( It.GetIndex() ) == 0 )
          {
          rgbPixel.SetRed( pixel );
          rgbPixel.SetGreen( pixel );
          rgbPixel.SetBlue( pixel );
          }
        else
          {
          rgbPixel.SetRed( static_cast<RgbComponentType>( ( 1.0 - alpha ) * pixel + alpha * rgbPixel.GetRed() ) );
          rgbPixel.SetGreen( static_cast<RgbComponentType>( ( 1.0 - alpha ) * pixel + alpha * rgbPixel.GetGreen() ) );
          rgbPixel.SetBlue( static_cast<RgbComponentType>( ( 1.0 - alpha ) * pixel + alpha * rgbPixel.GetBlue() ) );
          }

        ItRgb.Set( rgbPixel );
        }

      rgbTileFilter->SetInput( i, outputRgbSlice2 );
      }
    else
      {
      tileFilter->SetInput( i, outputSlice2 );
      }
    }

  itk::ants::CommandLineParser::OptionType::Pointer outputOption =
    parser->GetOption( "output" );
  if( outputOption && outputOption->GetNumberOfFunctions() )
    {
    std::string outputFile = outputOption->GetFunction( 0 )->GetName();
    if( rgbImage )
      {
      rgbTileFilter->Update();
      WriteImage<RgbSliceType>( rgbTileFilter->GetOutput(), outputFile.c_str() );
      }
    else
      {
      tileFilter->Update();
      WriteImage<SliceType>( tileFilter->GetOutput(), outputFile.c_str() );
      }
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
      std::string( "An optional Rgb image can be added as an overlay.  ")
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
      std::string( "Specifies the ROI of the RGB voxels used.  ");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "mask-image" );
    option->SetShortName( 'x' );
    option->SetUsageOption( 0, "maskImageFilename" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "If an Rgb image is provided, render the overlay using the specified " )
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
      + std::string( "user can also specify a constant pad value (default = 0). If a mask is " )
      + std::string( "specified, the user can use the mask to define the region, by using " )
      + std::string( "the keyword \"mask\" plus an offset, e.g. \"-p mask+3\"." );

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
      + std::string( "the number slices to increment with the optional specification of " )
      + std::string( "which slices to start and end the sequence.  A negative value " )
      + std::string( "for the numberOfSlicesToIncrement causes rendering in the reverse " )
      + std::string( "order.  For the third option, minSlice < maxSlice.  If a mask is " )
      + std::string( "specified, the user can use the mask to define the region, by using " )
      + std::string( "the keyword \"mask\" plus an offset, e.g. \"-s [1,mask-3,200]\"." )
      + std::string( "For the third option, minSlice < maxSlice." );

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "slices" );
    option->SetShortName( 's' );
    option->SetUsageOption( 0, "Slice1xSlice2xSlice3..." );
    option->SetUsageOption( 1, "numberOfSlicesToIncrement" );
    option->SetUsageOption( 2, "[numberOfSlicesToIncrement,<minSlice=0>,<maxSlice=lastSlice>]" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Flip individual slice images horizontally and/or vertically, specified " )
      + std::string( "e.g. as \'0x1\' or \'1x1\'.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "flip-slice" );
    option->SetShortName( 'f' );
    option->SetUsageOption( 0, "flipXxflipY" );
    option->SetDescription( description );
    parser->AddOption( option );
    }

    {
    std::string description =
      std::string( "Permute (or swap) the axes of the individual slice images.");

    OptionType::Pointer option = OptionType::New();
    option->SetLongName( "permute-axes" );
    option->SetShortName( 'g' );
    option->SetUsageOption( 0, "doPermute" );
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
    std::string( "Render a 3-D image volume with optional Rgb overlay." );

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
