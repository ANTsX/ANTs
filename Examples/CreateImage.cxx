#include <stdio.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"

#include <string>
#include <vector>

template <class TValue>
TValue Convert( std::string optionString )
{
  TValue             value;
  std::istringstream iss( optionString );

  iss >> value;
  return value;
}

template <class TValue>
std::vector<TValue> ConvertVector( std::string optionString )
{
  std::vector<TValue>    values;
  std::string::size_type crosspos = optionString.find( 'x', 0 );

  if( crosspos == std::string::npos )
    {
    values.push_back( Convert<TValue>( optionString ) );
    }
  else
    {
    std::string        element = optionString.substr( 0, crosspos );
    TValue             value;
    std::istringstream iss( element );
    iss >> value;
    values.push_back( value );
    while( crosspos != std::string::npos )
      {
      std::string::size_type crossposfrom = crosspos;
      crosspos = optionString.find( 'x', crossposfrom + 1 );
      if( crosspos == std::string::npos )
        {
        element = optionString.substr( crossposfrom + 1, optionString.length() );
        }
      else
        {
        element = optionString.substr( crossposfrom + 1, crosspos );
        }
      iss.str(element);
      iss >> value;
      values.push_back( value );
      }
    }
  return values;
}

template <int ImageDimension>
int CreateZeroImage( int argc, char *argv[] )
{
  typedef float                                 PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typename ImageType::Pointer image = ImageType::New();

  typedef typename
    itk::Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
  typename GeneratorType::Pointer generator = GeneratorType::New();
  generator->Initialize();
  generator->SetSeed();

  std::string which = std::string( argv[3] );
  typename std::string::size_type pos = which.find( "." );

  typename std::string::size_type pos3 = std::string::npos;
  if( argc > 6 )
    {
    std::string pixelValues = std::string( argv[6] );
    pos3 = pixelValues.find( "x" );
    }
  if( pos3 != std::string::npos )
    {
    std::vector<float> og = ConvertVector<float>( std::string( argv[3] ) );
    std::vector<float> sp = ConvertVector<float>( std::string( argv[4] ) );
    std::vector<int>   sz = ConvertVector<int>( std::string( argv[5] ) );
    std::vector<float> values = ConvertVector<float>( std::string( argv[6] ) );

    unsigned long numberOfPixels = 1;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      numberOfPixels *= sz[d];
      }
    if( values.size() > numberOfPixels )
      {
      std::cerr << "Number of specified pixel values is greater than "
                << "the size of the image." << std::endl;
      return EXIT_FAILURE;
      }

    if( og.size() != ImageDimension )
      {
      std::cerr << "Invalid origin size." << std::endl;
      return EXIT_FAILURE;
      }
    if( sp.size() != ImageDimension )
      {
      std::cerr << "Invalid spacing size." << std::endl;
      return EXIT_FAILURE;
      }
    if( sz.size() != ImageDimension )
      {
      std::cerr << "Invalid Size size." << std::endl;
      return EXIT_FAILURE;
      }

    typename ImageType::PointType origin;
    typename ImageType::SpacingType spacing;
    typename ImageType::SizeType size;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      origin[d] = og[d];
      spacing[d] = sp[d];
      size[d] = sz[d];
      }
    image->SetOrigin( origin );
    image->SetSpacing( spacing );
    image->SetRegions( size );
    image->Allocate();
    image->FillBuffer( 0.0 );

    unsigned long                       count = 0;
    itk::ImageRegionIterator<ImageType> It( image, image->GetRequestedRegion() );

    It.GoToBegin();
    while( !It.IsAtEnd() && count < values.size() )
      {
      It.Set( values[count] );

      ++It;
      ++count;
      }

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[2] );
    writer->SetInput( image );
    writer->Update();
    }
  else if( pos != std::string::npos )
    {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[2] );
    reader->Update();

    image->SetOrigin( reader->GetOutput()->GetOrigin() );
    image->SetSpacing( reader->GetOutput()->GetSpacing() );
    image->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
    image->Allocate();
    image->FillBuffer( atof( argv[4] ) );

    if( argc > 5  )
      {
      switch( atoi( argv[5] ) )
        {
        case 1: default:
          {
          itk::ImageRegionIterator<ImageType> It( image,
                                                  image->GetLargestPossibleRegion() );
          for( It.GoToBegin(); !It.IsAtEnd(); ++It )
            {
            It.Set( static_cast<PixelType>(
                      generator->GetIntegerVariate( static_cast<int>( It.Get() ) ) ) );
            }
          break;
          }
//        case 2:
//         {
//                                    itk::ImageRegionIteratorWithIndex<ImageType> ItI( image,
//                                            image->GetLargestPossibleRegion() );
//                                    for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
//                                            {
//                                            ItI.Set( constant - ItI.GetIndex()[d] );
//                                            }
//          break;
//          }
//        default:
//          std::cerr << "Incorrect choice" << std::endl;
//          return EXIT_FAILURE;
//          break;
        }
      }

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[3] );
    writer->SetInput( image );
    writer->Update();
    }
  else
    {
    std::vector<float> og = ConvertVector<float>( std::string( argv[3] ) );
    std::vector<float> sp = ConvertVector<float>( std::string( argv[4] ) );
    std::vector<int>   sz = ConvertVector<int>( std::string( argv[5] ) );

    if( og.size() != ImageDimension )
      {
      std::cerr << "Invalid origin size." << std::endl;
      return EXIT_FAILURE;
      }
    if( sp.size() != ImageDimension )
      {
      std::cerr << "Invalid spacing size." << std::endl;
      return EXIT_FAILURE;
      }
    if( sz.size() != ImageDimension )
      {
      std::cerr << "Invalid Size size." << std::endl;
      return EXIT_FAILURE;
      }

    typename ImageType::PointType origin;
    typename ImageType::SpacingType spacing;
    typename ImageType::SizeType size;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      origin[d] = og[d];
      spacing[d] = sp[d];
      size[d] = sz[d];
      }
    image->SetOrigin( origin );
    image->SetSpacing( spacing );
    image->SetRegions( size );
    image->Allocate();
    image->FillBuffer( atof( argv[6] ) );

    if( argc > 7  )
      {
      switch( atoi( argv[7] ) )
        {
        case 1: default:
          {
          itk::ImageRegionIterator<ImageType> It( image,
                                                  image->GetLargestPossibleRegion() );
          for( It.GoToBegin(); !It.IsAtEnd(); ++It )
            {
            It.Set( static_cast<PixelType>(
                      generator->GetIntegerVariate( static_cast<int>( It.Get() ) ) ) );
            }
          break;
          }
//        case 2:
//         {
//                                    itk::ImageRegionIteratorWithIndex<ImageType> ItI( image,
//                                            image->GetLargestPossibleRegion() );
//                                    for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
//                                            {
//                                            ItI.Set( constant - ItI.GetIndex()[d] );
//                                            }
//          break;
//          }
//        default:
//          {
//          std::cerr << "Incorrect choice" << std::endl;
//          return EXIT_FAILURE;
//          break;
//          }
        }
      }

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[2] );
    writer->SetInput( image );
    writer->Update();
    }

  return 0;
}

int main( int argc, char *argv[] )
{
  if( argc < 5 )
    {
    std::cout << "Usage 1: " << argv[0] << " imageDimension referenceImage outputImage constant [random?]" << std::endl;
    std::cout << "Usage 2: " << argv[0] << " imageDimension outputImage origin spacing size constant [random?]"
              << std::endl;
    std::cout << "Usage 3: " << argv[0] << " imageDimension outputImage origin spacing size pixelValues" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
    {
    case 1:
      {
      CreateZeroImage<1>( argc, argv );
      break;
      }
    case 2:
      {
      CreateZeroImage<2>( argc, argv );
      break;
      }
    case 3:
      {
      CreateZeroImage<3>( argc, argv );
      break;
      }
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
}
