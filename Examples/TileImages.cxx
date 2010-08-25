#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkTileImageFilter.h"

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
      std::istringstream iss( element );
      iss >> value;
      values.push_back( value );
      }
    }
  return values;
}

template <unsigned int ImageDimension>
int TileImages( unsigned int argc, char *argv[] )
{
  typedef float                                 PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::TileImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  typename FilterType::LayoutArrayType array;

  std::vector<unsigned int> layout = ConvertVector<unsigned int>( std::string( argv[3] ) );
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    array[d] = layout[d];
    }
  filter->SetLayout( array );
  for( unsigned int n = 4; n < argc; n++ )
    {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[n] );
    reader->Update();

    filter->SetInput( n - 4, reader->GetOutput() );
    }
  filter->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
    std::cout << argv[0] << " imageDimension outputImage layout inputImage1 ... inputImageN" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
    {
    case 2:
      TileImages<2>( argc, argv );
      break;
    case 3:
      TileImages<3>( argc, argv );
      break;
    case 4:
      TileImages<4>( argc, argv );
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
}
