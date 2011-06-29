#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageRegionIteratorWithIndex.h"

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
      std::istringstream iss2( element );
      iss2 >> value;
      values.push_back( value );
      }
    }
  return values;
}

template <unsigned int ImageDimension>
int PasteImageIntoImage( unsigned int argc, char *argv[] )
{
  typedef float                                 PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  std::vector<unsigned int> startIndex =
    ConvertVector<unsigned int>( std::string( argv[5] ) );

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[2] );
  reader1->Update();

  typename ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[3] );
  reader2->Update();

  PixelType backgroundValue = 0;
  if( argc > 6 )
    {
    backgroundValue = static_cast<PixelType>( atof( argv[6] ) );
    }
  bool writeOver = true;
  if( argc > 7 )
    {
    writeOver = static_cast<bool>( atoi( argv[7] ) );
    }

  itk::ImageRegionIteratorWithIndex<ImageType> It(
    reader2->GetOutput(), reader2->GetOutput()->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    typename ImageType::IndexType index = It.GetIndex();
    PixelType value = It.Get();
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      index[d] += startIndex[d];
      }
    if( value != backgroundValue )
      {
      if( reader1->GetOutput()->GetLargestPossibleRegion().IsInside( index ) )
        {
        PixelType canvasValue = reader1->GetOutput()->GetPixel( index );

        if( canvasValue == backgroundValue || writeOver == true )
          {
          reader1->GetOutput()->SetPixel( index, value );
          }
        }
      }
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[4] );
  writer->SetInput( reader1->GetOutput() );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if( argc < 6 )
    {
    std::cout << argv[0] << " imageDimension inputCanvasImage inputImage "
              << "outputImage startIndex [backgroundLabel=0] [writeOverNonBackgroundLabels=1]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
    {
    case 2:
      PasteImageIntoImage<2>( argc, argv );
      break;
    case 3:
      PasteImageIntoImage<3>( argc, argv );
      break;
    case 4:
      PasteImageIntoImage<4>( argc, argv );
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
}
