
#include "antsUtilities.h"
#include <algorithm>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageRegionIteratorWithIndex.h"

#include <string>
#include <vector>

namespace ants
{
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
  unsigned int writeOver = 1;
  if( argc > 7 )
    {
    writeOver = static_cast<unsigned int>( atoi( argv[7] ) );
    }
  PixelType conflictLabel = -1;
  if( argc > 8 )
    {
    conflictLabel = static_cast<PixelType>( atof( argv[8] ) );
    }

  itk::ImageRegionIteratorWithIndex<ImageType> It(
    reader2->GetOutput(), reader2->GetOutput()->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    typename ImageType::IndexType index = It.GetIndex();
    PixelType paintValue = It.Get();
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      index[d] += startIndex[d];
      }
    if( paintValue != backgroundValue )
      {
      if( reader1->GetOutput()->GetLargestPossibleRegion().IsInside( index ) )
        {
        PixelType canvasValue = reader1->GetOutput()->GetPixel( index );
        if( canvasValue == backgroundValue || writeOver == 1 )
          {
          reader1->GetOutput()->SetPixel( index, paintValue );
          }
        else if( canvasValue != backgroundValue && writeOver == 0 )
          {
          continue;
          }
        else
          {
          reader1->GetOutput()->SetPixel( index, conflictLabel );
          }
        }
      }
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[4] );
  writer->SetInput( reader1->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int PasteImageIntoImage( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "PasteImageIntoImage" );

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

  if( argc < 6 )
    {
    std::cout << argv[0] << " imageDimension inputCanvasImage inputImage "
             << "outputImage startIndex [backgroundLabel=0] [paintOverNonBackgroundVoxels=0] [conflictLabel=-1]"
             << std::endl;
    std::cout
      <<
      "   If the current painting image voxel is nonbackground and corresponds to a background voxel in the canvas image "
      << std::endl;
    std::cout << "     paintOverNonBackgroundVoxels = 0 -> leave the canvas voxel as is." << std::endl;
    std::cout
      << "     paintOverNonBackgroundVoxels = 1 -> replace canvas voxel value with painting image voxel value"
      << std::endl;
    std::cout << "     paintOverNonBackgroundVoxels = 2 -> replace canvas voxel walue with conflictLabel"  << std::endl;

    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) )
    {
    case 2:
      {
      PasteImageIntoImage<2>( argc, argv );
      }
      break;
    case 3:
      {
      PasteImageIntoImage<3>( argc, argv );
      }
      break;
    case 4:
      {
      PasteImageIntoImage<4>( argc, argv );
      }
      break;
    default:
      std::cout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
