
#include "antsUtilities.h"
#include <algorithm>

#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"

namespace ants
{
template <unsigned int ImageDimension>
int ExtractSliceFromImage( int itkNotUsed( argc ), char *argv[] )
{
  typedef float PixelType;

  typedef itk::Image<PixelType, ImageDimension>     ImageType;
  typedef itk::Image<PixelType, ImageDimension - 1> SliceType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typename ImageType::RegionType region;
  typename ImageType::RegionType::SizeType size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
  size[atoi( argv[4] )] = 0;
  typename ImageType::IndexType index;
  index.Fill( 0 );
  index[atoi( argv[4] )] = atoi( argv[5] );
  region.SetIndex( index );
  region.SetSize( size );

  typedef itk::ExtractImageFilter<ImageType, SliceType> ExtracterType;
  typename ExtracterType::Pointer extracter = ExtracterType::New();
  extracter->SetInput( reader->GetOutput() );
  extracter->SetExtractionRegion( region );
  extracter->SetDirectionCollapseToIdentity();
  extracter->Update();

  typedef itk::ImageFileWriter<SliceType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( extracter->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int ExtractSliceFromImage( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "ExtractSliceFromImage" );

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

  antscout->set_stream( out_stream );

  if( argc != 6 )
    {
    antscout << "Usage: " << argv[0]
             << " imageDimension inputImage outputSlice direction(e.g. 0, 1, 2) slice_number" << std::endl;
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
      ExtractSliceFromImage<2>( argc, argv );
      }
      break;
    case 3:
      {
      ExtractSliceFromImage<3>( argc, argv );
      }
      break;
    case 4:
      {
      ExtractSliceFromImage<4>( argc, argv );
      }
      break;
    default:
      antscout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
