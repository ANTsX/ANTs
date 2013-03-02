
#include "antsUtilities.h"
#include <algorithm>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkMatrixOffsetTransformBase.h"
#include "itkWarpImageMultiTransformFilter.h"
#include "itkGridImageSource.h"

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
int CreateWarpedGridImage( int argc, char *argv[] )
{
  typedef float                                  RealType;
  typedef itk::Image<RealType, ImageDimension>   RealImageType;
  typedef itk::Vector<RealType, ImageDimension>  VectorType;
  typedef itk::Image<VectorType, ImageDimension> VectorImageType;

  /**
   * Read in vector field
   */
  typedef itk::ImageFileReader<VectorImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::GridImageSource<RealImageType> GridSourceType;
  typename GridSourceType::Pointer gridder = GridSourceType::New();
  gridder->SetSpacing( reader->GetOutput()->GetSpacing() );
  gridder->SetOrigin( reader->GetOutput()->GetOrigin() );
  gridder->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );

  typename GridSourceType::ArrayType gridSpacing;
  typename GridSourceType::ArrayType gridSigma;
  typename GridSourceType::BoolArrayType which;
  which.Fill( false );
  for( unsigned int i = 0; i < 2; i++ )
    {
    which[i] = true;
    }

  if( argc > 4 )
    {
    std::vector<unsigned int> directions
      = ConvertVector<unsigned int>( std::string( argv[4] ) );
    if( directions.size() != ImageDimension )
      {
      antscout << "Incorrect direction size." << std::endl;
      return EXIT_FAILURE;
      }
    else
      {
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        which[i] = static_cast<bool>( directions[i] );
        }
      }
    }
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    gridSpacing[i] = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[i]
      * reader->GetOutput()->GetSpacing()[i] / 25.0;
    gridSigma[i] = gridSpacing[i] / 10.0;
    }
  if( argc > 5 )
    {
    std::vector<RealType> spacing
      = ConvertVector<RealType>( std::string( argv[5] ) );
    if( spacing.size() != ImageDimension )
      {
      antscout << "Incorrect spacing size." << std::endl;
      return EXIT_FAILURE;
      }
    else
      {
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        gridSpacing[i] = spacing[i];
        gridSigma[i] = gridSpacing[i] / 10.0;
        }
      }
    }
  if( argc > 6 )
    {
    std::vector<RealType> sigma
      = ConvertVector<RealType>( std::string( argv[6] ) );
    if( sigma.size() != ImageDimension )
      {
      antscout << "Incorrect sigma size." << std::endl;
      return EXIT_FAILURE;
      }
    else
      {
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        gridSigma[i] = sigma[i] / 10.0;
        }
      }
    }

  gridder->SetGridSpacing( gridSpacing );
  gridder->SetSigma( gridSigma );
  gridder->SetWhichDimensions( which );
  gridder->Update();
  typename RealImageType::Pointer grid = gridder->GetOutput();
  grid->SetDirection(reader->GetOutput()->GetDirection() );
  grid->SetOrigin(reader->GetOutput()->GetOrigin() );
  grid->SetSpacing(reader->GetOutput()->GetSpacing() );

  typedef itk::MatrixOffsetTransformBase<double, ImageDimension,
                                         ImageDimension>
    TransformType;
  typedef itk::WarpImageMultiTransformFilter<RealImageType, RealImageType, VectorImageType, TransformType> WarperType;
  typename WarperType::Pointer  warper = WarperType::New();
  warper->SetInput(grid);
  warper->SetEdgePaddingValue( 0);
  warper->SetSmoothScale(1);
  warper->PushBackDisplacementFieldTransform(reader->GetOutput() );
  warper->SetOutputParametersFromImage(  reader->GetOutput() );
  warper->Update();

  std::string file = std::string( argv[3] );
  typedef itk::ImageFileWriter<RealImageType> ImageWriterType;
  typename ImageWriterType::Pointer gridWriter = ImageWriterType::New();
  gridWriter->SetFileName( file.c_str() );
  gridWriter->SetInput( warper->GetOutput() );
  gridWriter->Update();

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int CreateWarpedGridImage( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "CreateWarpedGridImage" );

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

  if( argc < 4 )
    {
    antscout << "Usage: " << argv[0] << " ImageDimension deformationField "
             << "outputImage [directions,e.g. 1x0x0] [gridSpacing] [gridSigma]"
             << std::endl;
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
      CreateWarpedGridImage<2>( argc, argv );
      }
      break;
    case 3:
      {
      CreateWarpedGridImage<3>( argc, argv );
      }
      break;
    default:
      antscout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
