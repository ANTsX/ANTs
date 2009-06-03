#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorImageFileReader.h"
#include "itkVector.h"

#include "itkWarpImageFilter.h"
#include "itkGridImageSource.h"

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
int CreateWarpedGridImage( int argc, char *argv[] )
{
  typedef float                                  RealType;
  typedef itk::Image<RealType, ImageDimension>   RealImageType;
  typedef itk::Vector<RealType, ImageDimension>  VectorType;
  typedef itk::Image<VectorType, ImageDimension> VectorImageType;

  /**
   * Read in vector field
   */
  typedef itk::VectorImageFileReader<RealImageType, VectorImageType> ReaderType;
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
      std::cerr << "Incorrect direction size." << std::endl;
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
      std::cerr << "Incorrect spacing size." << std::endl;
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
      std::cerr << "Incorrect sigma size." << std::endl;
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

  typedef itk::WarpImageFilter<RealImageType, RealImageType, VectorImageType> WarperType;
  typename WarperType::Pointer warper = WarperType::New();
  warper->SetDeformationField( reader->GetOutput() );
  warper->SetInput( gridder->GetOutput() );
  warper->SetOutputOrigin( gridder->GetOutput()->GetOrigin() );
  warper->SetOutputSpacing( gridder->GetOutput()->GetSpacing() );
  warper->Update();

  std::string file = std::string( argv[3] );
  typedef itk::ImageFileWriter<RealImageType> ImageWriterType;
  typename ImageWriterType::Pointer gridWriter = ImageWriterType::New();
  gridWriter->SetFileName( file.c_str() );
  gridWriter->SetInput( warper->GetOutput() );
  gridWriter->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " ImageDimension deformationField "
              << "outputImage [directions,e.g. 1x0x0] [gridSpacing] [gridSigma]"
              << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
    {
    case 2:
      CreateWarpedGridImage<2>( argc, argv );
      break;
    case 3:
      CreateWarpedGridImage<3>( argc, argv );
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
}
