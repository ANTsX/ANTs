#include <stdio.h>

#include "itkCastImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkExtractImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"

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
      iss.str( element );
      iss >> value;
      values.push_back( value );
      }
    }
  return values;
}

template <unsigned int ImageDimension>
int ExtractRegionFromImage( int argc, char *argv[] )
{
  typedef float PixelType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType>       ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typename ImageType::RegionType region;
  typename ImageType::RegionType::SizeType size;
  typename ImageType::RegionType::IndexType index;

  if( argc == 6 )
    {
    std::vector<int> minIndex;
    std::vector<int> maxIndex;
    minIndex = ConvertVector<int>( std::string( argv[4] ) );
    maxIndex = ConvertVector<int>( std::string( argv[5] ) );
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      index[i] = minIndex[i];
      size[i] = maxIndex[i] - minIndex[i] + 1;
      }
    region.SetSize( size );
    region.SetIndex( index );
    }
  else
    {
    typedef itk::Image<unsigned short, ImageDimension>      ShortImageType;
    typedef itk::CastImageFilter<ImageType, ShortImageType> CasterType;
    typename CasterType::Pointer caster = CasterType::New();
    caster->SetInput( reader->GetOutput() );
    caster->Update();

    typedef itk::LabelStatisticsImageFilter<ShortImageType, ShortImageType>
      StatsFilterType;
    typename StatsFilterType::Pointer stats = StatsFilterType::New();
    stats->SetLabelInput( caster->GetOutput() );
    stats->SetInput( caster->GetOutput() );
    stats->Update();

    region = stats->GetRegion( atoi( argv[4] ) );
    }

  std::cout << region << std::endl;

  typedef itk::ExtractImageFilter<ImageType, ImageType> CropperType;
  typename CropperType::Pointer cropper = CropperType::New();
  cropper->SetInput( reader->GetOutput() );
  cropper->SetExtractionRegion( region );
  cropper->SetDirectionCollapseToSubmatrix();
  cropper->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( cropper->GetOutput() );
  writer->SetFileName( argv[3] );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if( argc < 5 || argc > 6 )
    {
    std::cout << "Usage 1: " << argv[0] << " ImageDimension "
              << "inputImage outputImage minIndex maxIndex " << std::endl;
    std::cout << "Usage 2: " << argv[0] << " ImageDimension "
              << "inputImage outputImage label " << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
    {
    case 2:
      ExtractRegionFromImage<2>( argc, argv );
      break;
    case 3:
      ExtractRegionFromImage<3>( argc, argv );
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
}
