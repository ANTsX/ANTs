#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"

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

  return 0;
}

int main( int argc, char *argv[] )
{
  if( argc != 6 )
    {
    std::cout << "Usage: " << argv[0]
              << " imageDimension inputImage outputSlice direction(e.g. 0, 1, 2) slice_number" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
    {
    case 2:
      ExtractSliceFromImage<2>( argc, argv );
      break;
    case 3:
      ExtractSliceFromImage<3>( argc, argv );
      break;
    case 4:
      ExtractSliceFromImage<4>( argc, argv );
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
}
