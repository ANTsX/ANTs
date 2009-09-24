#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkLabelOverlapMeasuresImageFilter.h"

template <unsigned int ImageDimension>
int LabelOverlapMeasures( int argc, char * argv[] )
{
  typedef unsigned int                          PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[2] );
  typename ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[3] );

  typedef itk::LabelOverlapMeasuresImageFilter<ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetSourceImage( reader1->GetOutput() );
  filter->SetTargetImage( reader2->GetOutput() );
  filter->Update();

  std::cout << "All Labels" << std::endl;
  std::cout << "Volumetric overlap measures." << std::endl;
  std::cout << "  Total overlap:            "
            << filter->GetVolumeTotalOverlap() << std::endl;
  std::cout << "  Union overlap (jaccard):  "
            << filter->GetVolumeUnionOverlap() << std::endl;
  std::cout << "  Mean overlap (dice):      "
            << filter->GetVolumeMeanOverlap() << std::endl;
  std::cout << "  False negative error:     "
            << filter->GetVolumeFalseNegativeError() << std::endl;
  std::cout << "  False positive error:     "
            << filter->GetVolumeFalsePositiveError() << std::endl;
  std::cout << "Surface overlap measures." << std::endl;
  std::cout << "  Total overlap:            "
            << filter->GetSurfaceTotalOverlap() << std::endl;
  std::cout << "  Union overlap (jaccard):  "
            << filter->GetSurfaceUnionOverlap() << std::endl;
  std::cout << "  Mean overlap (dice):      "
            << filter->GetSurfaceMeanOverlap() << std::endl;
  std::cout << "  False negative error:     "
            << filter->GetSurfaceFalseNegativeError() << std::endl;
  std::cout << "  False positive error:     "
            << filter->GetSurfaceFalsePositiveError() << std::endl;
  std::cout << "Other related measures." << std::endl;
  std::cout << "  Volume similarity:        "
            << filter->GetVolumeSimilarity() << std::endl;
  std::cout << std::endl;

  std::cout << "Overlap measures over individual labels." << std::endl;
  typename FilterType::MapType labelMap = filter->GetLabelSetMeasures();
  typename FilterType::MapType::const_iterator it = labelMap.begin();
  while( it != labelMap.end() )
    {
    if( (*it).first == 0 )
      {
      ++it;
      continue;
      }
    std::cout << "Label " << (*it).first << std::endl;
    std::cout << "  Volumetric overlap measures." << std::endl;
    std::cout << "    Total overlap:            "
              << filter->GetVolumeTotalOverlap( (*it).first ) << std::endl;
    std::cout << "    Union overlap (jaccard):  "
              << filter->GetVolumeUnionOverlap( (*it).first ) << std::endl;
    std::cout << "    Mean overlap (dice):      "
              << filter->GetVolumeMeanOverlap( (*it).first ) << std::endl;
    std::cout << "    False negative error:     "
              << filter->GetVolumeFalseNegativeError( (*it).first ) << std::endl;
    std::cout << "    False positive error:     "
              << filter->GetVolumeFalsePositiveError( (*it).first ) << std::endl;
    std::cout << "  Surface overlap measures." << std::endl;
    std::cout << "    Total overlap:            "
              << filter->GetSurfaceTotalOverlap( (*it).first ) << std::endl;
    std::cout << "    Union overlap (jaccard):  "
              << filter->GetSurfaceUnionOverlap( (*it).first ) << std::endl;
    std::cout << "    Mean overlap (dice):      "
              << filter->GetSurfaceMeanOverlap( (*it).first ) << std::endl;
    std::cout << "    False negative error:     "
              << filter->GetSurfaceFalseNegativeError( (*it).first ) << std::endl;
    std::cout << "    False positive error:     "
              << filter->GetSurfaceFalsePositiveError( (*it).first ) << std::endl;
    std::cout << "  Other related measures." << std::endl;
    std::cout << "    Volume similarity:        "
              << filter->GetVolumeSimilarity( (*it).first ) << std::endl;
    std::cout << "    Hausdorff distance:       "
              << filter->GetHausdorffDistance( (*it).first )
              << " (in pixels) " << std::endl;
    std::cout << "    Directed Hausdorff:       "
              << filter->GetDirectedHausdorffDistance( (*it).first )
              << " (in pixels) " << std::endl;
    std::cout << "    Contour mean distance:    "
              << filter->GetContourMeanDistance( (*it).first )
              << " (in pixel spacing)" << std::endl;
    std::cout << "    Directed contour mean:    "
              << filter->GetDirectedContourMeanDistance( (*it).first )
              << " (in pixel spacing)" << std::endl;

    ++it;
    }

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " ImageDimension image1 "
              << "image2" << std::endl;
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) )
    {
    case 2:
      LabelOverlapMeasures<2>( argc, argv );
      break;
    case 3:
      LabelOverlapMeasures<3>( argc, argv );
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
}
