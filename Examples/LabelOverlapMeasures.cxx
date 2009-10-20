#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkLabelOverlapMeasuresImageFilter.h"

#include <iomanip>

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

  std::cout << "                                  "
            << "************  All Labels *************" << std::endl;
  std::cout << std::setw( 17 ) << "   " << std::setw( 17 ) << "Total"
            << std::setw( 17 ) << "Union (jaccard)"
            << std::setw( 17 ) << "Mean (dice)" << std::setw( 17 ) << "False negative"
            << std::setw( 17 ) << "False positive" << std::endl;
  std::cout << std::setw( 17 ) << "Volumetric";
  std::cout << std::setw( 17 ) << filter->GetVolumeTotalOverlap();
  std::cout << std::setw( 17 ) << filter->GetVolumeUnionOverlap();
  std::cout << std::setw( 17 ) << filter->GetVolumeMeanOverlap();
  std::cout << std::setw( 17 ) << filter->GetVolumeFalseNegativeError();
  std::cout << std::setw( 17 ) << filter->GetVolumeFalsePositiveError();
  std::cout << std::endl;
  std::cout << std::setw( 17 ) << "Surface";
  std::cout << std::setw( 17 ) << filter->GetSurfaceTotalOverlap();
  std::cout << std::setw( 17 ) << filter->GetSurfaceUnionOverlap();
  std::cout << std::setw( 17 ) << filter->GetSurfaceMeanOverlap();
  std::cout << std::setw( 17 ) << filter->GetSurfaceFalseNegativeError();
  std::cout << std::setw( 17 ) << filter->GetSurfaceFalsePositiveError();
  std::cout << std::endl;
  std::cout << "Other related measures." << std::endl;
  std::cout << "  Volume similarity: "
            << filter->GetVolumeSimilarity() << std::endl;
  std::cout << std::endl;

  std::cout << "                               "
            << "************ Individual Labels *************" << std::endl;
  std::cout << "Volumetric overlap measures." << std::endl;
  std::cout << std::setw( 17 ) << "Label"
            << std::setw( 17 ) << "Total"
            << std::setw( 17 ) << "Union (jaccard)"
            << std::setw( 17 ) << "Mean (dice)"
            << std::setw( 17 ) << "False negative"
            << std::setw( 17 ) << "False positive" << std::endl;

  typename FilterType::MapType labelMap = filter->GetLabelSetMeasures();
  typename FilterType::MapType::const_iterator it;
  for( it = labelMap.begin(); it != labelMap.end(); ++it )
    {
    if( (*it).first == 0 )
      {
      continue;
      }

    int label = (*it).first;

    std::cout << std::setw( 17 ) << label;
    std::cout << std::setw( 17 ) << filter->GetVolumeTotalOverlap( label );
    std::cout << std::setw( 17 ) << filter->GetVolumeUnionOverlap( label );
    std::cout << std::setw( 17 ) << filter->GetVolumeMeanOverlap( label );
    std::cout << std::setw( 17 ) << filter->GetVolumeFalseNegativeError( label );
    std::cout << std::setw( 17 ) << filter->GetVolumeFalsePositiveError( label );
    std::cout << std::endl;
    }

  std::cout << "Surface overlap measures." << std::endl;
  std::cout << std::setw( 17 ) << "Label"
            << std::setw( 17 ) << "Total"
            << std::setw( 17 ) << "Union (jaccard)"
            << std::setw( 17 ) << "Mean (dice)"
            << std::setw( 17 ) << "False negative"
            << std::setw( 17 ) << "False positive" << std::endl;
  for( it = labelMap.begin(); it != labelMap.end(); ++it )
    {
    if( (*it).first == 0 )
      {
      continue;
      }
    int label = (*it).first;

    std::cout << std::setw( 17 ) << label;
    std::cout << std::setw( 17 ) << filter->GetSurfaceTotalOverlap( label );
    std::cout << std::setw( 17 ) << filter->GetSurfaceUnionOverlap( label );
    std::cout << std::setw( 17 ) << filter->GetSurfaceMeanOverlap( label );
    std::cout << std::setw( 17 ) << filter->GetSurfaceFalseNegativeError( label );
    std::cout << std::setw( 17 ) << filter->GetSurfaceFalsePositiveError( label );
    std::cout << std::endl;
    }

  std::cout << "Other related measures." << std::endl;
  std::cout << std::setw( 17 ) << "Label"
            << std::setw( 17 ) << "Volume sim."
            << std::setw( 17 ) << "Hausdorff"
            << std::setw( 17 ) << "Dir. Hausdorff"
            << std::setw( 17 ) << "Contour"
            << std::setw( 17 ) << "Dir. Contour" << std::endl;
  for( it = labelMap.begin(); it != labelMap.end(); ++it )
    {
    if( (*it).first == 0 )
      {
      continue;
      }
    int label = (*it).first;

    std::cout << std::setw( 17 ) << label;
    std::cout << std::setw( 17 ) << filter->GetVolumeSimilarity( label );
    std::cout << std::setw( 17 ) << filter->GetHausdorffDistance( label );
    std::cout << std::setw( 17 ) << filter->GetDirectedHausdorffDistance( label );
    std::cout << std::setw( 17 ) << filter->GetContourMeanDistance( label );
    std::cout << std::setw( 17 ) << filter->GetDirectedContourMeanDistance( label );
    std::cout << std::endl;
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
