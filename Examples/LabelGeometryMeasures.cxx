#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkLabelGeometryImageFilter.h"

#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>

template <unsigned int ImageDimension>
int LabelGeometryMeasures( int argc, char * argv[] )
{
  typedef int                                   LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;
  typedef itk::ImageFileReader<LabelImageType>  LabelReaderType;

  typedef float                                RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::ImageFileReader<RealImageType>  ReaderType;

  typename LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( argv[2] );
  labelReader->Update();

  typename ReaderType::Pointer reader = ReaderType::New();
  if( argc > 3 )
    {
    reader->SetFileName( argv[3] );
    reader->Update();
    }

  typedef itk::LabelGeometryImageFilter<LabelImageType, RealImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( labelReader->GetOutput() );
  if( argc > 3 )
    {
    filter->SetIntensityInput( reader->GetOutput() );
    }
  filter->CalculatePixelIndicesOff();
  filter->CalculateOrientedBoundingBoxOff();
  filter->CalculateOrientedLabelRegionsOff();
  // These generate optional outputs.
//   filter->CalculatePixelIndicesOn();
//   filter->CalculateOrientedBoundingBoxOn();;
//   filter->CalculateOrientedLabelRegionsOn();
  filter->Update();

  typename FilterType::LabelsType allLabels = filter->GetLabels();
  typename FilterType::LabelsType::iterator allLabelsIt;
//   std::cout << "Number of labels: " << labelGeometryFilter->GetNumberOfLabels() << std::endl;
//   std::cout << "Label geometry measures." << std::endl;
  std::cout << std::left << std::setw( 7 )  << "Label"
            << std::left << std::setw( 10 ) << "Volume"
            << std::left << std::setw( 15 ) << "Eccentricity"
            << std::left << std::setw( 15 ) << "Elongation"
            << std::left << std::setw( 15 ) << "Orientation"
            << std::left << std::setw( 30 ) << "Centroid"
            << std::left << std::setw( 30 ) << "Axes Length"
            << std::left << std::setw( 30 ) << "Bounding Box";
  if( filter->GetIntensityInput() )
    {
    std::cout << std::left << std::setw( 20 )  << "Integrated Int."
              << std::left << std::setw( 30 ) << "Weighted Centroid";
    }
  std::cout << std::endl;
  for( allLabelsIt = allLabels.begin(); allLabelsIt != allLabels.end(); allLabelsIt++ )
    {
    if( *allLabelsIt == 0 )
      {
      continue;
      }
    std::cout << std::setw( 7 ) << *allLabelsIt;
    std::cout << std::setw( 10 ) << filter->GetVolume( *allLabelsIt );
    std::cout << std::setw( 15 ) << filter->GetEccentricity( *allLabelsIt );
    std::cout << std::setw( 15 ) << filter->GetElongation( *allLabelsIt );
    std::cout << std::setw( 15 ) << filter->GetOrientation( *allLabelsIt );

    std::stringstream oss;
    oss << filter->GetCentroid( *allLabelsIt );
    std::cout << std::setw( 30 ) << ( oss.str() ).c_str();
    oss.str( "" );

    oss << filter->GetAxesLength( *allLabelsIt );
    std::cout << std::setw( 30 ) << ( oss.str() ).c_str();
    oss.str( "" );

    oss << filter->GetBoundingBox( *allLabelsIt );
    std::cout << std::setw( 30 ) << ( oss.str() ).c_str();
    oss.str( "" );

//     std::cout << filter->GetMajorAxisLength( *allLabelsIt ) << "\t";
//     std::cout << filter->GetMinorAxisLength( *allLabelsIt ) << "\t";
    if( filter->GetIntensityInput() )
      {
      oss << filter->GetIntegratedIntensity( *allLabelsIt );
      std::cout << std::setw( 20 ) << ( oss.str() ).c_str();
      oss.str( "" );

      oss << filter->GetWeightedCentroid( *allLabelsIt );
      std::cout << std::setw( 30 ) << ( oss.str() ).c_str();
      oss.str( "" );
      }
    std::cout << std::endl;
    }

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension labelImage [intensityImage]"
              << std::endl;
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) )
    {
    case 2:
      LabelGeometryMeasures<2>( argc, argv );
      break;
    case 3:
      LabelGeometryMeasures<3>( argc, argv );
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }
}
