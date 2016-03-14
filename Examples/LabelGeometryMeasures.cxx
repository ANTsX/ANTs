#include "antsUtilities.h"
#include "antsAllocImage.h"
#include <algorithm>
#include "ReadWriteData.h"

#include "itkAffineTransform.h"
#include "itkCSVArray2DDataObject.h"
#include "itkCSVArray2DFileReader.h"
#include "itkCSVNumericObjectFileWriter.h"
#include "itkImage.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkLabelPerimeterEstimationCalculator.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
#include "itkTransformFileWriter.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>

namespace ants
{
template <unsigned int ImageDimension>
int LabelGeometryMeasures( int argc, char * argv[] )
{
  typedef unsigned int                          LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;

  typedef float                                RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typename LabelImageType::Pointer labelImage = LabelImageType::New();
  ReadImage<LabelImageType>( labelImage, argv[2] );

  typename RealImageType::Pointer intensityImage = RealImageType::New();
  bool intensityImageUsed = false;
  if( argc > 3 )
    {
    try
      {
      ReadImage<RealImageType>( intensityImage, argv[3] );
      intensityImageUsed = true;
      }
    catch( ... )
      {
      }
    }

  bool outputCSVFormat = false;
  if( argc > 4  )
    {
    outputCSVFormat = true;
    }

  typedef itk::LabelGeometryImageFilter<LabelImageType, RealImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( labelImage );
  if( intensityImageUsed )
    {
    filter->SetIntensityInput( intensityImage );
    }
  filter->CalculatePixelIndicesOff();
  filter->CalculateOrientedBoundingBoxOff();
  filter->CalculateOrientedLabelRegionsOff();
  // These generate optional outputs.
//   filter->CalculatePixelIndicesOn();
//   filter->CalculateOrientedBoundingBoxOn();;
//   filter->CalculateOrientedLabelRegionsOn();
  filter->Update();

  typedef itk::LabelPerimeterEstimationCalculator<LabelImageType> AreaFilterType;
  typename AreaFilterType::Pointer areafilter = AreaFilterType::New();
  areafilter->SetImage( labelImage );
  areafilter->SetFullyConnected( false );
  areafilter->Compute();

  typename FilterType::LabelsType allLabels = filter->GetLabels();
  std::sort( allLabels.begin(), allLabels.end() );

  if( outputCSVFormat )
    {
    typename FilterType::LabelsType::iterator allLabelsIt;

    std::vector<std::string>   columnHeaders;

    columnHeaders.push_back( std::string( "Label" ) );
    columnHeaders.push_back( std::string( "VolumeInVoxels" ) );
    columnHeaders.push_back( std::string( "SurfaceAreaInMillimetersSquared" ) );
    columnHeaders.push_back( std::string( "Eccentricity" ) );
    columnHeaders.push_back( std::string( "Elongation" ) );
    columnHeaders.push_back( std::string( "Orientation" ) );
    columnHeaders.push_back( std::string( "Centroid_x" ) );
    columnHeaders.push_back( std::string( "Centroid_y" ) );
    if( ImageDimension == 3 )
      {
      columnHeaders.push_back( std::string( "Centroid_z" ) );
      }
    columnHeaders.push_back( std::string( "AxesLength_x" ) );
    columnHeaders.push_back( std::string( "AxesLength_y" ) );
    if( ImageDimension == 3 )
      {
      columnHeaders.push_back( std::string( "AxesLength_z" ) );
      }
    columnHeaders.push_back( std::string( "BoundingBoxLower_x" ) );
    columnHeaders.push_back( std::string( "BoundingBoxLower_y" ) );
    if( ImageDimension == 3 )
      {
      columnHeaders.push_back( std::string( "BoundingBoxLower_z" ) );
      }
    columnHeaders.push_back( std::string( "BoundingBoxUpper_x" ) );
    columnHeaders.push_back( std::string( "BoundingBoxUpper_y" ) );
    if( ImageDimension == 3 )
      {
      columnHeaders.push_back( std::string( "BoundingBoxUpper_z" ) );
      }

    if( filter->GetIntensityInput() )
      {
      columnHeaders.push_back( std::string( "IntegratedIntensity" ) );
      columnHeaders.push_back( std::string( "WeightedCentroid_x" ) );
      columnHeaders.push_back( std::string( "WeightedCentroid_y" ) );
      if( ImageDimension == 3 )
        {
        columnHeaders.push_back( std::string( "WeightedCentroid_z" ) );
        }
      }

    std::vector<std::string>   rowHeaders;
    std::ostringstream convert;// stream used for the conversion
    for( allLabelsIt = allLabels.begin(); allLabelsIt != allLabels.end(); allLabelsIt++ )
      {
      convert << *allLabelsIt;   // insert the textual representation of 'Number' in the characters in the stream
      rowHeaders.push_back( convert.str() ); // set 'Result' to the contents of the stream
      }

    vnl_matrix<double> measures( allLabels.size(), columnHeaders.size() - 1 );

    unsigned int rowIndex = 0;
    for( allLabelsIt = allLabels.begin(); allLabelsIt != allLabels.end(); allLabelsIt++ )
      {
      unsigned int columnIndex = 0;
      measures( rowIndex, columnIndex++ ) = filter->GetVolume( *allLabelsIt );
      measures( rowIndex, columnIndex++ ) = areafilter->GetPerimeter( *allLabelsIt );
      measures( rowIndex, columnIndex++ ) = filter->GetEccentricity( *allLabelsIt );
      measures( rowIndex, columnIndex++ ) = filter->GetElongation( *allLabelsIt );
      measures( rowIndex, columnIndex++ ) = filter->GetOrientation( *allLabelsIt );
      measures( rowIndex, columnIndex++ ) = filter->GetCentroid( *allLabelsIt )[0];
      measures( rowIndex, columnIndex++ ) = filter->GetCentroid( *allLabelsIt )[1];
      if( ImageDimension == 3 )
        {
        measures( rowIndex, columnIndex++ ) = filter->GetCentroid( *allLabelsIt )[2];
        }
      measures( rowIndex, columnIndex++ ) = filter->GetAxesLength( *allLabelsIt )[0];
      measures( rowIndex, columnIndex++ ) = filter->GetAxesLength( *allLabelsIt )[1];
      if( ImageDimension == 3 )
        {
        measures( rowIndex, columnIndex++ ) = filter->GetAxesLength( *allLabelsIt )[2];
        }

      unsigned int arrayIndex = 0;

      measures( rowIndex, columnIndex++ ) = filter->GetBoundingBox( *allLabelsIt )[arrayIndex++];
      measures( rowIndex, columnIndex++ ) = filter->GetBoundingBox( *allLabelsIt )[arrayIndex++];
      if( ImageDimension == 3 )
        {
        measures( rowIndex, columnIndex++ ) = filter->GetBoundingBox( *allLabelsIt )[arrayIndex++];
        }
      measures( rowIndex, columnIndex++ ) = filter->GetBoundingBox( *allLabelsIt )[arrayIndex++];
      measures( rowIndex, columnIndex++ ) = filter->GetBoundingBox( *allLabelsIt )[arrayIndex++];
      if( ImageDimension == 3 )
        {
        measures( rowIndex, columnIndex++ ) = filter->GetBoundingBox( *allLabelsIt )[arrayIndex++];
        }

      if( filter->GetIntensityInput() )
        {
        measures( rowIndex, columnIndex++ ) = filter->GetIntegratedIntensity( *allLabelsIt );
        measures( rowIndex, columnIndex++ ) = filter->GetWeightedCentroid( *allLabelsIt )[0];
        measures( rowIndex, columnIndex++ ) = filter->GetWeightedCentroid( *allLabelsIt )[1];
        if( ImageDimension == 3 )
          {
          measures( rowIndex, columnIndex++ ) = filter->GetWeightedCentroid( *allLabelsIt )[2];
          }
        }
      rowIndex++;
      }

    typedef itk::CSVNumericObjectFileWriter<double, 1, 1> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[4] );
    writer->SetColumnHeaders( columnHeaders );
    writer->SetRowHeaders( rowHeaders );
    writer->SetInput( &measures );
    try
      {
      writer->Write();
      }
    catch( itk::ExceptionObject& exp )
      {
      std::cerr << "Exception caught!" << std::endl;
      std::cerr << exp << std::endl;
      return EXIT_FAILURE;
      }
    }
  else
    {
    typename FilterType::LabelsType::iterator allLabelsIt;
  //   std::cout << "Number of labels: " << labelGeometryFilter->GetNumberOfLabels() << std::endl;
  //   std::cout << "Label geometry measures." << std::endl;
    std::cout << std::left << std::setw( 7 ) << "Label"
             << std::left << std::setw( 10 ) << "Volume"
             << std::left << std::setw( 15 ) << "SurfArea(mm^2)"
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
      std::cout << std::setw( 15 ) << areafilter->GetPerimeter( *allLabelsIt );
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
    }

  return EXIT_SUCCESS;
}

template <class TransformType, class LabelImageType>
unsigned int GetNumberOfLabelVoxelsInUpperRightQuadrant( TransformType *affineTransform, LabelImageType *labelImage )
{
  typedef itk::ResampleImageFilter<LabelImageType, LabelImageType> ResampleFilterType;
  typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetTransform( affineTransform->GetInverseTransform() );
  resampler->SetSize( labelImage->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputSpacing( labelImage->GetSpacing() );
  resampler->SetOutputOrigin( labelImage->GetOrigin() );

  typedef itk::NearestNeighborInterpolateImageFunction<LabelImageType, double> InterpolatorType;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  resampler->SetInterpolator( interpolator );
  resampler->SetInput( labelImage );
  resampler->Update();

  typedef typename LabelImageType::RegionType RegionType;

  typename RegionType::IndexType index;
  typename RegionType::SizeType size;
  for( unsigned int d = 0; d < LabelImageType::ImageDimension; d++ )
    {
    index[d] = std::floor( resampler->GetOutput()->GetLargestPossibleRegion().GetIndex()[d] + 0.5 * resampler->GetOutput()->GetLargestPossibleRegion().GetSize()[d] );
    size[d] = std::floor( 0.5 * resampler->GetOutput()->GetLargestPossibleRegion().GetSize()[d] );
    }

  index[1] = std::floor( resampler->GetOutput()->GetLargestPossibleRegion().GetIndex()[1] );
  size[1] = std::floor( resampler->GetOutput()->GetLargestPossibleRegion().GetSize()[1] );

  RegionType region;
  region.SetIndex( index );
  region.SetSize( size );

  unsigned int numberOfVoxels = 0;

  itk::ImageRegionIterator<LabelImageType> It( resampler->GetOutput(), region );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if( It.Get() == 1 )
      {
      numberOfVoxels++;
      }
    }

  return numberOfVoxels;
}

int GetReorientationRigidTransform( int argc, char * argv[] )
{
  const unsigned int ImageDimension = 2;

  if( argc < 4 )
    {
    std::cerr << "Not enough arguments.  See help." << std::endl;
    return EXIT_FAILURE;
    }

  typedef unsigned int                          LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;

  LabelImageType::Pointer labelImage = LabelImageType::New();
  ReadImage<LabelImageType>( labelImage, argv[2] );

  typedef itk::LabelGeometryImageFilter<LabelImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( labelImage );

//   filter->CalculatePixelIndicesOff();
//   filter->CalculateOrientedBoundingBoxOff();
//   filter->CalculateOrientedLabelRegionsOff();
//   // These generate optional outputs.
  filter->CalculatePixelIndicesOn();
  filter->CalculateOrientedBoundingBoxOn();;
  filter->CalculateOrientedLabelRegionsOn();
  filter->Update();

  typedef itk::AffineTransform<double, ImageDimension> TransformType;
  TransformType::Pointer transform = TransformType::New();
  TransformType::MatrixType reorientationMatrix( filter->GetRotationMatrix( 2 ) );
  TransformType::CenterType center;

  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    center[i] = filter->GetCentroid( 2 )[i] *
      filter->GetInput()->GetSpacing()[i];
    }

  TransformType::OutputVectorType translation;
  translation.Fill( 0 );
  transform->SetCenter( center );
  transform->SetTranslation( translation );
  transform->SetMatrix( reorientationMatrix );
  bool doReflections = true;
  float scale = 1.0;
  if( argc > 4 )
    {
    doReflections = static_cast<bool>( atoi( argv[4] ) );
    }
  if( argc > 5 )
    {
    scale = atof( argv[5] );
    }
  if( doReflections )
    {
    // check which has a higher concentraion

    TransformType::MatrixType flipXMatrix;
    flipXMatrix.SetIdentity();
    flipXMatrix( 0, 0 ) = -1.0 * scale;

    TransformType::MatrixType flipYMatrix;
    flipYMatrix.SetIdentity();
    flipYMatrix( 1, 1 ) = -1.0 * scale;

    TransformType::MatrixType testRotationMatrix;
    std::vector<int> numberOfVoxelsInUpperRightQuadrant( 2 );

    // check current rotation matrix
    testRotationMatrix = reorientationMatrix;
    transform->SetMatrix( testRotationMatrix );

    numberOfVoxelsInUpperRightQuadrant[0] =
      GetNumberOfLabelVoxelsInUpperRightQuadrant<TransformType, LabelImageType>( transform, labelImage );

    // check flip x
    testRotationMatrix = flipXMatrix * reorientationMatrix;
    transform->SetMatrix( testRotationMatrix );
    numberOfVoxelsInUpperRightQuadrant[1] =
      GetNumberOfLabelVoxelsInUpperRightQuadrant<TransformType, LabelImageType>( transform, labelImage );

      // check flip y
//     testRotationMatrix = flipYMatrix * reorientationMatrix;
//     transform->SetMatrix( testRotationMatrix );
//     numberOfVoxelsInUpperRightQuadrant[2] =
//       GetNumberOfLabelVoxelsInUpperRightQuadrant<TransformType, LabelImageType>( transform, labelImage );

    // check flip x and flip y
//     testRotationMatrix = flipXMatrix * flipYMatrix * reorientationMatrix;
//     transform->SetMatrix( testRotationMatrix );
//     numberOfVoxelsInUpperRightQuadrant[3] =
//       GetNumberOfLabelVoxelsInUpperRightQuadrant<TransformType, LabelImageType>( transform, labelImage );

    std::vector<int>::iterator result = std::max_element( numberOfVoxelsInUpperRightQuadrant.begin(), numberOfVoxelsInUpperRightQuadrant.end() );
    unsigned int index = std::distance( numberOfVoxelsInUpperRightQuadrant.begin(), result );

    switch( index )
      {
      case 0:  default:
        {
        testRotationMatrix = reorientationMatrix * scale;
        transform->SetMatrix( testRotationMatrix );
        break;
        }
      case 1:
        {
        testRotationMatrix = flipXMatrix * reorientationMatrix * scale;
        transform->SetMatrix( testRotationMatrix );
        break;
        }
//       case 2:
//         {
//         testRotationMatrix = flipYMatrix * reorientationMatrix * scale;
//         transform->SetMatrix( testRotationMatrix );
//         break;
//         }
//       case 3:
//         {
//         testRotationMatrix = flipXMatrix * flipYMatrix * reorientationMatrix * scale;
//         transform->SetMatrix( testRotationMatrix );
//         break;
//         }
      }
    }
  typedef itk::TransformFileWriter TransformWriterType;
  TransformWriterType::Pointer transformWriter = TransformWriterType::New();
  transformWriter->SetInput( transform );
  transformWriter->SetFileName( argv[3] );
  transformWriter->Update();
  return EXIT_SUCCESS;
}




// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int LabelGeometryMeasures( std::vector<std::string> args, std::ostream* itkNotUsed( out_stream ) )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "LabelGeometryMeasures" );

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
  argv[argc] = ITK_NULLPTR;
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

  if( argc < 3 )
    {
    std::cout << "Usage 1: " << argv[0] << " imageDimension labelImage [intensityImage=none] [csvFile]" << std::endl;
//     std::cout << "Usage 2: " << argv[0] << " X singleLabelImage outputTransform <doReflection=1> <scaleFactor=1>" << std::endl;

    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }
  std::string xstring=std::string("X");
  std::string xstringtest(argv[1]);
  if( xstring.compare(xstringtest) == 0 )
    {
    GetReorientationRigidTransform( argc, argv );
    }
  else
    {
    switch( atoi( argv[1] ) )
      {
      case 2:
        {
        LabelGeometryMeasures<2>( argc, argv );
        }
        break;
      case 3:
        {
        LabelGeometryMeasures<3>( argc, argv );
        }
        break;
      default:
        std::cout << "Unsupported dimension" << std::endl;
        return EXIT_FAILURE;
      }
    }
  return EXIT_SUCCESS;
}
} // namespace ants
