
#include "antsUtilities.h"
#include <algorithm>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkLabelGeometryImageFilter.h"

#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>

namespace ants
{
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
//   antscout << "Number of labels: " << labelGeometryFilter->GetNumberOfLabels() << std::endl;
//   antscout << "Label geometry measures." << std::endl;
  antscout << std::left << std::setw( 7 )  << "Label"
           << std::left << std::setw( 10 ) << "Volume"
           << std::left << std::setw( 15 ) << "Eccentricity"
           << std::left << std::setw( 15 ) << "Elongation"
           << std::left << std::setw( 15 ) << "Orientation"
           << std::left << std::setw( 30 ) << "Centroid"
           << std::left << std::setw( 30 ) << "Axes Length"
           << std::left << std::setw( 30 ) << "Bounding Box";
  if( filter->GetIntensityInput() )
    {
    antscout << std::left << std::setw( 20 )  << "Integrated Int."
             << std::left << std::setw( 30 ) << "Weighted Centroid";
    }
  antscout << std::endl;
  for( allLabelsIt = allLabels.begin(); allLabelsIt != allLabels.end(); allLabelsIt++ )
    {
    if( *allLabelsIt == 0 )
      {
      continue;
      }
    antscout << std::setw( 7 ) << *allLabelsIt;
    antscout << std::setw( 10 ) << filter->GetVolume( *allLabelsIt );
    antscout << std::setw( 15 ) << filter->GetEccentricity( *allLabelsIt );
    antscout << std::setw( 15 ) << filter->GetElongation( *allLabelsIt );
    antscout << std::setw( 15 ) << filter->GetOrientation( *allLabelsIt );

    std::stringstream oss;
    oss << filter->GetCentroid( *allLabelsIt );
    antscout << std::setw( 30 ) << ( oss.str() ).c_str();
    oss.str( "" );

    oss << filter->GetAxesLength( *allLabelsIt );
    antscout << std::setw( 30 ) << ( oss.str() ).c_str();
    oss.str( "" );

    oss << filter->GetBoundingBox( *allLabelsIt );
    antscout << std::setw( 30 ) << ( oss.str() ).c_str();
    oss.str( "" );

//     antscout << filter->GetMajorAxisLength( *allLabelsIt ) << "\t";
//     antscout << filter->GetMinorAxisLength( *allLabelsIt ) << "\t";
    if( filter->GetIntensityInput() )
      {
      oss << filter->GetIntegratedIntensity( *allLabelsIt );
      antscout << std::setw( 20 ) << ( oss.str() ).c_str();
      oss.str( "" );

      oss << filter->GetWeightedCentroid( *allLabelsIt );
      antscout << std::setw( 30 ) << ( oss.str() ).c_str();
      oss.str( "" );
      }
    antscout << std::endl;
    }

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int LabelGeometryMeasures( std::vector<std::string> args, std::ostream* out_stream = NULL )
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

  if( argc < 3 )
    {
    antscout << "Usage: " << argv[0] << " imageDimension labelImage [intensityImage]"
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
      LabelGeometryMeasures<2>( argc, argv );
      }
      break;
    case 3:
      {
      LabelGeometryMeasures<3>( argc, argv );
      }
      break;
    default:
      antscout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
