#include "antsUtilities.h"

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkPointSet.h"

#include <stdio.h>
#include <vector>
#include <fstream>
#include <string>

namespace ants
{

template <unsigned int PointDimension>
int FitBSplineCurveToPoints( unsigned int argc, char *argv[] )
{
  typedef float RealType;

  typedef itk::Vector<RealType, PointDimension> VectorType;
  typedef itk::Image<VectorType, 1> CurveImageType;

  typedef itk::PointSet<VectorType, 1> PointSetType;
  typename PointSetType::Pointer pointSet = PointSetType::New();
  pointSet->Initialize();

  typedef itk::BSplineScatteredDataPointSetToImageFilter
     <PointSetType, CurveImageType>  FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  typename FilterType::WeightsContainerType::Pointer weights = FilterType::WeightsContainerType::New();

  RealType totalDistance = 0.0;

  std::ifstream file( argv[2] );
  std::string line;

  unsigned int count = 0;
  if( file.is_open() )
    {
    while( std::getline( file, line ) )
      {
      VectorType vector( 0.0 );
      float weight = 1.0;

      std::string delimiter = ",";
      size_t pos = 0;
      std::string token;

      unsigned int dimensionCount = 0;
      while ( ( pos = line.find( delimiter ) ) != std::string::npos )
        {
        token = line.substr( 0, pos );

        std::istringstream stream( token );
        float element;
        stream >> element;

        if( dimensionCount == PointDimension )
          {
          weight = element;
          break;
          }
        else
          {
          vector[dimensionCount++] = element;
          }
        line.erase( 0, pos + delimiter.length() );
        }

      pointSet->SetPointData( count, vector );

      if ( count > 0 )
        {
        VectorType previous( 0.0 );
        pointSet->GetPointData( count-1, &previous );
        totalDistance += ( previous - vector ).GetNorm();
        }

      typename PointSetType::PointType point;
      point[0] = 0.0;
      pointSet->SetPoint( count, point );

      weights->InsertElement( count, weight );
      count++;
      }
    }

  RealType cumSum = 0.0;
  for ( unsigned int i = 1; i < pointSet->GetNumberOfPoints(); i++ )
    {
    VectorType vector( 0.0 ), previous( 0.0 );
    pointSet->GetPointData( i, &vector );
    pointSet->GetPointData( i-1, &previous );

    cumSum += ( vector - previous ).GetNorm();
    typename PointSetType::PointType point;
    point[0] = cumSum / totalDistance;

    pointSet->SetPoint( i, point );
    }

  filter->SetInput( pointSet );
  filter->SetGenerateOutputImage( true );

  typename CurveImageType::PointType origin;
  origin.Fill( 0.0 );
  filter->SetOrigin( origin );
  typename CurveImageType::SpacingType spacing;
  spacing[0] = 0.001;
  if ( argc > 6 )
    {
    spacing[0] = atof( argv[6] );
    }

  filter->SetSpacing( spacing );
  typename CurveImageType::SizeType size;
  size[0] = static_cast<unsigned int>( 1.0 / spacing[0] + 1 );
  filter->SetSize( size );
  typename FilterType::ArrayType order;
  order[0] = 3;
  if ( argc > 3 )
    {
    order[0] = atoi( argv[3] );
    }
  filter->SetSplineOrder( order );
  typename FilterType::ArrayType ncps;
  ncps[0] = order[0] + 1;
  if ( argc > 5 )
    {
    ncps[0] = atoi( argv[5] );
    }
  filter->SetNumberOfControlPoints( ncps );
  typename FilterType::ArrayType nlevels;
  nlevels[0] = 5;
  if ( argc > 4 )
    {
    nlevels[0] = atoi( argv[4] );
    }
  filter->SetNumberOfLevels( nlevels );
  typename FilterType::ArrayType close;
  close[0] = false;
  if ( argc > 7 )
    {
    close[0] = atoi( argv[7] );
    }
  filter->SetCloseDimension( close );

  filter->Update();

  itk::ImageRegionIterator<CurveImageType> It(
    filter->GetOutput(), filter->GetOutput()->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    VectorType vector = It.Get();
    for( unsigned int d = 0; d < PointDimension-1; d++ )
      {
      std::cout << vector[d] << ",";
      }
    std::cout << vector[PointDimension-1] << std::endl;
    }

//   {
//   std::string filename = std::string( argv[2] ) + std::string( "_cps.txt" );
//   std::ofstream ostr( filename.c_str() );
//   ostr << "0 0 0 0" << std::endl;
//
//   itk::ImageRegionIterator<CurveImageType> It(
//     filter->GetPhiLattice(), filter->GetPhiLattice()->GetLargestPossibleRegion() );
//   for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
//     {
//     ostr << It.Get()[0] << " " << It.Get()[1] << " " << It.Get()[2] << " 1" << std::endl;
//     }
//   ostr << "0 0 0 0" << std::endl;
//   ostr.close();
//   }

  return 0;
}


// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int FitBSplineCurveToPoints( std::vector<std::string> args, std::ostream* /*out_stream = NULL */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "FitBSplineCurveToPoints" );

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

  // antscout->set_stream( out_stream );
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " pointDimension inputLandmarksFile "
      << " [order=3] [nlevels=10] "
      << " [numberOfControlPoints=4] [sampleSpacing=0.001] [closed?=0]" << std::endl;
    std::cout << "  Note2:  1. Points are assumed to be parametrically ordered. " << std::endl
              << "          2. The last column (pointDimension+1) is used for weights." << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     FitBSplineCurveToPoints<2>( argc, argv );
     break;
   case 3:
     FitBSplineCurveToPoints<3>( argc, argv );
     break;
   case 4:
     FitBSplineCurveToPoints<4>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }

  return EXIT_SUCCESS;
}
} // namespace ants

