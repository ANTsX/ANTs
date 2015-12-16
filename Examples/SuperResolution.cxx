#include "antsUtilities.h"
#include <algorithm>

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkPointSet.h"
#include "itkTimeProbe.h"
#include "itkVector.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

#include "ReadWriteData.h"

namespace ants
{
template <unsigned int ImageDimension>
int SuperResolution( unsigned int argc, char *argv[] )
{

  typedef float                                   RealType;
  typedef itk::Image<RealType, ImageDimension>    ImageType;

  typedef itk::Vector<RealType, 1>                         ScalarType;
  typedef itk::Image<ScalarType, ImageDimension>           ScalarImageType;
  typedef itk::PointSet<ScalarType, ImageDimension>        PointSetType;
  typedef itk::BSplineScatteredDataPointSetToImageFilter
    <PointSetType, ScalarImageType>                        BSplineFilterType;

  typename ImageType::Pointer domainImage = ITK_NULLPTR;
  ReadImage<ImageType>( domainImage, argv[3] );

  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();

  typename PointSetType::Pointer bsplinePoints = PointSetType::New();
  bsplinePoints->Initialize();

  typename BSplineFilterType::WeightsContainerType::Pointer weights =
    BSplineFilterType::WeightsContainerType::New();
  weights->Initialize();

  unsigned int splineOrder = 3;
  typename BSplineFilterType::ArrayType numberOfLevels;
  typename BSplineFilterType::ArrayType ncps;

  std::vector<unsigned int> nlevels = ConvertVector<unsigned int>( std::string( argv[5] ) );
  if ( nlevels.size() == 1 )
    {
    numberOfLevels.Fill( nlevels[0] );
    }
  else if ( nlevels.size() == ImageDimension )
    {
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      numberOfLevels[d] = nlevels[d];
      }
    }
  else
    {
    std::cerr << "Invalid nlevels format." << std::endl;
    return EXIT_FAILURE;
    }

  std::vector<unsigned int> meshSize = ConvertVector<unsigned int>( std::string( argv[4] ) );
  if ( meshSize.size() == 1 )
    {
    ncps.Fill( meshSize[0] + splineOrder );
    }
  else if ( meshSize.size() == ImageDimension )
    {
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      ncps[d] = meshSize[d] + splineOrder;
      }
    }
  else
    {
    std::cerr << "Invalid ncps format." << std::endl;
    return EXIT_FAILURE;
    }

  unsigned int N = 0;
  for( unsigned int n = 6; n < argc; n++ )
    {
    typename ImageType::Pointer inputImage = ITK_NULLPTR;
    ReadImage<ImageType>( inputImage, argv[n] );

    itk::ImageRegionConstIteratorWithIndex<ImageType> It( inputImage, inputImage->GetRequestedRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      typename ImageType::PointType imagePoint;
      inputImage->TransformIndexToPhysicalPoint( It.GetIndex(), imagePoint );

      typename ImageType::IndexType index;
      bool isInside = domainImage->TransformPhysicalPointToIndex( imagePoint, index );

      if( !isInside )
        {
        continue;
        }

      ScalarType scalar;
      scalar[0] = It.Get();

      bsplinePoints->SetPointData( N, scalar );
      bsplinePoints->SetPoint( N, imagePoint );
      weights->InsertElement( N, 1 );

      N++;
      }
    }

  itk::TimeProbe timer;
  timer.Start();

  bspliner->SetOrigin( domainImage->GetOrigin() );
  bspliner->SetSpacing( domainImage->GetSpacing() );
  bspliner->SetSize( domainImage->GetRequestedRegion().GetSize() );
  bspliner->SetGenerateOutputImage( true );
  bspliner->SetNumberOfLevels( numberOfLevels );
  bspliner->SetSplineOrder( splineOrder );
  bspliner->SetNumberOfControlPoints( ncps );
  bspliner->SetInput( bsplinePoints );
  bspliner->SetPointWeights( weights );
  bspliner->Update();

  timer.Stop();

  std::cout << "Elapsed Time:  " << timer.GetMean() << std::endl;

  typedef itk::VectorIndexSelectionCastImageFilter<ScalarImageType, ImageType> SelectorType;
  typename SelectorType::Pointer selector = SelectorType::New();
  selector->SetInput( bspliner->GetOutput() );
  selector->SetIndex( 0 );
  selector->Update();

  WriteImage<ImageType>( selector->GetOutput(), argv[2] );

  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int SuperResolution( std::vector<std::string> args, std::ostream* /*out_stream = NULL */ )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "SuperResolution" );

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

  if ( argc < 7 )
    {
    std::cerr << argv[0] << " imageDimension outputImage domainImage meshSize numberOfLevels inputImage1 ... inputImageN" << std::endl;

    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

  const int ImageDimension = static_cast<int>( atoi( argv[1] ) );

  switch( ImageDimension )
     {
     case 2:
       return SuperResolution<2>( argc, argv );
       break;
     case 3:
       return SuperResolution<3>( argc, argv );
       break;
     case 4:
       return SuperResolution<4>( argc, argv );
       break;
     default:
        std::cerr << "Unsupported dimension" << std::endl;
        exit( EXIT_FAILURE );
     }
  return EXIT_SUCCESS;
}
} // namespace ants

