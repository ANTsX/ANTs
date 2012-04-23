/** ANTS Landmarks used to initialize an b-spline displacement field ... */

#include "antsUtilities.h"

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkPointSet.h"
#include "itkVector.h"

#include <string>
#include <vector>

namespace ants
{
template <class TValue>
TValue Convert( std::string optionString );

template <class TValue>
std::vector<TValue> ConvertVector( std::string optionString );

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
      std::istringstream iss2( element );
      iss2 >> value;
      values.push_back( value );
      }
    }
  return values;
}

template <unsigned int ImageDimension>
int LandmarkBasedDisplacementFieldTransformInitializer( int argc, char *argv[] )
{
  typedef float                                RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<RealImageType> ImageReaderType;
  typename ImageReaderType::Pointer fixedReader = ImageReaderType::New();
  fixedReader->SetFileName( argv[1] );
  fixedReader->Update();

  typename RealImageType::Pointer fixedImage = fixedReader->GetOutput();
  typename RealImageType::DirectionType fixedDirection = fixedImage->GetDirection();
  typename RealImageType::DirectionType fixedDirectionInverse( fixedDirection.GetInverse() );
  typename RealImageType::PointType fixedOrigin = fixedImage->GetOrigin();

  typename RealImageType::DirectionType identityDirection;
  identityDirection.SetIdentity();

  typename RealImageType::PointType zeroOrigin;
  zeroOrigin.Fill( 0.0 );

  fixedImage->SetDirection( identityDirection );
  fixedImage->SetOrigin( zeroOrigin );

  typename ImageReaderType::Pointer movingReader = ImageReaderType::New();
  movingReader->SetFileName( argv[2] );
  movingReader->Update();
  typename RealImageType::Pointer movingImage = movingReader->GetOutput();

  typedef itk::Vector<RealType, ImageDimension>  VectorType;
  typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;

  typedef itk::PointSet<long, ImageDimension> PointSetType;

  typename PointSetType::Pointer fixedPoints = PointSetType::New();
  fixedPoints->Initialize();

  std::vector<RealType> fixedLabels;

  itk::ImageRegionIteratorWithIndex<RealImageType> ItF( fixedImage, fixedImage->GetLargestPossibleRegion() );

  unsigned int fixedCount = 0;
  for( ItF.GoToBegin(); !ItF.IsAtEnd(); ++ItF )
    {
    if( ItF.Get() != 0 )
      {
      if( std::find( fixedLabels.begin(), fixedLabels.end(), ItF.Get() ) == fixedLabels.end() )
        {
        fixedLabels.push_back( ItF.Get() );
        }
      typename PointSetType::PointType fixedPoint;
      fixedImage->TransformIndexToPhysicalPoint( ItF.GetIndex(), fixedPoint );
      fixedPoints->SetPointData( fixedCount, ItF.Get() );
      fixedPoints->SetPoint( fixedCount++, fixedPoint );
      }
    }

  typename PointSetType::Pointer movingPoints = PointSetType::New();
  movingPoints->Initialize();

  std::vector<RealType> movingLabels;

  itk::ImageRegionIteratorWithIndex<RealImageType> ItM( movingImage, movingImage->GetLargestPossibleRegion() );
  unsigned int                                     movingCount = 0;
  for( ItM.GoToBegin(); !ItM.IsAtEnd(); ++ItM )
    {
    if( ItM.Get() != 0 )
      {
      if( std::find( movingLabels.begin(), movingLabels.end(), ItM.Get() ) == movingLabels.end() )
        {
        movingLabels.push_back( ItM.Get() );
        }
      typename PointSetType::PointType movingPoint;
      movingImage->TransformIndexToPhysicalPoint( ItM.GetIndex(), movingPoint );
      movingPoint -= ( fixedOrigin.GetVectorFromOrigin() );
      movingPoint = fixedDirectionInverse * movingPoint;
      movingPoints->SetPointData( movingCount, ItM.Get() );
      movingPoints->SetPoint( movingCount++, movingPoint );
      }
    }

  // Get moving center points
  typename PointSetType::Pointer movingCenters = PointSetType::New();
  movingCenters->Initialize();
  for( unsigned int n = 0; n < movingLabels.size(); n++ )
    {
    int currentLabel = movingLabels[n];
    typename PointSetType::PointType center;
    center.Fill( 0 );
    float N = 0;
    typename PointSetType::PointsContainerConstIterator ItP =
      movingPoints->GetPoints()->Begin();
    typename PointSetType::PointDataContainerIterator ItD =
      movingPoints->GetPointData()->Begin();
    while( ItP != movingPoints->GetPoints()->End() )
      {
      if( ItD.Value() == currentLabel )
        {
        typename PointSetType::PointType point = ItP.Value();
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          center[d] += point[d];
          }
        N += 1.0;
        }
      ++ItP;
      ++ItD;
      }
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      center[d] /= N;
      }
    movingCenters->SetPoint( n, center );
    movingCenters->SetPointData( n, currentLabel );
    }

  // Get fixed center points
  typename PointSetType::Pointer fixedCenters = PointSetType::New();
  fixedCenters->Initialize();
  for( unsigned int n = 0; n < fixedLabels.size(); n++ )
    {
    int currentLabel = fixedLabels[n];
    typename PointSetType::PointType center;
    center.Fill( 0 );
    float N = 0;
    typename PointSetType::PointsContainerConstIterator ItP =
      fixedPoints->GetPoints()->Begin();
    typename PointSetType::PointDataContainerIterator ItD =
      fixedPoints->GetPointData()->Begin();
    while( ItP != fixedPoints->GetPoints()->End() )
      {
      if( ItD.Value() == currentLabel )
        {
        typename PointSetType::PointType point = ItP.Value();
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          center[d] += point[d];
          }
        N += 1.0;
        }
      ++ItP;
      ++ItD;
      }
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      center[d] /= N;
      }
    fixedCenters->SetPoint( n, center );
    fixedCenters->SetPointData( n, currentLabel );
    }

  if( fixedCenters->GetNumberOfPoints() != movingCenters->GetNumberOfPoints() )
    {
    std::cerr << "The number of fixed points and moving points must be the same." << std::endl;
    return EXIT_FAILURE;
    }

  // Now match up the center points

  typedef itk::PointSet<VectorType, ImageDimension> DisplacementFieldPointSetType;
  typedef itk::BSplineScatteredDataPointSetToImageFilter
    <DisplacementFieldPointSetType, DisplacementFieldType> BSplineFilterType;
  typedef typename BSplineFilterType::WeightsContainerType WeightsContainerType;

  typename WeightsContainerType::Pointer weights = WeightsContainerType::New();
  weights->Initialize();
  const typename WeightsContainerType::Element boundaryWeight = 1.0e10;
  const typename WeightsContainerType::Element weight = 1.0;

  typename DisplacementFieldPointSetType::Pointer fieldPoints = DisplacementFieldPointSetType::New();
  fieldPoints->Initialize();
  unsigned long count = 0;

  typename PointSetType::PointsContainerConstIterator mIt =
    movingCenters->GetPoints()->Begin();
  typename PointSetType::PointDataContainerIterator mItD =
    movingCenters->GetPointData()->Begin();

  while( mItD != movingCenters->GetPointData()->End() )
    {
    typename PointSetType::PointsContainerConstIterator fIt =
      fixedCenters->GetPoints()->Begin();
    typename PointSetType::PointDataContainerIterator fItD =
      fixedCenters->GetPointData()->Begin();

    while( fItD != fixedCenters->GetPointData()->End() )
      {
      if( fItD.Value() == mItD.Value() )
        {
        typename PointSetType::PointType fpoint = fIt.Value();
        typename PointSetType::PointType mpoint = mIt.Value();

        typename DisplacementFieldType::PointType point;
        VectorType vector;
        typename DisplacementFieldPointSetType::PointType fieldPoint;
        for( unsigned int i = 0; i < ImageDimension; i++ )
          {
          fieldPoint[i] = mpoint[i];
          vector[i] = fpoint[i] - mpoint[i];
          }

        fieldPoints->SetPoint( count, fieldPoint );
        fieldPoints->SetPointData( count, vector );

        weights->InsertElement( count, weight );
        count++;

        break;
        }
      ++fItD;
      ++fIt;
      }

    ++mItD;
    ++mIt;
    }

  bool enforceStationaryBoundary = true;
  if( argc > 7 )
    {
    enforceStationaryBoundary = static_cast<bool>( atoi( argv[7] ) );
    }
  if( enforceStationaryBoundary )
    {
    typename RealImageType::IndexType startIndex = fixedImage->GetLargestPossibleRegion().GetIndex();

    typename RealImageType::SizeType inputSize = fixedImage->GetLargestPossibleRegion().GetSize();
    for( ItF.GoToBegin(); !ItF.IsAtEnd(); ++ItF )
      {
      typename RealImageType::IndexType index = ItF.GetIndex();

      bool isOnStationaryBoundary = false;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        if( index[d] == startIndex[d] || index[d] == startIndex[d] + static_cast<int>( inputSize[d] ) - 1 )
          {
          isOnStationaryBoundary = true;
          break;
          }
        }

      if( isOnStationaryBoundary )
        {
        VectorType vector;

        vector.Fill( 0.0 );

        typename PointSetType::PointType fixedPoint;
        fixedImage->TransformIndexToPhysicalPoint( index, fixedPoint );

        fieldPoints->SetPoint( count, fixedPoint );
        fieldPoints->SetPointData( count, vector );
        weights->InsertElement( count, boundaryWeight );
        count++;
        }
      }
    }

  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();

  unsigned int numberOfLevels = atoi( argv[5] );

  unsigned int splineOrder = 3;
  if( argc > 6 )
    {
    atoi( argv[6] );
    }

  std::vector<unsigned int> meshSize = ConvertVector<unsigned int>( std::string( argv[4] ) );
  typename BSplineFilterType::ArrayType ncps;

  if( meshSize.size() == 1 )
    {
    ncps.Fill( meshSize[0] + splineOrder );
    }
  else if( meshSize.size() == ImageDimension )
    {
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      ncps[d] = meshSize[d] + splineOrder;
      }
    }
  else
    {
    std::cerr << "Invalid meshSize format." << std::endl;
    }

//   std::cout << ncps << std::endl;
//
//   bspliner->DebugOn();
  bspliner->SetOrigin( fixedReader->GetOutput()->GetOrigin() );
  bspliner->SetSpacing( fixedReader->GetOutput()->GetSpacing() );
  bspliner->SetSize( fixedReader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  bspliner->SetGenerateOutputImage( true );
  bspliner->SetNumberOfLevels( numberOfLevels );
  bspliner->SetSplineOrder( splineOrder );
  bspliner->SetNumberOfControlPoints( ncps );
  bspliner->SetInput( fieldPoints );
  bspliner->SetPointWeights( weights );
  bspliner->Update();

  bspliner->GetOutput()->SetOrigin( fixedOrigin );
  bspliner->GetOutput()->SetDirection( fixedDirection );

  typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( bspliner->GetOutput() );
  writer->Update();
  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int ANTSUseLandmarkImagesToGetBSplineDisplacementField( std::vector<std::string> args, std::ostream* out_stream = NULL )
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "ANTSUseLandmarkImagesToGetBSplineDisplacementField" );
  std::remove( args.begin(), args.end(), std::string( "" ) );
  std::remove( args.begin(), args.end(), std::string( "" ) );
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

  if( argc < 4 )
    {
    antscout << "Usage:   " << argv[0]
             << " fixedImageWithLabeledLandmarks  movingImageWithLabeledLandmarks outputDisplacementField "
             << "meshSize[0]xmeshSize[1]x... numberOfLevels [order=3] [enforceStationaryBoundaries=1]"
             << std::endl;
    antscout
      << " we expect the input images to be (1) N-ary  (2) in the same physical space as the images you want to "
      << std::endl;
    antscout << " register and (3 ) to have the same landmark points defined within them ... " << std::endl;
    antscout << " landmarks will be defined from the center of mass of the labels in the input images . " << std::endl;
    antscout << " You can use ITK-snap to generate the label images. " << std::endl;
    return 1;
    }

  // Get the image dimension
  std::string               fn = std::string(argv[1]);
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(fn.c_str(), itk::ImageIOFactory::ReadMode);
  imageIO->SetFileName(fn.c_str() );
  imageIO->ReadImageInformation();

  switch( imageIO->GetNumberOfDimensions() )
    {
    case 2:
      {
      LandmarkBasedDisplacementFieldTransformInitializer<2>(argc, argv);
      }
      break;
    case 3:
      {
      LandmarkBasedDisplacementFieldTransformInitializer<3>(argc, argv);
      }
      break;
    default:
      antscout << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }

  return 0;
}
} // namespace ants
