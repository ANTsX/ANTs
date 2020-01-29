/** ANTS Landmarks used to initialize an b-spline displacement field ... */

#include "antsUtilities.h"

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkContinuousIndex.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImportImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkPointSet.h"

#include <string>
#include <vector>

namespace ants
{
template <unsigned int ImageDimension>
int LandmarkBasedDisplacementFieldTransformInitializer( int argc, char *argv[] )
{
  typedef float                                 RealType;
  typedef unsigned int                          LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<LabelImageType> ImageReaderType;
  typename ImageReaderType::Pointer fixedReader = ImageReaderType::New();
  fixedReader->SetFileName( argv[1] );
  fixedReader->Update();

  typename LabelImageType::Pointer fixedImage = fixedReader->GetOutput();
  typename LabelImageType::DirectionType fixedDirection = fixedImage->GetDirection();
  typename LabelImageType::DirectionType fixedDirectionInverse( fixedDirection.GetInverse() );

  typename LabelImageType::DirectionType identityDirection;
  identityDirection.SetIdentity();

  const typename LabelImageType::RegionType & bufferedRegion = fixedImage->GetBufferedRegion();
  const itk::SizeValueType numberOfPixels = bufferedRegion.GetNumberOfPixels();

  const bool filterHandlesMemory = false;

  typedef itk::ImportImageFilter<LabelType, ImageDimension> ImporterType;
  typename ImporterType::Pointer importer = ImporterType::New();
  importer->SetImportPointer( const_cast<LabelType *>( fixedImage->GetBufferPointer() ), numberOfPixels, filterHandlesMemory );
  importer->SetRegion( fixedImage->GetBufferedRegion() );
  importer->SetOrigin( fixedImage->GetOrigin() );
  importer->SetSpacing( fixedImage->GetSpacing() );
  importer->SetDirection( identityDirection );
  importer->Update();

  const typename ImporterType::OutputImageType * parametricInputImage = importer->GetOutput();

  typename ImageReaderType::Pointer movingReader = ImageReaderType::New();
  movingReader->SetFileName( argv[2] );
  movingReader->Update();
  typename LabelImageType::Pointer movingImage = movingReader->GetOutput();

///////////////////////////////////////////////////////////////////////////////////////////

  typedef itk::Vector<RealType, ImageDimension>  VectorType;
  typedef itk::Image<VectorType, ImageDimension> DisplacementFieldType;

  typedef itk::PointSet<LabelType, ImageDimension> PointSetType;

  typename PointSetType::Pointer fixedPoints = PointSetType::New();
  fixedPoints->Initialize();

  std::vector<LabelType> fixedLabels;

  itk::ImageRegionIteratorWithIndex<LabelImageType> ItF( fixedImage, fixedImage->GetLargestPossibleRegion() );

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
  std::sort( fixedLabels.begin(), fixedLabels.end() );

  typename PointSetType::Pointer movingPoints = PointSetType::New();
  movingPoints->Initialize();

  std::vector<LabelType> movingLabels;

  itk::ImageRegionIteratorWithIndex<LabelImageType> ItM( movingImage, movingImage->GetLargestPossibleRegion() );
  unsigned int                                      movingCount = 0;
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
      movingPoints->SetPointData( movingCount, ItM.Get() );
      movingPoints->SetPoint( movingCount++, movingPoint );
      }
    }
  std::sort( movingLabels.begin(), movingLabels.end() );

  // Get moving center points
  typename PointSetType::Pointer movingCenters = PointSetType::New();
  movingCenters->Initialize();
  for( unsigned int n = 0; n < movingLabels.size(); n++ )
    {
    LabelType currentLabel = movingLabels[n];
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
        N += 1.0f;
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
    LabelType currentLabel = fixedLabels[n];
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
        N += 1.0f;
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

  // Read in the optional label weights

  std::vector<float>     labelWeights;
  std::vector<LabelType> userLabels;

  bool useWeights = false;

  unsigned int labelCount = 0;
  if( argc > 8 )
    {
    useWeights = true;

    std::fstream labelStr( argv[8] );

    if( labelStr.is_open() )
      {
      while( !labelStr.eof() )
        {
        char line[256];
        labelStr.getline( line, 256 );

        std::string lineString = std::string( line );
        std::size_t pos = lineString.find( ',' );

        RealType value;
        if( pos == std::string::npos )
          {
          std::istringstream iss( lineString );
          iss >> value;
          labelWeights.push_back( value );
          userLabels.push_back( movingLabels[labelCount++] );
          }
        else
          {
          unsigned int localLabel;

          std::string        element = lineString.substr( 0, pos );
          std::istringstream iss( element );
          iss >> localLabel;
          userLabels.push_back( localLabel );

          element = lineString.substr( pos + 1, lineString.length() );
          std::istringstream iss2( element );
          iss2 >> value;
          labelWeights.push_back( value );
          }
        }

      labelStr.close();
      }
    else
      {
      std::cerr << "File " << argv[8] << " cannot be opened." << std::endl;
      return EXIT_FAILURE;
      }
    }

  // Now match up the center points

  typedef itk::PointSet<VectorType, ImageDimension> DisplacementFieldPointSetType;
  typedef itk::BSplineScatteredDataPointSetToImageFilter
    <DisplacementFieldPointSetType, DisplacementFieldType> BSplineFilterType;
  typedef typename BSplineFilterType::WeightsContainerType WeightsContainerType;

  typename WeightsContainerType::Pointer weights = WeightsContainerType::New();
  weights->Initialize();
  const typename WeightsContainerType::Element boundaryWeight = 1.0e10;
  typename WeightsContainerType::Element weight = 1.0;

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

        VectorType vector;

        typename LabelImageType::PointType fixedPhysicalPoint;
        for( unsigned int i = 0; i < ImageDimension; i++ )
          {
          fixedPhysicalPoint[i] = fpoint[i];
          vector[i] = mpoint[i] - fpoint[i];
          }

        itk::ContinuousIndex<double, ImageDimension> fixedCidx;
        fixedImage->TransformPhysicalPointToContinuousIndex( fixedPhysicalPoint, fixedCidx );

        typename DisplacementFieldType::PointType fieldPoint;
        parametricInputImage->TransformContinuousIndexToPhysicalPoint( fixedCidx, fieldPoint );

        fieldPoints->SetPoint( count, fieldPoint );
        fieldPoints->SetPointData( count, vector );

        if( useWeights )
          {
          auto it = std::find( userLabels.begin(), userLabels.end(), mItD.Value() );
          if( it != userLabels.end() )
            {
            weights->InsertElement( count, labelWeights[it - userLabels.begin()] );
            }
          else
            {
            std::cerr << "Unspecified label " << mItD.Value() << " in specified user label weights." << std::endl;
            return EXIT_FAILURE;
            }
          }
        else
          {
          weights->InsertElement( count, weight );
          }

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
    enforceStationaryBoundary = static_cast<bool>( std::stoi( argv[7] ) );
    }
  if( enforceStationaryBoundary )
    {
    typename LabelImageType::IndexType startIndex2 = fixedImage->GetLargestPossibleRegion().GetIndex();

    typename LabelImageType::SizeType inputSize2 = fixedImage->GetLargestPossibleRegion().GetSize();
    for( ItF.GoToBegin(); !ItF.IsAtEnd(); ++ItF )
      {
      typename LabelImageType::IndexType index = ItF.GetIndex();

      bool isOnStationaryBoundary = false;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        if( index[d] == startIndex2[d] || index[d] == startIndex2[d] + static_cast<int>( inputSize2[d] ) - 1 )
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
        parametricInputImage->TransformIndexToPhysicalPoint( index, fixedPoint );

        fieldPoints->SetPoint( count, fixedPoint );
        fieldPoints->SetPointData( count, vector );
        weights->InsertElement( count, boundaryWeight );
        count++;
        }
      }
    }

  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();

  unsigned int numberOfLevels = std::stoi( argv[5] );

  unsigned int splineOrder = 3;
  if( argc > 6 )
    {
    splineOrder = std::stoi( argv[6] );
    }

  std::vector<unsigned int> meshSize = ConvertVector<unsigned int>( std::string( argv[4] ) );
  typename BSplineFilterType::ArrayType ncps;
  ncps.Fill( 0 );

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
  bspliner->SetOrigin( fixedImage->GetOrigin() );
  bspliner->SetSpacing( fixedImage->GetSpacing() );
  bspliner->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
  bspliner->SetDirection( fixedImage->GetDirection() );
  bspliner->SetGenerateOutputImage( true );
  bspliner->SetNumberOfLevels( numberOfLevels );
  bspliner->SetSplineOrder( splineOrder );
  bspliner->SetNumberOfControlPoints( ncps );
  bspliner->SetInput( fieldPoints );
  bspliner->SetPointWeights( weights );
  bspliner->Update();

  typedef itk::VectorLinearInterpolateImageFunction<DisplacementFieldType, RealType> InterpolatorType;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage( bspliner->GetOutput() );

  std::cout << "Distance errors:" << std::endl;

  mIt = movingCenters->GetPoints()->Begin();
  mItD = movingCenters->GetPointData()->Begin();

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

        VectorType displacement = ( mpoint - fpoint );

        typename InterpolatorType::PointType point;
        for( unsigned int i = 0; i < ImageDimension; i++ )
          {
          point[i] = fpoint[i];
          }
        VectorType vector = interpolator->Evaluate( point );

        RealType error = ( vector - displacement ).GetNorm();
        std::cout << "  " << fItD.Value() << ": " << error << std::endl;

        break;
        }
      ++fItD;
      ++fIt;
      }

    ++mItD;
    ++mIt;
    }

  typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( bspliner->GetOutput() );
  writer->Update();
  return EXIT_SUCCESS;
}

// entry point for the library; parameter 'args' is equivalent to 'argv' in (argc,argv) of commandline parameters to
// 'main()'
int ANTSUseLandmarkImagesToGetBSplineDisplacementField( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */)
{
  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle
  args.insert( args.begin(), "ANTSUseLandmarkImagesToGetBSplineDisplacementField" );
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
  argv[argc] = nullptr;
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

  if( argc < 4 )
    {
    std::cerr << "Usage:   " << argv[0]
             << " fixedImageWithLabeledLandmarks  movingImageWithLabeledLandmarks outputDisplacementField "
             <<
      "meshSize[0]xmeshSize[1]x... numberOfLevels [order=3] [enforceStationaryBoundaries=1] [landmarkWeights]"
             << std::endl;
    std::cerr
      << " we expect the input images to be (1) N-ary  (2) in the same physical space as the images you want to "
      << std::endl;
    std::cerr << " register and (3 ) to have the same landmark points defined within them ... " << std::endl;
    std::cerr << " landmarks will be defined from the center of mass of the labels in the input images . " << std::endl;
    std::cerr << " You can use ITK-snap to generate the label images. " << std::endl;
    std::cerr << " The optional landmarks weights are read from a text file where each row is either:" << std::endl;
    std::cerr << " \"label,labelWeight\" or " << std::endl;
    std::cerr << " \"labelWeight\" or " << std::endl;
    std::cerr
      << " If the latter format is used, the label weights are assumed to be arranged in ascending order by label."
      << std::endl;
    if( argc >= 2 &&
        ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
      {
      return EXIT_SUCCESS;
      }
    return EXIT_FAILURE;
    }

  // Get the image dimension
  std::string               fn = std::string(argv[1]);
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(fn.c_str(), itk::ImageIOFactory::FileModeEnum::ReadMode);
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
      std::cerr << "Unsupported dimension" << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
} // namespace ants
