#include "antsAllocImage.h"
#include "antsCommandLineParser.h"
#include "antsUtilities.h"
#include "ReadWriteData.h"

#include "itkAffineTransform.h"
#include "itkConjugateGradientLineSearchOptimizerv4.h"
#include "itkEuler2DTransform.h"
#include "itkEuler3DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkCorrelationImageToImageMetricv4.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImageMomentsCalculator.h"
#include "itkImageToImageMetricv4.h"
#include "itkJointHistogramMutualInformationImageToImageMetricv4.h"
#include "itkLandmarkBasedTransformInitializer.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkMultiScaleLaplacianBlobDetectorImageFilter.h"
#include "itkMultiStartOptimizerv4.h"
#include "itkRegistrationParameterScalesFromPhysicalShift.h"
#include "itkRigid2DTransform.h"
#include "itkSimilarity2DTransform.h"
#include "itkSimilarity3DTransform.h"
#include "itkTranslationTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "itkTransformFileWriter.h"

#include "vnl/vnl_cross.h"
#include "vnl/vnl_inverse.h"

namespace ants
{

// ##########################################################################
//      Transform traits to generalize the different linear transforms
// ##########################################################################

template <typename TComputeType, unsigned int ImageDimension>
class RigidTransformTraits
{
  // Don't worry about the fact that the default option is the
  // affine Transform, that one will not actually be instantiated.
public:
typedef itk::AffineTransform<TComputeType, ImageDimension> TransformType;
};

template <typename TComputeType, unsigned int ImageDimension>
class LandmarkRigidTransformTraits
{
  // Don't worry about the fact that the default option is the
  // affine Transform, that one will not actually be instantiated.
public:
typedef itk::AffineTransform<TComputeType, ImageDimension> TransformType;
};

template <>
class RigidTransformTraits<double, 2>
{
public:
typedef itk::Euler2DTransform<double> TransformType;
};

template <>
class RigidTransformTraits<float, 2>
{
public:
typedef itk::Euler2DTransform<float> TransformType;
};

template <>
class RigidTransformTraits<double, 3>
{
public:
  // typedef itk::VersorRigid3DTransform<double>    TransformType;
  // typedef itk::QuaternionRigidTransform<double>  TransformType;
typedef itk::Euler3DTransform<double> TransformType;
};

template <>
class RigidTransformTraits<float, 3>
{
public:
  // typedef itk::VersorRigid3DTransform<float>    TransformType;
  // typedef itk::QuaternionRigidTransform<float>  TransformType;
typedef itk::Euler3DTransform<float> TransformType;
};

template <>
class LandmarkRigidTransformTraits<double, 2>
{
public:
typedef itk::Rigid2DTransform<double> TransformType;
};

template <>
class LandmarkRigidTransformTraits<float, 2>
{
public:
typedef itk::Rigid2DTransform<float> TransformType;
};

template <>
class LandmarkRigidTransformTraits<double, 3>
{
public:
typedef itk::VersorRigid3DTransform<double> TransformType;
};

template <>
class LandmarkRigidTransformTraits<float, 3>
{
public:
typedef itk::VersorRigid3DTransform<float> TransformType;
};

template <typename TComputeType, unsigned int ImageDimension>
class SimilarityTransformTraits
{
// Don't worry about the fact that the default option is the
// affine Transform, that one will not actually be instantiated.
public:
typedef itk::AffineTransform<TComputeType, ImageDimension> TransformType;
};

template <>
class SimilarityTransformTraits<double, 2>
{
public:
typedef itk::Similarity2DTransform<double> TransformType;
};

template <>
class SimilarityTransformTraits<float, 2>
{
public:
typedef itk::Similarity2DTransform<float> TransformType;
};

template <>
class SimilarityTransformTraits<double, 3>
{
public:
typedef itk::Similarity3DTransform<double> TransformType;
};

template <>
class SimilarityTransformTraits<float, 3>
{
public:
typedef itk::Similarity3DTransform<float> TransformType;
};

// ##########################################################################
// ##########################################################################
template<typename TImage, typename TGradientImage, typename TInterpolator, typename TReal>
TReal PatchCorrelation( itk::NeighborhoodIterator<TImage> fixedNeighborhood,
                        itk::NeighborhoodIterator<TImage> movingNeighborhood,
                        std::vector<unsigned int> activeIndex,
                        std::vector<TReal> weights,
                        typename TGradientImage::Pointer fixedGradientImage,
                        typename TGradientImage::Pointer movingGradientImage,
                        typename TInterpolator::Pointer movingInterpolator )
{
  typedef TReal                                          RealType;
  typedef TImage                                         ImageType;

  const unsigned int ImageDimension = ImageType::ImageDimension;

  typedef typename ImageType::PointType                  PointType;
  typedef itk::CovariantVector<RealType, ImageDimension> GradientPixelType;
  typedef vnl_vector<RealType>                           VectorType;
  typedef typename ImageType::IndexType                  IndexType;

  unsigned int           numberOfIndices = activeIndex.size();
  std::vector<PointType> movingImagePatchPoints;
  VectorType             fixedSamples( numberOfIndices, 0 );
  VectorType             movingSamples( numberOfIndices, 0 );

  vnl_matrix<RealType> fixedGradientMatrix( numberOfIndices, ImageDimension, 0.0 );
  vnl_matrix<RealType> movingGradientMatrix( numberOfIndices, ImageDimension, 0.0 );

  PointType movingPointCentroid( 0.0 );

  RealType weight = itk::NumericTraits<RealType>::OneValue() / static_cast<RealType>( numberOfIndices );
  for( unsigned int i = 0; i < numberOfIndices; i++ )
    {
    fixedSamples[i] = fixedNeighborhood.GetPixel( activeIndex[i] );
    movingSamples[i] = movingNeighborhood.GetPixel( activeIndex[i] );

    IndexType fixedIndex = fixedNeighborhood.GetIndex( activeIndex[i] );
    IndexType movingIndex = movingNeighborhood.GetIndex( activeIndex[i] );

    if( fixedGradientImage->GetRequestedRegion().IsInside( fixedIndex ) &&
        movingGradientImage->GetRequestedRegion().IsInside( movingIndex ) )
      {
      GradientPixelType fixedGradient = fixedGradientImage->GetPixel( fixedIndex ) * weights[i];
      GradientPixelType movingGradient = movingGradientImage->GetPixel( movingIndex ) * weights[i];
      for( unsigned int j = 0; j < ImageDimension; j++ )
        {
        fixedGradientMatrix(i, j) = fixedGradient[j];
        movingGradientMatrix(i, j) = movingGradient[j];
        }
      PointType movingPoint( 0.0 );
      movingGradientImage->TransformIndexToPhysicalPoint( movingIndex, movingPoint );
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        movingPointCentroid[d] += movingPoint[d] * static_cast<typename PointType::CoordRepType>( weight );
        }
      movingImagePatchPoints.push_back( movingPoint );
      }
    else
      {
      return 0.0;
      }
    }

  fixedSamples -= fixedSamples.mean();
  movingSamples -= movingSamples.mean();

  RealType fixedSd = std::sqrt( fixedSamples.squared_magnitude() );

  // compute patch orientation

  vnl_matrix<RealType> fixedCovariance = fixedGradientMatrix.transpose() * fixedGradientMatrix;
  vnl_matrix<RealType> movingCovariance = movingGradientMatrix.transpose() * movingGradientMatrix;

  vnl_symmetric_eigensystem<RealType> fixedImagePrincipalAxes( fixedCovariance );
  vnl_symmetric_eigensystem<RealType> movingImagePrincipalAxes( movingCovariance );

  vnl_vector<RealType> fixedPrimaryEigenVector;
  vnl_vector<RealType> fixedSecondaryEigenVector;
  vnl_vector<RealType> fixedTertiaryEigenVector;
  vnl_vector<RealType> movingPrimaryEigenVector;
  vnl_vector<RealType> movingSecondaryEigenVector;

  vnl_matrix<RealType> B;

  if( ImageDimension == 2 )
    {
    fixedPrimaryEigenVector = fixedImagePrincipalAxes.get_eigenvector( 1 );
    movingPrimaryEigenVector = movingImagePrincipalAxes.get_eigenvector( 1 );

    B = outer_product( movingPrimaryEigenVector, fixedPrimaryEigenVector );
    }
  else if( ImageDimension == 3 )
    {
    fixedPrimaryEigenVector = fixedImagePrincipalAxes.get_eigenvector( 2 );
    fixedSecondaryEigenVector = fixedImagePrincipalAxes.get_eigenvector( 1 );

    movingPrimaryEigenVector = movingImagePrincipalAxes.get_eigenvector( 2 );
    movingSecondaryEigenVector = movingImagePrincipalAxes.get_eigenvector( 1 );

    B = outer_product( movingPrimaryEigenVector, fixedPrimaryEigenVector ) +
      outer_product( movingSecondaryEigenVector, fixedSecondaryEigenVector );
    }

  vnl_svd<RealType> wahba( B );
  vnl_matrix<RealType> A = wahba.V() * wahba.U().transpose();

  bool patchIsInside = true;
  for( unsigned int i = 0; i < numberOfIndices; i++ )
    {
    PointType movingImagePoint = movingImagePatchPoints[i];
    vnl_vector<RealType> movingImageVector( ImageDimension, 0.0 );
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      movingImagePoint[d] -= movingPointCentroid[d];
      movingImageVector[d] = movingImagePoint[d];
      }

    vnl_vector<RealType> movingImagePointRotated = A * movingImageVector;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      movingImagePoint[d] = static_cast<typename PointType::CoordRepType>( movingImagePointRotated[d] ) + movingPointCentroid[d];
      }
    if( movingInterpolator->IsInsideBuffer( movingImagePoint ) )
      {
      movingSamples[i] = movingInterpolator->Evaluate( movingImagePoint );
      }
    else
      {
      patchIsInside = false;
      break;
      }
    }

  RealType correlation = 0.0;
  if( patchIsInside )
    {
    movingSamples -= movingSamples.mean();
    RealType movingSd = std::sqrt( movingSamples.squared_magnitude() );
    correlation = inner_product( fixedSamples, movingSamples ) / ( fixedSd * movingSd );

    if( std::isnan( correlation ) || std::isinf( correlation )  )
      {
      correlation = 0.0;
      }
    }

  return correlation;
}

template<typename TImage, typename TBlobFilter>
void GetBlobCorrespondenceMatrix( typename TImage::Pointer fixedImage,  typename TImage::Pointer movingImage,
                                  typename TBlobFilter::BlobsListType fixedBlobs, typename TBlobFilter::BlobsListType movingBlobs,
                                  float sigma, unsigned int radiusValue,
                                  vnl_matrix<float> & correspondenceMatrix )
{
  typedef TImage                                 ImageType;
  typedef float                                  RealType;
  typedef typename ImageType::IndexType          IndexType;
  typedef itk::NeighborhoodIterator<ImageType>   NeighborhoodIteratorType;

  const unsigned int ImageDimension = ImageType::ImageDimension;

  typedef itk::CovariantVector<RealType, ImageDimension>       GradientVectorType;
  typedef itk::Image<GradientVectorType, ImageDimension>       GradientImageType;
  typedef itk::GradientRecursiveGaussianImageFilter<ImageType,
                                            GradientImageType> GradientImageFilterType;
  typedef typename GradientImageFilterType::Pointer            GradientImageFilterPointer;

  typename GradientImageFilterType::Pointer fixedGradientFilter = GradientImageFilterType::New();
  fixedGradientFilter->SetInput( fixedImage );
  fixedGradientFilter->SetSigma( sigma );

  typename GradientImageType::Pointer fixedGradientImage = fixedGradientFilter->GetOutput();
  fixedGradientImage->Update();
  fixedGradientImage->DisconnectPipeline();

  GradientImageFilterPointer movingGradientFilter = GradientImageFilterType::New();
  movingGradientFilter->SetInput( movingImage );
  movingGradientFilter->SetSigma( sigma );

  typename GradientImageType::Pointer movingGradientImage = movingGradientFilter->GetOutput();
  movingGradientImage->Update();
  movingGradientImage->DisconnectPipeline();

  correspondenceMatrix.set_size( fixedBlobs.size(), movingBlobs.size() );
  correspondenceMatrix.fill( 0 );

  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( radiusValue );

  NeighborhoodIteratorType ItF( radius, fixedImage, fixedImage->GetLargestPossibleRegion() );
  NeighborhoodIteratorType ItM( radius, movingImage, movingImage->GetLargestPossibleRegion() );

  IndexType zeroIndex;
  zeroIndex.Fill( radiusValue );
  ItF.SetLocation( zeroIndex );

  std::vector<unsigned int> activeIndex;
  std::vector<RealType> weights;
  RealType weightSum = 0.0;
  for( unsigned int i = 0; i < ItF.Size(); i++ )
    {
    IndexType index = ItF.GetIndex( i );
    RealType distance = 0.0;
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      distance += itk::Math::sqr( index[j] - zeroIndex[j] );
      }
    distance = std::sqrt( distance );
    if( distance <= radiusValue )
      {
      activeIndex.push_back( i );
      RealType weight = std::exp( -distance / itk::Math::sqr( radiusValue ) );
      weights.push_back( weight );
      weightSum += ( weight );
      }
    }
  for(float & weight : weights)
    {
    weight /= weightSum;
    }

  typedef itk::LinearInterpolateImageFunction<ImageType, RealType> ScalarInterpolatorType;
  typename ScalarInterpolatorType::Pointer movingInterpolator =  ScalarInterpolatorType::New();
  movingInterpolator->SetInputImage( movingImage );

  for( unsigned int i = 0; i < fixedBlobs.size(); i++ )
    {
    IndexType fixedIndex = fixedBlobs[i]->GetCenter();
    if( fixedImage->GetPixel( fixedIndex ) > static_cast<typename ImageType::PixelType>( 1.0e-4 ) )
      {
      ItF.SetLocation( fixedIndex );
      for( unsigned int j = 0; j < movingBlobs.size(); j++ )
        {
        IndexType movingIndex = movingBlobs[j]->GetCenter();
        if( movingImage->GetPixel( movingIndex ) > static_cast<typename ImageType::PixelType>( 1.0e-4 ) )
          {
          ItM.SetLocation( movingIndex );

          RealType correlation = PatchCorrelation<ImageType, GradientImageType, ScalarInterpolatorType, RealType>
            ( ItF, ItM, activeIndex, weights, fixedGradientImage, movingGradientImage, movingInterpolator );

          if( correlation < itk::NumericTraits<RealType>::ZeroValue() )
            {
            correlation = itk::NumericTraits<RealType>::ZeroValue();
            }

          correspondenceMatrix(i, j) = correlation;
          }
        }
      }
    }

  return;
}

template <typename TImage, typename TTransform>
typename TTransform::Pointer GetTransformFromFeatureMatching( typename TImage::Pointer fixedImage,
  typename TImage::Pointer movingImage, unsigned int numberOfBlobsToExtract, unsigned int numberOfBlobsToMatch )
{
  typedef float                                     RealType;
  typedef TImage                                    ImageType;
  typedef TTransform                                TransformType;
  typedef typename ImageType::PointType             PointType;

  /////////////////////////////////////////////////////////////////
  //
  //         Values taken from Brian in ImageMath.cxx
  //
  /////////////////////////////////////////////////////////////////

  unsigned int radiusValue = 20;
  RealType     gradientSigma = 1.0;
  RealType     maxRadiusDifferenceAllowed = 0.25;
  RealType     distancePreservationThreshold = 0.1;
  RealType     minimumNumberOfNeighborhoodNodes = 3;

  RealType minScale = std::pow( 1.0, 1.0 );
  RealType maxScale = std::pow( 2.0, 10.0 );
  unsigned int stepsPerOctave = 10;

  ///////////////////////////////////////////////////////////////

  typedef itk::MultiScaleLaplacianBlobDetectorImageFilter<ImageType> BlobFilterType;
  typedef typename BlobFilterType::BlobPointer                       BlobPointer;
  typedef std::pair<BlobPointer, BlobPointer>                        BlobPairType;
  typedef typename BlobFilterType::BlobsListType                     BlobsListType;

  typename BlobFilterType::Pointer blobFixedImageFilter = BlobFilterType::New();
  blobFixedImageFilter->SetStartT( minScale );
  blobFixedImageFilter->SetEndT( maxScale );
  blobFixedImageFilter->SetStepsPerOctave( stepsPerOctave );
  blobFixedImageFilter->SetNumberOfBlobs( numberOfBlobsToExtract );
  blobFixedImageFilter->SetInput( fixedImage );
  blobFixedImageFilter->Update();

  BlobsListType fixedImageBlobs = blobFixedImageFilter->GetBlobs();

  typename BlobFilterType::Pointer blobMovingImageFilter = BlobFilterType::New();
  blobMovingImageFilter->SetStartT( minScale );
  blobMovingImageFilter->SetEndT( maxScale );
  blobMovingImageFilter->SetStepsPerOctave( stepsPerOctave );
  blobMovingImageFilter->SetNumberOfBlobs( numberOfBlobsToExtract );
  blobMovingImageFilter->SetInput( movingImage );
  blobMovingImageFilter->Update();

  BlobsListType movingImageBlobs = blobMovingImageFilter->GetBlobs();

  if( movingImageBlobs.empty() || fixedImageBlobs.empty() )
    {
    std::cerr << "The moving image or fixed image blobs list is empty." << std::endl;
    return nullptr;
    }

  vnl_matrix<RealType> correspondenceMatrix;

  /////////////////////////////////////////////////////////////////
  //
  //         Get the correspondence matrix based on local image patches
  //
  /////////////////////////////////////////////////////////////////

  GetBlobCorrespondenceMatrix<ImageType, BlobFilterType>
     ( fixedImage, movingImage, fixedImageBlobs, movingImageBlobs, gradientSigma, radiusValue, correspondenceMatrix );

  /////////////////////////////////////////////////////////////////
  //
  //         Pick the first numberOfBlobsToMatch pairs based on the correlation matrix
  //
  /////////////////////////////////////////////////////////////////

  std::vector<BlobPairType> blobPairs;

  BlobPointer bestBlob = nullptr;

  unsigned int matchPoint = 1;
  unsigned int fixedCount = 0;
  while( ( matchPoint <= numberOfBlobsToMatch ) && ( fixedCount < fixedImageBlobs.size() ) )
    {
    unsigned int maxPair = correspondenceMatrix.arg_max();
    auto maxRow = static_cast<unsigned int>( maxPair / correspondenceMatrix.cols() );
    unsigned int maxCol = maxPair - maxRow * correspondenceMatrix.cols();
    BlobPointer fixedBlob = fixedImageBlobs[maxRow];
    bestBlob = movingImageBlobs[maxCol];

    if( bestBlob && bestBlob->GetObjectRadius() > 1 )
      {
      if( static_cast<RealType>( std::fabs( bestBlob->GetObjectRadius() - fixedBlob->GetObjectRadius() ) ) < maxRadiusDifferenceAllowed )
        {
        if( fixedImage->GetPixel( fixedBlob->GetCenter() ) > static_cast<typename ImageType::PixelType>( 1.0e-4 ) &&
            movingImage->GetPixel( bestBlob->GetCenter() ) > static_cast<typename ImageType::PixelType>( 1.0e-4 ) )
          {
          BlobPairType blobPairing = std::make_pair( fixedBlob, bestBlob );
          blobPairs.push_back( blobPairing );
          matchPoint++;
          }
        }
      }
    correspondenceMatrix.set_row( maxRow, correspondenceMatrix.get_row( 0 ).fill( 0 ) );
    correspondenceMatrix.set_column( maxCol, correspondenceMatrix.get_column( 0 ).fill( 0 ) );
    fixedCount++;
    }

  /////////////////////////////////////////////////////////////////
  //
  //         For every blob pair compute the relative distance from its
  //         neighbors in both the fixed and moving images.
  //
  /////////////////////////////////////////////////////////////////

  vnl_matrix<RealType> distanceRatios( blobPairs.size(), blobPairs.size(), 0.0 );
  for( unsigned int i = 0; i < blobPairs.size(); i++ )
    {
    PointType fixedPoint( 0.0 );
    fixedImage->TransformIndexToPhysicalPoint( blobPairs[i].first->GetCenter(), fixedPoint );

    PointType movingPoint( 0.0 );
    movingImage->TransformIndexToPhysicalPoint( blobPairs[i].second->GetCenter(), movingPoint );

    for( unsigned int j = 0; j < blobPairs.size(); j++ )
      {
      PointType fixedNeighborPoint( 0.0 );
      fixedImage->TransformIndexToPhysicalPoint( blobPairs[j].first->GetCenter(), fixedNeighborPoint );

      PointType movingNeighborPoint( 0.0 );
      movingImage->TransformIndexToPhysicalPoint( blobPairs[j].second->GetCenter(), movingNeighborPoint );

      distanceRatios(i, j) = distanceRatios(j, i) =
        fixedNeighborPoint.SquaredEuclideanDistanceTo( movingNeighborPoint ) /
        fixedPoint.SquaredEuclideanDistanceTo( movingPoint );
      }
    }

  /////////////////////////////////////////////////////////////////
  //
  //         Keep those blob pairs which have close to the same
  //         distance in both the fixed and moving images.
  //
  /////////////////////////////////////////////////////////////////

  typename TransformType::Pointer transform = TransformType::New();

  typedef itk::LandmarkBasedTransformInitializer<TransformType, ImageType, ImageType> TransformInitializerType;

  typedef typename TransformInitializerType::LandmarkPointContainer PointsContainerType;
  PointsContainerType fixedLandmarks;
  PointsContainerType movingLandmarks;

  for( unsigned int i = 0; i < blobPairs.size(); i++ )
    {
    unsigned int neighborhoodCount = 0;
    for( unsigned int j = 0; j < blobPairs.size(); j++ )
      {
      if( ( j != i ) && ( std::fabs( distanceRatios( i, j ) - itk::NumericTraits<RealType>::OneValue() ) <  distancePreservationThreshold ) )
        {
        neighborhoodCount++;
        }
      }
    if( neighborhoodCount >= minimumNumberOfNeighborhoodNodes )
      {
      PointType fixedPoint( 0.0 );
      fixedImage->TransformIndexToPhysicalPoint( blobPairs[i].first->GetCenter(), fixedPoint );
      fixedLandmarks.push_back( fixedPoint );

      PointType movingPoint( 0.0 );
      movingImage->TransformIndexToPhysicalPoint( blobPairs[i].second->GetCenter(), movingPoint );
      movingLandmarks.push_back( movingPoint );
      }
    }

  /////////////////////////////////////////////////////////////////
  //
  //         Compute initial transform from the landmarks
  //
  /////////////////////////////////////////////////////////////////

  typename TransformInitializerType::Pointer transformInitializer = TransformInitializerType::New();
  transformInitializer->SetFixedLandmarks( fixedLandmarks );
  transformInitializer->SetMovingLandmarks( movingLandmarks );

  transformInitializer->SetTransform( transform );
  transformInitializer->InitializeTransform();

  return transform;
}

template <unsigned int ImageDimension>
int antsAI( itk::ants::CommandLineParser *parser )
{
  typedef double RealType;
  typedef float  PixelType;
  typedef float  FloatType;

  typedef itk::Vector<FloatType, ImageDimension>      VectorType;
  typedef itk::Image<PixelType, ImageDimension>       ImageType;

  typedef itk::AffineTransform<RealType, ImageDimension>                                 AffineTransformType;
  typedef typename RigidTransformTraits<RealType, ImageDimension>::TransformType         RigidTransformType;
  typedef typename SimilarityTransformTraits<RealType, ImageDimension>::TransformType    SimilarityTransformType;
  typedef typename LandmarkRigidTransformTraits<RealType, ImageDimension>::TransformType LandmarkRigidTransformType;

  enum SamplingStrategyType { NONE, REGULAR, RANDOM };

  bool verbose = false;
  itk::ants::CommandLineParser::OptionType::Pointer verboseOption = parser->GetOption( "verbose" );
  if( verboseOption && verboseOption->GetNumberOfFunctions() )
    {
    verbose = parser->Convert<bool>( verboseOption->GetFunction( 0 )->GetName() );
    }

  /////////////////////////////////////////////////////////////////
  //
  //         Read in the metric type and fixed/moving images
  //
  /////////////////////////////////////////////////////////////////

  typename ImageType::Pointer fixedImage;
  typename ImageType::Pointer movingImage;

  std::string metric;
  unsigned int numberOfBins = 32;
  SamplingStrategyType samplingStrategy = NONE;
  RealType samplingPercentage = 1.0;

  itk::ants::CommandLineParser::OptionType::Pointer metricOption = parser->GetOption( "metric" );
  if( metricOption && metricOption->GetNumberOfFunctions() )
    {
    metric = metricOption->GetFunction( 0 )->GetName();
    ConvertToLowerCase( metric );
    if( metricOption->GetFunction( 0 )->GetNumberOfParameters() < 2 )
      {
      if( verbose )
        {
        std::cerr << "Fixed and/or moving images not specified." << std::endl;
        }
      return EXIT_FAILURE;
      }
    ReadImage<ImageType>( fixedImage, ( metricOption->GetFunction( 0 )->GetParameter( 0 ) ).c_str() );
    ReadImage<ImageType>( movingImage, ( metricOption->GetFunction( 0 )->GetParameter( 1 ) ).c_str() );

    if( metricOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
      {
      numberOfBins = parser->Convert<unsigned int>( metricOption->GetFunction( 0 )->GetParameter( 2 ) );
      }
    if( metricOption->GetFunction( 0 )->GetNumberOfParameters() > 3 )
      {
      std::string samplingStrategyString = metricOption->GetFunction( 0 )->GetParameter( 3 );
      ConvertToLowerCase( samplingStrategyString );
      if( std::strcmp( samplingStrategyString.c_str(), "none" ) == 0 )
        {
        samplingStrategy = NONE;
        }
      else if( std::strcmp( samplingStrategyString.c_str(), "regular" ) == 0 )
        {
        samplingStrategy = REGULAR;
        }
      else if( std::strcmp( samplingStrategyString.c_str(), "random" ) == 0 )
        {
        samplingStrategy = RANDOM;
        }
      else
        {
        if( verbose )
          {
          std::cerr << "Unrecognized sampling strategy." << std::endl;
          }
        return EXIT_FAILURE;
        }
      }
    if( metricOption->GetFunction( 0 )->GetNumberOfParameters() > 4 )
      {
      samplingPercentage = parser->Convert<RealType>( metricOption->GetFunction( 0 )->GetParameter( 4 ) );
      }
    }
  else
    {
    std::cerr << "Input metric not specified." << std::endl;
    return EXIT_FAILURE;
    }

  /////////////////////////////////////////////////////////////////
  //
  //         Read in convergence options
  //
  /////////////////////////////////////////////////////////////////

  unsigned int numberOfIterations = 0;
  unsigned int convergenceWindowSize = 5;
  RealType convergenceThreshold = 1e-6;

  itk::ants::CommandLineParser::OptionType::Pointer convergenceOption = parser->GetOption( "convergence" );
  if( convergenceOption && convergenceOption->GetNumberOfFunctions() )
    {
    if( convergenceOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
      {
      numberOfIterations = parser->Convert<unsigned int>( convergenceOption->GetFunction( 0 )->GetName() );
      }
    else
      {
      numberOfIterations = parser->Convert<unsigned int>( convergenceOption->GetFunction( 0 )->GetParameter( 0 ) );
      if( convergenceOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        convergenceThreshold = parser->Convert<RealType>( convergenceOption->GetFunction( 0 )->GetParameter( 1 ) );
        }
      if( convergenceOption->GetFunction( 0 )->GetNumberOfParameters() > 2 )
        {
        convergenceWindowSize = parser->Convert<unsigned int>( convergenceOption->GetFunction( 0 )->GetParameter( 2 ) );
        }
      }
    }


  /////////////////////////////////////////////////////////////////
  //
  //         Get the transform
  //
  /////////////////////////////////////////////////////////////////

  std::string transform = "";
  // std::string outputTransformTypeName = "";
  RealType learningRate = 0.1;
  RealType searchFactor = 20.0 * itk::Math::pi / 180.0;
  RealType arcFraction = 1.0;

  itk::ants::CommandLineParser::OptionType::Pointer searchFactorOption = parser->GetOption( "search-factor" );
  if( searchFactorOption && searchFactorOption->GetNumberOfFunctions() )
    {
    if( searchFactorOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
      {
      searchFactor = parser->Convert<RealType>( searchFactorOption->GetFunction( 0 )->GetName() ) * itk::Math::pi / 180.0;
      }
    if( searchFactorOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      searchFactor = parser->Convert<RealType>( searchFactorOption->GetFunction( 0 )->GetParameter( 0 ) ) * itk::Math::pi / 180.0;
      }
    if( searchFactorOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
      {
      arcFraction = parser->Convert<RealType>( searchFactorOption->GetFunction( 0 )->GetParameter( 1 ) );
      }
    }

  std::vector<RealType> translationSearchGrid;

  RealType translationSearchStepSize = 25.0;

  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    translationSearchGrid.push_back(0.0);
    }

  itk::ants::CommandLineParser::OptionType::Pointer translationSearchGridOption = parser->GetOption( "translation-search-grid" );
  if( translationSearchGridOption && translationSearchGridOption->GetNumberOfFunctions() )
    {

    translationSearchStepSize = parser->Convert<RealType>( translationSearchGridOption->GetFunction( 0 )->GetParameter( 0 ) );

    if( translationSearchGridOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
      {
      translationSearchGrid = parser->ConvertVector<RealType>( translationSearchGridOption->GetFunction( 0 )->GetParameter( 1 ) );
      }
    }

  itk::ants::CommandLineParser::OptionType::Pointer transformOption = parser->GetOption( "transform" );
  if( transformOption && transformOption->GetNumberOfFunctions() )
    {
    transform = transformOption->GetFunction( 0 )->GetName();
    ConvertToLowerCase( transform );
    if( transformOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
      {
      learningRate = parser->Convert<RealType>( transformOption->GetFunction( 0 )->GetParameter( 0 ) );
      }
    }

  typename AffineTransformType::Pointer affineSearchTransform = AffineTransformType::New();
  typename RigidTransformType::Pointer rigidSearchTransform = RigidTransformType::New();
  typename SimilarityTransformType::Pointer similaritySearchTransform = SimilarityTransformType::New();

  typename AffineTransformType::Pointer initialTransform = AffineTransformType::New();
  initialTransform->SetIdentity();

  bool initialTransformInitializedWithImages = false;

  unsigned int numberOfTransformParameters = 0;
  if( strcmp( transform.c_str(), "affine" ) == 0 )
    {
    numberOfTransformParameters = AffineTransformType::ParametersDimension;
    // outputTransformTypeName = std::string( "Affine" );
    }
  else if( strcmp( transform.c_str(), "rigid" ) == 0 )
    {
    numberOfTransformParameters = RigidTransformType::ParametersDimension;
    // outputTransformTypeName = std::string( "Rigid" );
    }
  else if( strcmp( transform.c_str(), "similarity" ) == 0 )
    {
    numberOfTransformParameters = SimilarityTransformType::ParametersDimension;
    // outputTransformTypeName = std::string( "Similarity" );
    }
  else if( strcmp( transform.c_str(), "aligngeometriccenters" ) == 0 ||
    strcmp( transform.c_str(), "aligncentersofmass" ) == 0 )
    {
    typedef itk::CenteredTransformInitializer<AffineTransformType, ImageType, ImageType> TransformInitializerType;
    typename TransformInitializerType::Pointer initializer = TransformInitializerType::New();

    initializer->SetTransform( initialTransform );

    initializer->SetFixedImage( fixedImage );
    initializer->SetMovingImage( movingImage );

    if( strcmp( transform.c_str(), "aligngeometriccenters" ) == 0 )
      {
      initializer->GeometryOn();
      }
    else // strcmp( transform.c_str(), "aligncentersofmass" ) == 0
      {
      initializer->MomentsOn();
      }
    initializer->InitializeTransform();

    initialTransformInitializedWithImages = true;

    // outputTransformTypeName = std::string( "Translation" );
    }
  else
    {
    if( verbose )
      {
      std::cerr << "Unrecognized transform option." << std::endl;
      }
    return EXIT_FAILURE;
    }

  /////////////////////////////////////////////////////////////////
  //
  //         Align the images using center of mass and principal axes
  //         or feature blobs if the initial transform hasn't already
  //         been initialized.
  //
  /////////////////////////////////////////////////////////////////

  itk::Vector<RealType, ImageDimension> axis1( 0.0 );
  itk::Vector<RealType, ImageDimension> axis2( 0.0 );
  axis1[0] = 1.0;
  axis2[1] = 1.0;

  RealType bestScale = 1.0;

  if( initialTransformInitializedWithImages == false )
    {
    constexpr unsigned int minimumNumberOfBlobs = 3;  // should a different min number of blobs be expected?

    unsigned int numberOfBlobsToExtract = 0;
    unsigned int numberOfBlobsToMatch = 0;
    itk::ants::CommandLineParser::OptionType::Pointer blobsOption = parser->GetOption( "align-blobs" );
    if( blobsOption && blobsOption->GetNumberOfFunctions() )
      {
      if( blobsOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
        {
        numberOfBlobsToExtract = parser->Convert<unsigned int>( blobsOption->GetFunction( 0 )->GetName() );
        numberOfBlobsToMatch = numberOfBlobsToExtract;
        }
      if( blobsOption->GetFunction( 0 )->GetNumberOfParameters() > 0 )
        {
        numberOfBlobsToExtract = parser->Convert<unsigned int>( blobsOption->GetFunction( 0 )->GetParameter( 0 ) );
        numberOfBlobsToMatch = numberOfBlobsToExtract;
        }
      if( blobsOption->GetFunction( 0 )->GetNumberOfParameters() > 1 )
        {
        numberOfBlobsToMatch = parser->Convert<unsigned int>( blobsOption->GetFunction( 0 )->GetParameter( 1 ) );
        }
      if( numberOfBlobsToExtract < minimumNumberOfBlobs )
        {
        std::cerr << "Please specify a greater number of blobs (>=" << minimumNumberOfBlobs << ")." << std::endl;
        return EXIT_FAILURE;
        }
      }

    if( numberOfBlobsToExtract >= minimumNumberOfBlobs )
      {
      if( strcmp( transform.c_str(), "affine" ) == 0 )
        {
        typename AffineTransformType::Pointer initialAffineTransform =
          GetTransformFromFeatureMatching<ImageType, AffineTransformType>( fixedImage, movingImage, numberOfBlobsToExtract, numberOfBlobsToMatch );

        initialTransform->SetOffset( initialAffineTransform->GetOffset() );
        initialTransform->SetMatrix( initialAffineTransform->GetMatrix() );
        }
      else  // RigidTransform or SimilarityTransform
        {

        typename LandmarkRigidTransformType::Pointer initialRigidTransform =
          GetTransformFromFeatureMatching<ImageType, LandmarkRigidTransformType>( fixedImage, movingImage, numberOfBlobsToExtract, numberOfBlobsToMatch );

        initialTransform->SetOffset( initialRigidTransform->GetOffset() );
        initialTransform->SetMatrix( initialRigidTransform->GetMatrix() );
        }
      }
    else
      {
      bool doAlignPrincipalAxes = false;

      itk::ants::CommandLineParser::OptionType::Pointer axesOption = parser->GetOption( "align-principal-axes" );
      if( axesOption && axesOption->GetNumberOfFunctions() )
        {
        doAlignPrincipalAxes = parser->Convert<bool>( axesOption->GetFunction( 0 )->GetName() );
        }

      typedef typename itk::ImageMomentsCalculator<ImageType> ImageMomentsCalculatorType;
      typedef typename ImageMomentsCalculatorType::MatrixType MatrixType;

      typename ImageMomentsCalculatorType::Pointer fixedImageMomentsCalculator = ImageMomentsCalculatorType::New();
      typename ImageMomentsCalculatorType::Pointer movingImageMomentsCalculator = ImageMomentsCalculatorType::New();

      fixedImageMomentsCalculator->SetImage( fixedImage );
      fixedImageMomentsCalculator->Compute();
      VectorType fixedImageCenterOfGravity = fixedImageMomentsCalculator->GetCenterOfGravity();
      MatrixType fixedImagePrincipalAxes = fixedImageMomentsCalculator->GetPrincipalAxes();

      movingImageMomentsCalculator->SetImage( movingImage );
      movingImageMomentsCalculator->Compute();
      VectorType movingImageCenterOfGravity = movingImageMomentsCalculator->GetCenterOfGravity();
      MatrixType movingImagePrincipalAxes = movingImageMomentsCalculator->GetPrincipalAxes();

      // RealType bestScale = 1.0; // movingImageMomentsCalculator->GetTotalMass() / fixedImageMomentsCalculator->GetTotalMass();

      typename AffineTransformType::OutputVectorType translation;
      itk::Point<RealType, ImageDimension> center;
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        translation[i] = movingImageCenterOfGravity[i] - fixedImageCenterOfGravity[i];
        center[i] = fixedImageCenterOfGravity[i];
        }
      initialTransform->SetTranslation( translation );

      /** Solve Wahba's problem --- http://en.wikipedia.org/wiki/Wahba%27s_problem */

      vnl_vector<RealType> fixedPrimaryEigenVector;
      vnl_vector<RealType> fixedSecondaryEigenVector;
      vnl_vector<RealType> fixedTertiaryEigenVector;
      vnl_vector<RealType> movingPrimaryEigenVector;
      vnl_vector<RealType> movingSecondaryEigenVector;

      vnl_matrix<RealType> B;

      if( ImageDimension == 2 )
        {
        fixedPrimaryEigenVector = fixedImagePrincipalAxes.GetVnlMatrix().get_row( 1 ).as_ref();
        movingPrimaryEigenVector = movingImagePrincipalAxes.GetVnlMatrix().get_row( 1 ).as_ref();

        B = outer_product( movingPrimaryEigenVector, fixedPrimaryEigenVector );
        }
      else if( ImageDimension == 3 )
        {
        fixedPrimaryEigenVector = fixedImagePrincipalAxes.GetVnlMatrix().get_row( 2 ).as_ref();
        fixedSecondaryEigenVector = fixedImagePrincipalAxes.GetVnlMatrix().get_row( 1 ).as_ref();

        movingPrimaryEigenVector = movingImagePrincipalAxes.GetVnlMatrix().get_row( 2 ).as_ref();
        movingSecondaryEigenVector = movingImagePrincipalAxes.GetVnlMatrix().get_row( 1 ).as_ref();

        B = outer_product( movingPrimaryEigenVector, fixedPrimaryEigenVector ) +
          outer_product( movingSecondaryEigenVector, fixedSecondaryEigenVector );
        }

      if( doAlignPrincipalAxes )
        {
        vnl_svd<RealType> wahba( B );
        vnl_matrix<RealType> A = wahba.V() * wahba.U().transpose();
        A = vnl_inverse( A );
        RealType det = vnl_determinant( A );

        if( det < 0.0 )
          {
          if( verbose )
            {
            std::cout << "Bad determinant = " << det << std::endl;
            std::cout <<  "  det( V ) = " <<  vnl_determinant( wahba.V() ) << std::endl;
            std::cout <<  "  det( U ) = " << vnl_determinant( wahba.U() )  << std::endl;
            }
          vnl_matrix<RealType> I( A );
          I.set_identity();
          for( unsigned int i = 0; i < ImageDimension; i++ )
            {
            if( A( i, i ) < 0.0 )
              {
              I( i, i ) = -1.0;
              }
            }
          A = A * I.transpose();
          det = vnl_determinant( A );

          if( verbose )
            {
            std::cout << "New determinant = " << det << std::endl;
            }
          }
        initialTransform->SetMatrix( A );
        }
      initialTransform->SetCenter( center );

      if( ImageDimension == 2 )
        {
        fixedTertiaryEigenVector = fixedSecondaryEigenVector;
        fixedSecondaryEigenVector = fixedPrimaryEigenVector;
        }
      if( ImageDimension == 3 )
        {
        fixedTertiaryEigenVector = vnl_cross_3d( fixedPrimaryEigenVector, fixedSecondaryEigenVector );
        }

      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        axis1[d] = fixedTertiaryEigenVector[d];
        axis2[d] = fixedSecondaryEigenVector[d];
        }
      }
    }
  /////////////////////////////////////////////////////////////////
  //
  //         Write the output if the number of iterations == 0
  //
  /////////////////////////////////////////////////////////////////

  if( numberOfIterations == 0 )
    {
    itk::ants::CommandLineParser::OptionType::Pointer outputOption = parser->GetOption( "output" );
    if( outputOption && outputOption->GetNumberOfFunctions() )
      {
      std::string outputName = std::string( "" );
      if( outputOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
        {
        outputName = outputOption->GetFunction( 0 )->GetName();
        }
      else
        {
        outputName = outputOption->GetFunction( 0 )->GetParameter( 0 );
        }
      typedef itk::TransformFileWriter TransformWriterType;
      typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();

      if( strcmp( transform.c_str(), "affine" ) == 0 )
        {
        typename AffineTransformType::Pointer bestAffineTransform = AffineTransformType::New();
        bestAffineTransform->SetCenter( initialTransform->GetCenter() );
        bestAffineTransform->SetMatrix( initialTransform->GetMatrix() );
        bestAffineTransform->SetOffset( initialTransform->GetOffset() );
        transformWriter->SetInput( bestAffineTransform );
        }
      else if( strcmp( transform.c_str(), "rigid" ) == 0 )
        {
        typename RigidTransformType::Pointer bestRigidTransform = RigidTransformType::New();
        bestRigidTransform->SetCenter( initialTransform->GetCenter() );
        bestRigidTransform->SetMatrix( initialTransform->GetMatrix() );
        bestRigidTransform->SetOffset( initialTransform->GetOffset() );
        transformWriter->SetInput( bestRigidTransform );
        }
      else if( strcmp( transform.c_str(), "similarity" ) == 0 )
        {
        typename SimilarityTransformType::Pointer bestSimilarityTransform = SimilarityTransformType::New();
        bestSimilarityTransform->SetCenter( initialTransform->GetCenter() );
        bestSimilarityTransform->SetMatrix( initialTransform->GetMatrix() );
        bestSimilarityTransform->SetOffset( initialTransform->GetOffset() );
        transformWriter->SetInput( bestSimilarityTransform );
        }
      else if( initialTransformInitializedWithImages == true )
        {
        typedef itk::TranslationTransform<RealType, ImageDimension> TranslationTransformType;
        typename TranslationTransformType::Pointer bestTranslationTransform = TranslationTransformType::New();
        bestTranslationTransform->SetOffset( initialTransform->GetOffset() );
        transformWriter->SetInput( bestTranslationTransform );
        }

      transformWriter->SetFileName( outputName.c_str() );
#if ITK_VERSION_MAJOR >= 5
      transformWriter->SetUseCompression( true );
#endif
      transformWriter->Update();
      }

    return EXIT_SUCCESS;
    }

  /////////////////////////////////////////////////////////////////
  //
  //         Read in the masks
  //
  /////////////////////////////////////////////////////////////////

  typedef itk::ImageMaskSpatialObject<ImageDimension> ImageMaskSpatialObjectType;
  typedef typename ImageMaskSpatialObjectType::ImageType MaskImageType;

  typename MaskImageType::Pointer fixedMask = nullptr;
  typename MaskImageType::Pointer movingMask = nullptr;

  itk::ants::CommandLineParser::OptionType::Pointer maskOption = parser->GetOption( "masks" );
  if( maskOption && maskOption->GetNumberOfFunctions() )
    {
    if( maskOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
      {
      ReadImage<MaskImageType>( fixedMask, ( maskOption->GetFunction( 0 )->GetName() ).c_str() );
      }
    else if( maskOption->GetFunction( 0 )->GetNumberOfParameters() == 1 )
      {
      ReadImage<MaskImageType>( fixedMask, ( maskOption->GetFunction( 0 )->GetParameter( 0 ) ).c_str() );
      }
    else if( maskOption->GetFunction( 0 )->GetNumberOfParameters() == 2 )
      {
      ReadImage<MaskImageType>( fixedMask, ( maskOption->GetFunction( 0 )->GetParameter( 0 ) ).c_str() );
      ReadImage<MaskImageType>( movingMask, ( maskOption->GetFunction( 0 )->GetParameter( 1 ) ).c_str() );
      }
    }

  typename ImageMaskSpatialObjectType::Pointer fixedMaskSpatialObject = nullptr;
  if( fixedMask.IsNotNull() )
    {
    fixedMaskSpatialObject = ImageMaskSpatialObjectType::New();
    fixedMaskSpatialObject->SetImage( const_cast<MaskImageType *>( fixedMask.GetPointer() ) );
    }

  typename ImageMaskSpatialObjectType::Pointer movingMaskSpatialObject = nullptr;
  if( movingMask.IsNotNull() )
    {
    movingMaskSpatialObject = ImageMaskSpatialObjectType::New();
    movingMaskSpatialObject->SetImage( const_cast<MaskImageType *>( movingMask.GetPointer() ) );
    }

  /////////////////////////////////////////////////////////////////
  //
  //         Set up the image metric
  //
  /////////////////////////////////////////////////////////////////

  typedef itk::ImageToImageMetricv4<ImageType, ImageType, ImageType, RealType> ImageMetricType;
  typename ImageMetricType::Pointer imageMetric = nullptr;

  if( std::strcmp( metric.c_str(), "mattes" ) == 0 )
    {
    if( verbose )
      {
      std::cout << "Using the Mattes MI metric (number of bins = " << numberOfBins << ")" << std::endl;
      }
    typedef itk::MattesMutualInformationImageToImageMetricv4<ImageType, ImageType, ImageType, RealType> MutualInformationMetricType;
    typename MutualInformationMetricType::Pointer mutualInformationMetric = MutualInformationMetricType::New();
    //mutualInformationMetric = mutualInformationMetric;
    mutualInformationMetric->SetNumberOfHistogramBins( numberOfBins );
    mutualInformationMetric->SetUseMovingImageGradientFilter( true );
    mutualInformationMetric->SetUseFixedImageGradientFilter( true );

    imageMetric = mutualInformationMetric;
    }
  else if( std::strcmp( metric.c_str(), "mi" ) == 0 )
    {
    if( verbose )
      {
      std::cout << "Using the joint histogram MI metric (number of bins = " << numberOfBins << ")" << std::endl;
      }
    typedef itk::JointHistogramMutualInformationImageToImageMetricv4<ImageType, ImageType, ImageType,
                                                                     RealType> MutualInformationMetricType;
    typename MutualInformationMetricType::Pointer mutualInformationMetric = MutualInformationMetricType::New();
    //mutualInformationMetric = mutualInformationMetric;
    mutualInformationMetric->SetNumberOfHistogramBins( numberOfBins );
    mutualInformationMetric->SetUseMovingImageGradientFilter( true );
    mutualInformationMetric->SetUseFixedImageGradientFilter( true );
    mutualInformationMetric->SetVarianceForJointPDFSmoothing( 1.0 );

    imageMetric = mutualInformationMetric;
    }
  else if( std::strcmp( metric.c_str(), "gc" ) == 0 )
    {
    if( verbose )
      {
      std::cout << "Using the global correlation metric " << std::endl;
      }
    typedef itk::CorrelationImageToImageMetricv4<ImageType, ImageType, ImageType, RealType> corrMetricType;
    typename corrMetricType::Pointer corrMetric = corrMetricType::New();

    imageMetric = corrMetric;
    }
  else
    {
    if( verbose )
      {
      std::cerr << "ERROR: Unrecognized metric. " << std::endl;
      }
    return EXIT_FAILURE;
    }

  imageMetric->SetFixedImage( fixedImage );
  imageMetric->SetVirtualDomainFromImage( fixedImage );
  imageMetric->SetMovingImage( movingImage );
  imageMetric->SetFixedImageMask( fixedMaskSpatialObject );
  imageMetric->SetMovingImageMask( movingMaskSpatialObject );
  imageMetric->SetUseSampledPointSet( false );

  /** Sample the image domain **/

  if( samplingStrategy != NONE )
    {
    const typename ImageType::SpacingType oneThirdVirtualSpacing = fixedImage->GetSpacing() / 3.0;

    typedef typename ImageMetricType::FixedSampledPointSetType MetricSamplePointSetType;
    typename MetricSamplePointSetType::Pointer samplePointSet = MetricSamplePointSetType::New();
    samplePointSet->Initialize();

    typedef typename MetricSamplePointSetType::PointType SamplePointType;

    typedef typename itk::Statistics::MersenneTwisterRandomVariateGenerator RandomizerType;
    typename RandomizerType::Pointer randomizer = RandomizerType::New();

    int antsRandomSeed = 1234;

    itk::ants::CommandLineParser::OptionType::Pointer randomSeedOption = parser->GetOption( "random-seed" );
    if( randomSeedOption && randomSeedOption->GetNumberOfFunctions() )
      {
      antsRandomSeed = parser->Convert<int>( randomSeedOption->GetFunction(0)->GetName() );
      }
    else
      {
      char* envSeed = getenv( "ANTS_RANDOM_SEED" );

      if ( envSeed != nullptr )
	{
	antsRandomSeed = std::stoi( envSeed );
	}
      }

    if ( antsRandomSeed != 0 )
      {
      randomizer->SetSeed( antsRandomSeed );
      }

    unsigned long index = 0;

    switch( samplingStrategy )
      {
      case REGULAR:
        {
        const auto sampleCount = static_cast<unsigned long>( std::ceil( 1.0 / samplingPercentage ) );
        unsigned long count = sampleCount; //Start at sampleCount to keep behavior backwards identical, using first element.
        itk::ImageRegionConstIteratorWithIndex<ImageType> It( fixedImage, fixedImage->GetRequestedRegion() );
        for( It.GoToBegin(); !It.IsAtEnd(); ++It )
          {
          if( count == sampleCount )
            {
            count = 0; //Reset counter
            SamplePointType point;
            fixedImage->TransformIndexToPhysicalPoint( It.GetIndex(), point );

            // randomly perturb the point within a voxel (approximately)
            for( unsigned int d = 0; d < ImageDimension; d++ )
              {
              point[d] += static_cast<typename SamplePointType::CoordRepType>( randomizer->GetNormalVariate() ) * static_cast<typename SamplePointType::CoordRepType>( oneThirdVirtualSpacing[d] );
              }
            if( !fixedMaskSpatialObject || fixedMaskSpatialObject->IsInsideInWorldSpace( point ) )
              {
              samplePointSet->SetPoint( index, point );
              ++index;
              }
            }
          ++count;
          }
        break;
        }
      case RANDOM:
        {
        const unsigned long totalVirtualDomainVoxels = fixedImage->GetRequestedRegion().GetNumberOfPixels();
        const auto sampleCount = static_cast<unsigned long>( static_cast<float>( totalVirtualDomainVoxels ) * static_cast<float>( samplingPercentage ) );
        itk::ImageRandomConstIteratorWithIndex<ImageType> ItR( fixedImage, fixedImage->GetRequestedRegion() );
        ItR.SetNumberOfSamples( sampleCount );
        for( ItR.GoToBegin(); !ItR.IsAtEnd(); ++ItR )
          {
          SamplePointType point;
          fixedImage->TransformIndexToPhysicalPoint( ItR.GetIndex(), point );

          // randomly perturb the point within a voxel (approximately)
          for ( unsigned int d = 0; d < ImageDimension; d++ )
            {
            point[d] += static_cast<typename SamplePointType::CoordRepType>( randomizer->GetNormalVariate() * oneThirdVirtualSpacing[d] );
            }
          if( !fixedMaskSpatialObject || fixedMaskSpatialObject->IsInsideInWorldSpace( point ) )
            {
            samplePointSet->SetPoint( index, point );
            ++index;
            }
          }
        break;
        }
      case NONE:
        break;
      }
    imageMetric->SetFixedSampledPointSet( samplePointSet );
    imageMetric->SetUseSampledPointSet( true );
    }

  imageMetric->Initialize();

  if( strcmp( transform.c_str(), "affine" ) == 0 )
    {
    imageMetric->SetMovingTransform( affineSearchTransform );
    }
  else if( strcmp( transform.c_str(), "rigid" ) == 0 )
    {
    imageMetric->SetMovingTransform( rigidSearchTransform );
    }
  else if( strcmp( transform.c_str(), "similarity" ) == 0 )
    {
    imageMetric->SetMovingTransform( similaritySearchTransform );
    }


  /////////////////////////////////////////////////////////////////
  //
  //         Set up the optimizer
  //
  /////////////////////////////////////////////////////////////////

  typedef itk::RegistrationParameterScalesFromPhysicalShift<ImageMetricType> RegistrationParameterScalesFromPhysicalShiftType;
  typename RegistrationParameterScalesFromPhysicalShiftType::Pointer scalesEstimator = RegistrationParameterScalesFromPhysicalShiftType::New();
  scalesEstimator->SetMetric( imageMetric );
  scalesEstimator->SetTransformForward( true );

  typename RegistrationParameterScalesFromPhysicalShiftType::ScalesType movingScales( numberOfTransformParameters );
  scalesEstimator->EstimateScales( movingScales );

  typedef  itk::ConjugateGradientLineSearchOptimizerv4 LocalOptimizerType;
  typename LocalOptimizerType::Pointer localOptimizer = LocalOptimizerType::New();
  localOptimizer->SetLowerLimit( 0 );
  localOptimizer->SetUpperLimit( 2 );
  localOptimizer->SetEpsilon( 0.1 );
  localOptimizer->SetMaximumLineSearchIterations( 10 );
  localOptimizer->SetLearningRate( learningRate );
  localOptimizer->SetMaximumStepSizeInPhysicalUnits( learningRate );
  localOptimizer->SetNumberOfIterations( numberOfIterations );
  localOptimizer->SetMinimumConvergenceValue( convergenceThreshold );
  localOptimizer->SetConvergenceWindowSize( convergenceWindowSize );
  localOptimizer->SetDoEstimateLearningRateOnce( true );
  localOptimizer->SetScales( movingScales );
  localOptimizer->SetMetric( imageMetric );

  typedef itk::MultiStartOptimizerv4 MultiStartOptimizerType;
  typename MultiStartOptimizerType::Pointer multiStartOptimizer = MultiStartOptimizerType::New();
  multiStartOptimizer->SetScales( movingScales );
  multiStartOptimizer->SetMetric( imageMetric );

  unsigned int trialCounter = 0;

  typename MultiStartOptimizerType::ParametersListType parametersList = multiStartOptimizer->GetParametersList();
  for( RealType angle1 = ( itk::Math::pi * -arcFraction ); angle1 <= ( itk::Math::pi * arcFraction + 0.000001 ); angle1 += searchFactor )
    {
    if( ImageDimension == 2 )
      {
      for ( RealType translation1 = -1.0 * translationSearchGrid[0];
           translation1 <= translationSearchGrid[0] + 0.000001; translation1 += translationSearchStepSize )
        {
        for ( RealType translation2 = -1.0 * translationSearchGrid[1];
             translation2 <= translationSearchGrid[1] + 0.000001; translation2 += translationSearchStepSize )
          {
          typename AffineTransformType::OutputVectorType searchTranslation;
          searchTranslation[0] = translation1;
          searchTranslation[1] = translation2;

          affineSearchTransform->SetIdentity();
          affineSearchTransform->SetCenter( initialTransform->GetCenter() );
          affineSearchTransform->SetMatrix( initialTransform->GetMatrix() );
          affineSearchTransform->SetOffset( initialTransform->GetOffset() );
          affineSearchTransform->Translate( searchTranslation , 1 );
          affineSearchTransform->Rotate2D( angle1, 1 );

          if( strcmp( transform.c_str(), "affine" ) == 0 )
            {
            affineSearchTransform->Scale( bestScale );
            parametersList.push_back( affineSearchTransform->GetParameters() );
            }
          else if( strcmp( transform.c_str(), "rigid" ) == 0 )
            {
            rigidSearchTransform->SetIdentity();
            rigidSearchTransform->SetCenter( initialTransform->GetCenter() );
            rigidSearchTransform->SetMatrix( affineSearchTransform->GetMatrix() );
            rigidSearchTransform->SetOffset( initialTransform->GetOffset() );

            parametersList.push_back( rigidSearchTransform->GetParameters() );
            }
          else if( strcmp( transform.c_str(), "similarity" ) == 0 )
            {
            similaritySearchTransform->SetIdentity();
            similaritySearchTransform->SetCenter( initialTransform->GetCenter() );
            similaritySearchTransform->SetMatrix( affineSearchTransform->GetMatrix() );
            similaritySearchTransform->SetOffset( initialTransform->GetOffset() );
            similaritySearchTransform->SetScale( bestScale );

            similaritySearchTransform->SetScale( bestScale );

            parametersList.push_back( similaritySearchTransform->GetParameters() );
            }
          trialCounter++;
          }
        }
      }
    if( ImageDimension == 3 )
      {
      for( RealType angle2 = ( itk::Math::pi * -arcFraction ); angle2 <= ( itk::Math::pi * arcFraction + 0.000001 ); angle2 += searchFactor )
        {
        for( RealType angle3 = ( itk::Math::pi * -arcFraction ); angle3 <= ( itk::Math::pi * arcFraction + 0.000001 ); angle3 += searchFactor )
          {
          for ( RealType translation1 = -1.0 * translationSearchGrid[0];
                translation1 <= translationSearchGrid[0] + 0.000001; translation1 += translationSearchStepSize )
            {
            for ( RealType translation2 = -1.0 * translationSearchGrid[1];
                  translation2 <= translationSearchGrid[1] + 0.000001; translation2 += translationSearchStepSize )
              {
              for ( RealType translation3 = -1.0 * translationSearchGrid[2];
                    translation3 <= translationSearchGrid[2] + 0.000001; translation3 += translationSearchStepSize )
                {
                typename AffineTransformType::OutputVectorType searchTranslation;
                searchTranslation[0] = translation1;
                searchTranslation[1] = translation2;
                searchTranslation[2] = translation3;

                affineSearchTransform->SetIdentity();
                affineSearchTransform->SetCenter( initialTransform->GetCenter() );
                affineSearchTransform->SetOffset( initialTransform->GetOffset() );
                affineSearchTransform->SetMatrix( initialTransform->GetMatrix() );
                affineSearchTransform->Translate( searchTranslation, 0 );
                affineSearchTransform->Rotate3D( axis1, angle1, 1 );
                affineSearchTransform->Rotate3D( axis2, angle2, 1 );
                affineSearchTransform->Rotate3D( axis1, angle3, 1 );

                if( strcmp( transform.c_str(), "affine" ) == 0 )
                  {
                  affineSearchTransform->Scale( bestScale );
                  parametersList.push_back( affineSearchTransform->GetParameters() );
		  }
                else if( strcmp( transform.c_str(), "rigid" ) == 0 )
                  {
                  rigidSearchTransform->SetIdentity();
                  rigidSearchTransform->SetCenter( initialTransform->GetCenter() );
                  rigidSearchTransform->SetOffset( initialTransform->GetOffset() );
                  rigidSearchTransform->Translate( searchTranslation, 0 );
                  rigidSearchTransform->SetMatrix( affineSearchTransform->GetMatrix() );
                  parametersList.push_back( rigidSearchTransform->GetParameters() );
		  }
                else if( strcmp( transform.c_str(), "similarity" ) == 0 )
                  {
                  similaritySearchTransform->SetIdentity();
                  similaritySearchTransform->SetCenter( initialTransform->GetCenter() );
                  similaritySearchTransform->SetOffset( initialTransform->GetOffset() );
                  similaritySearchTransform->SetMatrix( affineSearchTransform->GetMatrix() );
                  similaritySearchTransform->SetScale( bestScale );

                  parametersList.push_back( similaritySearchTransform->GetParameters() );
                  }
		trialCounter++;
                }
              }
            }
          }
        }
      }
    }

  if( verbose )
    {
    std::cout << "Starting optimizer with " << trialCounter << " starting points" << std::endl;
    }

  multiStartOptimizer->SetParametersList( parametersList );
  multiStartOptimizer->SetLocalOptimizer( localOptimizer );
  multiStartOptimizer->StartOptimization();


  /////////////////////////////////////////////////////////////////
  //
  //         Write the output after convergence
  //
  /////////////////////////////////////////////////////////////////

  itk::ants::CommandLineParser::OptionType::Pointer outputOption = parser->GetOption( "output" );
  if( outputOption && outputOption->GetNumberOfFunctions() )
    {
    std::string outputName = std::string( "" );
    if( outputOption->GetFunction( 0 )->GetNumberOfParameters() == 0 )
      {
      outputName = outputOption->GetFunction( 0 )->GetName();
      }
    else
      {
      outputName = outputOption->GetFunction( 0 )->GetParameter( 0 );
      }

    typedef itk::TransformFileWriter TransformWriterType;
    typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();

    if( strcmp( transform.c_str(), "affine" ) == 0 )
      {
      typename AffineTransformType::Pointer bestAffineTransform = AffineTransformType::New();
      bestAffineTransform->SetCenter( initialTransform->GetCenter() );
      bestAffineTransform->SetParameters( multiStartOptimizer->GetBestParameters() );
      transformWriter->SetInput( bestAffineTransform );
      }
    else if( strcmp( transform.c_str(), "rigid" ) == 0 )
      {
      typename RigidTransformType::Pointer bestRigidTransform = RigidTransformType::New();
      bestRigidTransform->SetCenter( initialTransform->GetCenter() );
      bestRigidTransform->SetParameters( multiStartOptimizer->GetBestParameters() );
      transformWriter->SetInput( bestRigidTransform );
      }
    else if( strcmp( transform.c_str(), "similarity" ) == 0 )
      {
      typename SimilarityTransformType::Pointer bestSimilarityTransform = SimilarityTransformType::New();
      bestSimilarityTransform->SetCenter( initialTransform->GetCenter() );
      bestSimilarityTransform->SetParameters( multiStartOptimizer->GetBestParameters() );
      transformWriter->SetInput( bestSimilarityTransform );
      }

    transformWriter->SetFileName( outputName.c_str() );
    transformWriter->Update();
    }

  return EXIT_SUCCESS;
}

void InitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

  {
  std::string description =
    std::string( "This option forces the image to be treated as a specified-" )
    + std::string( "dimensional image.  If not specified, we try to " )
    + std::string( "infer the dimensionality from the input image." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "dimensionality" );
  option->SetShortName( 'd' );
  option->SetUsageOption( 0, "2/3" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "These image metrics are available:  " )
    + std::string( "MI:  joint histogram and Mattes: mutual information  and  GC:  global correlation." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "metric" );
  option->SetShortName( 'm' );
  option->SetUsageOption(
    0,
    "MI[fixedImage,movingImage,<numberOfBins=32>,<samplingStrategy={None,Regular,Random}>,<samplingPercentage=[0,1]>]" );
  option->SetUsageOption(
    1,
    "Mattes[fixedImage,movingImage,<numberOfBins=32>,<samplingStrategy={None,Regular,Random}>,<samplingPercentage=[0,1]>]" );
  option->SetUsageOption(
    2,
    "GC[fixedImage,movingImage,<radius=NA>,<samplingStrategy={None,Regular,Random}>,<samplingPercentage=[0,1]>]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Several transform options are available.  For the rigid, " )
    + std::string( "affine, and similarity transforms, the gradientStep characterizes the gradient " )
    + std::string( "descent optimization and is scaled appropriately for each transform using the " )
    + std::string( "shift scales estimator. The other two transform types finds the simple translation " )
    + std::string( "transform which aligns the specified image feature.  Note that the images are read " )
    + std::string( "from the similarity metric option although the metric isn't actually used." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "transform" );
  option->SetShortName( 't' );
  option->SetUsageOption(  0, "Rigid[gradientStep]" );
  option->SetUsageOption(  1, "Affine[gradientStep]" );
  option->SetUsageOption(  2, "Similarity[gradientStep]" );
  option->SetUsageOption(  3, "AlignGeometricCenters" );
  option->SetUsageOption(  4, "AlignCentersOfMass" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Boolean indicating alignment by principal axes.  " )
    + std::string( "Alternatively, one can align using blobs (see -b option)." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "align-principal-axes" );
  option->SetShortName( 'p' );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Boolean indicating alignment by a set of blobs.  " )
    + std::string( "Alternatively, one can align using blobs (see -p option)." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "align-blobs" );
  option->SetShortName( 'b' );
  option->SetUsageOption(  0, "numberOfBlobsToExtract" );
  option->SetUsageOption(  1, "[numberOfBlobsToExtract,<numberOfBlobsToMatch=numberOfBlobsToExtract>]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Incremental search factor (in degrees) which will " )
   + std::string( "sample the arc fraction around the principal axis or default axis." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "search-factor" );
  option->SetShortName( 's' );
  option->SetUsageOption( 0, "searchFactor" );
  option->SetUsageOption( 1, "[searchFactor=20,<arcFraction=1.0>]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Translation search grid in mm, which will " )
    + std::string( "translate the moving image in each dimension in increments of the step size.");

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "translation-search-grid" );
  option->SetShortName( 'g' );
  option->SetUsageOption( 0, "[stepSize=25, AxBxC=0x0x0]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Number of iterations." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "convergence" );
  option->SetShortName( 'c' );
  option->SetUsageOption( 0, "numberOfIterations" );
  option->SetUsageOption( 1, "[numberOfIterations,<convergenceThreshold=1e-6>,<convergenceWindowSize=10>]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Image masks to limit voxels considered by the metric." );
  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "masks" );
  option->SetShortName( 'x' );
  option->SetUsageOption( 0, "fixedImageMask" );
  option->SetUsageOption( 1, "[fixedImageMask,movingImageMask]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Specify the output transform (output format an ITK .mat file). " );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "output" );
  option->SetShortName( 'o' );
  option->SetUsageOption( 0, "outputFileName" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Use a fixed seed for random number generation. " )
    + std::string( "The default fixed seed is overwritten by this value. " )
    + std::string( "The fixed seed can be any nonzero int value. If the specified seed is zero, " )
    + std::string( "the system time will be used." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "random-seed" );
  option->SetUsageOption( 0, "seedValue" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Verbose output." );

  OptionType::Pointer option = OptionType::New();
  option->SetShortName( 'v' );
  option->SetLongName( "verbose" );
  option->SetUsageOption( 0, "(0)/1" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Print the help menu (short version)." );

  OptionType::Pointer option = OptionType::New();
  option->SetShortName( 'h' );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Print the help menu." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "help" );
  option->SetDescription( description );
  parser->AddOption( option );
  }
}

int antsAI( std::vector<std::string> args, std::ostream* /*out_stream = nullptr */ )
{

  // put the arguments coming in as 'args' into standard (argc,argv) format;
  // 'args' doesn't have the command name as first, argument, so add it manually;
  // 'args' may have adjacent arguments concatenated into one argument,
  // which the parser should handle

  args.insert( args.begin(), "antsAI" );

  int    argc = args.size();
  char** argv = new char *[args.size() + 1];
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
      Cleanup_argv( char* * argv_, int argc_plus_one_ ) : argv( argv_ ), argc_plus_one( argc_plus_one_ ){}

      ~Cleanup_argv()
        {
        for( unsigned int i = 0; i < argc_plus_one; ++i )
          {
          delete[] argv[i];
          }
        delete[] argv;
        }
    private:
    char**       argv;
    unsigned int argc_plus_one;
    };
  Cleanup_argv cleanup_argv( argv, argc + 1 );

  // antscout->set_stream( out_stream );
  typedef itk::ants::CommandLineParser ParserType;

  ParserType::Pointer parser = ParserType::New();
  parser->SetCommand( argv[0] );

  std::string commandDescription = std::string( "Program to calculate the optimal" )
    + std::string( "linear transform parameters for aligning two images." );

  parser->SetCommandDescription( commandDescription );
  InitializeCommandLineOptions( parser );

  if( parser->Parse( argc, argv ) == EXIT_FAILURE )
    {
    return EXIT_FAILURE;
    }

  bool verbose = false;
  itk::ants::CommandLineParser::OptionType::Pointer verboseOption =
    parser->GetOption( "verbose" );
  if( verboseOption && verboseOption->GetNumberOfFunctions() )
    {
    verbose = parser->Convert<bool>( verboseOption->GetFunction( 0 )->GetName() );
    }

  if( argc == 1 )
    {
    parser->PrintMenu( std::cout, 5, false );
    return EXIT_FAILURE;
    }
  else if( parser->GetOption( "help" )->GetFunction() && parser->Convert<bool>( parser->GetOption( "help" )->GetFunction()->GetName() ) )
    {
    parser->PrintMenu( std::cout, 5, false );
    return EXIT_SUCCESS;
    }
  else if( parser->GetOption( 'h' )->GetFunction() && parser->Convert<bool>( parser->GetOption( 'h' )->GetFunction()->GetName() ) )
    {
    parser->PrintMenu( std::cout, 5, true );
    return EXIT_SUCCESS;
    }

  unsigned int dimension = 3;

  ParserType::OptionType::Pointer dimOption = parser->GetOption( "dimensionality" );
  if( dimOption && dimOption->GetNumberOfFunctions() )
    {
    dimension = parser->Convert<unsigned int>( dimOption->GetFunction( 0 )->GetName() );
    }
  else
    {
    if( verbose )
      {
      std::cerr << "Image dimensionality not specified.  See command line option --dimensionality" << std::endl;
      }
    return EXIT_FAILURE;
    }

  switch( dimension )
    {
    case 2:
      {
      return antsAI<2>( parser );
      break;
      }
    case 3:
      {
      return antsAI<3>( parser );
      break;
      }
    default:
      {
      std::cerr << "Unrecognized dimensionality.  Please see help menu." << std::endl;
      return EXIT_FAILURE;
      }
    }
  return EXIT_SUCCESS;
}
} // namespace ants
