#include "antsAllocImage.h"
#include "antsCommandLineParser.h"
#include "antsUtilities.h"
#include "ReadWriteData.h"

#include "itkAffineTransform.h"
#include "itkConjugateGradientLineSearchOptimizerv4.h"
#include "itkEuler2DTransform.h"
#include "itkEuler3DTransform.h"
#include "itkCorrelationImageToImageMetricv4.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImageMomentsCalculator.h"
#include "itkImageToImageMetricv4.h"
#include "itkJointHistogramMutualInformationImageToImageMetricv4.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkMultiStartOptimizerv4.h"
#include "itkRegistrationParameterScalesFromPhysicalShift.h"
#include "itkSimilarity2DTransform.h"
#include "itkSimilarity3DTransform.h"
#include "itkTransformFileWriter.h"

#include "vnl/vnl_cross.h"
#include "vnl/vnl_inverse.h"

namespace ants
{

// ##########################################################################
//      Transform traits to generalize the different linear transforms
// ##########################################################################

template <class TComputeType, unsigned int ImageDimension>
class RigidTransformTraits
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

template <class TComputeType, unsigned int ImageDimension>
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

template <unsigned int ImageDimension>
int antsAI( itk::ants::CommandLineParser *parser )
{
  typedef double RealType;
  typedef float  PixelType;
  typedef float  FloatType;

  typedef itk::Vector<FloatType, ImageDimension>      VectorType;
  typedef itk::Image<PixelType, ImageDimension>       ImageType;

  typedef itk::AffineTransform<RealType, ImageDimension>                              AffineTransformType;
  typedef typename RigidTransformTraits<RealType, ImageDimension>::TransformType      RigidTransformType;
  typedef typename SimilarityTransformTraits<RealType, ImageDimension>::TransformType SimilarityTransformType;

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
  //         Align the images
  //           1. align the image centers of mass
  //           2. if desired, align the principal axes
  //
  /////////////////////////////////////////////////////////////////

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

  RealType bestScale = 1.0; // movingImageMomentsCalculator->GetTotalMass() / fixedImageMomentsCalculator->GetTotalMass();

  /** Align centers of mass **/

  typename AffineTransformType::Pointer initialTransform = AffineTransformType::New();
  initialTransform->SetIdentity();

  typename AffineTransformType::OffsetType offset;
  itk::Point<RealType, ImageDimension> center;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    offset[i] = movingImageCenterOfGravity[i] - fixedImageCenterOfGravity[i];
    center[i] = fixedImageCenterOfGravity[i];
    }
  initialTransform->SetOffset( offset );

  /** Solve Wahba's problem --- http://en.wikipedia.org/wiki/Wahba%27s_problem */

  vnl_vector<RealType> fixedPrimaryEigenVector;
  vnl_vector<RealType> fixedSecondaryEigenVector;
  vnl_vector<RealType> fixedTertiaryEigenVector;
  vnl_vector<RealType> movingPrimaryEigenVector;
  vnl_vector<RealType> movingSecondaryEigenVector;

  vnl_matrix<RealType> B;

  if( ImageDimension == 2 )
    {
    fixedPrimaryEigenVector = fixedImagePrincipalAxes.GetVnlMatrix().get_row( 1 );
    movingPrimaryEigenVector = movingImagePrincipalAxes.GetVnlMatrix().get_row( 1 );

    B = outer_product( movingPrimaryEigenVector, fixedPrimaryEigenVector );
    }
  else if( ImageDimension == 3 )
    {
    fixedPrimaryEigenVector = fixedImagePrincipalAxes.GetVnlMatrix().get_row( 2 );
    fixedSecondaryEigenVector = fixedImagePrincipalAxes.GetVnlMatrix().get_row( 1 );

    movingPrimaryEigenVector = movingImagePrincipalAxes.GetVnlMatrix().get_row( 2 );
    movingSecondaryEigenVector = movingImagePrincipalAxes.GetVnlMatrix().get_row( 1 );

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

  /////////////////////////////////////////////////////////////////
  //
  //         Read in the masks
  //
  /////////////////////////////////////////////////////////////////

  typedef itk::ImageMaskSpatialObject<ImageDimension> ImageMaskSpatialObjectType;
  typedef typename ImageMaskSpatialObjectType::ImageType MaskImageType;

  typename MaskImageType::Pointer fixedMask = ITK_NULLPTR;
  typename MaskImageType::Pointer movingMask = ITK_NULLPTR;

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

  typename ImageMaskSpatialObjectType::Pointer fixedMaskSpatialObject = ITK_NULLPTR;
  if( fixedMask.IsNotNull() )
    {
    fixedMaskSpatialObject = ImageMaskSpatialObjectType::New();
    fixedMaskSpatialObject->SetImage( const_cast<MaskImageType *>( fixedMask.GetPointer() ) );
    }

  typename ImageMaskSpatialObjectType::Pointer movingMaskSpatialObject = ITK_NULLPTR;
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
  typename ImageMetricType::Pointer imageMetric = ITK_NULLPTR;

  if( std::strcmp( metric.c_str(), "mattes" ) == 0 )
    {
    if( verbose )
      {
      std::cout << "Using the Mattes MI metric (number of bins = " << numberOfBins << ")" << std::endl;
      }
    typedef itk::MattesMutualInformationImageToImageMetricv4<ImageType, ImageType, ImageType, RealType> MutualInformationMetricType;
    typename MutualInformationMetricType::Pointer mutualInformationMetric = MutualInformationMetricType::New();
    mutualInformationMetric = mutualInformationMetric;
    mutualInformationMetric->SetNumberOfHistogramBins( numberOfBins );
    mutualInformationMetric->SetUseMovingImageGradientFilter( true );
    mutualInformationMetric->SetUseFixedImageGradientFilter( true );

    imageMetric = mutualInformationMetric;
    }
  else if( std::strcmp( metric.c_str(), "mi" ) == 0 )
    {
    if( verbose )
      {
      std::cout << "Using the Mattes MI metric (number of bins = " << numberOfBins << ")" << std::endl;
      }
    typedef itk::JointHistogramMutualInformationImageToImageMetricv4<ImageType, ImageType, ImageType,
                                                                     RealType> MutualInformationMetricType;
    typename MutualInformationMetricType::Pointer mutualInformationMetric = MutualInformationMetricType::New();
    mutualInformationMetric = mutualInformationMetric;
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
  imageMetric->SetUseFixedSampledPointSet( false );

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
    randomizer->SetSeed( 1234 );

    unsigned long index = 0;

    switch( samplingStrategy )
      {
      case REGULAR:
        {
        const unsigned long sampleCount = static_cast<unsigned long>( std::ceil( 1.0 / samplingPercentage ) );
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
              point[d] += randomizer->GetNormalVariate() * oneThirdVirtualSpacing[d];
              }
            typename ImageType::IndexType imageIndex;
            fixedImage->TransformPhysicalPointToIndex( point, imageIndex );
            if( !fixedMask || fixedMask->GetPixel( imageIndex ) != 0 )
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
        const unsigned long sampleCount = static_cast<unsigned long>( static_cast<float>( totalVirtualDomainVoxels ) * samplingPercentage );
        itk::ImageRandomConstIteratorWithIndex<ImageType> ItR( fixedImage, fixedImage->GetRequestedRegion() );
        ItR.SetNumberOfSamples( sampleCount );
        for( ItR.GoToBegin(); !ItR.IsAtEnd(); ++ItR )
          {
          SamplePointType point;
          fixedImage->TransformIndexToPhysicalPoint( ItR.GetIndex(), point );

          // randomly perturb the point within a voxel (approximately)
          for ( unsigned int d = 0; d < ImageDimension; d++ )
            {
            point[d] += randomizer->GetNormalVariate() * oneThirdVirtualSpacing[d];
            }
          typename ImageType::IndexType imageIndex;
          fixedImage->TransformPhysicalPointToIndex( point, imageIndex );
          if( !fixedMask || fixedMask->GetPixel( imageIndex ) != 0 )
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
    imageMetric->SetUseFixedSampledPointSet( true );
    }

  imageMetric->Initialize();

  /////////////////////////////////////////////////////////////////
  //
  //         Get the transform and set up the optimizers
  //
  /////////////////////////////////////////////////////////////////

  std::string transform = "";
  std::string outputTransformTypeName = "";
  RealType learningRate = 0.1;
  RealType searchFactor = 10.0 * vnl_math::pi / 180.0;

  itk::ants::CommandLineParser::OptionType::Pointer searchFactorOption = parser->GetOption( "search-factor" );
  if( searchFactorOption && searchFactorOption->GetNumberOfFunctions() )
    {
    searchFactor = parser->Convert<RealType>( searchFactorOption->GetFunction( 0 )->GetName() ) * vnl_math::pi / 180.0;
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

  unsigned int numberOfTransformParameters = 0;
  if( strcmp( transform.c_str(), "affine" ) == 0 )
    {
    numberOfTransformParameters = AffineTransformType::ParametersDimension;
    outputTransformTypeName = std::string( "Affine" );

    imageMetric->SetMovingTransform( affineSearchTransform );
    }
  else if( strcmp( transform.c_str(), "rigid" ) == 0 )
    {
    numberOfTransformParameters = RigidTransformType::ParametersDimension;
    outputTransformTypeName = std::string( "Rigid" );

    imageMetric->SetMovingTransform( rigidSearchTransform );
    }
  else if( strcmp( transform.c_str(), "similarity" ) == 0 )
    {
    numberOfTransformParameters = SimilarityTransformType::ParametersDimension;
    outputTransformTypeName = std::string( "Similarity" );

    imageMetric->SetMovingTransform( similaritySearchTransform );
    }
  else
    {
    if( verbose )
      {
      std::cerr << "Unrecognized transform option." << std::endl;
      }
    return EXIT_FAILURE;
    }

  unsigned int numberOfIterations = 20;
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

  if( ImageDimension == 2 )
    {
    fixedTertiaryEigenVector = fixedSecondaryEigenVector;
    fixedSecondaryEigenVector = fixedPrimaryEigenVector;
    }
  if( ImageDimension == 3 )
    {
    fixedTertiaryEigenVector = vnl_cross_3d( fixedPrimaryEigenVector, fixedSecondaryEigenVector );
    }

  itk::Vector<RealType, ImageDimension> axis1;
  itk::Vector<RealType, ImageDimension> axis2;

  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    axis1[d] = fixedTertiaryEigenVector[d];
    axis2[d] = fixedSecondaryEigenVector[d];
    }

  typename MultiStartOptimizerType::ParametersListType parametersList = multiStartOptimizer->GetParametersList();
  for( RealType angle1 = ( vnl_math::pi_over_4 * -1.0 ); angle1 < ( vnl_math::pi_over_4 + searchFactor ); angle1 += searchFactor )
    {
    if( ImageDimension == 2 )
      {
      affineSearchTransform->SetIdentity();
      affineSearchTransform->SetCenter( initialTransform->GetCenter() );
      affineSearchTransform->SetOffset( initialTransform->GetOffset() );
      affineSearchTransform->SetMatrix( initialTransform->GetMatrix() );

      affineSearchTransform->Rotate2D( angle1, 1 );

      if( strcmp( transform.c_str(), "affine" ) == 0 )
        {
        affineSearchTransform->Scale( bestScale );
        parametersList.push_back( affineSearchTransform->GetParameters() );
        }
      else if( strcmp( transform.c_str(), "rigid" ) == 0 )
        {
        rigidSearchTransform->SetIdentity();
        rigidSearchTransform->SetCenter( affineSearchTransform->GetCenter() );
        rigidSearchTransform->SetOffset( affineSearchTransform->GetOffset() );
        rigidSearchTransform->SetMatrix( affineSearchTransform->GetMatrix() );

        parametersList.push_back( rigidSearchTransform->GetParameters() );
        }
      else if( strcmp( transform.c_str(), "similarity" ) == 0 )
        {
        similaritySearchTransform->SetIdentity();
        similaritySearchTransform->SetCenter( affineSearchTransform->GetCenter() );
        similaritySearchTransform->SetOffset( affineSearchTransform->GetOffset() );
        similaritySearchTransform->SetMatrix( affineSearchTransform->GetMatrix() );

        similaritySearchTransform->SetScale( bestScale );

        parametersList.push_back( similaritySearchTransform->GetParameters() );
        }
      }
    if( ImageDimension == 3 )
      {
      for( RealType angle2 = ( vnl_math::pi_over_4 * -1.0 ); angle2 < ( vnl_math::pi_over_4 + searchFactor ); angle2 += searchFactor )
        {
        affineSearchTransform->SetIdentity();
        affineSearchTransform->SetCenter( initialTransform->GetCenter() );
        affineSearchTransform->SetOffset( initialTransform->GetOffset() );
        affineSearchTransform->SetMatrix( initialTransform->GetMatrix() );

        affineSearchTransform->Rotate3D( axis1, angle1, 1 );
        affineSearchTransform->Rotate3D( axis2, angle2, 1 );

        if( strcmp( transform.c_str(), "affine" ) == 0 )
          {
          affineSearchTransform->Scale( bestScale );
          parametersList.push_back( affineSearchTransform->GetParameters() );
          }
        else if( strcmp( transform.c_str(), "rigid" ) == 0 )
          {
          rigidSearchTransform->SetIdentity();
          rigidSearchTransform->SetCenter( affineSearchTransform->GetCenter() );
          rigidSearchTransform->SetOffset( affineSearchTransform->GetOffset() );
          rigidSearchTransform->SetMatrix( affineSearchTransform->GetMatrix() );

          parametersList.push_back( rigidSearchTransform->GetParameters() );
          }
        else if( strcmp( transform.c_str(), "similarity" ) == 0 )
          {
          similaritySearchTransform->SetIdentity();
          similaritySearchTransform->SetCenter( affineSearchTransform->GetCenter() );
          similaritySearchTransform->SetOffset( affineSearchTransform->GetOffset() );
          similaritySearchTransform->SetMatrix( affineSearchTransform->GetMatrix() );
          similaritySearchTransform->SetScale( bestScale );

          parametersList.push_back( similaritySearchTransform->GetParameters() );
          }
        }
      }
    }
  multiStartOptimizer->SetParametersList( parametersList );

  if( numberOfIterations > 0 )
    {
    multiStartOptimizer->SetLocalOptimizer( localOptimizer );
    }
  multiStartOptimizer->StartOptimization();

  /////////////////////////////////////////////////////////////////
  //
  //         Write the output
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
      bestAffineTransform->SetParameters( multiStartOptimizer->GetBestParameters() );
      transformWriter->SetInput( bestAffineTransform );
      }
    else if( strcmp( transform.c_str(), "rigid" ) == 0 )
      {
      typename RigidTransformType::Pointer bestRigidTransform = RigidTransformType::New();
      bestRigidTransform->SetParameters( multiStartOptimizer->GetBestParameters() );
      transformWriter->SetInput( bestRigidTransform );
      }
    else if( strcmp( transform.c_str(), "similarity" ) == 0 )
      {
      typename SimilarityTransformType::Pointer bestSimilarityTransform = SimilarityTransformType::New();
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
  std::string description = std::string( "Several transform options are available.  The gradientStep or " )
    + std::string( "learningRate characterizes the gradient descent optimization and is scaled appropriately " )
    + std::string( "for each transform using the shift scales estimator. " );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "transform" );
  option->SetShortName( 't' );
  option->SetUsageOption(  0, "Rigid[gradientStep]" );
  option->SetUsageOption(  1, "Affine[gradientStep]" );
  option->SetUsageOption(  2, "Similarity[gradientStep]" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Boolean indicating alignment by principal axes." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "align-principal-axes" );
  option->SetShortName( 'p' );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Incremental search factor (in degrees)." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "search-factor" );
  option->SetShortName( 's' );
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

int antsAI( std::vector<std::string> args, std::ostream* /*out_stream = NULL */ )
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
  argv[argc] = ITK_NULLPTR;

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
