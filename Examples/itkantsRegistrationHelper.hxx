#ifndef __itkantsRegistrationHelper_hxx
#define __itkantsRegistrationHelper_hxx

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_copy.h>

namespace ants
{

/**
* GetShrinkImageOutputInformation provides a consistent way to compute the
* outputImage space for each level of a registration in a consistent way.
* By always using the same reference image, we can ensure that the same
* shrink results always are produced.
*/
template <typename TComputeType, unsigned VImageDimension>
typename itk::ImageBase<VImageDimension>::Pointer
RegistrationHelper<TComputeType, VImageDimension>::GetShrinkImageOutputInformation(const itk::ImageBase<VImageDimension> * inputImageInformation,
                                const typename RegistrationHelper<TComputeType, VImageDimension>::ShrinkFactorsPerDimensionContainerType &shrinkFactorsPerDimensionForCurrentLevel) const
{
  typedef itk::Image<unsigned char, VImageDimension> DummyImageType;

  typename DummyImageType::Pointer dummyImage = AllocImage<DummyImageType>( inputImageInformation, 0 );

  // We use the shrink image filter to calculate the fixed parameters of the virtual
  // domain at each level.  To speed up calculation and avoid unnecessary memory
  // usage, we could calculate these fixed parameters directly.

  typedef itk::ShrinkImageFilter<DummyImageType, DummyImageType> ShrinkFilterType;
  typename ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
  shrinkFilter->SetShrinkFactors( shrinkFactorsPerDimensionForCurrentLevel );
  shrinkFilter->SetInput( dummyImage );
  shrinkFilter->GenerateOutputInformation(); //Don't need to allocate space or run the filter, just create output information
  typename itk::ImageBase<VImageDimension>::Pointer returnImageBase=shrinkFilter->GetOutput();
  return returnImageBase;
}


template <typename TComputeType, unsigned VImageDimension>
RegistrationHelper<TComputeType, VImageDimension>
::RegistrationHelper() :
  m_CompositeTransform( nullptr ),
  m_RegistrationState( nullptr ),
  m_FixedInitialTransform( nullptr ),
  m_NumberOfStages( 0 ),
  m_Metrics(),
  m_TransformMethods(),
  m_Iterations(),
  m_SmoothingSigmas(),
  m_RestrictDeformationOptimizerWeights(),
  m_ShrinkFactors(),
  m_UseHistogramMatching( true ),
  m_WinsorizeImageIntensities( false ),
  m_DoEstimateLearningRateAtEachIteration( true ),
  m_LowerQuantile( 0.0 ),
  m_UpperQuantile( 1.0 ),
  m_LogStream( &std::cout ),
  m_PrintSimilarityMeasureInterval( 0 ),
  m_WriteIntervalVolumes( 0 ),
  m_InitializeTransformsPerStage( false ),
  m_AllPreviousTransformsAreLinear( true )
{
  typedef itk::LinearInterpolateImageFunction<ImageType, RealType> LinearInterpolatorType;
  typename LinearInterpolatorType::Pointer linearInterpolator = LinearInterpolatorType::New();
  this->m_Interpolator = linearInterpolator;
}

template <typename TComputeType, unsigned VImageDimension>
RegistrationHelper<TComputeType, VImageDimension>
::~RegistrationHelper()
= default;

template <typename ImageType>
typename ImageType::Pointer PreprocessImage( typename ImageType::ConstPointer  inputImage,
                                             typename ImageType::PixelType lowerScaleValue,
                                             typename ImageType::PixelType upperScaleValue,
                                             float winsorizeLowerQuantile, float winsorizeUpperQuantile,
                                             typename ImageType::ConstPointer histogramMatchSourceImage = nullptr )
{
  typedef itk::Statistics::ImageToHistogramFilter<ImageType>   HistogramFilterType;
  typedef typename HistogramFilterType::InputBooleanObjectType InputBooleanObjectType;
  typedef typename HistogramFilterType::HistogramSizeType      HistogramSizeType;

  HistogramSizeType histogramSize( 1 );
  histogramSize[0] = 256;

  typename InputBooleanObjectType::Pointer autoMinMaxInputObject = InputBooleanObjectType::New();
  autoMinMaxInputObject->Set( true );

  typename HistogramFilterType::Pointer histogramFilter = HistogramFilterType::New();
  histogramFilter->SetInput( inputImage );
  histogramFilter->SetAutoMinimumMaximumInput( autoMinMaxInputObject );
  histogramFilter->SetHistogramSize( histogramSize );
  histogramFilter->SetMarginalScale( 10.0 );
  histogramFilter->Update();

  float lowerValue = histogramFilter->GetOutput()->Quantile( 0, winsorizeLowerQuantile );
  float upperValue = histogramFilter->GetOutput()->Quantile( 0, winsorizeUpperQuantile );

  typedef itk::IntensityWindowingImageFilter<ImageType, ImageType> IntensityWindowingImageFilterType;

  typename IntensityWindowingImageFilterType::Pointer windowingFilter = IntensityWindowingImageFilterType::New();
  windowingFilter->SetInput( inputImage );
  windowingFilter->SetWindowMinimum( lowerValue );
  windowingFilter->SetWindowMaximum( upperValue );
  windowingFilter->SetOutputMinimum( lowerScaleValue );
  windowingFilter->SetOutputMaximum( upperScaleValue );
  windowingFilter->Update();

  typename ImageType::Pointer outputImage = nullptr;
  if( histogramMatchSourceImage )
    {
    typedef itk::HistogramMatchingImageFilter<ImageType, ImageType> HistogramMatchingFilterType;
    typename HistogramMatchingFilterType::Pointer matchingFilter = HistogramMatchingFilterType::New();
    matchingFilter->SetSourceImage( windowingFilter->GetOutput() );
    matchingFilter->SetReferenceImage( histogramMatchSourceImage );
    matchingFilter->SetNumberOfHistogramLevels( 256 );
    matchingFilter->SetNumberOfMatchPoints( 12 );
    matchingFilter->ThresholdAtMeanIntensityOn();
    matchingFilter->Update();

    outputImage = matchingFilter->GetOutput();
    outputImage->Update();
    outputImage->DisconnectPipeline();
    }
  else
    {
    outputImage = windowingFilter->GetOutput();
    outputImage->Update();
    outputImage->DisconnectPipeline();
    }
  return outputImage;
}

template <typename TComputeType, unsigned VImageDimension>
typename RegistrationHelper<TComputeType, VImageDimension>::MetricEnumeration
RegistrationHelper<TComputeType, VImageDimension>
::StringToMetricType( const std::string & str ) const
{
  if( str == "cc" )
    {
    return CC;
    }
  else if( str == "mi2" )
    {
    return MI;
    }
  else if( str == "mattes" || str == "mi" )
    {
    return Mattes;
    }
  else if( str == "meansquares" || str == "msq" || str == "ssd" )
    {
    return MeanSquares;
    }
  else if( str == "demons" )
    {
    return Demons;
    }
  else if( str == "gc" )
    {
    return GC;
    }
  else if( str == "icp" )
    {
    return ICP;
    }
  else if( str == "pse" )
    {
    return PSE;
    }
  else if( str == "jhct" )
    {
    return JHCT;
    }
  else if( str == "igdm" )
    {
    return IGDM;
    }
  return IllegalMetric;
}

template <typename TComputeType, unsigned VImageDimension>
typename RegistrationHelper<TComputeType, VImageDimension>::XfrmMethod
RegistrationHelper<TComputeType, VImageDimension>
::StringToXfrmMethod(const std::string & str) const
{
  if( str == "rigid" )
    {
    return Rigid;
    }
  else if( str == "affine" )
    {
    return Affine;
    }
  if( str == "compositeaffine" || str == "compaff" )
    {
    return CompositeAffine;
    }
  if( str == "similarity" )
    {
    return Similarity;
    }
  if( str == "translation" )
    {
    return Translation;
    }
  if( str == "bspline" ||
      str == "ffd" )
    {
    return BSpline;
    }
  if( str == "gaussiandisplacementfield" ||
      str == "gdf" )
    {
    return GaussianDisplacementField;
    }
  if( str == "bsplinedisplacementfield" ||
      str == "dmffd" )
    {
    return BSplineDisplacementField;
    }
  if( str == "timevaryingvelocityfield" ||
      str == "tvf" )
    {
    return TimeVaryingVelocityField;
    }
  if( str == "timevaryingbsplinevelocityfield" ||
      str == "tvdmffd" )
    {
    return TimeVaryingBSplineVelocityField;
    }
  if( str == "syn" ||
      str == "symmetricnormalization" )
    {
    return SyN;
    }
  if( str == "bsplinesyn" )
    {
    return BSplineSyN;
    }
  if( str == "exp" ||
      str == "exponential" )
    {
    return Exponential;
    }
  if( str == "bsplineexponential" )
    {
    return BSplineExponential;
    }
  return UnknownXfrm;
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddMetric( MetricEnumeration metricType,
             ImageType *fixedImage,
             ImageType *movingImage,
             LabeledPointSetType *fixedLabeledPointSet,
             LabeledPointSetType *movingLabeledPointSet,
             IntensityPointSetType *fixedIntensityPointSet,
             IntensityPointSetType *movingIntensityPointSet,
             unsigned int stageID,
             RealType weighting,
             SamplingStrategy samplingStrategy,
             int numberOfBins,
             unsigned int  radius,
             bool useBoundaryPointsOnly,
             RealType pointSetSigma,
             unsigned int evaluationKNeighborhood,
             RealType alpha,
             bool useAnisotropicCovariances,
             RealType samplingPercentage,
             RealType intensityDistanceSigma,
             RealType euclideanDistanceSigma )
{
  Metric init( metricType, fixedImage, movingImage,
               fixedLabeledPointSet, movingLabeledPointSet,
               fixedIntensityPointSet, movingIntensityPointSet,
               stageID, weighting, samplingStrategy, numberOfBins, radius,
               useBoundaryPointsOnly, pointSetSigma, evaluationKNeighborhood,
               alpha, useAnisotropicCovariances, samplingPercentage,
               intensityDistanceSigma, euclideanDistanceSigma );

  this->m_Metrics.push_back( init );
}

template <typename TComputeType, unsigned VImageDimension>
typename RegistrationHelper<TComputeType, VImageDimension>::MetricListType
RegistrationHelper<TComputeType, VImageDimension>
::GetMetricListPerStage( unsigned int stageID )
{
  MetricListType stageMetricList;

  typename MetricListType::const_iterator it;
  for( it = this->m_Metrics.begin(); it != this->m_Metrics.end(); ++it )
    {
    if( ( *it ).m_StageID == stageID )
      {
      stageMetricList.push_back( *it );
      }
    }

  return stageMetricList;
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddRigidTransform(RealType GradientStep)
{
  TransformMethod init;

  init.m_XfrmMethod = Rigid;
  init.m_GradientStep = GradientStep;
  this->m_TransformMethods.push_back( init );
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddAffineTransform(RealType GradientStep)
{
  TransformMethod init;

  init.m_XfrmMethod = Affine;
  init.m_GradientStep = GradientStep;
  this->m_TransformMethods.push_back( init );
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddCompositeAffineTransform(RealType GradientStep)
{
  TransformMethod init;

  init.m_XfrmMethod = CompositeAffine;
  init.m_GradientStep = GradientStep;
  this->m_TransformMethods.push_back( init );
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddSimilarityTransform(RealType GradientStep)
{
  TransformMethod init;

  init.m_XfrmMethod = Similarity;
  init.m_GradientStep = GradientStep;
  this->m_TransformMethods.push_back( init );
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddTranslationTransform(RealType GradientStep)
{
  TransformMethod init;

  init.m_XfrmMethod = Translation;
  init.m_GradientStep = GradientStep;
  this->m_TransformMethods.push_back( init );
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddBSplineTransform(RealType GradientStep, std::vector<unsigned int> & MeshSizeAtBaseLevel)
{
  TransformMethod init;

  init.m_XfrmMethod = BSpline;
  init.m_GradientStep = GradientStep;
  init.m_MeshSizeAtBaseLevel = MeshSizeAtBaseLevel;
  this->m_TransformMethods.push_back( init );
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddGaussianDisplacementFieldTransform(RealType GradientStep, RealType UpdateFieldVarianceInVarianceSpace,
                                        RealType TotalFieldVarianceInVarianceSpace)
{
  TransformMethod init;

  init.m_XfrmMethod = GaussianDisplacementField;
  init.m_GradientStep = GradientStep;
  init.m_UpdateFieldVarianceInVarianceSpace = UpdateFieldVarianceInVarianceSpace;
  init.m_TotalFieldVarianceInVarianceSpace = TotalFieldVarianceInVarianceSpace;
  this->m_TransformMethods.push_back( init );
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddBSplineDisplacementFieldTransform(RealType GradientStep,
                                       std::vector<unsigned int> & UpdateFieldMeshSizeAtBaseLevel,
                                       std::vector<unsigned int> & TotalFieldMeshSizeAtBaseLevel,
                                       unsigned int SplineOrder)
{
  TransformMethod init;

  init.m_XfrmMethod = BSplineDisplacementField;
  init.m_GradientStep = GradientStep;
  init.m_UpdateFieldMeshSizeAtBaseLevel = UpdateFieldMeshSizeAtBaseLevel;
  init.m_TotalFieldMeshSizeAtBaseLevel = TotalFieldMeshSizeAtBaseLevel;
  init.m_SplineOrder = SplineOrder;
  this->m_TransformMethods.push_back( init );
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddTimeVaryingVelocityFieldTransform( RealType GradientStep,
                                        unsigned int NumberOfTimeIndices,
                                        RealType UpdateFieldVarianceInVarianceSpace,
                                        RealType UpdateFieldTimeSigma,
                                        RealType TotalFieldVarianceInVarianceSpace,
                                        RealType TotalFieldTimeSigma )
{
  TransformMethod init;

  init.m_XfrmMethod = TimeVaryingVelocityField;
  init.m_GradientStep = GradientStep;
  init.m_NumberOfTimeIndices = NumberOfTimeIndices;
  init.m_UpdateFieldVarianceInVarianceSpace = UpdateFieldVarianceInVarianceSpace;
  init.m_UpdateFieldTimeSigma = UpdateFieldTimeSigma;
  init.m_TotalFieldVarianceInVarianceSpace = TotalFieldVarianceInVarianceSpace;
  init.m_TotalFieldTimeSigma = TotalFieldTimeSigma;
  this->m_TransformMethods.push_back( init );
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddTimeVaryingBSplineVelocityFieldTransform( RealType GradientStep, std::vector<unsigned int> VelocityFieldMeshSize,
                                               unsigned int NumberOfTimePointSamples, unsigned int SplineOrder )
{
  TransformMethod init;

  init.m_XfrmMethod = TimeVaryingBSplineVelocityField;;
  init.m_GradientStep = GradientStep;
  init.m_VelocityFieldMeshSize = VelocityFieldMeshSize;
  init.m_NumberOfTimePointSamples = NumberOfTimePointSamples;
  init.m_SplineOrder = SplineOrder;
  this->m_TransformMethods.push_back( init );
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddSyNTransform( RealType GradientStep, RealType UpdateFieldVarianceInVarianceSpace,
                   RealType TotalFieldVarianceInVarianceSpace )
{
  TransformMethod init;

  init.m_XfrmMethod = SyN;
  init.m_GradientStep = GradientStep;
  init.m_UpdateFieldVarianceInVarianceSpace = UpdateFieldVarianceInVarianceSpace;
  init.m_TotalFieldVarianceInVarianceSpace = TotalFieldVarianceInVarianceSpace;
  this->m_TransformMethods.push_back( init );
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddBSplineSyNTransform( RealType GradientStep, std::vector<unsigned int> &  UpdateFieldMeshSizeAtBaseLevel,
                          std::vector<unsigned int> &  TotalFieldMeshSizeAtBaseLevel,
                          unsigned int SplineOrder )
{
  TransformMethod init;

  init.m_XfrmMethod = BSplineSyN;
  init.m_GradientStep = GradientStep;
  init.m_UpdateFieldMeshSizeAtBaseLevel = UpdateFieldMeshSizeAtBaseLevel;
  init.m_TotalFieldMeshSizeAtBaseLevel = TotalFieldMeshSizeAtBaseLevel;
  init.m_SplineOrder = SplineOrder;
  this->m_TransformMethods.push_back( init );
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddExponentialTransform( RealType GradientStep, RealType UpdateFieldVarianceInVarianceSpace,
                           RealType VelocityFieldVarianceInVarianceSpace, unsigned int NumberOfIntegrationSteps )
{
  TransformMethod init;

  init.m_XfrmMethod = Exponential;
  init.m_GradientStep = GradientStep;
  init.m_UpdateFieldVarianceInVarianceSpace = UpdateFieldVarianceInVarianceSpace;
  init.m_VelocityFieldVarianceInVarianceSpace = VelocityFieldVarianceInVarianceSpace;
  init.m_NumberOfTimeIndices = NumberOfIntegrationSteps;

  this->m_TransformMethods.push_back( init );
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddBSplineExponentialTransform( RealType GradientStep, std::vector<unsigned int> &  UpdateFieldMeshSizeAtBaseLevel,
                                  std::vector<unsigned int> & VelocityFieldMeshSizeAtBaseLevel,
                                  unsigned int NumberOfIntegrationSteps,
                                  unsigned int SplineOrder )
{
  TransformMethod init;

  init.m_XfrmMethod = BSplineExponential;
  init.m_GradientStep = GradientStep;
  init.m_UpdateFieldMeshSizeAtBaseLevel = UpdateFieldMeshSizeAtBaseLevel;
  init.m_VelocityFieldMeshSizeAtBaseLevel = VelocityFieldMeshSizeAtBaseLevel;
  init.m_SplineOrder = SplineOrder;
  init.m_NumberOfTimeIndices = NumberOfIntegrationSteps;

  this->m_TransformMethods.push_back( init );
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetIterations( const std::vector<std::vector<unsigned int> > & Iterations )
{
  this->m_Iterations = Iterations;
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetConvergenceThresholds( const std::vector<RealType> & thresholds )
{
  this->m_ConvergenceThresholds = thresholds;
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetConvergenceWindowSizes( const std::vector<unsigned int> & windowSizes )
{
  this->m_ConvergenceWindowSizes = windowSizes;
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetSmoothingSigmas( const std::vector<std::vector<float> > & SmoothingSigmas )
{
  this->m_SmoothingSigmas = SmoothingSigmas;
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetRestrictDeformationOptimizerWeights( const std::vector<std::vector<RealType> > & restrictDeformationWeights )
{
  this->m_RestrictDeformationOptimizerWeights = restrictDeformationWeights;
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetSmoothingSigmasAreInPhysicalUnits( const std::vector<bool> & SmoothingSigmasAreInPhysicalUnits )
{
  this->m_SmoothingSigmasAreInPhysicalUnits = SmoothingSigmasAreInPhysicalUnits;
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetShrinkFactors( const std::vector<std::vector<unsigned int> > & ShrinkFactors )
{
  this->m_ShrinkFactors = ShrinkFactors;
}

template <typename TComputeType, unsigned VImageDimension>
typename RegistrationHelper<TComputeType, VImageDimension>::ShrinkFactorsPerDimensionContainerType
RegistrationHelper<TComputeType, VImageDimension>
::CalculateShrinkFactorsPerDimension( unsigned int factor, ImageSpacingType spacing )
{
  using SpacingValueType = typename ImageSpacingType::ComponentType;

  SpacingValueType minSpacing = spacing[0];
  unsigned int minIndex = 0;
  for( unsigned int n = 1; n < VImageDimension; n++ )
    {
    if( minSpacing > static_cast<SpacingValueType>( spacing[n] ) )
      {
      minSpacing = spacing[n];
      minIndex = n;
      }
    }

  ShrinkFactorsPerDimensionContainerType shrinkFactorsPerDimension;
  shrinkFactorsPerDimension.Fill( 0 );
  shrinkFactorsPerDimension[minIndex] = factor;

  ImageSpacingType newSpacing;
  newSpacing[minIndex] = spacing[minIndex] * factor;

  for( unsigned int n = 0; n < VImageDimension; n++ )
    {
    if( shrinkFactorsPerDimension[n] == 0 )
      {
      SpacingValueType newMinSpacing = static_cast<SpacingValueType>( spacing[n] ) *
        static_cast<SpacingValueType>( factor );
      RealType minDifferenceFromMinSpacing = static_cast<RealType>( std::fabs( newMinSpacing - newSpacing[minIndex] ) );
      unsigned int minFactor = factor;
      for( unsigned int f = factor - 1; f > 0; f-- )
        {
        newMinSpacing = static_cast<SpacingValueType>( spacing[n] ) * static_cast<SpacingValueType>( f );

        // We use <= such that the smaller factor is preferred if distances are the same
        if( static_cast<RealType>( std::fabs( newMinSpacing - newSpacing[minIndex] ) ) <= minDifferenceFromMinSpacing )
          {
          minDifferenceFromMinSpacing = itk::Math::abs ( newMinSpacing - newSpacing[minIndex] );
          minFactor = f;
          }
        }
      shrinkFactorsPerDimension[n] = minFactor;
      }
    }
  return shrinkFactorsPerDimension;
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetWinsorizeImageIntensities( bool Winsorize, float LowerQuantile, float UpperQuantile )
{
  this->m_WinsorizeImageIntensities = Winsorize;
  this->m_LowerQuantile = LowerQuantile;
  this->m_UpperQuantile = UpperQuantile;
}

template <typename TComputeType, unsigned VImageDimension>
int
RegistrationHelper<TComputeType, VImageDimension>
::ValidateParameters()
{
  if( this->m_NumberOfStages == 0 )
    {
    this->Logger() << "No transformations are specified." << std::endl;
    return EXIT_FAILURE;
    }
  if( this->m_Iterations.size() != this->m_NumberOfStages )
    {
    this->Logger() << "The number of iteration sets specified does not match the number of stages." << std::endl;
    return EXIT_FAILURE;
    }
  if( this->m_ShrinkFactors.size() != this->m_NumberOfStages )
    {
    this->Logger() << "The number of shrinkFactors specified does not match the number of stages." << std::endl;
    return EXIT_FAILURE;
    }
  if( this->m_SmoothingSigmas.size() != this->m_NumberOfStages )
    {
    this->Logger() << "The number of smoothing sigma sets specified does not match the number of stages."
                     << std::endl;
    return EXIT_FAILURE;
    }
  if( this->m_SmoothingSigmasAreInPhysicalUnits.size() != this->m_NumberOfStages )
    {
    this->Logger()
      << "The number of smoothing sigma in physical units bool values does not match the number of stages."
      << std::endl;
    return EXIT_FAILURE;
    }
  for( unsigned int i = 0; i < this->m_Metrics.size(); i++ )
    {
    if( !this->IsPointSetMetric( this->m_Metrics[i].m_MetricType ) )
      {
      if( this->m_Metrics[i].m_FixedImage.IsNull() ||
        this->m_Metrics[i].m_MovingImage.IsNull() )
        {
        this->Logger() << "The image metric has no fixed and/or moving image." << std::endl;
        return EXIT_FAILURE;
        }
      }
    }

  // Check the number of masks.  We are going to allow the user 2 options w.r.t.
  // mask specification:
  //   1. Either the user specifies a single mask to be used for all stages or
  //   2. the user specifies a mask for each stage.
  // Note that we handle the fixed and moving masks separately to enforce this constraint.

  if( this->m_FixedImageMasks.size() > 1 && this->m_FixedImageMasks.size() != this->m_NumberOfStages )
    {
    this->Logger() << "The number of fixed masks must be equal to 1 (use the mask for all "
                   << "stages) or the number of fixed masks must be equal to the number of stages." << std::endl;
    return EXIT_FAILURE;
    }

  if( this->m_MovingImageMasks.size() > 1 && this->m_MovingImageMasks.size() != this->m_NumberOfStages )
    {
    this->Logger() << "The number of moving masks must be equal to 1 (i.e., use the mask for all "
                   << "stages) or the number of moving masks must be equal to the number of stages." << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

template <typename TComputeType, unsigned VImageDimension>
typename RegistrationHelper<TComputeType, VImageDimension>::ImageType::Pointer
RegistrationHelper<TComputeType, VImageDimension>
::GetWarpedImage() const
{
  typename ImageType::Pointer fixedImage = this->m_Metrics[0].m_FixedImage;
  typename ImageType::Pointer movingImage = this->m_Metrics[0].m_MovingImage;

  typedef itk::ResampleImageFilter<ImageType, ImageType, RealType> ResampleFilterType;
  typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetTransform( this->m_CompositeTransform );
  resampler->SetInput( movingImage );
  resampler->SetOutputParametersFromImage( fixedImage );
  resampler->SetInterpolator( this->m_Interpolator );
  resampler->SetDefaultPixelValue( 0 );
  resampler->Update();

  typename ImageType::Pointer  WarpedImage;
  WarpedImage = resampler->GetOutput();
  return WarpedImage.GetPointer();
}

template <typename TComputeType, unsigned VImageDimension>
typename RegistrationHelper<TComputeType, VImageDimension>::ImageType::Pointer
RegistrationHelper<TComputeType, VImageDimension>
::GetInverseWarpedImage() const
{
  typename ImageType::Pointer fixedImage = this->m_Metrics[0].m_FixedImage;
  typename ImageType::Pointer movingImage = this->m_Metrics[0].m_MovingImage;

  if( this->m_CompositeTransform->GetInverseTransform().IsNull() )
    {
    return nullptr;
    }
  typedef itk::ResampleImageFilter<ImageType, ImageType, RealType> ResampleFilterType;
  typename ResampleFilterType::Pointer inverseResampler = ResampleFilterType::New();
  inverseResampler->SetTransform( this->m_CompositeTransform->GetInverseTransform() );
  inverseResampler->SetInput( fixedImage );
  inverseResampler->SetOutputParametersFromImage( movingImage );
  inverseResampler->SetInterpolator( this->m_Interpolator );
  inverseResampler->SetDefaultPixelValue( 0 );
  inverseResampler->Update();

  typename ImageType::Pointer InverseWarpedImage;
  InverseWarpedImage = inverseResampler->GetOutput();
  return InverseWarpedImage.GetPointer();
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddFixedImageMask( typename MaskImageType::Pointer & fixedImageMask )
{
  typename ImageMaskSpatialObjectType::Pointer so = nullptr;

  if( fixedImageMask.IsNotNull() )
    {
    so = ImageMaskSpatialObjectType::New();
    so->SetImage( fixedImageMask.GetPointer() );
    }
  this->AddFixedImageMask( so );
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddMovingImageMask( typename MaskImageType::Pointer & movingImageMask )
{
  typename ImageMaskSpatialObjectType::Pointer so = nullptr;

  if( movingImageMask.IsNotNull() )
    {
    so = ImageMaskSpatialObjectType::New();
    so->SetImage( movingImageMask.GetPointer() );
    }
  this->AddMovingImageMask( so );
}

template <typename TComputeType, unsigned VImageDimension>
int
RegistrationHelper<TComputeType, VImageDimension>
::DoRegistration()
{
  /** Can really impact performance */
  const bool     gradientfilter = false;
  itk::TimeProbe totalTimer;

  totalTimer.Start();

  this->m_NumberOfStages = this->m_TransformMethods.size();

  if( this->ValidateParameters() != EXIT_SUCCESS )
    {
    return EXIT_FAILURE;
    }
  this->PrintState();

  this->Logger() << "Registration using " << this->m_NumberOfStages << " total stages." << std::endl;

  // NOTE:  the -1 is to ignore the initial identity identity transform
  if( this->m_CompositeTransform.IsNull() )
    {
    this->m_CompositeTransform = CompositeTransformType::New();
    }
  if( this->m_FixedInitialTransform.IsNull() )
    {
    this->m_FixedInitialTransform = CompositeTransformType::New();
    }

  // ########################################################################################
  // ########################################################################################
  // ##The main loop for exstimating the total composite transform
  // ########################################################################################
  // ########################################################################################
  for( unsigned int currentStageNumber = 0; currentStageNumber < this->m_NumberOfStages; currentStageNumber++ )
    {
    itk::TimeProbe timer;
    timer.Start();

    this->Logger() << std::endl << "Stage " << currentStageNumber << std::endl;
    std::stringstream currentStageString;
    currentStageString << currentStageNumber;

    // Get the number of iterations and use that information to specify the number of levels

    const std::vector<unsigned int> & currentStageIterations = this->m_Iterations[currentStageNumber];
    this->Logger() << "  iterations = ";
    for( unsigned int m = 0; m < currentStageIterations.size(); m++ )
      {
      this->Logger() << currentStageIterations[m];
      if( m < currentStageIterations.size() - 1 )
        {
        this->Logger() << 'x';
        }
      }
    this->Logger() << std::endl;

    const RealType convergenceThreshold = this->m_ConvergenceThresholds[currentStageNumber];
    this->Logger() << "  convergence threshold = " << convergenceThreshold << std::endl;
    const unsigned int convergenceWindowSize = this->m_ConvergenceWindowSizes[currentStageNumber];
    this->Logger() << "  convergence window size = " << convergenceWindowSize << std::endl;

    const unsigned int numberOfLevels = currentStageIterations.size();
    this->Logger() << "  number of levels = " << numberOfLevels << std::endl;

    unsigned int fixedMaskIndex = itk::NumericTraits<unsigned int>::max();
    unsigned int movingMaskIndex = itk::NumericTraits<unsigned int>::max();
    bool useFixedImageMaskForThisStage = false;
    bool useMovingImageMaskForThisStage = false;

    // We already checked that number of masks = 1 or = number of stages
    if( this->m_FixedImageMasks.size() > 0 )
      {
      useFixedImageMaskForThisStage = true;

      if( this->m_FixedImageMasks.size() == 1 )
        {
        fixedMaskIndex = 0;
        }
      else
        {
        fixedMaskIndex = currentStageNumber;
        }
      }
    if( this->m_MovingImageMasks.size() > 0 )
      {
      useMovingImageMaskForThisStage = true;

      if( this->m_MovingImageMasks.size() == 1 )
        {
        movingMaskIndex = 0;
        }
      else
        {
        movingMaskIndex = currentStageNumber;
        }
      }

    // Get the number of metrics at the current stage.  If more than one metric
    // then we need to use the MultiMetricType.  Due to the way the metrics are
    // pulled off the command line stack, we need to iterate from the top down.

    MetricListType stageMetricList = this->GetMetricListPerStage( this->m_NumberOfStages - currentStageNumber - 1 );

    typename ObjectMetricType::Pointer singleMetric;
    typename MultiMetricType::Pointer multiMetric;

    typename MultiMetricType::WeightsArrayType metricWeights( stageMetricList.size() );
    metricWeights.Fill( 1.0 );

    bool useMultiMetric = false;
    if( stageMetricList.size() > 1 )
      {
      useMultiMetric = true;
      multiMetric = MultiMetricType::New();
      }

    std::vector<typename ImageType::Pointer> preprocessedFixedImagesPerStage;
    std::vector<typename ImageType::Pointer> preprocessedMovingImagesPerStage;

    typename ImageBaseType::Pointer virtualDomainImage = nullptr;

    for( unsigned int currentMetricNumber = 0; currentMetricNumber < stageMetricList.size(); currentMetricNumber++ )
      {
      MetricEnumeration currentMetricType = stageMetricList[currentMetricNumber].m_MetricType;

      typename ImageMetricType::Pointer imageMetric = nullptr;

      typedef itk::LabeledPointSetToPointSetMetricv4<LabeledPointSetType, LabeledPointSetType, RealType> LabeledPointSetMetricType;
      typename LabeledPointSetMetricType::Pointer labeledPointSetMetric = LabeledPointSetMetricType::New();

      typedef itk::MeanSquaresPointSetToPointSetIntensityMetricv4<IntensityPointSetType, IntensityPointSetType, RealType> IntensityPointSetMetricType;
      typename IntensityPointSetMetricType::Pointer intensityPointSetMetric = nullptr;

      switch( currentMetricType )
        {
        case CC:
          {
          const unsigned int radiusOption = stageMetricList[currentMetricNumber].m_Radius;
          this->Logger() << "  using the CC metric (radius = "
                         << radiusOption << ", weight = "
                         << stageMetricList[currentMetricNumber].m_Weighting << ")" << std::endl;
          typedef itk::ANTSNeighborhoodCorrelationImageToImageMetricv4<ImageType, ImageType, ImageType, TComputeType> CorrelationMetricType;
          typename CorrelationMetricType::Pointer correlationMetric = CorrelationMetricType::New();
            {
            typename CorrelationMetricType::RadiusType radius;
            radius.Fill( radiusOption );
            correlationMetric->SetRadius( radius );
            }
          correlationMetric->SetUseMovingImageGradientFilter( gradientfilter );
          correlationMetric->SetUseFixedImageGradientFilter( gradientfilter );

          imageMetric = correlationMetric;
          }
          break;
        case Mattes:
          {
          const unsigned int binOption = stageMetricList[currentMetricNumber].m_NumberOfBins;
          this->Logger() << "  using the Mattes MI metric (number of bins = "
                         << binOption << ", weight = "
                         << stageMetricList[currentMetricNumber].m_Weighting << ")" << std::endl;
          typedef itk::MattesMutualInformationImageToImageMetricv4<ImageType, ImageType, ImageType, TComputeType> MutualInformationMetricType;
          typename MutualInformationMetricType::Pointer mutualInformationMetric = MutualInformationMetricType::New();
          // mutualInformationMetric = mutualInformationMetric;
          mutualInformationMetric->SetNumberOfHistogramBins( binOption );
          mutualInformationMetric->SetUseMovingImageGradientFilter( gradientfilter );
          mutualInformationMetric->SetUseFixedImageGradientFilter( gradientfilter );
          mutualInformationMetric->SetUseSampledPointSet( false );

          imageMetric = mutualInformationMetric;
          }
          break;
        case MI:
          {
          const unsigned int binOption = stageMetricList[currentMetricNumber].m_NumberOfBins;
          this->Logger() << "  using the joint histogram MI metric (number of bins = "
                         << binOption << ", weight = "
                         << stageMetricList[currentMetricNumber].m_Weighting << ")" << std::endl;
          typedef itk::JointHistogramMutualInformationImageToImageMetricv4<ImageType, ImageType, ImageType,
                                                                           TComputeType> MutualInformationMetricType;
          typename MutualInformationMetricType::Pointer mutualInformationMetric = MutualInformationMetricType::New();
          //mutualInformationMetric = mutualInformationMetric;
          mutualInformationMetric->SetNumberOfHistogramBins( binOption );
          mutualInformationMetric->SetUseMovingImageGradientFilter( gradientfilter );
          mutualInformationMetric->SetUseFixedImageGradientFilter( gradientfilter );
          mutualInformationMetric->SetUseSampledPointSet( false );
          mutualInformationMetric->SetVarianceForJointPDFSmoothing( 1.0 );

          imageMetric = mutualInformationMetric;
          }
          break;
        case MeanSquares:
          {
          this->Logger() << "  using the MeanSquares metric (weight = "
                         << stageMetricList[currentMetricNumber].m_Weighting << ")" << std::endl;

          typedef itk::MeanSquaresImageToImageMetricv4<ImageType, ImageType, ImageType, TComputeType> MeanSquaresMetricType;
          typename MeanSquaresMetricType::Pointer meanSquaresMetric = MeanSquaresMetricType::New();
          //meanSquaresMetric = meanSquaresMetric;

          imageMetric = meanSquaresMetric;
          }
          break;
        case Demons:
          {
          this->Logger() << "  using the Demons metric (weight = "
                         << stageMetricList[currentMetricNumber].m_Weighting << ")" << std::endl;

          typedef itk::DemonsImageToImageMetricv4<ImageType, ImageType, ImageType, TComputeType> DemonsMetricType;
          typename DemonsMetricType::Pointer demonsMetric = DemonsMetricType::New();

          imageMetric = demonsMetric;
          }
          break;
        case GC:
          {
          this->Logger() << "  using the global correlation metric (weight = "
                         << stageMetricList[currentMetricNumber].m_Weighting << ")" << std::endl;
          typedef itk::CorrelationImageToImageMetricv4<ImageType, ImageType, ImageType, TComputeType> corrMetricType;
          typename corrMetricType::Pointer corrMetric = corrMetricType::New();

          imageMetric = corrMetric;
          }
          break;
        case ICP:
          {
          this->Logger() << "  using the ICP metric (weight = "
                         << stageMetricList[currentMetricNumber].m_Weighting << ")" << std::endl;
          typedef itk::EuclideanDistancePointSetToPointSetMetricv4<LabeledPointSetType, LabeledPointSetType, RealType> IcpPointSetMetricType;
          typename IcpPointSetMetricType::Pointer icpMetric = IcpPointSetMetricType::New();

          labeledPointSetMetric->SetPointSetMetric( icpMetric.GetPointer() );
          }
          break;
        case PSE:
          {
          this->Logger() << "  using the PSE metric (weight = "
                         << stageMetricList[currentMetricNumber].m_Weighting << ")" << std::endl;
          typedef itk::ExpectationBasedPointSetToPointSetMetricv4<LabeledPointSetType, LabeledPointSetType, RealType> PsePointSetMetricType;
          typename PsePointSetMetricType::Pointer pseMetric = PsePointSetMetricType::New();
          pseMetric->SetPointSetSigma( stageMetricList[currentMetricNumber].m_PointSetSigma );
          pseMetric->SetEvaluationKNeighborhood( stageMetricList[currentMetricNumber].m_EvaluationKNeighborhood );

          labeledPointSetMetric->SetPointSetMetric( pseMetric.GetPointer() );
          }
          break;
        case JHCT:
          {
          this->Logger() << "  using the JHCT metric (weight = "
                         << stageMetricList[currentMetricNumber].m_Weighting << ")" << std::endl;
          typedef itk::JensenHavrdaCharvatTsallisPointSetToPointSetMetricv4<LabeledPointSetType, RealType> JhctPointSetMetricType;
          typename JhctPointSetMetricType::Pointer jhctMetric = JhctPointSetMetricType::New();
          jhctMetric->SetPointSetSigma( stageMetricList[currentMetricNumber].m_PointSetSigma );
          jhctMetric->SetKernelSigma( 10.0 );
          jhctMetric->SetUseAnisotropicCovariances( stageMetricList[currentMetricNumber].m_UseAnisotropicCovariances );
          jhctMetric->SetCovarianceKNeighborhood( 5 );
          jhctMetric->SetEvaluationKNeighborhood( stageMetricList[currentMetricNumber].m_EvaluationKNeighborhood );
          jhctMetric->SetAlpha( stageMetricList[currentMetricNumber].m_Alpha );

          labeledPointSetMetric->SetPointSetMetric( jhctMetric.GetPointer() );
          }
          break;
        case IGDM:
          {
          this->Logger() << "  using the IGDM metric (weight = "
                         << stageMetricList[currentMetricNumber].m_Weighting << ")" << std::endl;
          typedef itk::MeanSquaresPointSetToPointSetIntensityMetricv4<IntensityPointSetType, IntensityPointSetType, RealType> MsqPointSetMetricType;
          typename MsqPointSetMetricType::Pointer msqMetric = MsqPointSetMetricType::New();

          msqMetric->SetIntensityDistanceSigma( stageMetricList[currentMetricNumber].m_IntensityDistanceSigma );
          msqMetric->SetEuclideanDistanceSigma( stageMetricList[currentMetricNumber].m_EuclideanDistanceSigma );
          if( msqMetric->GetEuclideanDistanceSigma() <= itk::NumericTraits<RealType>::ZeroValue() )
            {
            msqMetric->EstimateEuclideanDistanceSigmaAutomaticallyOn();
            }
          else
            {
            msqMetric->EstimateEuclideanDistanceSigmaAutomaticallyOff();
            }
          if( msqMetric->GetIntensityDistanceSigma() <= itk::NumericTraits<RealType>::ZeroValue() )
            {
            msqMetric->EstimateIntensityDistanceSigmaAutomaticallyOn();
            }
          else
            {
            msqMetric->EstimateIntensityDistanceSigmaAutomaticallyOff();
            }

          intensityPointSetMetric = msqMetric;
          }
          break;
        default:
          this->Logger() << "ERROR: Unrecognized metric. " << std::endl;
          return EXIT_FAILURE;
        }

      if( !this->IsPointSetMetric( currentMetricType ) )
        {
        // Get the fixed and moving images
        const typename ImageType::ConstPointer fixedImage =
          stageMetricList[currentMetricNumber].m_FixedImage.GetPointer();
        const typename ImageType::ConstPointer movingImage =
          stageMetricList[currentMetricNumber].m_MovingImage.GetPointer();

        // Preprocess images

        std::string outputPreprocessingString = "";

        PixelType lowerScaleValue = 0.0;
        PixelType upperScaleValue = 1.0;
        if( this->m_WinsorizeImageIntensities )
          {
          outputPreprocessingString += "  preprocessing:  winsorizing the image intensities\n";
          }

        typename ImageType::Pointer preprocessFixedImage =
          PreprocessImage<ImageType>( fixedImage.GetPointer(), lowerScaleValue,
                                      upperScaleValue, this->m_LowerQuantile, this->m_UpperQuantile,
                                      nullptr );

        preprocessedFixedImagesPerStage.push_back( preprocessFixedImage.GetPointer() );

        typename ImageType::Pointer preprocessMovingImage =
          PreprocessImage<ImageType>( movingImage.GetPointer(), lowerScaleValue,
                                      upperScaleValue, this->m_LowerQuantile, this->m_UpperQuantile,
                                      nullptr );

        if( this->m_UseHistogramMatching )
          {
          outputPreprocessingString += "  preprocessing:  histogram matching the images\n";
          preprocessMovingImage =
            PreprocessImage<ImageType>( movingImage.GetPointer(),
                                        lowerScaleValue, upperScaleValue,
                                        this->m_LowerQuantile, this->m_UpperQuantile,
                                        preprocessFixedImage.GetPointer() );
          }
        preprocessedMovingImagesPerStage.push_back( preprocessMovingImage.GetPointer() );

        this->Logger() << outputPreprocessingString << std::flush;

        // Set up the image metric and scales estimator

        imageMetric->SetVirtualDomainFromImage( fixedImage );
        imageMetric->SetUseMovingImageGradientFilter( gradientfilter );
        imageMetric->SetUseFixedImageGradientFilter( gradientfilter );
        metricWeights[currentMetricNumber] = stageMetricList[currentMetricNumber].m_Weighting;
        if( useFixedImageMaskForThisStage )
          {
          imageMetric->SetFixedImageMask( this->m_FixedImageMasks[fixedMaskIndex] );
          }
        if( useMovingImageMaskForThisStage )
          {
          imageMetric->SetMovingImageMask( this->m_MovingImageMasks[movingMaskIndex] );
          }
        if( virtualDomainImage.IsNull() )
          {
          virtualDomainImage = imageMetric->GetModifiableVirtualImage();
          }

        if( useMultiMetric )
          {
          multiMetric->AddMetric( imageMetric );
          }
        if( !useMultiMetric || currentMetricNumber == 0 )
          {
          singleMetric = static_cast<ObjectMetricType *>( imageMetric );
          }
        }
      else
        {
        preprocessedFixedImagesPerStage.push_back( nullptr );
        preprocessedMovingImagesPerStage.push_back( nullptr );

        metricWeights[currentMetricNumber] = stageMetricList[currentMetricNumber].m_Weighting;

        if( currentMetricType == IGDM  )
          {
          if( useFixedImageMaskForThisStage )
            {
            typedef itk::CastImageFilter<MaskImageType, typename LabeledPointSetMetricType::VirtualImageType> CasterType;
            typename CasterType::Pointer caster = CasterType::New();
            caster->SetInput( this->m_FixedImageMasks[fixedMaskIndex]->GetImage() );
            caster->Update();

            intensityPointSetMetric->SetVirtualDomainFromImage( caster->GetOutput() );
            if( virtualDomainImage.IsNull() )
              {
              virtualDomainImage = intensityPointSetMetric->GetModifiableVirtualImage();
              }
            }

          if( useMultiMetric )
            {
            multiMetric->AddMetric( intensityPointSetMetric );
            }
          if( !useMultiMetric || currentMetricNumber == 0 )
            {
            intensityPointSetMetric->SetFixedPointSet( stageMetricList[currentMetricNumber].m_FixedIntensityPointSet );
            intensityPointSetMetric->SetMovingPointSet( stageMetricList[currentMetricNumber].m_MovingIntensityPointSet );
            singleMetric = static_cast<ObjectMetricType *>( intensityPointSetMetric );
            }
          }
        else
          {
          if( useFixedImageMaskForThisStage )
            {
            typedef itk::CastImageFilter<MaskImageType, typename LabeledPointSetMetricType::VirtualImageType> CasterType;
            typename CasterType::Pointer caster = CasterType::New();
            caster->SetInput( this->m_FixedImageMasks[fixedMaskIndex]->GetImage() );
            caster->Update();

            labeledPointSetMetric->SetVirtualDomainFromImage( caster->GetOutput() );
            if( virtualDomainImage.IsNull() )
              {
              virtualDomainImage = labeledPointSetMetric->GetModifiableVirtualImage();
              }
            }
          if( useMultiMetric )
            {
            multiMetric->AddMetric( labeledPointSetMetric );
            }
          if( !useMultiMetric || currentMetricNumber == 0 )
            {
            labeledPointSetMetric->SetFixedPointSet( stageMetricList[currentMetricNumber].m_FixedLabeledPointSet );
            labeledPointSetMetric->SetMovingPointSet( stageMetricList[currentMetricNumber].m_MovingLabeledPointSet );
            singleMetric = static_cast<ObjectMetricType *>( labeledPointSetMetric );
            }
          }
        }
      }
    if( useMultiMetric )
      {
      multiMetric->SetMetricWeights( metricWeights );
      }

    // These two variables are specified in setting up the registration method.
    // However, for point set metrics, they are not required.
    std::vector<ShrinkFactorsPerDimensionContainerType> shrinkFactorsPerDimensionForAllLevels;
    typename AffineRegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;

    // Get shrink factors and adjust according to the current image
    const std::vector<unsigned int> factors( this->m_ShrinkFactors[currentStageNumber] );
    if( factors.size() != numberOfLevels )
      {
      this->Logger() << "\n\n\n"
                       << "ERROR:  The number of shrink factors does not match the number of levels."
                       << "\nShrink Factors: " << factors.size()
                       << "\nNumber Of Levels: " << numberOfLevels
                       << "\n\n\n"
                       << std::endl;
      return EXIT_FAILURE;
      }

    for( unsigned int n = 0; n < numberOfLevels; n++ )
      {
      ShrinkFactorsPerDimensionContainerType shrinkFactorsPerDimension =
        this->CalculateShrinkFactorsPerDimension( factors[n], virtualDomainImage->GetSpacing() );
      shrinkFactorsPerDimensionForAllLevels.push_back( shrinkFactorsPerDimension );
      this->Logger() << "  Shrink factors (level " << n+1 << " out of " << numberOfLevels << "): " << shrinkFactorsPerDimension << std::endl;
      }

    // Get smoothing sigmas
    const std::vector<float> sigmas( this->m_SmoothingSigmas[currentStageNumber] );
    smoothingSigmasPerLevel.SetSize( sigmas.size() );

    if( sigmas.size() != numberOfLevels )
      {
      this->Logger() << "ERROR:  The number of smoothing sigmas "
                       << "does not match the number of levels." << std::endl;
      return EXIT_FAILURE;
      }
    for( unsigned int n = 0; n < smoothingSigmasPerLevel.Size(); n++ )
      {
      smoothingSigmasPerLevel[n] = sigmas[n];
      }
    this->Logger() << "  smoothing sigmas per level: " << smoothingSigmasPerLevel << std::endl;

    // The sampling strategy/percentage is only specified once for the image registration
    // method.  We might need to change this in the future.

    const float samplingPercentage = stageMetricList[0].m_SamplingPercentage;

    const SamplingStrategy samplingStrategy = stageMetricList[0].m_SamplingStrategy;
    typename AffineRegistrationType::MetricSamplingStrategyType metricSamplingStrategy = AffineRegistrationType::NONE;
    if( samplingStrategy == random )
      {
      this->Logger() << "  random sampling (percentage = " << samplingPercentage << ")" << std::endl;
      metricSamplingStrategy = AffineRegistrationType::RANDOM;
      }
    else if( samplingStrategy == regular )
      {
      this->Logger() << "  regular sampling (percentage = " << samplingPercentage << ")" << std::endl;
      metricSamplingStrategy = AffineRegistrationType::REGULAR;
      }
    else if( samplingStrategy == none )
      {
      this->Logger() << "  Using default NONE metricSamplingStrategy " << std::endl;
      }
    else
      {
      this->Logger() << "ERROR: samplingStrategy is incorrectly specified" << std::endl;
      return EXIT_FAILURE;
      }

    // Set up the optimizers.  To change the iteration number for each level we rely
    // on the command observer.

    const RealType learningRate = this->m_TransformMethods[currentStageNumber].m_GradientStep;

    // There's a scale issue here.  Currently we are using the first metric to estimate the
    // scales but we might need to change this.
    typedef itk::RegistrationParameterScalesFromPhysicalShift<ObjectMetricType> ScalesEstimatorType;

    typename ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();
    scalesEstimator->SetMetric( singleMetric );
    scalesEstimator->SetTransformForward( true );

    typedef itk::LabeledPointSetToPointSetMetricv4<LabeledPointSetType, LabeledPointSetType, RealType> LabeledPointSetMetricType;
    typename LabeledPointSetMetricType::Pointer labeledPointSetMetric2 = dynamic_cast<LabeledPointSetMetricType *>( singleMetric.GetPointer() );
    if( labeledPointSetMetric2.IsNotNull() )
      {
      typedef typename ScalesEstimatorType::VirtualPointSetType VirtualPointSetType;
      typename VirtualPointSetType::Pointer virtualPointSet = VirtualPointSetType::New();
      virtualPointSet->Initialize();
      virtualPointSet->SetPoints(
        const_cast<typename LabeledPointSetType::PointsContainer *>( labeledPointSetMetric2->GetFixedPointSet()->GetPoints() ) );
      scalesEstimator->SetVirtualDomainPointSet( virtualPointSet );
      }
    else
      {
      typedef itk::MeanSquaresPointSetToPointSetIntensityMetricv4<IntensityPointSetType, IntensityPointSetType, RealType> IntensityPointSetMetricType;
      typename IntensityPointSetMetricType::Pointer intensityPointSetMetric2 = dynamic_cast<IntensityPointSetMetricType *>( singleMetric.GetPointer() );
      if( intensityPointSetMetric2.IsNotNull() )
        {
        typedef typename ScalesEstimatorType::VirtualPointSetType VirtualPointSetType;
        typename VirtualPointSetType::Pointer virtualPointSet = VirtualPointSetType::New();
        virtualPointSet->Initialize();
        virtualPointSet->SetPoints(
          const_cast<typename IntensityPointSetType::PointsContainer *>( intensityPointSetMetric2->GetFixedPointSet()->GetPoints() ) );
        scalesEstimator->SetVirtualDomainPointSet( virtualPointSet );
        }
      }

    typename ConjugateGradientDescentOptimizerType::Pointer optimizer = ConjugateGradientDescentOptimizerType::New();
    optimizer->SetLowerLimit( 0 );
    optimizer->SetUpperLimit( 2 );
    optimizer->SetEpsilon( 0.2 );
    //    optimizer->SetMaximumLineSearchIterations( 20 );
    optimizer->SetLearningRate( learningRate );
    optimizer->SetMaximumStepSizeInPhysicalUnits( learningRate );
    optimizer->SetNumberOfIterations( currentStageIterations[0] );
    optimizer->SetScalesEstimator( scalesEstimator );
    optimizer->SetMinimumConvergenceValue( convergenceThreshold );
    optimizer->SetConvergenceWindowSize( convergenceWindowSize );
    optimizer->SetDoEstimateLearningRateAtEachIteration( this->m_DoEstimateLearningRateAtEachIteration );
    optimizer->SetDoEstimateLearningRateOnce( !this->m_DoEstimateLearningRateAtEachIteration );

    typedef antsRegistrationOptimizerCommandIterationUpdate<TComputeType, VImageDimension,
                                                            ConjugateGradientDescentOptimizerType> OptimizerCommandType;
    typename OptimizerCommandType::Pointer optimizerObserver = OptimizerCommandType::New();
    optimizerObserver->SetLogStream( *this->m_LogStream );
    optimizerObserver->SetNumberOfIterations( currentStageIterations );
    optimizerObserver->SetOptimizer( optimizer );

    if( !this->IsPointSetMetric( this->m_Metrics[0].m_MetricType ) )
      {
      optimizerObserver->SetOrigFixedImage( this->m_Metrics[0].m_FixedImage );
      optimizerObserver->SetOrigMovingImage( this->m_Metrics[0].m_MovingImage );
      }

    if( this->m_PrintSimilarityMeasureInterval != 0 )
      {
      optimizerObserver->SetComputeFullScaleCCInterval( this->m_PrintSimilarityMeasureInterval );
      }
    if( this->m_WriteIntervalVolumes != 0 )
      {
      optimizerObserver->SetWriteIterationsOutputsInIntervals( this->m_WriteIntervalVolumes );
      optimizerObserver->SetCurrentStageNumber( currentStageNumber );
      }

    typename GradientDescentOptimizerType::Pointer optimizer2 = GradientDescentOptimizerType::New();
    //    optimizer2->SetLowerLimit( 0 );
    //    optimizer2->SetUpperLimit( 2 );
    //    optimizer2->SetEpsilon( 0.2 );
    //    optimizer->SetMaximumLineSearchIterations( 20 );
    optimizer2->SetLearningRate( learningRate );
    optimizer2->SetMaximumStepSizeInPhysicalUnits( learningRate );
    optimizer2->SetNumberOfIterations( currentStageIterations[0] );
    optimizer2->SetScalesEstimator( nullptr );
    optimizer2->SetMinimumConvergenceValue( convergenceThreshold );
    optimizer2->SetConvergenceWindowSize( convergenceWindowSize );
    optimizer2->SetDoEstimateLearningRateAtEachIteration( this->m_DoEstimateLearningRateAtEachIteration );
    optimizer2->SetDoEstimateLearningRateOnce( !this->m_DoEstimateLearningRateAtEachIteration );

    typedef antsRegistrationOptimizerCommandIterationUpdate<TComputeType, VImageDimension,
                                                            GradientDescentOptimizerType> OptimizerCommandType2;
    typename OptimizerCommandType2::Pointer optimizerObserver2 = OptimizerCommandType2::New();
    optimizerObserver2->SetLogStream( *this->m_LogStream );
    optimizerObserver2->SetNumberOfIterations( currentStageIterations );
    optimizerObserver2->SetOptimizer( optimizer2 );
    if( !this->IsPointSetMetric( this->m_Metrics[0].m_MetricType ) )
      {
      optimizerObserver2->SetOrigFixedImage( this->m_Metrics[0].m_FixedImage );
      optimizerObserver2->SetOrigMovingImage( this->m_Metrics[0].m_MovingImage );
      }

    if( this->m_PrintSimilarityMeasureInterval != 0 )
      {
      optimizerObserver2->SetComputeFullScaleCCInterval( this->m_PrintSimilarityMeasureInterval );
      }
    if( this->m_WriteIntervalVolumes != 0 )
      {
      optimizerObserver2->SetWriteIterationsOutputsInIntervals( this->m_WriteIntervalVolumes );
      optimizerObserver2->SetCurrentStageNumber( currentStageNumber );
      }

    std::vector<typename LabeledPointSetType::Pointer> fixedLabeledPointSetsPerStage;
    std::vector<typename LabeledPointSetType::Pointer> movingLabeledPointSetsPerStage;
    for( unsigned int n = 0; n < stageMetricList.size(); n++ )
      {
      fixedLabeledPointSetsPerStage.push_back( stageMetricList[n].m_FixedLabeledPointSet.GetPointer() );
      movingLabeledPointSetsPerStage.push_back( stageMetricList[n].m_MovingLabeledPointSet.GetPointer() );
      }

    std::vector<typename IntensityPointSetType::Pointer> fixedIntensityPointSetsPerStage;
    std::vector<typename IntensityPointSetType::Pointer> movingIntensityPointSetsPerStage;
    for( unsigned int n = 0; n < stageMetricList.size(); n++ )
      {
      fixedIntensityPointSetsPerStage.push_back( stageMetricList[n].m_FixedIntensityPointSet.GetPointer() );
      movingIntensityPointSetsPerStage.push_back( stageMetricList[n].m_MovingIntensityPointSet.GetPointer() );
      }

    // Set up the image registration methods along with the transforms
    const XfrmMethod whichTransform( this->m_TransformMethods[currentStageNumber].m_XfrmMethod );

    switch( whichTransform )
      {
      case Affine:
        {
        if( stageMetricList[0].m_MetricType != IGDM )
          {
          this->AddLinearTransformToCompositeTransform<AffineRegistrationType>(
            this->m_CompositeTransform, currentStageNumber, AffineTransformType::ParametersDimension,
            preprocessedFixedImagesPerStage, preprocessedMovingImagesPerStage,
            fixedLabeledPointSetsPerStage, movingLabeledPointSetsPerStage, stageMetricList, singleMetric,
            multiMetric, optimizer, numberOfLevels, shrinkFactorsPerDimensionForAllLevels,
            smoothingSigmasPerLevel, metricSamplingStrategy, samplingPercentage );
          }
        else
          {
          typedef itk::ImageRegistrationMethodv4<ImageType, ImageType,
            AffineTransformType, ImageType, IntensityPointSetType>  AffineRegistrationType2;

          this->AddLinearTransformToCompositeTransform<AffineRegistrationType2>(
            this->m_CompositeTransform, currentStageNumber, AffineTransformType::ParametersDimension,
            preprocessedFixedImagesPerStage, preprocessedMovingImagesPerStage,
            fixedIntensityPointSetsPerStage, movingIntensityPointSetsPerStage, stageMetricList, singleMetric,
            multiMetric, optimizer, numberOfLevels, shrinkFactorsPerDimensionForAllLevels,
            smoothingSigmasPerLevel, metricSamplingStrategy, samplingPercentage );
          }
        }
        break;
      case Rigid:
        {
        typedef typename RigidTransformTraits<TComputeType, VImageDimension>::TransformType RigidTransformType;

        if( stageMetricList[0].m_MetricType != IGDM )
          {
          typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, RigidTransformType,
            ImageType, LabeledPointSetType> RigidRegistrationType;

          this->AddLinearTransformToCompositeTransform<RigidRegistrationType>(
            this->m_CompositeTransform, currentStageNumber, RigidTransformType::ParametersDimension,
            preprocessedFixedImagesPerStage, preprocessedMovingImagesPerStage,
            fixedLabeledPointSetsPerStage, movingLabeledPointSetsPerStage, stageMetricList, singleMetric,
            multiMetric, optimizer, numberOfLevels, shrinkFactorsPerDimensionForAllLevels,
            smoothingSigmasPerLevel, metricSamplingStrategy, samplingPercentage );
          }
        else
          {
          typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, RigidTransformType,
            ImageType, IntensityPointSetType> RigidRegistrationType;

          this->AddLinearTransformToCompositeTransform<RigidRegistrationType>(
            this->m_CompositeTransform, currentStageNumber, RigidTransformType::ParametersDimension,
            preprocessedFixedImagesPerStage, preprocessedMovingImagesPerStage,
            fixedIntensityPointSetsPerStage, movingIntensityPointSetsPerStage, stageMetricList, singleMetric,
            multiMetric, optimizer, numberOfLevels, shrinkFactorsPerDimensionForAllLevels,
            smoothingSigmasPerLevel, metricSamplingStrategy, samplingPercentage );
          }
        }
        break;
      case CompositeAffine:
        {
        typedef typename CompositeAffineTransformTraits<TComputeType, VImageDimension>::TransformType CompositeAffineTransformType;

        if( stageMetricList[0].m_MetricType != IGDM )
          {
          typedef itk::ImageRegistrationMethodv4<ImageType, ImageType,
            CompositeAffineTransformType, ImageType, LabeledPointSetType> CompositeAffineRegistrationType;

          this->AddLinearTransformToCompositeTransform<CompositeAffineRegistrationType>(
            this->m_CompositeTransform, currentStageNumber, CompositeAffineTransformType::ParametersDimension,
            preprocessedFixedImagesPerStage, preprocessedMovingImagesPerStage,
            fixedLabeledPointSetsPerStage, movingLabeledPointSetsPerStage, stageMetricList, singleMetric,
            multiMetric, optimizer, numberOfLevels, shrinkFactorsPerDimensionForAllLevels,
            smoothingSigmasPerLevel, metricSamplingStrategy, samplingPercentage );
          }
        else
          {
          typedef itk::ImageRegistrationMethodv4<ImageType, ImageType,
            CompositeAffineTransformType, ImageType, IntensityPointSetType> CompositeAffineRegistrationType;

          this->AddLinearTransformToCompositeTransform<CompositeAffineRegistrationType>(
            this->m_CompositeTransform, currentStageNumber, CompositeAffineTransformType::ParametersDimension,
            preprocessedFixedImagesPerStage, preprocessedMovingImagesPerStage,
            fixedIntensityPointSetsPerStage, movingIntensityPointSetsPerStage, stageMetricList, singleMetric,
            multiMetric, optimizer, numberOfLevels, shrinkFactorsPerDimensionForAllLevels,
            smoothingSigmasPerLevel, metricSamplingStrategy, samplingPercentage );
          }
        }
        break;
      case Similarity:
        {
        typedef typename SimilarityTransformTraits<TComputeType, VImageDimension>::TransformType SimilarityTransformType;

        if( stageMetricList[0].m_MetricType != IGDM )
          {
          typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, SimilarityTransformType,
            ImageType, LabeledPointSetType> SimilarityRegistrationType;

          this->AddLinearTransformToCompositeTransform<SimilarityRegistrationType>(
            this->m_CompositeTransform, currentStageNumber, SimilarityTransformType::ParametersDimension,
            preprocessedFixedImagesPerStage, preprocessedMovingImagesPerStage,
            fixedLabeledPointSetsPerStage, movingLabeledPointSetsPerStage, stageMetricList, singleMetric,
            multiMetric, optimizer, numberOfLevels, shrinkFactorsPerDimensionForAllLevels,
            smoothingSigmasPerLevel, metricSamplingStrategy, samplingPercentage );
          }
        else
          {
          typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, SimilarityTransformType,
            ImageType, IntensityPointSetType> SimilarityRegistrationType;

          this->AddLinearTransformToCompositeTransform<SimilarityRegistrationType>(
            this->m_CompositeTransform, currentStageNumber, SimilarityTransformType::ParametersDimension,
            preprocessedFixedImagesPerStage, preprocessedMovingImagesPerStage,
            fixedIntensityPointSetsPerStage, movingIntensityPointSetsPerStage, stageMetricList, singleMetric,
            multiMetric, optimizer, numberOfLevels, shrinkFactorsPerDimensionForAllLevels,
            smoothingSigmasPerLevel, metricSamplingStrategy, samplingPercentage );
          }
        }
        break;
      case Translation:
        {
        typedef itk::TranslationTransform<RealType, VImageDimension> TranslationTransformType;

        if( stageMetricList[0].m_MetricType != IGDM )
          {
          typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, TranslationTransformType,
            ImageType, LabeledPointSetType> TranslationRegistrationType;

          this->AddLinearTransformToCompositeTransform<TranslationRegistrationType>(
            this->m_CompositeTransform, currentStageNumber, TranslationTransformType::ParametersDimension,
            preprocessedFixedImagesPerStage, preprocessedMovingImagesPerStage,
            fixedLabeledPointSetsPerStage, movingLabeledPointSetsPerStage, stageMetricList, singleMetric,
            multiMetric, optimizer, numberOfLevels, shrinkFactorsPerDimensionForAllLevels,
            smoothingSigmasPerLevel, metricSamplingStrategy, samplingPercentage );
          }
        else
          {
          typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, TranslationTransformType,
            ImageType, IntensityPointSetType> TranslationRegistrationType;

          this->AddLinearTransformToCompositeTransform<TranslationRegistrationType>(
            this->m_CompositeTransform, currentStageNumber, TranslationTransformType::ParametersDimension,
            preprocessedFixedImagesPerStage, preprocessedMovingImagesPerStage,
            fixedIntensityPointSetsPerStage, movingIntensityPointSetsPerStage, stageMetricList, singleMetric,
            multiMetric, optimizer, numberOfLevels, shrinkFactorsPerDimensionForAllLevels,
            smoothingSigmasPerLevel, metricSamplingStrategy, samplingPercentage );
          }
        }
        break;
      case GaussianDisplacementField:
        {
        if( stageMetricList[0].m_MetricType == IGDM )
          {
          this->Logger() << "Intensity point set metric is not implemented yet for the specified transform." << std::endl;
          return EXIT_FAILURE;
          }

        typedef itk::GaussianSmoothingOnUpdateDisplacementFieldTransform<RealType, VImageDimension>
          GaussianDisplacementFieldTransformType;
        typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, GaussianDisplacementFieldTransformType,
          ImageType, LabeledPointSetType> DisplacementFieldRegistrationType;

        typename DisplacementFieldRegistrationType::Pointer registrationMethod =
          this->PrepareRegistrationMethod<DisplacementFieldRegistrationType>(
                this->m_CompositeTransform, currentStageNumber, VImageDimension,
                preprocessedFixedImagesPerStage, preprocessedMovingImagesPerStage,
                fixedLabeledPointSetsPerStage, movingLabeledPointSetsPerStage, stageMetricList, singleMetric,
                multiMetric, optimizer, numberOfLevels, shrinkFactorsPerDimensionForAllLevels,
                smoothingSigmasPerLevel, metricSamplingStrategy, samplingPercentage );

        typename GaussianDisplacementFieldTransformType::Pointer outputDisplacementFieldTransform =
          registrationMethod->GetModifiableTransform();

        typedef itk::Vector<RealType, VImageDimension> VectorType;
        VectorType zeroVector( 0.0 );
        typename DisplacementFieldType::Pointer displacementField = AllocImage<DisplacementFieldType>(
            preprocessedFixedImagesPerStage[0], zeroVector );
        outputDisplacementFieldTransform->SetDisplacementField( displacementField );

        // Create the transform adaptors

        typedef itk::GaussianSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor<
          GaussianDisplacementFieldTransformType> DisplacementFieldTransformAdaptorType;
        typename DisplacementFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;

        RealType varianceForUpdateField =
          this->m_TransformMethods[currentStageNumber].m_UpdateFieldVarianceInVarianceSpace;
        RealType varianceForTotalField  =
          this->m_TransformMethods[currentStageNumber].m_TotalFieldVarianceInVarianceSpace;

        outputDisplacementFieldTransform->SetGaussianSmoothingVarianceForTheUpdateField( varianceForUpdateField );
        outputDisplacementFieldTransform->SetGaussianSmoothingVarianceForTheTotalField( varianceForTotalField );
        // Create the transform adaptors
        // For the gaussian displacement field, the specified variances are in image spacing terms
        // and, in normal practice, we typically don't change these values at each level.  However,
        // if the user wishes to add that option, they can use the class
        // GaussianSmoothingOnUpdateDisplacementFieldTransformAdaptor
        for( unsigned int level = 0; level < numberOfLevels; level++ )
          {
          typename itk::ImageBase<VImageDimension>::Pointer shrunkSpace=
                     this->GetShrinkImageOutputInformation(
                          virtualDomainImage.GetPointer(),
                          shrinkFactorsPerDimensionForAllLevels[level]  );

          typename DisplacementFieldTransformAdaptorType::Pointer fieldTransformAdaptor =
            DisplacementFieldTransformAdaptorType::New();
          fieldTransformAdaptor->SetRequiredSpacing( shrunkSpace->GetSpacing() );
          fieldTransformAdaptor->SetRequiredSize( shrunkSpace->GetLargestPossibleRegion().GetSize() );
          fieldTransformAdaptor->SetRequiredDirection( shrunkSpace->GetDirection() );
          fieldTransformAdaptor->SetRequiredOrigin( shrunkSpace->GetOrigin() );
          fieldTransformAdaptor->SetTransform( outputDisplacementFieldTransform );

          adaptors.push_back( fieldTransformAdaptor.GetPointer() );
          }
        registrationMethod->SetOptimizer( optimizer );
        registrationMethod->SetTransformParametersAdaptorsPerLevel( adaptors );

        typedef antsRegistrationCommandIterationUpdate<DisplacementFieldRegistrationType> DisplacementFieldCommandType;
        typename DisplacementFieldCommandType::Pointer displacementFieldRegistrationObserver =
          DisplacementFieldCommandType::New();
        displacementFieldRegistrationObserver->SetLogStream( *this->m_LogStream );
        displacementFieldRegistrationObserver->SetNumberOfIterations( currentStageIterations );

        registrationMethod->AddObserver( itk::IterationEvent(), displacementFieldRegistrationObserver );
        registrationMethod->AddObserver( itk::InitializeEvent(), displacementFieldRegistrationObserver );

        try
          {
          this->Logger() << std::endl
                         << "*** Running gaussian displacement field registration (varianceForUpdateField = "
                         << varianceForUpdateField << ", varianceForTotalField = " << varianceForTotalField << ") ***"
                         << std::endl << std::endl;
          displacementFieldRegistrationObserver->Execute( registrationMethod, itk::StartEvent() );
          registrationMethod->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          this->Logger() << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }

        // Add calculated transform to the composite transform
        this->m_CompositeTransform->AddTransform( outputDisplacementFieldTransform );

        this->m_AllPreviousTransformsAreLinear = false;
        }
        break;
      case BSplineDisplacementField:
        {
        if( stageMetricList[0].m_MetricType == IGDM )
          {
          this->Logger() << "Intensity point set metric is not implemented yet for the specified transform." << std::endl;
          return EXIT_FAILURE;
          }

        typedef itk::BSplineSmoothingOnUpdateDisplacementFieldTransform<RealType, VImageDimension>
          BSplineDisplacementFieldTransformType;
        typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, BSplineDisplacementFieldTransformType,
          ImageType, LabeledPointSetType> DisplacementFieldRegistrationType;

        typename DisplacementFieldRegistrationType::Pointer registrationMethod =
          this->PrepareRegistrationMethod<DisplacementFieldRegistrationType>(
                this->m_CompositeTransform, currentStageNumber, VImageDimension,
                preprocessedFixedImagesPerStage, preprocessedMovingImagesPerStage,
                fixedLabeledPointSetsPerStage, movingLabeledPointSetsPerStage, stageMetricList, singleMetric,
                multiMetric, optimizer, numberOfLevels, shrinkFactorsPerDimensionForAllLevels,
                smoothingSigmasPerLevel, metricSamplingStrategy, samplingPercentage );

        typename BSplineDisplacementFieldTransformType::Pointer outputDisplacementFieldTransform =
          registrationMethod->GetModifiableTransform();

        typedef itk::Vector<RealType, VImageDimension> VectorType;
        VectorType zeroVector( 0.0 );
        typename DisplacementFieldType::Pointer displacementField = AllocImage<DisplacementFieldType>(
            preprocessedFixedImagesPerStage[0], zeroVector );
        outputDisplacementFieldTransform->SetDisplacementField( displacementField );

        // Create the transform adaptors

        typename DisplacementFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;

        const std::vector<unsigned int> & meshSizeForTheUpdateField =
          this->m_TransformMethods[currentStageNumber].m_UpdateFieldMeshSizeAtBaseLevel;
        std::vector<unsigned int> meshSizeForTheTotalField =
          this->m_TransformMethods[currentStageNumber].m_TotalFieldMeshSizeAtBaseLevel;

        outputDisplacementFieldTransform->SetSplineOrder( this->m_TransformMethods[currentStageNumber].m_SplineOrder );

        if( meshSizeForTheUpdateField.size() != VImageDimension || meshSizeForTheTotalField.size() != VImageDimension )
          {
          this->Logger() << "ERROR:  The mesh size(s) don't match the ImageDimension." << std::endl;
          return EXIT_FAILURE;
          }

        typename BSplineDisplacementFieldTransformType::ArrayType updateMeshSize;
        typename BSplineDisplacementFieldTransformType::ArrayType totalMeshSize;
        for( unsigned int d = 0; d < VImageDimension; d++ )
          {
          updateMeshSize[d] = meshSizeForTheUpdateField[d];
          totalMeshSize[d] = meshSizeForTheTotalField[d];
          }
        // Create the transform adaptors specific to B-splines
        for( unsigned int level = 0; level < numberOfLevels; level++ )
          {
          typename itk::ImageBase<VImageDimension>::Pointer shrunkSpace=
                     this->GetShrinkImageOutputInformation(
                          virtualDomainImage.GetPointer(),
                          shrinkFactorsPerDimensionForAllLevels[level]  );

          typedef itk::BSplineSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor<
              BSplineDisplacementFieldTransformType> BSplineDisplacementFieldTransformAdaptorType;
          typename BSplineDisplacementFieldTransformAdaptorType::Pointer bsplineFieldTransformAdaptor =
            BSplineDisplacementFieldTransformAdaptorType::New();
          bsplineFieldTransformAdaptor->SetRequiredSpacing( shrunkSpace->GetSpacing() );
          bsplineFieldTransformAdaptor->SetRequiredSize( shrunkSpace->GetLargestPossibleRegion().GetSize() );
          bsplineFieldTransformAdaptor->SetRequiredDirection( shrunkSpace->GetDirection() );
          bsplineFieldTransformAdaptor->SetRequiredOrigin( shrunkSpace->GetOrigin() );
          bsplineFieldTransformAdaptor->SetTransform( outputDisplacementFieldTransform );

          // A good heuristic is to double the b-spline mesh resolution at each level
          typename BSplineDisplacementFieldTransformType::ArrayType newUpdateMeshSize = updateMeshSize;
          typename BSplineDisplacementFieldTransformType::ArrayType newTotalMeshSize = totalMeshSize;
          for( unsigned int d = 0; d < VImageDimension; d++ )
            {
            newUpdateMeshSize[d] = newUpdateMeshSize[d] << ( level + 1 );
            newTotalMeshSize[d] = newTotalMeshSize[d] << ( level + 1 );
            }
          bsplineFieldTransformAdaptor->SetMeshSizeForTheUpdateField( newUpdateMeshSize );
          bsplineFieldTransformAdaptor->SetMeshSizeForTheTotalField( newTotalMeshSize );

          adaptors.push_back( bsplineFieldTransformAdaptor.GetPointer() );
          }

        registrationMethod->SetOptimizer( optimizer );
        registrationMethod->SetTransformParametersAdaptorsPerLevel( adaptors );

        typedef antsRegistrationCommandIterationUpdate<DisplacementFieldRegistrationType> DisplacementFieldCommandType;
        typename DisplacementFieldCommandType::Pointer displacementFieldRegistrationObserver =
          DisplacementFieldCommandType::New();
        displacementFieldRegistrationObserver->SetLogStream( *this->m_LogStream );
        displacementFieldRegistrationObserver->SetNumberOfIterations( currentStageIterations );

        registrationMethod->AddObserver( itk::IterationEvent(), displacementFieldRegistrationObserver );
        registrationMethod->AddObserver( itk::InitializeEvent(), displacementFieldRegistrationObserver );

        try
          {
          this->Logger() << std::endl
                         << "*** Running bspline displacement field registration (updateMeshSizeAtBaseLevel = "
                         << updateMeshSize << ", totalMeshSizeAtBaseLevel = " << totalMeshSize << ") ***" << std::endl
                         << std::endl;
          displacementFieldRegistrationObserver->Execute( registrationMethod, itk::StartEvent() );
          registrationMethod->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          this->Logger() << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }

        // Add calculated transform to the composite transform
        this->m_CompositeTransform->AddTransform( outputDisplacementFieldTransform );

        this->m_AllPreviousTransformsAreLinear = false;
        }
        break;
      case SyN:
        {
        if( stageMetricList[0].m_MetricType == IGDM )
          {
          this->Logger() << "Intensity point set metric is not implemented yet for the specified transform." << std::endl;
          return EXIT_FAILURE;
          }

        typedef itk::Vector<RealType, VImageDimension> VectorType;
        VectorType zeroVector( 0.0 );
        //typedef itk::Image<VectorType, VImageDimension> DisplacementFieldType;

        typename DisplacementFieldType::Pointer displacementField =
          AllocImage<DisplacementFieldType>( virtualDomainImage, zeroVector );
        typename DisplacementFieldType::Pointer inverseDisplacementField =
          AllocImage<DisplacementFieldType>( virtualDomainImage, zeroVector );

        typedef itk::SyNImageRegistrationMethod<ImageType, ImageType,
          DisplacementFieldTransformType, ImageType, LabeledPointSetType> DisplacementFieldRegistrationType;
        typename DisplacementFieldRegistrationType::Pointer displacementFieldRegistration =
          DisplacementFieldRegistrationType::New();

	if ( this->m_RegistrationRandomSeed != 0 )
	  {
	  displacementFieldRegistration->MetricSamplingReinitializeSeed( this->m_RegistrationRandomSeed );
	  }

        if( this->m_RestrictDeformationOptimizerWeights.size() > currentStageNumber )
          {
          if( this->m_RestrictDeformationOptimizerWeights[currentStageNumber].size() == VImageDimension )
            {
            typename DisplacementFieldRegistrationType::OptimizerWeightsType optimizerWeights( VImageDimension );
            for( unsigned int d = 0; d < VImageDimension; d++ )
              {
              optimizerWeights[d] = this->m_RestrictDeformationOptimizerWeights[currentStageNumber][d];
              }
            displacementFieldRegistration->SetOptimizerWeights( optimizerWeights );
            }
          }

        typename DisplacementFieldTransformType::Pointer outputDisplacementFieldTransform = displacementFieldRegistration->GetModifiableTransform();

        // Create the transform adaptors

        typedef itk::DisplacementFieldTransformParametersAdaptor<DisplacementFieldTransformType>
          DisplacementFieldTransformAdaptorType;
        typename DisplacementFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;
        // Create the transform adaptors
        // For the gaussian displacement field, the specified variances are in image spacing terms
        // and, in normal practice, we typically don't change these values at each level.  However,
        // if the user wishes to add that option, they can use the class
        // GaussianSmoothingOnUpdateDisplacementFieldTransformAdaptor
        for( unsigned int level = 0; level < numberOfLevels; level++ )
          {
          typename itk::ImageBase<VImageDimension>::Pointer shrunkSpace=
                     this->GetShrinkImageOutputInformation(
                          virtualDomainImage.GetPointer(),
                          shrinkFactorsPerDimensionForAllLevels[level]  );

          typename DisplacementFieldTransformAdaptorType::Pointer fieldTransformAdaptor =
            DisplacementFieldTransformAdaptorType::New();
          fieldTransformAdaptor->SetRequiredSpacing( shrunkSpace->GetSpacing() );
          fieldTransformAdaptor->SetRequiredSize( shrunkSpace->GetLargestPossibleRegion().GetSize() );
          fieldTransformAdaptor->SetRequiredDirection( shrunkSpace->GetDirection() );
          fieldTransformAdaptor->SetRequiredOrigin( shrunkSpace->GetOrigin() );
          fieldTransformAdaptor->SetTransform( outputDisplacementFieldTransform );

          adaptors.push_back( fieldTransformAdaptor.GetPointer() );
          }

        // Extract parameters
        typename DisplacementFieldRegistrationType::NumberOfIterationsArrayType numberOfIterationsPerLevel;
        numberOfIterationsPerLevel.SetSize( numberOfLevels );
        for( unsigned int d = 0; d < numberOfLevels; d++ )
          {
          numberOfIterationsPerLevel[d] = currentStageIterations[d];
          }

        const RealType varianceForUpdateField =
          this->m_TransformMethods[currentStageNumber].m_UpdateFieldVarianceInVarianceSpace;
        const RealType varianceForTotalField =
          this->m_TransformMethods[currentStageNumber].m_TotalFieldVarianceInVarianceSpace;
        for( unsigned int n = 0; n < stageMetricList.size(); n++ )
          {
          if( !this->IsPointSetMetric( stageMetricList[n].m_MetricType ) )
            {
            displacementFieldRegistration->SetFixedImage( n, preprocessedFixedImagesPerStage[n] );
            displacementFieldRegistration->SetMovingImage( n, preprocessedMovingImagesPerStage[n] );
            }
          else
            {
            displacementFieldRegistration->SetFixedPointSet( n, stageMetricList[n].m_FixedLabeledPointSet.GetPointer() );
            displacementFieldRegistration->SetMovingPointSet( n, stageMetricList[n].m_MovingLabeledPointSet.GetPointer() );
            }
          }
        if( useMultiMetric )
          {
          displacementFieldRegistration->SetMetric( multiMetric );
          }
        else
          {
          displacementFieldRegistration->SetMetric( singleMetric );
          }

        bool synIsInitialized = false;
        if( this->m_InitializeTransformsPerStage )
          {
          if( this->m_RegistrationState.IsNotNull() )
            {
            const unsigned int numOfTransforms = this->m_RegistrationState->GetNumberOfTransforms();
            typename TransformType::Pointer oneToEndTransform = this->m_RegistrationState->GetNthTransform( numOfTransforms-2 );
            typename TransformType::Pointer endTransform = this->m_RegistrationState->GetNthTransform( numOfTransforms-1 );

            typename DisplacementFieldTransformType::Pointer fixedToMiddle =
              dynamic_cast<DisplacementFieldTransformType *>( oneToEndTransform.GetPointer() );
            typename DisplacementFieldTransformType::Pointer movingToMiddle =
              dynamic_cast<DisplacementFieldTransformType *>( endTransform.GetPointer() );

            if( fixedToMiddle.IsNotNull() && movingToMiddle.IsNotNull()
               && fixedToMiddle->GetInverseDisplacementField() && movingToMiddle->GetInverseDisplacementField() )
              {
              this->Logger() << "Current SyN transform is directly initialized from the previous stage." << std::endl;
              displacementFieldRegistration->SetFixedToMiddleTransform( fixedToMiddle );
              displacementFieldRegistration->SetMovingToMiddleTransform( movingToMiddle );

              this->m_RegistrationState->RemoveTransform();
              this->m_RegistrationState->RemoveTransform();
              }

            // If there are components other than SyN state
            if( this->m_RegistrationState->GetNumberOfTransforms() > 0 )
              {
              displacementFieldRegistration->SetMovingInitialTransform( this->m_RegistrationState );
              }
            synIsInitialized = true;
            this->m_CompositeTransform->RemoveTransform();
            }
          }

        if( this->m_CompositeTransform->GetNumberOfTransforms() > 0 && !synIsInitialized )
          {
          displacementFieldRegistration->SetMovingInitialTransform( this->m_CompositeTransform );
          }

        if( this->m_FixedInitialTransform->GetNumberOfTransforms() > 0 )
          {
          displacementFieldRegistration->SetFixedInitialTransform( this->m_FixedInitialTransform );
          }
        displacementFieldRegistration->SetDownsampleImagesForMetricDerivatives( true );
        displacementFieldRegistration->SetAverageMidPointGradients( false );

        displacementFieldRegistration->SetNumberOfLevels( numberOfLevels );

        for( unsigned int level = 0; level < numberOfLevels; ++level )
          {
          displacementFieldRegistration->SetShrinkFactorsPerDimension( level, shrinkFactorsPerDimensionForAllLevels[level] );
          }
        displacementFieldRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
        displacementFieldRegistration->SetSmoothingSigmasAreSpecifiedInPhysicalUnits(
          this->m_SmoothingSigmasAreInPhysicalUnits[currentStageNumber] );

        displacementFieldRegistration->SetLearningRate( learningRate );
        displacementFieldRegistration->SetConvergenceThreshold( convergenceThreshold );
        displacementFieldRegistration->SetConvergenceWindowSize( convergenceWindowSize );
        displacementFieldRegistration->SetNumberOfIterationsPerLevel( numberOfIterationsPerLevel );
        displacementFieldRegistration->SetTransformParametersAdaptorsPerLevel( adaptors );
        displacementFieldRegistration->SetGaussianSmoothingVarianceForTheUpdateField( varianceForUpdateField );
        displacementFieldRegistration->SetGaussianSmoothingVarianceForTheTotalField( varianceForTotalField );
        outputDisplacementFieldTransform->SetDisplacementField( displacementField );
        outputDisplacementFieldTransform->SetInverseDisplacementField( inverseDisplacementField );

        // For all Velocity field and Displacement field registration types that are not using generic
        // itkImageRegistrationMethodv4 we use following type of observer:
        typedef antsDisplacementAndVelocityFieldRegistrationCommandIterationUpdate<DisplacementFieldRegistrationType>
          DisplacementFieldCommandType2;
        typename DisplacementFieldCommandType2::Pointer displacementFieldRegistrationObserver2 =
          DisplacementFieldCommandType2::New();
        displacementFieldRegistrationObserver2->SetLogStream(*this->m_LogStream);
        displacementFieldRegistrationObserver2->SetNumberOfIterations( currentStageIterations );
        displacementFieldRegistrationObserver2->SetOrigFixedImage( this->m_Metrics[0].m_FixedImage );
        displacementFieldRegistrationObserver2->SetOrigMovingImage( this->m_Metrics[0].m_MovingImage );
        if( this->m_PrintSimilarityMeasureInterval != 0 )
          {
          displacementFieldRegistrationObserver2->SetComputeFullScaleCCInterval( this->m_PrintSimilarityMeasureInterval );
          }
        if( this->m_WriteIntervalVolumes != 0 )
          {
          displacementFieldRegistrationObserver2->SetWriteIterationsOutputsInIntervals( this->m_WriteIntervalVolumes );
          displacementFieldRegistrationObserver2->SetCurrentStageNumber( currentStageNumber );
          }
        displacementFieldRegistration->AddObserver( itk::InitializeEvent(), displacementFieldRegistrationObserver2 );
        displacementFieldRegistration->AddObserver( itk::IterationEvent(), displacementFieldRegistrationObserver2 );

        try
          {
          this->Logger() << std::endl << "*** Running SyN registration (varianceForUpdateField = "
                         << varianceForUpdateField << ", varianceForTotalField = " << varianceForTotalField << ") ***"
                         << std::endl << std::endl;
          displacementFieldRegistrationObserver2->Execute( displacementFieldRegistration, itk::StartEvent() );
          displacementFieldRegistration->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          this->Logger() << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }

        // Add calculated internal transforms to the registration state
        if( this->m_RegistrationState.IsNull() )
          {
          this->m_RegistrationState = CompositeTransformType::New();
          }
        this->m_RegistrationState->ClearTransformQueue();
        this->m_RegistrationState->AddTransform( this->m_CompositeTransform );
        this->m_RegistrationState->AddTransform( displacementFieldRegistration->GetModifiableFixedToMiddleTransform() );
        this->m_RegistrationState->AddTransform( displacementFieldRegistration->GetModifiableMovingToMiddleTransform() );
        this->m_RegistrationState->FlattenTransformQueue();

        // Add calculated transform to the composite transform
        this->m_CompositeTransform->AddTransform( outputDisplacementFieldTransform );
        this->m_AllPreviousTransformsAreLinear = false;
        }
        break;
      case BSplineSyN:
        {
        typedef itk::BSplineSmoothingOnUpdateDisplacementFieldTransform<RealType, VImageDimension>
          BSplineDisplacementFieldTransformType;

        typedef itk::Vector<RealType, VImageDimension> VectorType;
        VectorType zeroVector( 0.0 );

        typename DisplacementFieldType::Pointer displacementField =
          AllocImage<DisplacementFieldType>( virtualDomainImage, zeroVector );
        typename DisplacementFieldType::Pointer inverseDisplacementField =
          AllocImage<DisplacementFieldType>( virtualDomainImage, zeroVector );

        const std::vector<unsigned int> & meshSizeForTheUpdateField =
          this->m_TransformMethods[currentStageNumber].m_UpdateFieldMeshSizeAtBaseLevel;
        std::vector<unsigned int> meshSizeForTheTotalField =
          this->m_TransformMethods[currentStageNumber].m_TotalFieldMeshSizeAtBaseLevel;

        if( meshSizeForTheUpdateField.size() != VImageDimension || meshSizeForTheTotalField.size() != VImageDimension )
          {
          this->Logger() << "ERROR:  The mesh size(s) don't match the ImageDimension." << std::endl;
          return EXIT_FAILURE;
          }

        typename BSplineDisplacementFieldTransformType::ArrayType updateMeshSize;
        typename BSplineDisplacementFieldTransformType::ArrayType totalMeshSize;
        for( unsigned int d = 0; d < VImageDimension; d++ )
          {
          updateMeshSize[d] = meshSizeForTheUpdateField[d];
          totalMeshSize[d] = meshSizeForTheTotalField[d];
          }

        if( stageMetricList[0].m_MetricType != IGDM )
          {
          typedef itk::BSplineSyNImageRegistrationMethod<ImageType, ImageType,
            BSplineDisplacementFieldTransformType, ImageType, LabeledPointSetType>
            DisplacementFieldRegistrationType;

          typename DisplacementFieldRegistrationType::Pointer registrationMethod =
            this->PrepareRegistrationMethod<DisplacementFieldRegistrationType>(
                  this->m_CompositeTransform, currentStageNumber, VImageDimension,
                  preprocessedFixedImagesPerStage, preprocessedMovingImagesPerStage,
                  fixedLabeledPointSetsPerStage, movingLabeledPointSetsPerStage, stageMetricList, singleMetric,
                  multiMetric, optimizer, numberOfLevels, shrinkFactorsPerDimensionForAllLevels,
                  smoothingSigmasPerLevel, metricSamplingStrategy, samplingPercentage );

          typename BSplineDisplacementFieldTransformType::Pointer outputDisplacementFieldTransform =
            registrationMethod->GetModifiableTransform();

          // Create the transform adaptors

          typename DisplacementFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;

          outputDisplacementFieldTransform->SetSplineOrder( this->m_TransformMethods[currentStageNumber].m_SplineOrder );

          // Create the transform adaptors
          for( unsigned int level = 0; level < numberOfLevels; level++ )
            {
            typename itk::ImageBase<VImageDimension>::Pointer shrunkSpace=
                     this->GetShrinkImageOutputInformation(
                          virtualDomainImage.GetPointer(),
                          shrinkFactorsPerDimensionForAllLevels[level]  );

            typedef itk::BSplineSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor<
                BSplineDisplacementFieldTransformType>
              BSplineDisplacementFieldTransformAdaptorType;
            typename BSplineDisplacementFieldTransformAdaptorType::Pointer bsplineFieldTransformAdaptor =
              BSplineDisplacementFieldTransformAdaptorType::New();
            bsplineFieldTransformAdaptor->SetRequiredSpacing( shrunkSpace->GetSpacing() );
            bsplineFieldTransformAdaptor->SetRequiredSize( shrunkSpace->GetLargestPossibleRegion().GetSize() );
            bsplineFieldTransformAdaptor->SetRequiredDirection( shrunkSpace->GetDirection() );
            bsplineFieldTransformAdaptor->SetRequiredOrigin( shrunkSpace->GetOrigin() );
            bsplineFieldTransformAdaptor->SetTransform( outputDisplacementFieldTransform );

            // A good heuristic is to RealType the b-spline mesh resolution at each level
            typename BSplineDisplacementFieldTransformType::ArrayType newUpdateMeshSize = updateMeshSize;
            typename BSplineDisplacementFieldTransformType::ArrayType newTotalMeshSize = totalMeshSize;
            for( unsigned int d = 0; d < VImageDimension; d++ )
              {
              newUpdateMeshSize[d] = newUpdateMeshSize[d] << ( level );
              newTotalMeshSize[d] = newTotalMeshSize[d] << ( level );
              }

            bsplineFieldTransformAdaptor->SetMeshSizeForTheUpdateField( newUpdateMeshSize );
            bsplineFieldTransformAdaptor->SetMeshSizeForTheTotalField( newTotalMeshSize );

            adaptors.push_back( bsplineFieldTransformAdaptor.GetPointer() );
            }

          registrationMethod->SetDownsampleImagesForMetricDerivatives( true );
          registrationMethod->SetAverageMidPointGradients( false );
          registrationMethod->SetNumberOfLevels( numberOfLevels );

          typename DisplacementFieldRegistrationType::NumberOfIterationsArrayType numberOfIterationsPerLevel;
          numberOfIterationsPerLevel.SetSize( numberOfLevels );
          for( unsigned int d = 0; d < numberOfLevels; d++ )
            {
            numberOfIterationsPerLevel[d] = currentStageIterations[d];
            }

          registrationMethod->SetLearningRate( learningRate );
          registrationMethod->SetConvergenceThreshold( convergenceThreshold );
          registrationMethod->SetConvergenceWindowSize( convergenceWindowSize );
          registrationMethod->SetNumberOfIterationsPerLevel( numberOfIterationsPerLevel );
          registrationMethod->SetTransformParametersAdaptorsPerLevel( adaptors );
          outputDisplacementFieldTransform->SetDisplacementField( displacementField );
          outputDisplacementFieldTransform->SetInverseDisplacementField( inverseDisplacementField );

          typedef antsRegistrationCommandIterationUpdate<DisplacementFieldRegistrationType> DisplacementFieldCommandType;
          typename DisplacementFieldCommandType::Pointer displacementFieldRegistrationObserver =
            DisplacementFieldCommandType::New();
          displacementFieldRegistrationObserver->SetLogStream( *this->m_LogStream );
          displacementFieldRegistrationObserver->SetNumberOfIterations( currentStageIterations );

          registrationMethod->AddObserver( itk::IterationEvent(), displacementFieldRegistrationObserver );
          registrationMethod->AddObserver( itk::InitializeEvent(), displacementFieldRegistrationObserver );

          try
            {
            this->Logger() << std::endl << "*** Running B-spline SyN registration (updateMeshSizeAtBaseLevel = "
                           << updateMeshSize << ", totalMeshSizeAtBaseLevel = " << totalMeshSize << ") ***" << std::endl
                           << std::endl;
            displacementFieldRegistrationObserver->Execute( registrationMethod, itk::StartEvent() );
            registrationMethod->Update();
            }
          catch( itk::ExceptionObject & e )
            {
            this->Logger() << "Exception caught: " << e << std::endl;
            return EXIT_FAILURE;
            }

          // Add calculated transform to the composite transform
          this->m_CompositeTransform->AddTransform( outputDisplacementFieldTransform );
          }
        else
          {
          typedef itk::BSplineSyNImageRegistrationMethod<ImageType, ImageType,
            BSplineDisplacementFieldTransformType, ImageType, IntensityPointSetType>
            DisplacementFieldRegistrationType;

          typename DisplacementFieldRegistrationType::Pointer registrationMethod =
            this->PrepareRegistrationMethod<DisplacementFieldRegistrationType>(
                  this->m_CompositeTransform, currentStageNumber, VImageDimension,
                  preprocessedFixedImagesPerStage, preprocessedMovingImagesPerStage,
                  fixedIntensityPointSetsPerStage, movingIntensityPointSetsPerStage, stageMetricList, singleMetric,
                  multiMetric, optimizer, numberOfLevels, shrinkFactorsPerDimensionForAllLevels,
                  smoothingSigmasPerLevel, metricSamplingStrategy, samplingPercentage );

          typename BSplineDisplacementFieldTransformType::Pointer outputDisplacementFieldTransform =
            registrationMethod->GetModifiableTransform();

          // Create the transform adaptors

          typename DisplacementFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;

          outputDisplacementFieldTransform->SetSplineOrder( this->m_TransformMethods[currentStageNumber].m_SplineOrder );

          // Create the transform adaptors
          for( unsigned int level = 0; level < numberOfLevels; level++ )
            {
            typename itk::ImageBase<VImageDimension>::Pointer shrunkSpace=
                     this->GetShrinkImageOutputInformation(
                          virtualDomainImage.GetPointer(),
                          shrinkFactorsPerDimensionForAllLevels[level]  );

            typedef itk::BSplineSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor<
                BSplineDisplacementFieldTransformType>
              BSplineDisplacementFieldTransformAdaptorType;
            typename BSplineDisplacementFieldTransformAdaptorType::Pointer bsplineFieldTransformAdaptor =
              BSplineDisplacementFieldTransformAdaptorType::New();
            bsplineFieldTransformAdaptor->SetRequiredSpacing( shrunkSpace->GetSpacing() );
            bsplineFieldTransformAdaptor->SetRequiredSize( shrunkSpace->GetLargestPossibleRegion().GetSize() );
            bsplineFieldTransformAdaptor->SetRequiredDirection( shrunkSpace->GetDirection() );
            bsplineFieldTransformAdaptor->SetRequiredOrigin( shrunkSpace->GetOrigin() );
            bsplineFieldTransformAdaptor->SetTransform( outputDisplacementFieldTransform );

            // A good heuristic is to RealType the b-spline mesh resolution at each level
            typename BSplineDisplacementFieldTransformType::ArrayType newUpdateMeshSize = updateMeshSize;
            typename BSplineDisplacementFieldTransformType::ArrayType newTotalMeshSize = totalMeshSize;
            for( unsigned int d = 0; d < VImageDimension; d++ )
              {
              newUpdateMeshSize[d] = newUpdateMeshSize[d] << ( level );
              newTotalMeshSize[d] = newTotalMeshSize[d] << ( level );
              }

            bsplineFieldTransformAdaptor->SetMeshSizeForTheUpdateField( newUpdateMeshSize );
            bsplineFieldTransformAdaptor->SetMeshSizeForTheTotalField( newTotalMeshSize );

            adaptors.push_back( bsplineFieldTransformAdaptor.GetPointer() );
            }

          registrationMethod->SetDownsampleImagesForMetricDerivatives( true );
          registrationMethod->SetAverageMidPointGradients( false );
          registrationMethod->SetNumberOfLevels( numberOfLevels );

          typename DisplacementFieldRegistrationType::NumberOfIterationsArrayType numberOfIterationsPerLevel;
          numberOfIterationsPerLevel.SetSize( numberOfLevels );
          for( unsigned int d = 0; d < numberOfLevels; d++ )
            {
            numberOfIterationsPerLevel[d] = currentStageIterations[d];
            }

          registrationMethod->SetLearningRate( learningRate );
          registrationMethod->SetConvergenceThreshold( convergenceThreshold );
          registrationMethod->SetConvergenceWindowSize( convergenceWindowSize );
          registrationMethod->SetNumberOfIterationsPerLevel( numberOfIterationsPerLevel );
          registrationMethod->SetTransformParametersAdaptorsPerLevel( adaptors );
          outputDisplacementFieldTransform->SetDisplacementField( displacementField );
          outputDisplacementFieldTransform->SetInverseDisplacementField( inverseDisplacementField );

          typedef antsRegistrationCommandIterationUpdate<DisplacementFieldRegistrationType> DisplacementFieldCommandType;
          typename DisplacementFieldCommandType::Pointer displacementFieldRegistrationObserver =
            DisplacementFieldCommandType::New();
          displacementFieldRegistrationObserver->SetLogStream( *this->m_LogStream );
          displacementFieldRegistrationObserver->SetNumberOfIterations( currentStageIterations );

          registrationMethod->AddObserver( itk::IterationEvent(), displacementFieldRegistrationObserver );
          registrationMethod->AddObserver( itk::InitializeEvent(), displacementFieldRegistrationObserver );

          try
            {
            this->Logger() << std::endl << "*** Running B-spline SyN registration (updateMeshSizeAtBaseLevel = "
                           << updateMeshSize << ", totalMeshSizeAtBaseLevel = " << totalMeshSize << ") ***" << std::endl
                           << std::endl;
            displacementFieldRegistrationObserver->Execute( registrationMethod, itk::StartEvent() );
            registrationMethod->Update();
            }
          catch( itk::ExceptionObject & e )
            {
            this->Logger() << "Exception caught: " << e << std::endl;
            return EXIT_FAILURE;
            }

          // Add calculated transform to the composite transform
          this->m_CompositeTransform->AddTransform( outputDisplacementFieldTransform );
          }

        this->m_AllPreviousTransformsAreLinear = false;
        }
        break;
      case TimeVaryingVelocityField:
        {
        typedef itk::Vector<RealType, VImageDimension> VectorType;
        VectorType zeroVector( 0.0 );

        // Determine the parameters (size, spacing, etc) for the time-varying velocity field

        typedef itk::Image<VectorType, VImageDimension + 1>  TimeVaryingVelocityFieldType;

        typename TimeVaryingVelocityFieldType::IndexType velocityFieldIndex;
        typename TimeVaryingVelocityFieldType::SizeType velocityFieldSize;
        typename TimeVaryingVelocityFieldType::PointType velocityFieldOrigin;
        typename TimeVaryingVelocityFieldType::SpacingType velocityFieldSpacing;
        typename TimeVaryingVelocityFieldType::DirectionType velocityFieldDirection;
        typename TimeVaryingVelocityFieldType::RegionType velocityFieldRegion;

        typename ImageType::IndexType fixedImageIndex =
          preprocessedFixedImagesPerStage[0]->GetBufferedRegion().GetIndex();
        typename ImageType::SizeType fixedImageSize = preprocessedFixedImagesPerStage[0]->GetBufferedRegion().GetSize();
        typename ImageType::PointType fixedImageOrigin = preprocessedFixedImagesPerStage[0]->GetOrigin();
        typename ImageType::SpacingType fixedImageSpacing = preprocessedFixedImagesPerStage[0]->GetSpacing();
        typename ImageType::DirectionType fixedImageDirection = preprocessedFixedImagesPerStage[0]->GetDirection();

        unsigned int numberOfTimeIndices = this->m_TransformMethods[currentStageNumber].m_NumberOfTimeIndices;

        velocityFieldIndex.Fill( 0 );
        velocityFieldSize.Fill( numberOfTimeIndices );
        velocityFieldOrigin.Fill( 0.0 );
        velocityFieldSpacing.Fill( 1.0 );
        velocityFieldDirection.SetIdentity();
        for( unsigned int i = 0; i < VImageDimension; i++ )
          {
          velocityFieldIndex[i] = fixedImageIndex[i];
          velocityFieldSize[i] = fixedImageSize[i];
          velocityFieldOrigin[i] = fixedImageOrigin[i];
          velocityFieldSpacing[i] = fixedImageSpacing[i];
          for( unsigned int j = 0; j < VImageDimension; j++ )
            {
            velocityFieldDirection[i][j] = fixedImageDirection[i][j];
            }
          }

        velocityFieldRegion.SetSize( velocityFieldSize );
        velocityFieldRegion.SetIndex( velocityFieldIndex );
        typename TimeVaryingVelocityFieldType::Pointer velocityField =
          AllocImage<TimeVaryingVelocityFieldType>(velocityFieldRegion,
                                                   velocityFieldSpacing,
                                                   velocityFieldOrigin,
                                                   velocityFieldDirection,
                                                   zeroVector);

        typename DisplacementFieldType::Pointer displacementField =
          AllocImage<DisplacementFieldType>( preprocessedFixedImagesPerStage[0]->GetBufferedRegion(),
                                                    fixedImageSpacing,
                                                    fixedImageOrigin,
                                                    fixedImageDirection,
                                                    zeroVector );
        typename DisplacementFieldType::Pointer inverseDisplacementField =
          AllocImage<DisplacementFieldType>( preprocessedFixedImagesPerStage[0]->GetBufferedRegion(),
                                                    fixedImageSpacing,
                                                    fixedImageOrigin,
                                                    fixedImageDirection,
                                                    zeroVector );

        // Extract parameters

        RealType varianceForUpdateField =
          this->m_TransformMethods[currentStageNumber].m_UpdateFieldVarianceInVarianceSpace;
        RealType varianceForUpdateFieldTime = this->m_TransformMethods[currentStageNumber].m_UpdateFieldTimeSigma;
        RealType varianceForTotalField =
          this->m_TransformMethods[currentStageNumber].m_TotalFieldVarianceInVarianceSpace;
        RealType varianceForTotalFieldTime = this->m_TransformMethods[currentStageNumber].m_TotalFieldTimeSigma;

        typedef itk::GaussianSmoothingOnUpdateTimeVaryingVelocityFieldTransform <TComputeType, ImageType::ImageDimension>
           TimeVaryingVelocityFieldOutputTransformType;

        typedef itk::TimeVaryingVelocityFieldImageRegistrationMethodv4<ImageType, ImageType,
          TimeVaryingVelocityFieldOutputTransformType, ImageType, LabeledPointSetType>  VelocityFieldRegistrationType;
        typename VelocityFieldRegistrationType::Pointer velocityFieldRegistration =
          VelocityFieldRegistrationType::New();

        if( this->m_RestrictDeformationOptimizerWeights.size() > currentStageNumber )
          {
          if( this->m_RestrictDeformationOptimizerWeights[currentStageNumber].size() == VImageDimension )
            {
            typename VelocityFieldRegistrationType::OptimizerWeightsType optimizerWeights( VImageDimension );
            for( unsigned int d = 0; d < VImageDimension; d++ )
              {
              optimizerWeights[d] = this->m_RestrictDeformationOptimizerWeights[currentStageNumber][d];
              }
            velocityFieldRegistration->SetOptimizerWeights( optimizerWeights );
            }
          }

        typedef typename VelocityFieldRegistrationType::OutputTransformType OutputTransformType;
        typename OutputTransformType::Pointer outputTransform = velocityFieldRegistration->GetModifiableTransform();
        for( unsigned int n = 0; n < stageMetricList.size(); n++ )
          {
          if( !this->IsPointSetMetric( stageMetricList[n].m_MetricType ) )
            {
            velocityFieldRegistration->SetFixedImage( n, preprocessedFixedImagesPerStage[n] );
            velocityFieldRegistration->SetMovingImage( n, preprocessedMovingImagesPerStage[n] );
            }
          else
            {
            velocityFieldRegistration->SetFixedPointSet( n, stageMetricList[n].m_FixedLabeledPointSet.GetPointer() );
            velocityFieldRegistration->SetMovingPointSet( n, stageMetricList[n].m_MovingLabeledPointSet.GetPointer() );
            }
          }
        if( useMultiMetric )
          {
          velocityFieldRegistration->SetMetric( multiMetric );
          }
        else
          {
          velocityFieldRegistration->SetMetric( singleMetric );
          }

        if( this->m_CompositeTransform->GetNumberOfTransforms() > 0 )
          {
          velocityFieldRegistration->SetMovingInitialTransform( this->m_CompositeTransform );
          }
        if( this->m_FixedInitialTransform->GetNumberOfTransforms() > 0 )
          {
          velocityFieldRegistration->SetFixedInitialTransform( this->m_FixedInitialTransform );
          }
        velocityFieldRegistration->SetNumberOfLevels( numberOfLevels );
        velocityFieldRegistration->SetMetricSamplingStrategy(
          static_cast<typename VelocityFieldRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
        velocityFieldRegistration->SetMetricSamplingPercentage( samplingPercentage );
        velocityFieldRegistration->SetLearningRate( learningRate );
        velocityFieldRegistration->SetConvergenceThreshold( convergenceThreshold );
        velocityFieldRegistration->SetConvergenceWindowSize( convergenceWindowSize );
        outputTransform->SetGaussianSpatialSmoothingVarianceForTheTotalField( varianceForTotalField );
        outputTransform->SetGaussianSpatialSmoothingVarianceForTheUpdateField( varianceForUpdateField );
        outputTransform->SetGaussianTemporalSmoothingVarianceForTheTotalField( varianceForTotalFieldTime );
        outputTransform->SetGaussianTemporalSmoothingVarianceForTheUpdateField( varianceForUpdateFieldTime );

        outputTransform->SetTimeVaryingVelocityField( velocityField );
        outputTransform->SetLowerTimeBound( 0.0 );
        outputTransform->SetUpperTimeBound( 1.0 );
        outputTransform->SetDisplacementField( displacementField );
        outputTransform->SetInverseDisplacementField( inverseDisplacementField );

        typename VelocityFieldRegistrationType::NumberOfIterationsArrayType numberOfIterationsPerLevel;
        numberOfIterationsPerLevel.SetSize( numberOfLevels );
        for( unsigned int d = 0; d < numberOfLevels; d++ )
          {
          numberOfIterationsPerLevel[d] = currentStageIterations[d];
          }
        velocityFieldRegistration->SetNumberOfIterationsPerLevel( numberOfIterationsPerLevel );

        for( unsigned int level = 0; level < numberOfLevels; ++level )
          {
          velocityFieldRegistration->SetShrinkFactorsPerDimension( level, shrinkFactorsPerDimensionForAllLevels[level] );
          }
        velocityFieldRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
        velocityFieldRegistration->SetSmoothingSigmasAreSpecifiedInPhysicalUnits(
          this->m_SmoothingSigmasAreInPhysicalUnits[currentStageNumber] );

        typedef itk::TimeVaryingVelocityFieldTransformParametersAdaptor<OutputTransformType>
          VelocityFieldTransformAdaptorType;

        typename VelocityFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;
        for( unsigned int level = 0; level < numberOfLevels; level++ )
          {
          typename itk::ImageBase<VImageDimension>::Pointer shrunkSpace=
                     this->GetShrinkImageOutputInformation(
                          virtualDomainImage.GetPointer(),
                          shrinkFactorsPerDimensionForAllLevels[level]  );

          // Although we shrink the images for the given levels,
          // we keep the size in time the same

          velocityFieldSize.Fill( numberOfTimeIndices );
          velocityFieldOrigin.Fill( 0.0 );
          velocityFieldSpacing.Fill( 1.0 );
          velocityFieldDirection.SetIdentity();

          fixedImageSize = shrunkSpace->GetLargestPossibleRegion().GetSize();
          fixedImageOrigin = shrunkSpace->GetOrigin();
          fixedImageSpacing = shrunkSpace->GetSpacing();
          fixedImageDirection = shrunkSpace->GetDirection();
          for( unsigned int i = 0; i < VImageDimension; i++ )
            {
            velocityFieldSize[i] = fixedImageSize[i];
            velocityFieldOrigin[i] = fixedImageOrigin[i];
            velocityFieldSpacing[i] = fixedImageSpacing[i];
            for( unsigned int j = 0; j < VImageDimension; j++ )
              {
              velocityFieldDirection[i][j] = fixedImageDirection[i][j];
              }
            }

          typename VelocityFieldTransformAdaptorType::Pointer fieldTransformAdaptor =
            VelocityFieldTransformAdaptorType::New();
          fieldTransformAdaptor->SetRequiredSpacing( velocityFieldSpacing );
          fieldTransformAdaptor->SetRequiredSize( velocityFieldSize );
          fieldTransformAdaptor->SetRequiredDirection( velocityFieldDirection );
          fieldTransformAdaptor->SetRequiredOrigin( velocityFieldOrigin );

          adaptors.push_back( fieldTransformAdaptor.GetPointer() );
          }

        velocityFieldRegistration->SetTransformParametersAdaptorsPerLevel( adaptors );

        typedef antsRegistrationCommandIterationUpdate<VelocityFieldRegistrationType> VelocityFieldCommandType;
        typename VelocityFieldCommandType::Pointer velocityFieldRegistrationObserver = VelocityFieldCommandType::New();
        velocityFieldRegistrationObserver->SetLogStream( *this->m_LogStream );
        velocityFieldRegistrationObserver->SetNumberOfIterations( currentStageIterations );

        velocityFieldRegistration->AddObserver( itk::IterationEvent(), velocityFieldRegistrationObserver );
        velocityFieldRegistration->AddObserver( itk::InitializeEvent(), velocityFieldRegistrationObserver );

        try
          {
          this->Logger() << std::endl
                         << "*** Running time-varying velocity field registration (varianceForUpdateField = "
                         << varianceForUpdateField << ", varianceForTotalField = " << varianceForTotalField
                         << ", varianceForUpdateFieldTime = "
                         << varianceForUpdateFieldTime << ", varianceForTotalFieldTime = " << varianceForTotalFieldTime
                         << ") ***" << std::endl << std::endl;
          velocityFieldRegistrationObserver->Execute( velocityFieldRegistration, itk::StartEvent() );
          velocityFieldRegistration->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          this->Logger() << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }
        // Add calculated transform to the composite transform
        this->m_CompositeTransform->AddTransform( outputTransform );

        this->m_AllPreviousTransformsAreLinear = false;
        }
        break;
      case TimeVaryingBSplineVelocityField:
        {
        typedef itk::Vector<RealType, VImageDimension> VectorType;
        VectorType zeroVector( 0.0 );

        // Determine the parameters (size, spacing, etc) for the time-varying velocity field control point lattice

        const std::vector<unsigned int> & meshSize =
          this->m_TransformMethods[currentStageNumber].m_VelocityFieldMeshSize;
        if( meshSize.size() != VImageDimension + 1 )
          {
          this->Logger() << "The transform domain mesh size does not have the correct number of elements."
                           << "For image dimension = " << VImageDimension << ", you need " << VImageDimension + 1
                           << "elements. " << std::endl;
          return EXIT_FAILURE;
          }

        unsigned int numberOfTimePointSamples =
          this->m_TransformMethods[currentStageNumber].m_NumberOfTimePointSamples;
        unsigned int splineOrder = this->m_TransformMethods[currentStageNumber].m_SplineOrder;

        typedef itk::Image<VectorType, VImageDimension + 1> TimeVaryingVelocityFieldControlPointLatticeType;

        typename ImageType::SizeType fixedImageSize = virtualDomainImage->GetBufferedRegion().GetSize();
        typename ImageType::PointType fixedImageOrigin = virtualDomainImage->GetOrigin();
        typename ImageType::SpacingType fixedImageSpacing = virtualDomainImage->GetSpacing();
        typename ImageType::DirectionType fixedImageDirection = virtualDomainImage->GetDirection();

        typename TimeVaryingVelocityFieldControlPointLatticeType::SizeType transformDomainMeshSize;
        typename TimeVaryingVelocityFieldControlPointLatticeType::PointType transformDomainOrigin;
        typename TimeVaryingVelocityFieldControlPointLatticeType::SpacingType transformDomainSpacing;
        typename TimeVaryingVelocityFieldControlPointLatticeType::SizeType transformDomainSize;
        typename TimeVaryingVelocityFieldControlPointLatticeType::DirectionType transformDomainDirection;

        transformDomainDirection.SetIdentity();
        transformDomainOrigin.Fill( 0.0 );
        transformDomainSpacing.Fill( 1.0 );
        transformDomainSize.Fill( 2 );
        for( unsigned int i = 0; i < VImageDimension; i++ )
          {
          transformDomainOrigin[i] = fixedImageOrigin[i];
          transformDomainMeshSize[i] = 3;
          transformDomainSpacing[i] = fixedImageSpacing[i];
          transformDomainSize[i] = fixedImageSize[i];
          for( unsigned int j = 0; j < VImageDimension; j++ )
            {
            transformDomainDirection[i][j] = fixedImageDirection[i][j];
            }
          }
        for( unsigned int i = 0; i < meshSize.size(); i++ )
          {
          transformDomainMeshSize[i] = meshSize[i];
          }
        typename TimeVaryingVelocityFieldControlPointLatticeType::SizeType initialTransformDomainMeshSize =
          transformDomainMeshSize;

        typedef itk::TimeVaryingBSplineVelocityFieldTransform <TComputeType, ImageType::ImageDimension>
          TimeVaryingBSplineVelocityFieldOutputTransformType;

        if( stageMetricList[0].m_MetricType != IGDM )
          {
          typedef itk::TimeVaryingBSplineVelocityFieldImageRegistrationMethod<ImageType, ImageType,
            TimeVaryingBSplineVelocityFieldOutputTransformType, ImageType, LabeledPointSetType>
            VelocityFieldRegistrationType;

          typename VelocityFieldRegistrationType::Pointer velocityFieldRegistration =
            this->PrepareRegistrationMethod<VelocityFieldRegistrationType>(
                  this->m_CompositeTransform, currentStageNumber, VImageDimension,
                  preprocessedFixedImagesPerStage, preprocessedMovingImagesPerStage,
                  fixedLabeledPointSetsPerStage, movingLabeledPointSetsPerStage, stageMetricList, singleMetric,
                  multiMetric, optimizer, numberOfLevels, shrinkFactorsPerDimensionForAllLevels,
                  smoothingSigmasPerLevel, metricSamplingStrategy, samplingPercentage );

          typename TimeVaryingBSplineVelocityFieldOutputTransformType::Pointer
            outputTransform = velocityFieldRegistration->GetModifiableTransform();

          if( useMultiMetric )
            {
            velocityFieldRegistration->SetMetric( multiMetric );
            }
          else
            {
            velocityFieldRegistration->SetMetric( singleMetric );
            }

          velocityFieldRegistration->SetNumberOfTimePointSamples( numberOfTimePointSamples );
          velocityFieldRegistration->SetLearningRate( learningRate );
          velocityFieldRegistration->SetConvergenceThreshold( convergenceThreshold );
          velocityFieldRegistration->SetConvergenceWindowSize( convergenceWindowSize );
          outputTransform->SetSplineOrder( splineOrder );
          outputTransform->SetLowerTimeBound( 0.0 );
          outputTransform->SetUpperTimeBound( 1.0 );

          typedef itk::TimeVaryingBSplineVelocityFieldTransformParametersAdaptor<TimeVaryingBSplineVelocityFieldOutputTransformType>
            VelocityFieldTransformAdaptorType;
          typename VelocityFieldTransformAdaptorType::Pointer initialFieldTransformAdaptor =
            VelocityFieldTransformAdaptorType::New();
          initialFieldTransformAdaptor->SetTransform( outputTransform );
          initialFieldTransformAdaptor->SetRequiredTransformDomainOrigin( transformDomainOrigin );
          initialFieldTransformAdaptor->SetRequiredTransformDomainSpacing( transformDomainSpacing );
          initialFieldTransformAdaptor->SetRequiredTransformDomainSize( transformDomainSize );
          initialFieldTransformAdaptor->SetRequiredTransformDomainMeshSize( transformDomainMeshSize );
          initialFieldTransformAdaptor->SetRequiredTransformDomainDirection( transformDomainDirection );

          typename TimeVaryingVelocityFieldControlPointLatticeType::Pointer
          velocityFieldLattice = AllocImage<TimeVaryingVelocityFieldControlPointLatticeType>
              ( initialFieldTransformAdaptor->GetRequiredControlPointLatticeSize(),
              initialFieldTransformAdaptor->GetRequiredControlPointLatticeSpacing(),
              initialFieldTransformAdaptor->GetRequiredControlPointLatticeOrigin(),
              initialFieldTransformAdaptor->GetRequiredControlPointLatticeDirection(),
              zeroVector );

          typename TimeVaryingBSplineVelocityFieldOutputTransformType::VelocityFieldPointType        sampledVelocityFieldOrigin;
          typename TimeVaryingBSplineVelocityFieldOutputTransformType::VelocityFieldSpacingType      sampledVelocityFieldSpacing;
          typename TimeVaryingBSplineVelocityFieldOutputTransformType::VelocityFieldSizeType         sampledVelocityFieldSize;
          typename TimeVaryingBSplineVelocityFieldOutputTransformType::VelocityFieldDirectionType    sampledVelocityFieldDirection;

          sampledVelocityFieldOrigin.Fill( 0.0 );
          sampledVelocityFieldSpacing.Fill( 1.0 );
          sampledVelocityFieldSize.Fill( numberOfTimePointSamples );
          sampledVelocityFieldDirection.SetIdentity();
          for( unsigned int i = 0; i < VImageDimension; i++ )
            {
            sampledVelocityFieldOrigin[i] = virtualDomainImage->GetOrigin()[i];
            sampledVelocityFieldSpacing[i] = virtualDomainImage->GetSpacing()[i];
            sampledVelocityFieldSize[i] = virtualDomainImage->GetRequestedRegion().GetSize()[i];
            for( unsigned int j = 0; j < VImageDimension; j++ )
              {
              sampledVelocityFieldDirection[i][j] = virtualDomainImage->GetDirection()[i][j];
              }
            }

          outputTransform->SetTimeVaryingVelocityFieldControlPointLattice( velocityFieldLattice );
          outputTransform->SetVelocityFieldOrigin( sampledVelocityFieldOrigin );
          outputTransform->SetVelocityFieldDirection( sampledVelocityFieldDirection );
          outputTransform->SetVelocityFieldSpacing( sampledVelocityFieldSpacing );
          outputTransform->SetVelocityFieldSize( sampledVelocityFieldSize );

          typename VelocityFieldRegistrationType::NumberOfIterationsArrayType numberOfIterationsPerLevel;
          numberOfIterationsPerLevel.SetSize( numberOfLevels );
          for( unsigned int d = 0; d < numberOfLevels; d++ )
            {
            numberOfIterationsPerLevel[d] = currentStageIterations[d];
            }
          velocityFieldRegistration->SetNumberOfIterationsPerLevel( numberOfIterationsPerLevel );

          for( unsigned int level = 0; level < numberOfLevels; ++level )
            {
            velocityFieldRegistration->SetShrinkFactorsPerDimension( level, shrinkFactorsPerDimensionForAllLevels[level] );
            }
          velocityFieldRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
          velocityFieldRegistration->SetSmoothingSigmasAreSpecifiedInPhysicalUnits(
            this->m_SmoothingSigmasAreInPhysicalUnits[currentStageNumber] );

          typename VelocityFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;
          for( unsigned int level = 0; level < numberOfLevels; level++ )
            {
            typename VelocityFieldTransformAdaptorType::Pointer fieldTransformAdaptor =
              VelocityFieldTransformAdaptorType::New();
            fieldTransformAdaptor->SetTransform( outputTransform );
            fieldTransformAdaptor->SetRequiredTransformDomainOrigin( transformDomainOrigin );
            fieldTransformAdaptor->SetRequiredTransformDomainMeshSize( transformDomainMeshSize );
            fieldTransformAdaptor->SetRequiredTransformDomainSpacing( transformDomainSpacing );
            fieldTransformAdaptor->SetRequiredTransformDomainSize( transformDomainSize );

            adaptors.push_back( fieldTransformAdaptor.GetPointer() );
            for( unsigned int i = 0; i <= VImageDimension; i++ )
              {
              transformDomainMeshSize[i] <<= 1;
              }
            }
          velocityFieldRegistration->SetTransformParametersAdaptorsPerLevel( adaptors );

          typedef antsRegistrationCommandIterationUpdate<VelocityFieldRegistrationType> VelocityFieldCommandType;
          typename VelocityFieldCommandType::Pointer velocityFieldRegistrationObserver = VelocityFieldCommandType::New();
          velocityFieldRegistrationObserver->SetLogStream( *this->m_LogStream );
          velocityFieldRegistrationObserver->SetNumberOfIterations( currentStageIterations );

          velocityFieldRegistration->AddObserver( itk::IterationEvent(), velocityFieldRegistrationObserver );
          velocityFieldRegistration->AddObserver( itk::InitializeEvent(), velocityFieldRegistrationObserver );

          try
            {
            this->Logger() << std::endl
                           << "*** Running time-varying b-spline velocity field registration (initial mesh size = "
                           << initialTransformDomainMeshSize << ") ***" << std::endl << std::endl;
            velocityFieldRegistrationObserver->Execute( velocityFieldRegistration, itk::StartEvent() );
            velocityFieldRegistration->Update();
            }
          catch( itk::ExceptionObject & e )
            {
            this->Logger() << "Exception caught: " << e << std::endl;
            return EXIT_FAILURE;
            }
          // Add calculated transform to the composite transform
          this->m_CompositeTransform->AddTransform( outputTransform );
          }
        else
          {
          typedef itk::TimeVaryingBSplineVelocityFieldImageRegistrationMethod<ImageType, ImageType,
            TimeVaryingBSplineVelocityFieldOutputTransformType, ImageType, IntensityPointSetType>
            VelocityFieldRegistrationType;

          typename VelocityFieldRegistrationType::Pointer velocityFieldRegistration =
            this->PrepareRegistrationMethod<VelocityFieldRegistrationType>(
                  this->m_CompositeTransform, currentStageNumber, VImageDimension,
                  preprocessedFixedImagesPerStage, preprocessedMovingImagesPerStage,
                  fixedIntensityPointSetsPerStage, movingIntensityPointSetsPerStage, stageMetricList, singleMetric,
                  multiMetric, optimizer, numberOfLevels, shrinkFactorsPerDimensionForAllLevels,
                  smoothingSigmasPerLevel, metricSamplingStrategy, samplingPercentage );

          typename TimeVaryingBSplineVelocityFieldOutputTransformType::Pointer
            outputTransform = velocityFieldRegistration->GetModifiableTransform();

          if( useMultiMetric )
            {
            velocityFieldRegistration->SetMetric( multiMetric );
            }
          else
            {
            velocityFieldRegistration->SetMetric( singleMetric );
            }

          velocityFieldRegistration->SetNumberOfTimePointSamples( numberOfTimePointSamples );
          velocityFieldRegistration->SetLearningRate( learningRate );
          velocityFieldRegistration->SetConvergenceThreshold( convergenceThreshold );
          velocityFieldRegistration->SetConvergenceWindowSize( convergenceWindowSize );
          outputTransform->SetSplineOrder( splineOrder );
          outputTransform->SetLowerTimeBound( 0.0 );
          outputTransform->SetUpperTimeBound( 1.0 );

          typedef itk::TimeVaryingBSplineVelocityFieldTransformParametersAdaptor<TimeVaryingBSplineVelocityFieldOutputTransformType>
            VelocityFieldTransformAdaptorType;
          typename VelocityFieldTransformAdaptorType::Pointer initialFieldTransformAdaptor =
            VelocityFieldTransformAdaptorType::New();
          initialFieldTransformAdaptor->SetTransform( outputTransform );
          initialFieldTransformAdaptor->SetRequiredTransformDomainOrigin( transformDomainOrigin );
          initialFieldTransformAdaptor->SetRequiredTransformDomainSpacing( transformDomainSpacing );
          initialFieldTransformAdaptor->SetRequiredTransformDomainSize( transformDomainSize );
          initialFieldTransformAdaptor->SetRequiredTransformDomainMeshSize( transformDomainMeshSize );
          initialFieldTransformAdaptor->SetRequiredTransformDomainDirection( transformDomainDirection );

          typename TimeVaryingVelocityFieldControlPointLatticeType::Pointer
          velocityFieldLattice = AllocImage<TimeVaryingVelocityFieldControlPointLatticeType>
              ( initialFieldTransformAdaptor->GetRequiredControlPointLatticeSize(),
              initialFieldTransformAdaptor->GetRequiredControlPointLatticeSpacing(),
              initialFieldTransformAdaptor->GetRequiredControlPointLatticeOrigin(),
              initialFieldTransformAdaptor->GetRequiredControlPointLatticeDirection(),
              zeroVector );

          typename TimeVaryingBSplineVelocityFieldOutputTransformType::VelocityFieldPointType        sampledVelocityFieldOrigin;
          typename TimeVaryingBSplineVelocityFieldOutputTransformType::VelocityFieldSpacingType      sampledVelocityFieldSpacing;
          typename TimeVaryingBSplineVelocityFieldOutputTransformType::VelocityFieldSizeType         sampledVelocityFieldSize;
          typename TimeVaryingBSplineVelocityFieldOutputTransformType::VelocityFieldDirectionType    sampledVelocityFieldDirection;

          sampledVelocityFieldOrigin.Fill( 0.0 );
          sampledVelocityFieldSpacing.Fill( 1.0 );
          sampledVelocityFieldSize.Fill( numberOfTimePointSamples );
          sampledVelocityFieldDirection.SetIdentity();
          for( unsigned int i = 0; i < VImageDimension; i++ )
            {
            sampledVelocityFieldOrigin[i] = virtualDomainImage->GetOrigin()[i];
            sampledVelocityFieldSpacing[i] = virtualDomainImage->GetSpacing()[i];
            sampledVelocityFieldSize[i] = virtualDomainImage->GetRequestedRegion().GetSize()[i];
            for( unsigned int j = 0; j < VImageDimension; j++ )
              {
              sampledVelocityFieldDirection[i][j] = virtualDomainImage->GetDirection()[i][j];
              }
            }

          outputTransform->SetTimeVaryingVelocityFieldControlPointLattice( velocityFieldLattice );
          outputTransform->SetVelocityFieldOrigin( sampledVelocityFieldOrigin );
          outputTransform->SetVelocityFieldDirection( sampledVelocityFieldDirection );
          outputTransform->SetVelocityFieldSpacing( sampledVelocityFieldSpacing );
          outputTransform->SetVelocityFieldSize( sampledVelocityFieldSize );

          typename VelocityFieldRegistrationType::NumberOfIterationsArrayType numberOfIterationsPerLevel;
          numberOfIterationsPerLevel.SetSize( numberOfLevels );
          for( unsigned int d = 0; d < numberOfLevels; d++ )
            {
            numberOfIterationsPerLevel[d] = currentStageIterations[d];
            }
          velocityFieldRegistration->SetNumberOfIterationsPerLevel( numberOfIterationsPerLevel );

          for( unsigned int level = 0; level < numberOfLevels; ++level )
            {
            velocityFieldRegistration->SetShrinkFactorsPerDimension( level, shrinkFactorsPerDimensionForAllLevels[level] );
            }
          velocityFieldRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
          velocityFieldRegistration->SetSmoothingSigmasAreSpecifiedInPhysicalUnits(
            this->m_SmoothingSigmasAreInPhysicalUnits[currentStageNumber] );

          typename VelocityFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;
          for( unsigned int level = 0; level < numberOfLevels; level++ )
            {
            typename VelocityFieldTransformAdaptorType::Pointer fieldTransformAdaptor =
              VelocityFieldTransformAdaptorType::New();
            fieldTransformAdaptor->SetTransform( outputTransform );
            fieldTransformAdaptor->SetRequiredTransformDomainOrigin( transformDomainOrigin );
            fieldTransformAdaptor->SetRequiredTransformDomainMeshSize( transformDomainMeshSize );
            fieldTransformAdaptor->SetRequiredTransformDomainSpacing( transformDomainSpacing );
            fieldTransformAdaptor->SetRequiredTransformDomainSize( transformDomainSize );

            adaptors.push_back( fieldTransformAdaptor.GetPointer() );
            for( unsigned int i = 0; i <= VImageDimension; i++ )
              {
              transformDomainMeshSize[i] <<= 1;
              }
            }
          velocityFieldRegistration->SetTransformParametersAdaptorsPerLevel( adaptors );

          typedef antsRegistrationCommandIterationUpdate<VelocityFieldRegistrationType> VelocityFieldCommandType;
          typename VelocityFieldCommandType::Pointer velocityFieldRegistrationObserver = VelocityFieldCommandType::New();
          velocityFieldRegistrationObserver->SetLogStream( *this->m_LogStream );
          velocityFieldRegistrationObserver->SetNumberOfIterations( currentStageIterations );

          velocityFieldRegistration->AddObserver( itk::IterationEvent(), velocityFieldRegistrationObserver );
          velocityFieldRegistration->AddObserver( itk::InitializeEvent(), velocityFieldRegistrationObserver );

          try
            {
            this->Logger() << std::endl
                           << "*** Running time-varying b-spline velocity field registration (initial mesh size = "
                           << initialTransformDomainMeshSize << ") ***" << std::endl << std::endl;
            velocityFieldRegistrationObserver->Execute( velocityFieldRegistration, itk::StartEvent() );
            velocityFieldRegistration->Update();
            }
          catch( itk::ExceptionObject & e )
            {
            this->Logger() << "Exception caught: " << e << std::endl;
            return EXIT_FAILURE;
            }
          // Add calculated transform to the composite transform
          this->m_CompositeTransform->AddTransform( outputTransform );
          }

        this->m_AllPreviousTransformsAreLinear = false;
        }
        break;
      case Exponential:
        {
        typedef itk::Vector<RealType, VImageDimension> VectorType;
        VectorType zeroVector( 0.0 );

        typedef itk::Image<VectorType, VImageDimension> ConstantVelocityFieldType;

        typename ConstantVelocityFieldType::Pointer constantVelocityField = AllocImage<ConstantVelocityFieldType>(
            preprocessedFixedImagesPerStage[0], zeroVector );

        typedef itk::GaussianExponentialDiffeomorphicTransform<RealType,
                                                               VImageDimension> GaussianDisplacementFieldTransformType;

        typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, GaussianDisplacementFieldTransformType,
          ImageType, LabeledPointSetType> DisplacementFieldRegistrationType;
        typename DisplacementFieldRegistrationType::Pointer displacementFieldRegistration =
          DisplacementFieldRegistrationType::New();

        if( this->m_RestrictDeformationOptimizerWeights.size() > currentStageNumber )
          {
          if( this->m_RestrictDeformationOptimizerWeights[currentStageNumber].size() == VImageDimension )
            {
            typename DisplacementFieldRegistrationType::OptimizerWeightsType optimizerWeights( VImageDimension );
            for( unsigned int d = 0; d < VImageDimension; d++ )
              {
              optimizerWeights[d] = this->m_RestrictDeformationOptimizerWeights[currentStageNumber][d];
              }
            displacementFieldRegistration->SetOptimizerWeights( optimizerWeights );
            }
          }

        typename GaussianDisplacementFieldTransformType::Pointer outputDisplacementFieldTransform = displacementFieldRegistration->GetModifiableTransform();

        // Create the transform adaptors

        typedef itk::GaussianExponentialDiffeomorphicTransformParametersAdaptor<GaussianDisplacementFieldTransformType>
          DisplacementFieldTransformAdaptorType;
        typename DisplacementFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;

        // Extract parameters
        RealType varianceForUpdateField =
          this->m_TransformMethods[currentStageNumber].m_UpdateFieldVarianceInVarianceSpace;
        RealType varianceForVelocityField  =
          this->m_TransformMethods[currentStageNumber].m_VelocityFieldVarianceInVarianceSpace;
        unsigned int numberOfIntegrationSteps = this->m_TransformMethods[currentStageNumber].m_NumberOfTimeIndices;

        outputDisplacementFieldTransform->SetGaussianSmoothingVarianceForTheUpdateField( varianceForUpdateField );
        outputDisplacementFieldTransform->SetGaussianSmoothingVarianceForTheConstantVelocityField(
          varianceForVelocityField );
        if( numberOfIntegrationSteps == 0 )
          {
          outputDisplacementFieldTransform->SetCalculateNumberOfIntegrationStepsAutomatically( true );
          }
        else
          {
          outputDisplacementFieldTransform->SetNumberOfIntegrationSteps( numberOfIntegrationSteps );
          }
        outputDisplacementFieldTransform->SetConstantVelocityField( constantVelocityField );
        outputDisplacementFieldTransform->SetDisplacementField( constantVelocityField );
        // Create the transform adaptors
        // For the gaussian displacement field, the specified variances are in image spacing terms
        // and, in normal practice, we typically don't change these values at each level.  However,
        // if the user wishes to add that option, they can use the class
        // GaussianSmoothingOnUpdateDisplacementFieldTransformAdaptor
        for( unsigned int level = 0; level < numberOfLevels; level++ )
          {
          typename itk::ImageBase<VImageDimension>::Pointer shrunkSpace=
                     this->GetShrinkImageOutputInformation(
                          virtualDomainImage.GetPointer(),
                          shrinkFactorsPerDimensionForAllLevels[level]  );

          typename DisplacementFieldTransformAdaptorType::Pointer fieldTransformAdaptor =
            DisplacementFieldTransformAdaptorType::New();
          fieldTransformAdaptor->SetRequiredSpacing( shrunkSpace->GetSpacing() );
          fieldTransformAdaptor->SetRequiredSize( shrunkSpace->GetLargestPossibleRegion().GetSize() );
          fieldTransformAdaptor->SetRequiredDirection( shrunkSpace->GetDirection() );
          fieldTransformAdaptor->SetRequiredOrigin( shrunkSpace->GetOrigin() );
          fieldTransformAdaptor->SetTransform( outputDisplacementFieldTransform );

          fieldTransformAdaptor->SetGaussianSmoothingVarianceForTheUpdateField( varianceForUpdateField );
          fieldTransformAdaptor->SetGaussianSmoothingVarianceForTheConstantVelocityField( varianceForVelocityField );

          adaptors.push_back( fieldTransformAdaptor.GetPointer() );
          }
        for( unsigned int n = 0; n < stageMetricList.size(); n++ )
          {
          if( !this->IsPointSetMetric( stageMetricList[n].m_MetricType ) )
            {
            displacementFieldRegistration->SetFixedImage( n, preprocessedFixedImagesPerStage[n] );
            displacementFieldRegistration->SetMovingImage( n, preprocessedMovingImagesPerStage[n] );
            }
          else
            {
            displacementFieldRegistration->SetFixedPointSet( n, stageMetricList[n].m_FixedLabeledPointSet.GetPointer() );
            displacementFieldRegistration->SetMovingPointSet( n, stageMetricList[n].m_MovingLabeledPointSet.GetPointer() );
            }
          }
        if( useMultiMetric )
          {
          displacementFieldRegistration->SetMetric( multiMetric );
          }
        else
          {
          displacementFieldRegistration->SetMetric( singleMetric );
          }

        displacementFieldRegistration->SetNumberOfLevels( numberOfLevels );

        for( unsigned int level = 0; level < numberOfLevels; ++level )
          {
          displacementFieldRegistration->SetShrinkFactorsPerDimension( level, shrinkFactorsPerDimensionForAllLevels[level] );
          }
        displacementFieldRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
        displacementFieldRegistration->SetMetricSamplingStrategy(
          static_cast<typename DisplacementFieldRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
        displacementFieldRegistration->SetMetricSamplingPercentage( samplingPercentage );
        displacementFieldRegistration->SetOptimizer( optimizer2 );

        displacementFieldRegistration->SetTransformParametersAdaptorsPerLevel( adaptors );
        if( this->m_CompositeTransform->GetNumberOfTransforms() > 0 )
          {
          displacementFieldRegistration->SetMovingInitialTransform( this->m_CompositeTransform );
          }
        if( this->m_FixedInitialTransform->GetNumberOfTransforms() > 0 )
          {
          displacementFieldRegistration->SetFixedInitialTransform( this->m_FixedInitialTransform );
          }

        typedef antsRegistrationCommandIterationUpdate<DisplacementFieldRegistrationType> DisplacementFieldCommandType;
        typename DisplacementFieldCommandType::Pointer displacementFieldRegistrationObserver =
          DisplacementFieldCommandType::New();
        displacementFieldRegistrationObserver->SetLogStream(*this->m_LogStream );
        displacementFieldRegistrationObserver->SetNumberOfIterations( currentStageIterations );

        displacementFieldRegistration->AddObserver( itk::IterationEvent(), displacementFieldRegistrationObserver );
        displacementFieldRegistration->AddObserver( itk::InitializeEvent(), displacementFieldRegistrationObserver );

        try
          {
          this->Logger() << std::endl
                         << "*** Running gaussian exponential field registration (varianceForUpdateField = "
                         << varianceForUpdateField << ", varianceForVelocityField = " << varianceForVelocityField
                         << ") ***"
                         << std::endl << std::endl;
          displacementFieldRegistrationObserver->Execute( displacementFieldRegistration, itk::StartEvent() );
          displacementFieldRegistration->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          this->Logger() << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }

        // Add calculated transform to the composite transform
        this->m_CompositeTransform->AddTransform( outputDisplacementFieldTransform );

        this->m_AllPreviousTransformsAreLinear = false;
        }
        break;
      case BSplineExponential:
        {
        typedef itk::Vector<RealType, VImageDimension> VectorType;
        VectorType zeroVector( 0.0 );
        typedef itk::Image<VectorType, VImageDimension> ConstantVelocityFieldType;

        typename ConstantVelocityFieldType::Pointer constantVelocityField = AllocImage<ConstantVelocityFieldType>(
            preprocessedFixedImagesPerStage[0], zeroVector );

        typedef itk::BSplineExponentialDiffeomorphicTransform<RealType, VImageDimension>
          BSplineDisplacementFieldTransformType;

        typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, BSplineDisplacementFieldTransformType,
          ImageType, LabeledPointSetType> DisplacementFieldRegistrationType;
        typename DisplacementFieldRegistrationType::Pointer displacementFieldRegistration =
          DisplacementFieldRegistrationType::New();

        if( this->m_RestrictDeformationOptimizerWeights.size() > currentStageNumber )
          {
          if( this->m_RestrictDeformationOptimizerWeights[currentStageNumber].size() == VImageDimension )
            {
            typename DisplacementFieldRegistrationType::OptimizerWeightsType optimizerWeights( VImageDimension );
            for( unsigned int d = 0; d < VImageDimension; d++ )
              {
              optimizerWeights[d] = this->m_RestrictDeformationOptimizerWeights[currentStageNumber][d];
              }
            displacementFieldRegistration->SetOptimizerWeights( optimizerWeights );
            }
          }

        typename BSplineDisplacementFieldTransformType::Pointer outputDisplacementFieldTransform = displacementFieldRegistration->GetModifiableTransform();

        // Create the transform adaptors

        typename DisplacementFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;

        // Extract parameters

        const std::vector<unsigned int> & meshSizeForTheUpdateField =
          this->m_TransformMethods[currentStageNumber].m_UpdateFieldMeshSizeAtBaseLevel;
        std::vector<unsigned int> meshSizeForTheVelocityField =
          this->m_TransformMethods[currentStageNumber].m_VelocityFieldMeshSizeAtBaseLevel;
        unsigned int numberOfIntegrationSteps = this->m_TransformMethods[currentStageNumber].m_NumberOfTimeIndices;

        if( numberOfIntegrationSteps == 0 )
          {
          outputDisplacementFieldTransform->SetCalculateNumberOfIntegrationStepsAutomatically( true );
          }
        else
          {
          outputDisplacementFieldTransform->SetNumberOfIntegrationSteps( numberOfIntegrationSteps );
          }
        outputDisplacementFieldTransform->SetSplineOrder( this->m_TransformMethods[currentStageNumber].m_SplineOrder );
        outputDisplacementFieldTransform->SetConstantVelocityField( constantVelocityField );
        outputDisplacementFieldTransform->SetDisplacementField( constantVelocityField );

        if( meshSizeForTheUpdateField.size() != VImageDimension || meshSizeForTheVelocityField.size() !=
            VImageDimension )
          {
          this->Logger() << "ERROR:  The mesh size(s) don't match the ImageDimension." << std::endl;
          return EXIT_FAILURE;
          }

        typename BSplineDisplacementFieldTransformType::ArrayType updateMeshSize;
        typename BSplineDisplacementFieldTransformType::ArrayType velocityMeshSize;
        for( unsigned int d = 0; d < VImageDimension; d++ )
          {
          updateMeshSize[d] = meshSizeForTheUpdateField[d];
          velocityMeshSize[d] = meshSizeForTheVelocityField[d];
          }
        // Create the transform adaptors specific to B-splines
        for( unsigned int level = 0; level < numberOfLevels; level++ )
          {
          typename itk::ImageBase<VImageDimension>::Pointer shrunkSpace=
                     this->GetShrinkImageOutputInformation(
                          virtualDomainImage.GetPointer(),
                          shrinkFactorsPerDimensionForAllLevels[level]  );


          typedef itk::BSplineExponentialDiffeomorphicTransformParametersAdaptor<BSplineDisplacementFieldTransformType>
            BSplineDisplacementFieldTransformAdaptorType;
          typename BSplineDisplacementFieldTransformAdaptorType::Pointer bsplineFieldTransformAdaptor =
            BSplineDisplacementFieldTransformAdaptorType::New();
          bsplineFieldTransformAdaptor->SetRequiredSpacing( shrunkSpace->GetSpacing() );
          bsplineFieldTransformAdaptor->SetRequiredSize( shrunkSpace->GetLargestPossibleRegion().GetSize() );
          bsplineFieldTransformAdaptor->SetRequiredDirection( shrunkSpace->GetDirection() );
          bsplineFieldTransformAdaptor->SetRequiredOrigin( shrunkSpace->GetOrigin() );
          bsplineFieldTransformAdaptor->SetTransform( outputDisplacementFieldTransform );

          // A good heuristic is to RealType the b-spline mesh resolution at each level
          typename BSplineDisplacementFieldTransformType::ArrayType newUpdateMeshSize = updateMeshSize;
          typename BSplineDisplacementFieldTransformType::ArrayType newVelocityMeshSize = velocityMeshSize;
          for( unsigned int d = 0; d < VImageDimension; d++ )
            {
            newUpdateMeshSize[d] = newUpdateMeshSize[d] << ( level + 1 );
            newVelocityMeshSize[d] = newVelocityMeshSize[d] << ( level + 1 );
            }
          bsplineFieldTransformAdaptor->SetMeshSizeForTheUpdateField( newUpdateMeshSize );
          bsplineFieldTransformAdaptor->SetMeshSizeForTheConstantVelocityField( newVelocityMeshSize );

          adaptors.push_back( bsplineFieldTransformAdaptor.GetPointer() );
          }
        for( unsigned int n = 0; n < stageMetricList.size(); n++ )
          {
          if( !this->IsPointSetMetric( stageMetricList[n].m_MetricType ) )
            {
            displacementFieldRegistration->SetFixedImage( n, preprocessedFixedImagesPerStage[n] );
            displacementFieldRegistration->SetMovingImage( n, preprocessedMovingImagesPerStage[n] );
            }
          else
            {
            displacementFieldRegistration->SetFixedPointSet( n, stageMetricList[n].m_FixedLabeledPointSet.GetPointer() );
            displacementFieldRegistration->SetMovingPointSet( n, stageMetricList[n].m_MovingLabeledPointSet.GetPointer() );
            }
          }
        if( useMultiMetric )
          {
          displacementFieldRegistration->SetMetric( multiMetric );
          }
        else
          {
          displacementFieldRegistration->SetMetric( singleMetric );
          }

        displacementFieldRegistration->SetNumberOfLevels( numberOfLevels );

        for( unsigned int level = 0; level < numberOfLevels; ++level )
          {
          displacementFieldRegistration->SetShrinkFactorsPerDimension( level, shrinkFactorsPerDimensionForAllLevels[level] );
          }
        displacementFieldRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );

        if( this->m_CompositeTransform->GetNumberOfTransforms() > 0 )
          {
          displacementFieldRegistration->SetMovingInitialTransform( this->m_CompositeTransform );
          }
        if( this->m_FixedInitialTransform->GetNumberOfTransforms() > 0 )
          {
          displacementFieldRegistration->SetFixedInitialTransform( this->m_FixedInitialTransform );
          }
        displacementFieldRegistration->SetMetricSamplingStrategy(
          static_cast<typename DisplacementFieldRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
        displacementFieldRegistration->SetMetricSamplingPercentage( samplingPercentage );
        displacementFieldRegistration->SetOptimizer( optimizer2 );
        displacementFieldRegistration->SetTransformParametersAdaptorsPerLevel( adaptors );

        typedef antsRegistrationCommandIterationUpdate<DisplacementFieldRegistrationType> DisplacementFieldCommandType;
        typename DisplacementFieldCommandType::Pointer displacementFieldRegistrationObserver =
          DisplacementFieldCommandType::New();
        displacementFieldRegistrationObserver->SetLogStream(*this->m_LogStream);
        displacementFieldRegistrationObserver->SetNumberOfIterations( currentStageIterations );

        displacementFieldRegistration->AddObserver( itk::IterationEvent(), displacementFieldRegistrationObserver );
        displacementFieldRegistration->AddObserver( itk::InitializeEvent(), displacementFieldRegistrationObserver );

        try
          {
          this->Logger() << std::endl
                         << "*** Running bspline exponential field registration (updateMeshSizeAtBaseLevel = "
                         << updateMeshSize << ", velocityMeshSizeAtBaseLevel = " << velocityMeshSize << ") ***"
                         << std::endl
                         << std::endl;
          displacementFieldRegistrationObserver->Execute( displacementFieldRegistration, itk::StartEvent() );
          displacementFieldRegistration->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          this->Logger() << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }

        // Add calculated transform to the composite transform
        this->m_CompositeTransform->AddTransform( outputDisplacementFieldTransform );

        this->m_AllPreviousTransformsAreLinear = false;
        }
        break;
      case BSpline:
        {
        constexpr unsigned int SplineOrder = 3;
        typedef itk::BSplineTransform<RealType, VImageDimension, SplineOrder> BSplineTransformType;
        typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, BSplineTransformType,
          ImageType, LabeledPointSetType> BSplineRegistrationType;

        typename BSplineRegistrationType::Pointer registrationMethod =
          this->PrepareRegistrationMethod<BSplineRegistrationType>(
                this->m_CompositeTransform, currentStageNumber, VImageDimension,
                preprocessedFixedImagesPerStage, preprocessedMovingImagesPerStage,
                fixedLabeledPointSetsPerStage, movingLabeledPointSetsPerStage, stageMetricList, singleMetric,
                multiMetric, optimizer, numberOfLevels, shrinkFactorsPerDimensionForAllLevels,
                smoothingSigmasPerLevel, metricSamplingStrategy, samplingPercentage );

        typename BSplineTransformType::Pointer outputBSplineTransform = registrationMethod->GetModifiableTransform();

        const std::vector<unsigned int> & size = this->m_TransformMethods[currentStageNumber].m_MeshSizeAtBaseLevel;

        typename BSplineTransformType::PhysicalDimensionsType physicalDimensions;
        typename BSplineTransformType::MeshSizeType meshSize;
        for( unsigned int d = 0; d < VImageDimension; d++ )
          {
          physicalDimensions[d] = static_cast<RealType>( preprocessedFixedImagesPerStage[0]->GetSpacing()[d] )
            * static_cast<RealType>( preprocessedFixedImagesPerStage[0]->GetLargestPossibleRegion().GetSize()[d] - 1 );
          meshSize[d] = size[d];
          }

        // Create the transform adaptors

        typename BSplineRegistrationType::TransformParametersAdaptorsContainerType adaptors;
        // Create the transform adaptors specific to B-splines
        for( unsigned int level = 0; level < numberOfLevels; level++ )
          {
          typename itk::ImageBase<VImageDimension>::Pointer shrunkSpace=
                     this->GetShrinkImageOutputInformation(
                          virtualDomainImage.GetPointer(),
                          shrinkFactorsPerDimensionForAllLevels[level]  );


          // A good heuristic is to RealType the b-spline mesh resolution at each level

          typename BSplineTransformType::MeshSizeType requiredMeshSize;
          for( unsigned int d = 0; d < VImageDimension; d++ )
            {
            requiredMeshSize[d] = meshSize[d] << level;
            }

          typedef itk::BSplineTransformParametersAdaptor<BSplineTransformType> BSplineAdaptorType;
          typename BSplineAdaptorType::Pointer bsplineAdaptor = BSplineAdaptorType::New();
          bsplineAdaptor->SetTransform( outputBSplineTransform );
          bsplineAdaptor->SetRequiredTransformDomainMeshSize( requiredMeshSize );
          bsplineAdaptor->SetRequiredTransformDomainOrigin( shrunkSpace->GetOrigin() );
          bsplineAdaptor->SetRequiredTransformDomainDirection( shrunkSpace->GetDirection() );
          bsplineAdaptor->SetRequiredTransformDomainPhysicalDimensions( physicalDimensions );

          adaptors.push_back( bsplineAdaptor.GetPointer() );
          }

        registrationMethod->SetTransformParametersAdaptorsPerLevel( adaptors );
        outputBSplineTransform->SetTransformDomainOrigin( preprocessedFixedImagesPerStage[0]->GetOrigin() );
        outputBSplineTransform->SetTransformDomainPhysicalDimensions( physicalDimensions );
        outputBSplineTransform->SetTransformDomainMeshSize( meshSize );
        outputBSplineTransform->SetTransformDomainDirection( preprocessedFixedImagesPerStage[0]->GetDirection() );
        outputBSplineTransform->SetIdentity();

        typedef antsRegistrationCommandIterationUpdate<BSplineRegistrationType> BSplineCommandType;
        typename BSplineCommandType::Pointer bsplineObserver = BSplineCommandType::New();
        bsplineObserver->SetLogStream( *this->m_LogStream );
        bsplineObserver->SetNumberOfIterations( currentStageIterations );

        registrationMethod->AddObserver( itk::IterationEvent(), bsplineObserver );
        registrationMethod->AddObserver( itk::InitializeEvent(), bsplineObserver );

        try
          {
          this->Logger() << std::endl << "*** Running bspline registration (meshSizeAtBaseLevel = " << meshSize
                         << ") ***"
                         << std::endl << std::endl;
          bsplineObserver->Execute( registrationMethod, itk::StartEvent() );
          registrationMethod->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          this->Logger() << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }
        // Add calculated transform to the composite transform
        this->m_CompositeTransform->AddTransform( outputBSplineTransform );

        this->m_AllPreviousTransformsAreLinear = false;
        }
        break;
      default:
        this->Logger() << "ERROR:  Unrecognized transform option - " << whichTransform << std::endl;
        return EXIT_FAILURE;
      }
    timer.Stop();
    this->Logger() << "  Elapsed time (stage " << currentStageNumber << "): " << timer.GetMean() << std::endl
                   << std::endl;
    }

  totalTimer.Stop();
  this->Logger() << std::endl << "Total elapsed time: " << totalTimer.GetMean() << std::endl;

  return EXIT_SUCCESS;
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetMovingInitialTransform( const TransformType *initialTransform )
{
  // Since the initial transform might be linear (or a composition of
  // linear transforms), we might want to add those initial transforms
  // to the moving image header for faster processing.

  typename CompositeTransformType::Pointer compToAdd;

  typename CompositeTransformType::ConstPointer compXfrm =
    dynamic_cast<const CompositeTransformType *>( initialTransform );
  if( compXfrm.IsNotNull() )
    {
    compToAdd = compXfrm->Clone();
    this->m_CompositeTransform = compToAdd;
    }
  else
    {
    compToAdd = CompositeTransformType::New();
    typename TransformType::Pointer xfrm = initialTransform->Clone();
    compToAdd->AddTransform( xfrm );
    this->m_CompositeTransform = compToAdd;
    }
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetFixedInitialTransform( const TransformType *initialTransform  )
{
  typename CompositeTransformType::Pointer compToAdd;

  typename CompositeTransformType::ConstPointer compXfrm =
    dynamic_cast<const CompositeTransformType *>( initialTransform );
  if( compXfrm.IsNotNull() )
    {
    compToAdd = compXfrm->Clone();

    this->m_FixedInitialTransform = compToAdd;
    this->m_AllPreviousTransformsAreLinear = false;
    }
  else
    {
    compToAdd = CompositeTransformType::New();
    typename TransformType::Pointer xfrm = initialTransform->Clone();
    compToAdd->AddTransform( xfrm );

    this->m_FixedInitialTransform = compToAdd;
    this->m_AllPreviousTransformsAreLinear = false;
    }
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetRestoreStateTransform( const TransformType *initialTransform )
{
  typename CompositeTransformType::Pointer compToRestore;
  typename CompositeTransformType::Pointer compToAdd;

  typename CompositeTransformType::ConstPointer compXfrm =
  dynamic_cast<const CompositeTransformType *>( initialTransform );
  if( compXfrm.IsNotNull() )
    {
    compToRestore = compXfrm->Clone();

    // If the last four transforms are displacementFieldType, we assume that they are
    // forward and inverse displacement fields of the FixedToMiddle and MovingToMiddle
    // transforms for a SyN registration.
    //
    unsigned int numTransforms = compToRestore->GetNumberOfTransforms();
    if( (compToRestore->GetNthTransform( numTransforms-1 )->GetTransformCategory() == TransformType::TransformCategoryEnum::DisplacementField)
       && (compToRestore->GetNthTransform( numTransforms-2 )->GetTransformCategory() == TransformType::TransformCategoryEnum::DisplacementField)
       && (compToRestore->GetNthTransform( numTransforms-3 )->GetTransformCategory() == TransformType::TransformCategoryEnum::DisplacementField)
       && (compToRestore->GetNthTransform( numTransforms-4 )->GetTransformCategory() == TransformType::TransformCategoryEnum::DisplacementField) )
      {
      typename DisplacementFieldTransformType::Pointer fixedToMiddleForwardTx =
        dynamic_cast<DisplacementFieldTransformType *>( compToRestore->GetNthTransform( numTransforms-4 ).GetPointer() );
      typename DisplacementFieldTransformType::Pointer fixedToMiddleInverseTx =
        dynamic_cast<DisplacementFieldTransformType *>( compToRestore->GetNthTransform( numTransforms-3 ).GetPointer() );
      typename DisplacementFieldTransformType::Pointer movingToMiddleForwardTx =
        dynamic_cast<DisplacementFieldTransformType *>( compToRestore->GetNthTransform( numTransforms-2 ).GetPointer() );
      typename DisplacementFieldTransformType::Pointer movingToMiddleInverseTx =
        dynamic_cast<DisplacementFieldTransformType *>( compToRestore->GetNthTransform( numTransforms-1 ).GetPointer() );

      typename DisplacementFieldTransformType::Pointer fixedToMiddleTransform = DisplacementFieldTransformType::New();
      fixedToMiddleTransform->SetDisplacementField( fixedToMiddleForwardTx->GetModifiableDisplacementField() );
      fixedToMiddleTransform->SetInverseDisplacementField( fixedToMiddleInverseTx->GetModifiableDisplacementField() );

      typename DisplacementFieldTransformType::Pointer movingToMiddleTransform = DisplacementFieldTransformType::New();
      movingToMiddleTransform->SetDisplacementField( movingToMiddleForwardTx->GetModifiableDisplacementField() );
      movingToMiddleTransform->SetInverseDisplacementField( movingToMiddleInverseTx->GetModifiableDisplacementField() );

      this->Logger() << "Initial FixedToMiddle and MovingToMiddle transforms are restored from the registration state file."
                << std::endl;

      compToRestore->RemoveTransform();
      compToRestore->RemoveTransform();
      compToRestore->RemoveTransform();
      compToRestore->RemoveTransform();
      compToRestore->AddTransform( fixedToMiddleTransform );
      compToRestore->AddTransform( movingToMiddleTransform );

      // m_RegistrationState has initial linear transforms + fixedToMiddle + movingToMiddle
      this->m_RegistrationState = compToRestore;

      // Now we restore the SyN transform from FixedToMiddle and MovingToMiddle transforms
      compToAdd = compToRestore->Clone();

      typename DisplacementFieldTransformType::Pointer initialSyNTransform = DisplacementFieldTransformType::New();

      typedef itk::ComposeDisplacementFieldsImageFilter<DisplacementFieldType, DisplacementFieldType> ComposerType;

      typename ComposerType::Pointer composer = ComposerType::New();
      composer->SetDisplacementField( movingToMiddleTransform->GetInverseDisplacementField() );
      composer->SetWarpingField( fixedToMiddleTransform->GetDisplacementField() );
      composer->Update();

      typename ComposerType::Pointer inverseComposer = ComposerType::New();
      inverseComposer->SetDisplacementField( fixedToMiddleTransform->GetInverseDisplacementField() );
      inverseComposer->SetWarpingField( movingToMiddleTransform->GetDisplacementField() );
      inverseComposer->Update();

      initialSyNTransform->SetDisplacementField( composer->GetOutput() );
      initialSyNTransform->SetInverseDisplacementField( inverseComposer->GetOutput() );

      compToAdd->RemoveTransform();
      compToAdd->RemoveTransform();
      compToAdd->AddTransform( initialSyNTransform );
      }
    else
      {
      this->m_RegistrationState = nullptr;
      }

    if( compToAdd.IsNull() )
      {
      compToAdd = compToRestore->Clone();
      }
    // m_CompositeTransform has initial linear transforms + initial SyN transform
    this->m_CompositeTransform = compToAdd;
    }
  else
    {
    this->m_CompositeTransform = nullptr;
    }
}

template <typename TComputeType, unsigned VImageDimension>
std::vector<unsigned int>
RegistrationHelper<TComputeType, VImageDimension>
::CalculateMeshSizeForSpecifiedKnotSpacing( ImageBaseType * const inputImage,
                                            const RealType knotSpacing,
                                            const unsigned int itkNotUsed( splineOrder ) )
{
  // The commented code is for use with itk::ConstantPadImageFilter.  Right now
  // the mesh size is simply an approximation.

  std::vector<unsigned int> meshSize;

//   unsigned long lowerBound[VImageDimension];
//   unsigned long upperBound[VImageDimension];

  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    RealType domain = static_cast<RealType>(
      inputImage->GetLargestPossibleRegion().GetSize()[d] - 1 ) *
      static_cast<RealType>( inputImage->GetSpacing()[d] );
    meshSize.push_back( static_cast<unsigned int>( std::ceil( domain / knotSpacing ) ) );
//     unsigned long extraPadding = static_cast<unsigned long>(
//       ( numberOfSpans * splineDistance - domain ) / inputImage->GetSpacing()[d] + 0.5 );
//     lowerBound[d] = static_cast<unsigned long>( 0.5 * extraPadding );
//     upperBound[d] = extraPadding - lowerBound[d];
//     numberOfControlPoints[d] = meshSize[d] + splineOrder;
    }

  return meshSize;
}

template <typename TComputeType, unsigned VImageDimension>
typename RegistrationHelper<TComputeType, VImageDimension>::AffineTransformType::Pointer
RegistrationHelper<TComputeType, VImageDimension>
::CollapseLinearTransforms( const CompositeTransformType * compositeTransform )
{
  if( !compositeTransform->IsLinear() )
    {
    itkExceptionMacro( "The composite transform is not linear." );
    }

  typename AffineTransformType::Pointer totalTransform = AffineTransformType::New();

  const unsigned int numberOfTransforms = compositeTransform->GetNumberOfTransforms();

  // Find the last transform that has a center, and set that as the fixed parameters of the total transform.
  // It should be set only once.
  for( unsigned int n = numberOfTransforms; n > 0; n--)
    {
    typename TransformType::Pointer transform = compositeTransform->GetNthTransform( n-1 );
    typename MatrixOffsetTransformBaseType::ConstPointer matrixOffsetTransform =
      dynamic_cast<MatrixOffsetTransformBaseType *>( transform.GetPointer() );
    if( matrixOffsetTransform.IsNotNull() )
     {
     totalTransform->SetCenter( matrixOffsetTransform->GetCenter() );
     break;
     }
    }

  typedef itk::TranslationTransform<RealType, VImageDimension> TranslationTransformType;

  for( unsigned int n = 0; n < numberOfTransforms; n++ )
    {
    typename TransformType::Pointer transform = compositeTransform->GetNthTransform( n );

    typename AffineTransformType::Pointer nthTransform = AffineTransformType::New();

    typename TranslationTransformType::Pointer translationTransform =
      dynamic_cast<TranslationTransformType *>( transform.GetPointer() );
    if( translationTransform.IsNotNull() )
      {
      nthTransform->SetOffset( translationTransform->GetOffset() );
      }
    else
      {
      typename MatrixOffsetTransformBaseType::ConstPointer matrixOffsetTransform =
        dynamic_cast<MatrixOffsetTransformBaseType *>( transform.GetPointer() );
      nthTransform->SetCenter( matrixOffsetTransform->GetCenter() );
      nthTransform->SetMatrix( matrixOffsetTransform->GetMatrix() );
      nthTransform->SetTranslation( matrixOffsetTransform->GetTranslation() );
      }
    totalTransform->Compose( nthTransform, true );
    }
  return totalTransform;
}

template <typename TComputeType, unsigned VImageDimension>
typename RegistrationHelper<TComputeType, VImageDimension>::CompositeTransformType::Pointer
RegistrationHelper<TComputeType, VImageDimension>
::CollapseDisplacementFieldTransforms( const CompositeTransformType * compositeTransform )
{
  typename CompositeTransformType::Pointer combinedCompositeTransform = CompositeTransformType::New();

  if( compositeTransform->GetTransformCategory() != TransformType::TransformCategoryEnum::DisplacementField  )
    {
    itkExceptionMacro( "The composite transform is not composed strictly of displacement fields." );
    }

  if( compositeTransform->GetNumberOfTransforms() == 0 )
    {
    itkWarningMacro( "The composite transform is empty.  Returning empty displacement field transform." );
    return combinedCompositeTransform;
    }

  typename TransformType::Pointer transform = compositeTransform->GetNthTransform( 0 );

  typename DisplacementFieldTransformType::Pointer currentTransform =
    dynamic_cast<DisplacementFieldTransformType *>( transform.GetPointer() );

  bool isCurrentTransformInvertible = false;
  if( currentTransform->GetInverseDisplacementField() )
    {
    isCurrentTransformInvertible = true;
    }

  for( unsigned int n = 1; n < compositeTransform->GetNumberOfTransforms(); n++ )
    {
    transform = compositeTransform->GetNthTransform( n );
    typename DisplacementFieldTransformType::Pointer nthTransform =
      dynamic_cast<DisplacementFieldTransformType *>( transform.GetPointer() );

    if( ( isCurrentTransformInvertible && nthTransform->GetInverseDisplacementField() ) ||
        ! ( isCurrentTransformInvertible || nthTransform->GetInverseDisplacementField() ) )
      {
      // Adjacent transforms are the same so we can combine
      typedef itk::ComposeDisplacementFieldsImageFilter<DisplacementFieldType> ComposerType;

      typename ComposerType::Pointer composer = ComposerType::New();
      composer->SetWarpingField( nthTransform->GetDisplacementField() );
      composer->SetDisplacementField( currentTransform->GetDisplacementField() );

      typename DisplacementFieldType::Pointer totalField = composer->GetOutput();
      totalField->Update();
      totalField->DisconnectPipeline();

      typename DisplacementFieldType::Pointer totalInverseField = nullptr;

      if( isCurrentTransformInvertible )
        {
        typename ComposerType::Pointer inverseComposer = ComposerType::New();
        inverseComposer->SetWarpingField( currentTransform->GetInverseDisplacementField() );
        inverseComposer->SetDisplacementField( nthTransform->GetInverseDisplacementField() );

        totalInverseField = inverseComposer->GetOutput();
        totalInverseField->Update();
        totalInverseField->DisconnectPipeline();
        }
      currentTransform->SetDisplacementField( totalField );
      currentTransform->SetInverseDisplacementField( totalInverseField );
      }
    else
      {
      DisplacementFieldTransformPointer displacementFieldTransform = DisplacementFieldTransformType::New();
      displacementFieldTransform->SetDisplacementField( currentTransform->GetModifiableDisplacementField() );
      if( isCurrentTransformInvertible )
        {
        displacementFieldTransform->SetInverseDisplacementField( currentTransform->GetModifiableInverseDisplacementField() );
        }

      combinedCompositeTransform->AddTransform( displacementFieldTransform );

      currentTransform->SetDisplacementField( nthTransform->GetModifiableDisplacementField() );
      currentTransform->SetInverseDisplacementField( nthTransform->GetModifiableInverseDisplacementField() );
      if( currentTransform->GetInverseDisplacementField() )
        {
        isCurrentTransformInvertible = true;
        }
      else
        {
        isCurrentTransformInvertible = false;
        }
      }
    }
  combinedCompositeTransform->AddTransform( currentTransform );

  return combinedCompositeTransform;
}

template <typename TComputeType, unsigned VImageDimension>
typename RegistrationHelper<TComputeType, VImageDimension>::CompositeTransformPointer
RegistrationHelper<TComputeType, VImageDimension>
::CollapseCompositeTransform( const CompositeTransformType * compositeTransform )
{
  CompositeTransformPointer collapsedCompositeTransform = CompositeTransformType::New();

  // Check for the simple cases where the composite transform is composed entirely
  // of linear transforms or displacement field transforms.
  if( compositeTransform->IsLinear() )
    {
    collapsedCompositeTransform->AddTransform( this->CollapseLinearTransforms( compositeTransform ) );
    return collapsedCompositeTransform;
    }
  else if( compositeTransform->GetTransformCategory() == TransformType::TransformCategoryEnum::DisplacementField )
    {
    collapsedCompositeTransform->AddTransform( this->CollapseDisplacementFieldTransforms( compositeTransform ) );
    collapsedCompositeTransform->FlattenTransformQueue();
    return collapsedCompositeTransform;
    }

  // Find the first linear or displacement field transform
  typename TransformType::TransformCategoryEnum currentTransformCategory = TransformType::TransformCategoryEnum::UnknownTransformCategory;
  unsigned int startIndex = 0;
  for( unsigned int n = 0; n < compositeTransform->GetNumberOfTransforms(); n++ )
    {
    typename TransformType::TransformCategoryEnum transformCategory =
      compositeTransform->GetNthTransform( n )->GetTransformCategory();
    if( transformCategory == TransformType::TransformCategoryEnum::Linear || transformCategory == TransformType::TransformCategoryEnum::DisplacementField )
      {
      currentTransformCategory = transformCategory;
      startIndex = n;
      break;
      }
    else
      {
      collapsedCompositeTransform->AddTransform( compositeTransform->GetNthTransform( n ) );
      }
    }

  // If a linear or displacement field transform is found then we can break down the
  // composite transform into neighboring sets of like transform types.
  if( currentTransformCategory != TransformType::TransformCategoryEnum::UnknownTransformCategory )
    {
    CompositeTransformPointer currentCompositeTransform = CompositeTransformType::New();
    currentCompositeTransform->AddTransform( compositeTransform->GetNthTransform( startIndex ) );
    for( unsigned int n = startIndex + 1; n < compositeTransform->GetNumberOfTransforms(); n++ )
      {
      typename TransformType::TransformCategoryEnum transformCategory =
        compositeTransform->GetNthTransform( n )->GetTransformCategory();
      if( transformCategory == currentTransformCategory )
        {
        currentCompositeTransform->AddTransform( compositeTransform->GetNthTransform( n ) );
        if( n == compositeTransform->GetNumberOfTransforms() - 1 )
          {
          if( currentTransformCategory == TransformType::TransformCategoryEnum::Linear )
            {
            collapsedCompositeTransform->AddTransform( this->CollapseLinearTransforms( currentCompositeTransform ) );
            }
          else if( currentTransformCategory == TransformType::TransformCategoryEnum::DisplacementField )
            {
            collapsedCompositeTransform->AddTransform( this->CollapseDisplacementFieldTransforms(
                                                         currentCompositeTransform ) );
            }
          }
        }
      else
        {
        if( currentTransformCategory == TransformType::TransformCategoryEnum::Linear )
          {
          collapsedCompositeTransform->AddTransform( this->CollapseLinearTransforms( currentCompositeTransform ) );
          currentCompositeTransform->ClearTransformQueue();
          }
        else if( currentTransformCategory == TransformType::TransformCategoryEnum::DisplacementField )
          {
          collapsedCompositeTransform->AddTransform( this->CollapseDisplacementFieldTransforms(
                                                       currentCompositeTransform ) );
          currentCompositeTransform->ClearTransformQueue();
          }
        currentTransformCategory = transformCategory;

        if( ( transformCategory == TransformType::TransformCategoryEnum::Linear || transformCategory == TransformType::TransformCategoryEnum::DisplacementField ) &&
            n < compositeTransform->GetNumberOfTransforms() - 1 )
          {
          currentCompositeTransform->AddTransform( compositeTransform->GetNthTransform( n ) );
          }
        else
          {
          collapsedCompositeTransform->AddTransform( compositeTransform->GetNthTransform( n ) );
          }
        }
      }
    }

  collapsedCompositeTransform->FlattenTransformQueue();
  return collapsedCompositeTransform;
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::ApplyCompositeLinearTransformToImageHeader( const CompositeTransformType * compositeTransform,
                                              ImageBaseType * const image,
                                              const bool applyInverse )
{
  if( !compositeTransform->IsLinear() )
    {
    itkExceptionMacro( "The composite transform is not linear.  Cannot collapse it to the image header." );
    }

  typename AffineTransformType::Pointer totalTransform = this->CollapseLinearTransforms( compositeTransform );

  typename ImageType::PointType origin = image->GetOrigin();
  typename ImageType::DirectionType direction = image->GetDirection();

  // Image direction matrix is type of RealType.
  // It should be converted to the current InternalComputationType before it is used to set transfrom parameters.
  vnl_matrix<typename ImageType::DirectionType::ValueType> DoubleLocalDirection( VImageDimension, VImageDimension );
  vnl_matrix<TComputeType> localDirection( VImageDimension, VImageDimension );
  DoubleLocalDirection = direction.GetVnlMatrix();
  vnl_copy( DoubleLocalDirection, localDirection );

  // Image origin is an itk point of type RealType.
  // It should be converted to the current InternalComputationType before it is used to set the offset parameters of transform.
  typename itk::Point<TComputeType, VImageDimension> localOrigin;
  localOrigin.CastFrom(origin);

  typename AffineTransformType::Pointer imageTransform = AffineTransformType::New();
  imageTransform->SetMatrix( localDirection );
  imageTransform->SetOffset( localOrigin.GetVectorFromOrigin() );

  if( applyInverse )
    {
    typename AffineTransformType::Pointer inverseImageTransform = AffineTransformType::New();
    inverseImageTransform->SetMatrix( dynamic_cast<MatrixOffsetTransformBaseType *>( imageTransform->
                                                                                     GetInverseTransform().GetPointer() )
                                      ->GetMatrix() );
    inverseImageTransform->SetOffset( -( inverseImageTransform->GetMatrix() * imageTransform->GetOffset() ) );

    totalTransform->Compose( inverseImageTransform.GetPointer(), false );

    typename AffineTransformType::MatrixType inverseMatrix =
      dynamic_cast<MatrixOffsetTransformBaseType *>( totalTransform->GetInverseTransform().GetPointer() )->GetMatrix();
    typename AffineTransformType::OffsetType inverseOffset = -( inverseMatrix * totalTransform->GetOffset() );
    for( unsigned int d = 0; d < VImageDimension; d++ )
      {
      origin[d] = inverseOffset[d];
      }
    //direction = inverseMatrix; // Does not work because they probably have different types!
    vnl_matrix<TComputeType> localInverseMatrix( VImageDimension, VImageDimension );
    localInverseMatrix = inverseMatrix.GetVnlMatrix();
    vnl_copy(localInverseMatrix, DoubleLocalDirection);
    direction = DoubleLocalDirection;
    }
  else
    {
    totalTransform->Compose( imageTransform, true );

    typename AffineTransformType::MatrixType matrix = totalTransform->GetMatrix();
    typename AffineTransformType::OffsetType offset = totalTransform->GetOffset();
    for( unsigned int d = 0; d < VImageDimension; d++ )
      {
      origin[d] = offset[d];
      }
    //direction = matrix; // Does not work because they probably have different types!
    vnl_matrix<TComputeType> localMatrix( VImageDimension, VImageDimension );
    localMatrix = matrix.GetVnlMatrix();
    vnl_copy(localMatrix, DoubleLocalDirection);
    direction = DoubleLocalDirection;
    }

  image->SetDirection( direction );
  image->SetOrigin( origin );
}

template <typename TComputeType, unsigned VImageDimension>
template <typename TTransformType>
bool
RegistrationHelper<TComputeType, VImageDimension>
::InitializeWithPreviousLinearTransform( const CompositeTransformType * compositeTransform,
                                         const std::string transformTypeName,
                                         typename TTransformType::Pointer & resultTransform )
{
  typedef itk::TranslationTransform<RealType, VImageDimension> TranslationTransformType;
  typedef typename RigidTransformTraits<TComputeType, VImageDimension>::TransformType RigidTransformType;

  std::string previousTxFileType = "";
  const typename TransformType::ConstPointer preTransform = compositeTransform->GetBackTransform();
  if( preTransform.IsNotNull() )
    {
    previousTxFileType = preTransform->GetNameOfClass();
    }
  else
    {
    this->Logger() << "ERROR: INITIALIZATION RETURNS FALSE. Previous Linear Transform is Null" << std::endl;
    return false;
    }
  this->Logger() << "Try to initialize the current " << transformTypeName
            << " from previous " << previousTxFileType << "." << std::endl;
/////
  if( transformTypeName == "Translation" )
    {
    typename TranslationTransformType::Pointer initialTransform =
      dynamic_cast<TranslationTransformType *>(resultTransform.GetPointer());
    initialTransform->SetIdentity();
    if( previousTxFileType == "TranslationTransform" )
      {
      typename TranslationTransformType::ConstPointer tempInitializerTransform =
        dynamic_cast<TranslationTransformType const *>( preTransform.GetPointer() );
      if( tempInitializerTransform.IsNull() )
        {
        this->Logger() << "WARNING: Initialization Failed" << std::endl;
        return false;
        }
      //Translation to Translation
      initialTransform->SetFixedParameters( tempInitializerTransform->GetFixedParameters() );
      initialTransform->SetParameters( tempInitializerTransform->GetParameters() );
      }
    else
      {
      this->Logger() << "WARNING: Initialization Failed" << std::endl;
      return false;
      }
    }
/////
  else if( transformTypeName == "Euler2D" || transformTypeName == "Euler3D" )
    {
    typename RigidTransformType::Pointer initialTransform =
      dynamic_cast<RigidTransformType *>(resultTransform.GetPointer());
    initialTransform->SetIdentity();
    if( previousTxFileType == "TranslationTransform" )
      {
      typename TranslationTransformType::ConstPointer tempInitializerTransform =
        dynamic_cast<TranslationTransformType const *>( preTransform.GetPointer() );
      if( tempInitializerTransform.IsNull() )
        {
        this->Logger() << "WARNING: Initialization Failed" << std::endl;
        return false;
        }
      //Translation to Rigid
      initialTransform->SetOffset( tempInitializerTransform->GetOffset() );
      }
    else if( previousTxFileType == "Euler3DTransform" || previousTxFileType == "Euler2DTransform" )
      {
      typename RigidTransformType::ConstPointer tempInitializerTransform =
        dynamic_cast<RigidTransformType const *>( preTransform.GetPointer() );
      if( tempInitializerTransform.IsNull() )
        {
        this->Logger() << "WARNING: Initialization Failed" << std::endl;
        return false;
        }
      //Rigid to Rigid
      initialTransform->SetFixedParameters( tempInitializerTransform->GetFixedParameters() );
      initialTransform->SetParameters( tempInitializerTransform->GetParameters() );
      }
    else
      {
      this->Logger() << "WARNING: Initialization Failed" << std::endl;
      return false;
      }
    }
/////
  else if( transformTypeName == "Affine" )
    {
    typename AffineTransformType::Pointer initialTransform =
      dynamic_cast<AffineTransformType *>(resultTransform.GetPointer());
    initialTransform->SetIdentity();

    if( previousTxFileType == "TranslationTransform" )
      {
      typename TranslationTransformType::ConstPointer tempInitializerTransform =
        dynamic_cast<TranslationTransformType const *>( preTransform.GetPointer() );
      if( tempInitializerTransform.IsNull() )
        {
        this->Logger() << "WARNING: Initialization Failed" << std::endl;
        return false;
        }
      //Translation to Affine
      initialTransform->SetOffset( tempInitializerTransform->GetOffset() );
      }
    else if( previousTxFileType == "Euler3DTransform" || previousTxFileType == "Euler2DTransform" )
      {
      typename RigidTransformType::ConstPointer tempInitializerTransform =
        dynamic_cast<RigidTransformType const *>( preTransform.GetPointer() );
      if( tempInitializerTransform.IsNull() )
        {
        this->Logger() << "WARNING: Initialization Failed" << std::endl;
        return false;
        }
      //Rigid to Affine
      initialTransform->SetCenter( tempInitializerTransform->GetCenter() );
      initialTransform->SetMatrix( tempInitializerTransform->GetMatrix() );
      initialTransform->SetTranslation( tempInitializerTransform->GetTranslation() );
      }
    else if( previousTxFileType == "AffineTransform" )
      {
      typename AffineTransformType::ConstPointer tempInitializerTransform =
        dynamic_cast<AffineTransformType const *>( preTransform.GetPointer() );
      if( tempInitializerTransform.IsNull() )
        {
        this->Logger() << "WARNING: Initialization Failed" << std::endl;
        return false;
        }
      //Affine to Affine
      initialTransform->SetFixedParameters( tempInitializerTransform->GetFixedParameters() );
      initialTransform->SetParameters( tempInitializerTransform->GetParameters() );
      }
    else
      {
      this->Logger() << "WARNING: Initialization Failed" << std::endl;
      return false;
      }
    }
  else
    {
    this->Logger() << "WARNING: Initialization Failed" << std::endl;
    return false;
    }
/////
  return true; // This function only returns flase or true (NOT FAILURE or SUCCESS).
               // If direct intialization fails, the program should NOT be stopped,
               // because the initial transform will be kept in the composite transform,
               // and the final results will be still correct.
}

template <typename TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::PrintState() const
{
  this->Logger() << "Dimension = " << Self::ImageDimension << std::endl
                 << "Number of stages = " << this->m_NumberOfStages << std::endl
                 << "Use Histogram Matching " << ( this->m_UseHistogramMatching ? "true" : "false" )
                 << std::endl
                 << "Winsorize image intensities "
                 << ( this->m_WinsorizeImageIntensities ? "true" : "false" ) << std::endl
                 << "Lower quantile = " << this->m_LowerQuantile << std::endl
                 << "Upper quantile = " << this->m_UpperQuantile << std::endl;

  for( unsigned i = 0; i < this->m_NumberOfStages; i++ )
    {
    this->Logger() << "Stage " << i + 1 << " State" << std::endl; // NOTE: + 1 for consistency.
    const Metric &          curMetric = this->m_Metrics[i];
    const TransformMethod & curTransform = this->m_TransformMethods[i];

    if( !this->IsPointSetMetric( curMetric.m_MetricType ) )
      {
      this->Logger() << "   Image metric = " << curMetric.GetMetricAsString() << std::endl
                     << "     Fixed image = " << curMetric.m_FixedImage << std::endl
                     << "     Moving image = " << curMetric.m_MovingImage << std::endl
                     << "     Weighting = " << curMetric.m_Weighting << std::endl
                     << "     Sampling strategy = "
                     << (curMetric.m_SamplingStrategy ==
          random ? "random" : (curMetric.m_SamplingStrategy == regular ) ? "regular" : (curMetric.m_SamplingStrategy ==
                                                                                        none ) ? "none" :
          "WARNING: UNKNOWN")
                     << std::endl
                     << "     Number of bins = " << curMetric.m_NumberOfBins << std::endl
                     << "     Radius = " << curMetric.m_Radius << std::endl
                     << "     Sampling percentage  = " << curMetric.m_SamplingPercentage << std::endl;
      }
    else
      {
      if( curMetric.m_MetricType == IGDM )
        {
        this->Logger() << "   Point Set Metric = " << curMetric.GetMetricAsString() << std::endl
                       << "     Fixed intensity point set = " << curMetric.m_FixedIntensityPointSet << std::endl
                       << "     Moving intensity point set = " << curMetric.m_MovingIntensityPointSet << std::endl
                       << "     Weighting = " << curMetric.m_Weighting << std::endl
                       << "     Intensity distance sigma = " << curMetric.m_IntensityDistanceSigma << std::endl
                       << "     Euclidean distance sigma = " << curMetric.m_EuclideanDistanceSigma << std::endl
                       << "     Evaluation K neighborhood = " << curMetric.m_EvaluationKNeighborhood << std::endl;
        }
      else
        {
        this->Logger() << "   Point Set Metric = " << curMetric.GetMetricAsString() << std::endl
                       << "     Fixed labeled point set = " << curMetric.m_FixedLabeledPointSet << std::endl
                       << "     Moving labeled point set = " << curMetric.m_MovingLabeledPointSet << std::endl
                       << "     Weighting = " << curMetric.m_Weighting << std::endl
                       << "     Use only boundary points = " << ( curMetric.m_UseBoundaryPointsOnly ? "true" : "false" ) << std::endl
                       << "     Point set sigma = " << curMetric.m_PointSetSigma << std::endl
                       << "     Evaluation K neighborhood = " << curMetric.m_EvaluationKNeighborhood << std::endl
                       << "     Alpha = " << curMetric.m_Alpha << std::endl
                       << "     Use anisotropic covariances = " << ( curMetric.m_UseAnisotropicCovariances ? "true" : "false" ) << std::endl
                       << "     Sampling percentage = " << curMetric.m_SamplingPercentage << std::endl;
        }
      }
    this->Logger() << "   Transform = " << curTransform.XfrmMethodAsString() << std::endl
                   << "     Gradient step = " << curTransform.m_GradientStep << std::endl
                   << "     Update field sigma (voxel space) = "
                   << curTransform.m_UpdateFieldVarianceInVarianceSpace << std::endl
                   << "     Total field sigma (voxel space) = "
                   << curTransform.m_TotalFieldVarianceInVarianceSpace << std::endl
                   << "     Update field time sigma = " << curTransform.m_UpdateFieldTimeSigma << std::endl
                   << "     Total field time sigma  = " << curTransform.m_TotalFieldTimeSigma << std::endl
                   << "     Number of time indices = " << curTransform.m_NumberOfTimeIndices << std::endl
                   << "     Number of time point samples = " << curTransform.m_NumberOfTimeIndices << std::endl;
    }

}
} // namespace ants

#endif // __itkantsRegistrationHelper_hxx
