#ifndef __itkantsRegistrationHelper_hxx
#define __itkantsRegistrationHelper_hxx

#include <iomanip>
#include <itkAffineTransform.h>

#include "antsRegistrationCommandIterationUpdate.h"
#include "antsRegistrationOptimizerCommandIterationUpdate.h"
#include "antsDisplacementAndVelocityFieldRegistrationCommandIterationUpdate.h"

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_copy.h>

namespace ants
{
/**
 * Transform traits to generalize the rigid transform
 */
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

template <class TComputeType, unsigned int ImageDimension>
class CompositeAffineTransformTraits
{
// Don't worry about the fact that the default option is the
// affine Transform, that one will not actually be instantiated.
public:
  typedef itk::AffineTransform<TComputeType, ImageDimension> TransformType;
};

template <>
class CompositeAffineTransformTraits<double, 2>
{
public:
  typedef itk::ANTSCenteredAffine2DTransform<double> TransformType;
};

template <>
class CompositeAffineTransformTraits<float, 2>
{
public:
typedef itk::ANTSCenteredAffine2DTransform<float> TransformType;
};

template <>
class CompositeAffineTransformTraits<double, 3>
{
public:
  typedef itk::ANTSAffine3DTransform<double> TransformType;
};

template <>
class CompositeAffineTransformTraits<float, 3>
{
public:
typedef itk::ANTSAffine3DTransform<float> TransformType;
};

template <class TComputeType, unsigned VImageDimension>
RegistrationHelper<TComputeType, VImageDimension>
::RegistrationHelper() :
  m_CompositeTransform( NULL ),
  m_FixedInitialTransform( NULL ),
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
  m_ApplyLinearTransformsToFixedImageHeader( true ),
  m_PrintSimilarityMeasureInterval( 0 ),
  m_WriteIntervalVolumes( 0 ),
  m_AllPreviousTransformsAreLinear( true ),
  m_CompositeLinearTransformForFixedImageHeader( NULL )
{
  typedef itk::LinearInterpolateImageFunction<ImageType, RealType> LinearInterpolatorType;
  typename LinearInterpolatorType::Pointer linearInterpolator = LinearInterpolatorType::New();
  this->m_Interpolator = linearInterpolator;
}

template <class TComputeType, unsigned VImageDimension>
RegistrationHelper<TComputeType, VImageDimension>
::~RegistrationHelper()
{
}

template <class ImageType>
typename ImageType::Pointer PreprocessImage( typename ImageType::ConstPointer  inputImage,
                                             typename ImageType::PixelType lowerScaleValue,
                                             typename ImageType::PixelType upperScaleValue,
                                             float winsorizeLowerQuantile, float winsorizeUpperQuantile,
                                             typename ImageType::ConstPointer histogramMatchSourceImage = NULL )
{
  typedef itk::Statistics::ImageToHistogramFilter<ImageType>   HistogramFilterType;
  typedef typename HistogramFilterType::InputBooleanObjectType InputBooleanObjectType;
  typedef typename HistogramFilterType::HistogramSizeType      HistogramSizeType;
  typedef typename HistogramFilterType::HistogramType          HistogramType;

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

  typename ImageType::Pointer outputImage = NULL;
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

template <class TComputeType, unsigned VImageDimension>
typename RegistrationHelper<TComputeType, VImageDimension>::MetricEnumeration
RegistrationHelper<TComputeType, VImageDimension>
::StringToMetricType(const std::string & str) const
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
  return IllegalMetric;
}

template <class TComputeType, unsigned VImageDimension>
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

template <class TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddMetric( MetricEnumeration metricType,
             typename ImageType::Pointer & fixedImage,
             typename ImageType::Pointer & movingImage,
             unsigned int stageID,
             RealType weighting,
             SamplingStrategy samplingStrategy,
             int numberOfBins,
             unsigned int  radius,
             RealType samplingPercentage )
{
  Metric init( metricType, fixedImage, movingImage, stageID,
               weighting, samplingStrategy, numberOfBins,
               radius,
               samplingPercentage );

  this->m_Metrics.push_back( init );
}

template <class TComputeType, unsigned VImageDimension>
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

template <class TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddRigidTransform(RealType GradientStep)
{
  TransformMethod init;

  init.m_XfrmMethod = Rigid;
  init.m_GradientStep = GradientStep;
  this->m_TransformMethods.push_back( init );
}

template <class TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddAffineTransform(RealType GradientStep)
{
  TransformMethod init;

  init.m_XfrmMethod = Affine;
  init.m_GradientStep = GradientStep;
  this->m_TransformMethods.push_back( init );
}

template <class TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddCompositeAffineTransform(RealType GradientStep)
{
  TransformMethod init;

  init.m_XfrmMethod = CompositeAffine;
  init.m_GradientStep = GradientStep;
  this->m_TransformMethods.push_back( init );
}

template <class TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddSimilarityTransform(RealType GradientStep)
{
  TransformMethod init;

  init.m_XfrmMethod = Similarity;
  init.m_GradientStep = GradientStep;
  this->m_TransformMethods.push_back( init );
}

template <class TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::AddTranslationTransform(RealType GradientStep)
{
  TransformMethod init;

  init.m_XfrmMethod = Translation;
  init.m_GradientStep = GradientStep;
  this->m_TransformMethods.push_back( init );
}

template <class TComputeType, unsigned VImageDimension>
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

template <class TComputeType, unsigned VImageDimension>
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

template <class TComputeType, unsigned VImageDimension>
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

template <class TComputeType, unsigned VImageDimension>
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

template <class TComputeType, unsigned VImageDimension>
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

template <class TComputeType, unsigned VImageDimension>
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

template <class TComputeType, unsigned VImageDimension>
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

template <class TComputeType, unsigned VImageDimension>
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

template <class TComputeType, unsigned VImageDimension>
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

template <class TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetIterations( const std::vector<std::vector<unsigned int> > & Iterations )
{
  this->m_Iterations = Iterations;
}

template <class TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetConvergenceThresholds( const std::vector<RealType> & thresholds )
{
  this->m_ConvergenceThresholds = thresholds;
}

template <class TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetConvergenceWindowSizes( const std::vector<unsigned int> & windowSizes )
{
  this->m_ConvergenceWindowSizes = windowSizes;
}

template <class TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetSmoothingSigmas( const std::vector<std::vector<float> > & SmoothingSigmas )
{
  this->m_SmoothingSigmas = SmoothingSigmas;
}

template <class TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetRestrictDeformationOptimizerWeights( const std::vector<RealType> & restrictDeformationWeights )
{
  this->m_RestrictDeformationOptimizerWeights = restrictDeformationWeights;
}

template <class TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetSmoothingSigmasAreInPhysicalUnits( const std::vector<bool> & SmoothingSigmasAreInPhysicalUnits )
{
  this->m_SmoothingSigmasAreInPhysicalUnits = SmoothingSigmasAreInPhysicalUnits;
}

template <class TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetShrinkFactors( const std::vector<std::vector<unsigned int> > & ShrinkFactors )
{
  this->m_ShrinkFactors = ShrinkFactors;
}

template <class TComputeType, unsigned VImageDimension>
typename RegistrationHelper<TComputeType, VImageDimension>::ShrinkFactorsPerDimensionContainerType
RegistrationHelper<TComputeType, VImageDimension>
::CalculateShrinkFactorsPerDimension( unsigned int factor, ImagePointer image )
{
  typedef typename ImageType::SpacingType SpacingType;
  typedef RealType                          SpacingValueType;

  SpacingType spacing = image->GetSpacing();

  SpacingValueType minSpacing = spacing[0];
  unsigned int minIndex = 0;
  for( unsigned int n = 1; n < VImageDimension; n++ )
    {
    if( minSpacing > spacing[n] )
      {
      minSpacing = spacing[n];
      minIndex = n;
      }
    }

  ShrinkFactorsPerDimensionContainerType shrinkFactorsPerDimension;
  shrinkFactorsPerDimension.Fill( 0 );
  shrinkFactorsPerDimension[minIndex] = factor;

  SpacingType newSpacing;
  newSpacing[minIndex] = spacing[minIndex] * factor;

  for( unsigned int n = 0; n < VImageDimension; n++ )
    {
    if( shrinkFactorsPerDimension[n] == 0 )
      {
      SpacingValueType newMinSpacing = spacing[n] * static_cast<SpacingValueType>( factor );
      RealType minDifferenceFromMinSpacing = vnl_math_abs( newMinSpacing - newSpacing[minIndex] );
      unsigned int minFactor = factor;
      for( unsigned int f = factor - 1; f > 0; f-- )
        {
        newMinSpacing = spacing[n] * static_cast<SpacingValueType>( f );

        // We use <= such that the smaller factor is preferred if distances are the same
        if( vnl_math_abs( newMinSpacing - newSpacing[minIndex] ) <= minDifferenceFromMinSpacing )
          {
          minDifferenceFromMinSpacing = vnl_math_abs( newMinSpacing - newSpacing[minIndex] );
          minFactor = f;
          }
        }
      shrinkFactorsPerDimension[n] = minFactor;
      }
    }
  return shrinkFactorsPerDimension;
}

template <class TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetWinsorizeImageIntensities( bool Winsorize, float LowerQuantile, float UpperQuantile )
{
  this->m_WinsorizeImageIntensities = Winsorize;
  this->m_LowerQuantile = LowerQuantile;
  this->m_UpperQuantile = UpperQuantile;
}

template <class TComputeType, unsigned VImageDimension>
int
RegistrationHelper<TComputeType, VImageDimension>
::ValidateParameters()
{
  if( this->m_NumberOfStages == 0 )
    {
    std::cout << "No transformations are specified." << std::endl;
    return EXIT_FAILURE;
    }
  if( this->m_Iterations.size() != this->m_NumberOfStages )
    {
    std::cout << "The number of iteration sets specified does not match the number of stages." << std::endl;
    return EXIT_FAILURE;
    }
  if( this->m_ShrinkFactors.size() != this->m_NumberOfStages )
    {
    std::cout << "The number of shrinkFactors specified does not match the number of stages." << std::endl;
    return EXIT_FAILURE;
    }
  if( this->m_SmoothingSigmas.size() != this->m_NumberOfStages )
    {
    std::cout << "The number of smoothing sigma sets specified does not match the number of stages."
                     << std::endl;
    return EXIT_FAILURE;
    }
  if( this->m_SmoothingSigmasAreInPhysicalUnits.size() != this->m_NumberOfStages )
    {
    std::cout
      << "The number of smoothing sigma in physical units bool values does not match the number of stages."
      << std::endl;
    return EXIT_FAILURE;
    }
  for( unsigned int i = 0; i < this->m_Metrics.size(); i++ )
    {
    if( this->m_Metrics[i].m_FixedImage.IsNull() ||
        this->m_Metrics[i].m_MovingImage.IsNull() )
      {
      std::cout << "Must either add Metrics with filenames, or pointers to images" << std::endl;
      return EXIT_FAILURE;
      }
    }
  return EXIT_SUCCESS;
}

template <class TComputeType, unsigned VImageDimension>
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

template <class TComputeType, unsigned VImageDimension>
typename RegistrationHelper<TComputeType, VImageDimension>::ImageType::Pointer
RegistrationHelper<TComputeType, VImageDimension>
::GetInverseWarpedImage() const
{
  typename ImageType::Pointer fixedImage = this->m_Metrics[0].m_FixedImage;
  typename ImageType::Pointer movingImage = this->m_Metrics[0].m_MovingImage;

  if( this->m_CompositeTransform->GetInverseTransform().IsNull() )
    {
    return 0;
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

template <class TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetFixedImageMask(typename MaskImageType::Pointer & fixedImageMask)
{
  typename ImageMaskSpatialObjectType::Pointer so =
    ImageMaskSpatialObjectType::New();
  so->SetImage( fixedImageMask.GetPointer() );
  this->SetFixedImageMask(so);
}

template <class TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::SetMovingImageMask(typename MaskImageType::Pointer & movingImageMask)
{
  typename ImageMaskSpatialObjectType::Pointer so =
    ImageMaskSpatialObjectType::New();
  so->SetImage( movingImageMask.GetPointer() );
  this->SetMovingImageMask(so);
}

template <class TComputeType, unsigned VImageDimension>
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
  if( this->m_CompositeLinearTransformForFixedImageHeader.IsNull() )
    {
    this->m_CompositeLinearTransformForFixedImageHeader = CompositeTransformType::New();
    }
  if( this->m_FixedInitialTransform.IsNull() )
    {
    this->m_FixedInitialTransform = CompositeTransformType::New();
    }

  // ########################################################################################
  // ########################################################################################
  // ##The main loop for exstimating the total composite tranform
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

    // Get the number of metrics at the current stage.  If more than one metric
    // then we need to use the MultiMetricType.  Due to the way the metrics are
    // pulled off the command line stack, we need to iterate from the top down.

    MetricListType stageMetricList = this->GetMetricListPerStage( this->m_NumberOfStages - currentStageNumber - 1 );

    typename MetricType::Pointer      singleMetric;
    typename MultiMetricType::Pointer multiMetric;

    typename MultiMetricType::WeightsArrayType metricWeights( stageMetricList.size() );
    metricWeights.Fill( 1.0 );

    bool useMultiMetric = false;
    if( stageMetricList.size() > 1 )
      {
      useMultiMetric = true;
      multiMetric = MultiMetricType::New();
      }

    // Get shrink factors and adjust according to the current image
    const std::vector<unsigned int> factors( this->m_ShrinkFactors[currentStageNumber] );
    if( factors.size() != numberOfLevels )
      {
      std::cout << "\n\n\n"
                       << "ERROR:  The number of shrink factors does not match the number of levels."
                       << "\nShrink Factors: " << factors.size()
                       << "\nNumber Of Levels: " << numberOfLevels
                       << "\n\n\n"
                       << std::endl;
      return EXIT_FAILURE;
      }

    std::vector<ShrinkFactorsPerDimensionContainerType> shrinkFactorsPerDimensionForAllLevels;
    for( unsigned int n = 0; n < numberOfLevels; n++ )
      {
      ShrinkFactorsPerDimensionContainerType shrinkFactorsPerDimension =
        this->CalculateShrinkFactorsPerDimension( factors[n], stageMetricList[0].m_FixedImage.GetPointer() );
      shrinkFactorsPerDimensionForAllLevels.push_back( shrinkFactorsPerDimension );
      this->Logger() << "  Shrink factors (level " << n+1 << " out of " << numberOfLevels << "): " << shrinkFactorsPerDimension << std::endl;
      }

    // Get smoothing sigmas
    const std::vector<float> sigmas( this->m_SmoothingSigmas[currentStageNumber] );
    typename AffineRegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
    smoothingSigmasPerLevel.SetSize( sigmas.size() );

    if( sigmas.size() != numberOfLevels )
      {
      std::cout << "ERROR:  The number of smoothing sigmas "
                       << "does not match the number of levels." << std::endl;
      return EXIT_FAILURE;
      }
    for( unsigned int n = 0; n < smoothingSigmasPerLevel.Size(); n++ )
      {
      smoothingSigmasPerLevel[n] = sigmas[n];
      }
    this->Logger() << "  smoothing sigmas per level: " << smoothingSigmasPerLevel << std::endl;

    std::vector<typename ImageType::Pointer> preprocessedFixedImagesPerStage;
    std::vector<typename ImageType::Pointer> preprocessedMovingImagesPerStage;

    for( unsigned int currentMetricNumber = 0; currentMetricNumber < stageMetricList.size(); currentMetricNumber++ )
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
                                    NULL );

      preprocessedFixedImagesPerStage.push_back( preprocessFixedImage.GetPointer() );

      typename ImageType::Pointer preprocessMovingImage =
        PreprocessImage<ImageType>( movingImage.GetPointer(), lowerScaleValue,
                                    upperScaleValue, this->m_LowerQuantile, this->m_UpperQuantile,
                                    NULL );

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

      if( this->m_ApplyLinearTransformsToFixedImageHeader )
        {
        this->ApplyCompositeLinearTransformToImageHeader( this->m_CompositeLinearTransformForFixedImageHeader,
                                                          dynamic_cast<ImageBaseType *>( preprocessFixedImage.
                                                                                         GetPointer() ), false );

        if( this->m_FixedImageMask.IsNotNull() )
          {
          this->ApplyCompositeLinearTransformToImageHeader( this->m_CompositeLinearTransformForFixedImageHeader,
                                                            dynamic_cast<ImageBaseType *>( const_cast<MaskImageType *>(
                                                                                             this->m_FixedImageMask->
                                                                                             GetImage() ) ), false );
          }
        }

      this->Logger() << outputPreprocessingString << std::flush;

      // Set up the image metric and scales estimator

      typename MetricType::Pointer metric;

      switch( stageMetricList[currentMetricNumber].m_MetricType )
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

          metric = correlationMetric;
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
          mutualInformationMetric = mutualInformationMetric;
          mutualInformationMetric->SetNumberOfHistogramBins( binOption );
          mutualInformationMetric->SetUseMovingImageGradientFilter( gradientfilter );
          mutualInformationMetric->SetUseFixedImageGradientFilter( gradientfilter );
          mutualInformationMetric->SetUseFixedSampledPointSet( false );
          metric = mutualInformationMetric;
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
          mutualInformationMetric = mutualInformationMetric;
          mutualInformationMetric->SetNumberOfHistogramBins( binOption );
          mutualInformationMetric->SetUseMovingImageGradientFilter( gradientfilter );
          mutualInformationMetric->SetUseFixedImageGradientFilter( gradientfilter );
          mutualInformationMetric->SetUseFixedSampledPointSet( false );
          mutualInformationMetric->SetVarianceForJointPDFSmoothing( 1.0 );
          metric = mutualInformationMetric;
          }
          break;
        case MeanSquares:
          {
          this->Logger() << "  using the MeanSquares metric (weight = "
                         << stageMetricList[currentMetricNumber].m_Weighting << ")" << std::endl;

          typedef itk::MeanSquaresImageToImageMetricv4<ImageType, ImageType, ImageType, TComputeType> MeanSquaresMetricType;
          typename MeanSquaresMetricType::Pointer meanSquaresMetric = MeanSquaresMetricType::New();
          meanSquaresMetric = meanSquaresMetric;
          metric = meanSquaresMetric;
          }
          break;
        case Demons:
          {
          this->Logger() << "  using the Demons metric (weight = "
                         << stageMetricList[currentMetricNumber].m_Weighting << ")" << std::endl;

          typedef itk::DemonsImageToImageMetricv4<ImageType, ImageType, ImageType, TComputeType> DemonsMetricType;
          typename DemonsMetricType::Pointer demonsMetric = DemonsMetricType::New();
          demonsMetric = demonsMetric;
          metric = demonsMetric;
          }
          break;
        case GC:
          {
          this->Logger() << "  using the global correlation metric (weight = "
                         << stageMetricList[currentMetricNumber].m_Weighting << ")" << std::endl;
          typedef itk::CorrelationImageToImageMetricv4<ImageType, ImageType, ImageType, TComputeType> corrMetricType;
          typename corrMetricType::Pointer corrMetric = corrMetricType::New();
          metric = corrMetric;
          }
          break;
        default:
          std::cout << "ERROR: Unrecognized image metric: " << std::endl;
        }
      metric->SetVirtualDomainFromImage( fixedImage );
      metric->SetUseMovingImageGradientFilter( gradientfilter );
      metric->SetUseFixedImageGradientFilter( gradientfilter );
      metricWeights[currentMetricNumber] = stageMetricList[currentMetricNumber].m_Weighting;
      if( this->m_FixedImageMask.IsNotNull() )
        {
        metric->SetFixedImageMask( this->m_FixedImageMask );
        }
      if( this->m_MovingImageMask.IsNotNull() )
        {
        metric->SetMovingImageMask( this->m_MovingImageMask );
        }

      if( useMultiMetric )
        {
        multiMetric->AddMetric( metric );
        }
      if( !useMultiMetric || currentMetricNumber == 0 )
        {
        singleMetric = metric;
        }
      }
    if( useMultiMetric )
      {
      multiMetric->SetMetricWeights( metricWeights );
      }

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
      std::cout << "  Using default NONE metricSamplingStrategy " << std::endl;
      }
    else
      {
      std::cout << "ERROR: samplingStrategy is incorrectly specified" << std::endl;
      exit( -1 );
      }

    // Set up the optimizers.  To change the iteration number for each level we rely
    // on the command observer.

    const RealType learningRate = this->m_TransformMethods[currentStageNumber].m_GradientStep;

    // There's a scale issue here.  Currently we are using the first metric to estimate the
    // scales but we might need to change this.
    typedef itk::RegistrationParameterScalesFromPhysicalShift<MetricType> ScalesEstimatorType;

    typename ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();
    scalesEstimator->SetMetric( singleMetric );
    scalesEstimator->SetTransformForward( true );

    typedef itk::ConjugateGradientLineSearchOptimizerv4Template<TComputeType> ConjugateGradientDescentOptimizerType;
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
    optimizerObserver->SetOrigFixedImage( this->m_Metrics[0].m_FixedImage );
    optimizerObserver->SetOrigMovingImage( this->m_Metrics[0].m_MovingImage );
    if( this->m_PrintSimilarityMeasureInterval != 0 )
      {
      optimizerObserver->SetComputeFullScaleCCInterval( this->m_PrintSimilarityMeasureInterval );
      }
    if( this->m_WriteIntervalVolumes != 0 )
      {
      optimizerObserver->SetWriteInterationsOutputsInIntervals( this->m_WriteIntervalVolumes );
      optimizerObserver->SetCurrentStageNumber( currentStageNumber );
      }

    typedef itk::GradientDescentLineSearchOptimizerv4Template<TComputeType> GradientDescentLSOptimizerType;
    typedef itk::GradientDescentOptimizerv4Template<TComputeType>           GradientDescentOptimizerType;
    typename GradientDescentOptimizerType::Pointer optimizer2 = GradientDescentOptimizerType::New();
    //    optimizer2->SetLowerLimit( 0 );
    //    optimizer2->SetUpperLimit( 2 );
    //    optimizer2->SetEpsilon( 0.2 );
    //    optimizer->SetMaximumLineSearchIterations( 20 );
    optimizer2->SetLearningRate( learningRate );
    optimizer2->SetMaximumStepSizeInPhysicalUnits( learningRate );
    optimizer2->SetNumberOfIterations( currentStageIterations[0] );
    optimizer2->SetScalesEstimator( NULL );
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
    optimizerObserver2->SetOrigFixedImage( this->m_Metrics[0].m_FixedImage );
    optimizerObserver2->SetOrigMovingImage( this->m_Metrics[0].m_MovingImage );
    if( this->m_PrintSimilarityMeasureInterval != 0 )
      {
      optimizerObserver2->SetComputeFullScaleCCInterval( this->m_PrintSimilarityMeasureInterval );
      }
    if( this->m_WriteIntervalVolumes != 0 )
      {
      optimizerObserver2->SetWriteInterationsOutputsInIntervals( this->m_WriteIntervalVolumes );
      optimizerObserver2->SetCurrentStageNumber( currentStageNumber );
      }

    // Set up the image registration methods along with the transforms
    const XfrmMethod whichTransform( this->m_TransformMethods[currentStageNumber].m_XfrmMethod );

    switch( whichTransform )
      {
      case Affine:
        {
        typename AffineRegistrationType::Pointer affineRegistration = AffineRegistrationType::New();
        for( unsigned int n = 0; n < stageMetricList.size(); n++ )
          {
          affineRegistration->SetFixedImage( n, preprocessedFixedImagesPerStage[n] );
          affineRegistration->SetMovingImage( n, preprocessedMovingImagesPerStage[n] );
          }
        if( useMultiMetric )
          {
          affineRegistration->SetMetric( multiMetric );
          }
        else
          {
          affineRegistration->SetMetric( singleMetric );
          }

        affineRegistration->SetNumberOfLevels( numberOfLevels );
        for( unsigned int level = 0; level < numberOfLevels; ++level )
          {
          affineRegistration->SetShrinkFactorsPerDimension( level, shrinkFactorsPerDimensionForAllLevels[level] );
          }
        affineRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
        affineRegistration->SetSmoothingSigmasAreSpecifiedInPhysicalUnits( this->m_SmoothingSigmasAreInPhysicalUnits[
                                                                             currentStageNumber] );
        affineRegistration->SetMetricSamplingStrategy( metricSamplingStrategy );
        affineRegistration->SetMetricSamplingPercentage( samplingPercentage );
        affineRegistration->SetOptimizer( optimizer );
        if( this->m_CompositeTransform->GetNumberOfTransforms() > 0 )
          {
          affineRegistration->SetMovingInitialTransform( this->m_CompositeTransform );
          }
        if( this->m_FixedInitialTransform->GetNumberOfTransforms() > 0 )
          {
          affineRegistration->SetFixedInitialTransform( this->m_FixedInitialTransform );
          }

        typedef antsRegistrationCommandIterationUpdate<AffineRegistrationType> AffineCommandType;
        typename AffineCommandType::Pointer affineObserver = AffineCommandType::New();
        affineObserver->SetLogStream( *this->m_LogStream );
        affineObserver->SetNumberOfIterations( currentStageIterations );
        affineRegistration->AddObserver( itk::IterationEvent(), affineObserver );
        affineRegistration->AddObserver( itk::InitializeEvent(), affineObserver );
        try
          {
          this->Logger() << std::endl << "*** Running affine registration ***" << std::endl << std::endl;
          affineObserver->Execute( affineRegistration, itk::StartEvent() );
          affineRegistration->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cout << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }

        // Add calculated transform to the composite transform or add it to the composite transform
        // which is incorporated into the fixed image header.
        if( this->m_ApplyLinearTransformsToFixedImageHeader && this->m_AllPreviousTransformsAreLinear )
          {
          this->m_CompositeLinearTransformForFixedImageHeader->AddTransform( const_cast<AffineTransformType *>(
                                                                               affineRegistration->GetOutput()->Get() ) );
          }
        else
          {
          this->m_CompositeTransform->AddTransform( const_cast<AffineTransformType *>(
                                                      affineRegistration->GetOutput()->Get() ) );
          }
        }
        break;
      case Rigid:
        {
        typedef typename RigidTransformTraits<TComputeType, VImageDimension>::TransformType RigidTransformType;

        typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, RigidTransformType> RigidRegistrationType;
        typename RigidRegistrationType::Pointer rigidRegistration = RigidRegistrationType::New();
        for( unsigned int n = 0; n < stageMetricList.size(); n++ )
          {
          rigidRegistration->SetFixedImage( n, preprocessedFixedImagesPerStage[n] );
          rigidRegistration->SetMovingImage( n, preprocessedMovingImagesPerStage[n] );
          }
        if( useMultiMetric )
          {
          rigidRegistration->SetMetric( multiMetric );
          }
        else
          {
          rigidRegistration->SetMetric( singleMetric );
          }

        rigidRegistration->SetNumberOfLevels( numberOfLevels );
        for( unsigned int level = 0; level < numberOfLevels; ++level )
          {
          rigidRegistration->SetShrinkFactorsPerDimension( level, shrinkFactorsPerDimensionForAllLevels[level] );
          }
        rigidRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
        rigidRegistration->SetSmoothingSigmasAreSpecifiedInPhysicalUnits( this->m_SmoothingSigmasAreInPhysicalUnits[
                                                                            currentStageNumber] );
        rigidRegistration->SetMetricSamplingStrategy(
          static_cast<typename RigidRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
        rigidRegistration->SetMetricSamplingPercentage( samplingPercentage );
        rigidRegistration->SetOptimizer( optimizer );
        if( this->m_CompositeTransform->GetNumberOfTransforms() > 0 )
          {
          rigidRegistration->SetMovingInitialTransform( this->m_CompositeTransform );
          }
        if( this->m_FixedInitialTransform->GetNumberOfTransforms() > 0 )
          {
          rigidRegistration->SetFixedInitialTransform( this->m_FixedInitialTransform );
          }

        typedef antsRegistrationCommandIterationUpdate<RigidRegistrationType> RigidCommandType;
        typename RigidCommandType::Pointer rigidObserver = RigidCommandType::New();
        rigidObserver->SetLogStream(*this->m_LogStream);
        rigidObserver->SetNumberOfIterations( currentStageIterations );

        rigidRegistration->AddObserver( itk::IterationEvent(), rigidObserver );
        rigidRegistration->AddObserver( itk::InitializeEvent(), rigidObserver );

        try
          {
          this->Logger() << std::endl << "*** Running rigid registration ***" << std::endl << std::endl;
          rigidObserver->Execute( rigidRegistration, itk::StartEvent() );
          rigidRegistration->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cout << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }

        // Add calculated transform to the composite transform or add it to the composite transform
        // which is incorporated into the fixed image header.
        if( this->m_ApplyLinearTransformsToFixedImageHeader && this->m_AllPreviousTransformsAreLinear )
          {
          this->m_CompositeLinearTransformForFixedImageHeader->AddTransform( const_cast<RigidTransformType *>(
                                                                               rigidRegistration->GetOutput()->Get() ) );
          }
        else
          {
          this->m_CompositeTransform->AddTransform( const_cast<RigidTransformType *>( rigidRegistration->GetOutput()->
                                                                                      Get() ) );
          }
        }
        break;
      case CompositeAffine:
        {
        typedef typename CompositeAffineTransformTraits<TComputeType, VImageDimension>::TransformType CompositeAffineTransformType;

        typedef itk::ImageRegistrationMethodv4<ImageType, ImageType,
                                               CompositeAffineTransformType> CompositeAffineRegistrationType;
        typename CompositeAffineRegistrationType::Pointer affineRegistration = CompositeAffineRegistrationType::New();
        for( unsigned int n = 0; n < stageMetricList.size(); n++ )
          {
          affineRegistration->SetFixedImage( n, preprocessedFixedImagesPerStage[n] );
          affineRegistration->SetMovingImage( n, preprocessedMovingImagesPerStage[n] );
          }
        if( useMultiMetric )
          {
          affineRegistration->SetMetric( multiMetric );
          }
        else
          {
          affineRegistration->SetMetric( singleMetric );
          }

        affineRegistration->SetNumberOfLevels( numberOfLevels );
        for( unsigned int level = 0; level < numberOfLevels; ++level )
          {
          affineRegistration->SetShrinkFactorsPerDimension( level, shrinkFactorsPerDimensionForAllLevels[level] );
          }
        affineRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
        affineRegistration->SetSmoothingSigmasAreSpecifiedInPhysicalUnits( this->m_SmoothingSigmasAreInPhysicalUnits[
                                                                             currentStageNumber] );
        affineRegistration->SetMetricSamplingStrategy(
          static_cast<typename CompositeAffineRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
        affineRegistration->SetMetricSamplingPercentage( samplingPercentage );
        affineRegistration->SetOptimizer( optimizer );
        if( this->m_CompositeTransform->GetNumberOfTransforms() > 0 )
          {
          affineRegistration->SetMovingInitialTransform( this->m_CompositeTransform );
          }
        if( this->m_FixedInitialTransform->GetNumberOfTransforms() > 0 )
          {
          affineRegistration->SetFixedInitialTransform( this->m_FixedInitialTransform );
          }

        typedef antsRegistrationCommandIterationUpdate<AffineRegistrationType> AffineCommandType;
        typename AffineCommandType::Pointer affineObserver = AffineCommandType::New();
        affineObserver->SetLogStream(*this->m_LogStream);
        affineObserver->SetNumberOfIterations( currentStageIterations );

        affineRegistration->AddObserver( itk::IterationEvent(), affineObserver );
        affineRegistration->AddObserver( itk::InitializeEvent(), affineObserver );

        try
          {
          this->Logger() << std::endl << "*** Running composite affine registration ***" << std::endl << std::endl;
          affineObserver->Execute( affineRegistration, itk::StartEvent() );
          affineRegistration->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cout << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }

        // Add calculated transform to the composite transform or add it to the composite transform
        // which is incorporated into the fixed image header.
        if( this->m_ApplyLinearTransformsToFixedImageHeader && this->m_AllPreviousTransformsAreLinear )
          {
          this->m_CompositeLinearTransformForFixedImageHeader->AddTransform( const_cast<CompositeAffineTransformType *>(
                                                                               affineRegistration->GetOutput()->Get() ) );
          }
        else
          {
          this->m_CompositeTransform->AddTransform( const_cast<CompositeAffineTransformType *>( affineRegistration->
                                                                                                GetOutput()->Get() ) );
          }
        }
        break;
      case Similarity:
        {
        typedef typename SimilarityTransformTraits<TComputeType, VImageDimension>::TransformType SimilarityTransformType;

        typedef itk::ImageRegistrationMethodv4<ImageType, ImageType,
                                               SimilarityTransformType> SimilarityRegistrationType;
        typename SimilarityRegistrationType::Pointer similarityRegistration = SimilarityRegistrationType::New();
        for( unsigned int n = 0; n < stageMetricList.size(); n++ )
          {
          similarityRegistration->SetFixedImage( n, preprocessedFixedImagesPerStage[n] );
          similarityRegistration->SetMovingImage( n, preprocessedMovingImagesPerStage[n] );
          }
        if( useMultiMetric )
          {
          similarityRegistration->SetMetric( multiMetric );
          }
        else
          {
          similarityRegistration->SetMetric( singleMetric );
          }

        similarityRegistration->SetNumberOfLevels( numberOfLevels );
        for( unsigned int level = 0; level < numberOfLevels; ++level )
          {
          similarityRegistration->SetShrinkFactorsPerDimension( level, shrinkFactorsPerDimensionForAllLevels[level] );
          }
        similarityRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
        similarityRegistration->SetSmoothingSigmasAreSpecifiedInPhysicalUnits(
          this->m_SmoothingSigmasAreInPhysicalUnits[currentStageNumber] );
        similarityRegistration->SetMetricSamplingStrategy(
          static_cast<typename SimilarityRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
        similarityRegistration->SetMetricSamplingPercentage( samplingPercentage );
        similarityRegistration->SetOptimizer( optimizer );
        if( this->m_CompositeTransform->GetNumberOfTransforms() > 0 )
          {
          similarityRegistration->SetMovingInitialTransform( this->m_CompositeTransform );
          }
        if( this->m_FixedInitialTransform->GetNumberOfTransforms() > 0 )
          {
          similarityRegistration->SetFixedInitialTransform( this->m_FixedInitialTransform );
          }

        typedef antsRegistrationCommandIterationUpdate<SimilarityRegistrationType> SimilarityCommandType;
        typename SimilarityCommandType::Pointer similarityObserver = SimilarityCommandType::New();
        similarityObserver->SetLogStream(*this->m_LogStream);
        similarityObserver->SetNumberOfIterations( currentStageIterations );

        similarityRegistration->AddObserver( itk::IterationEvent(), similarityObserver );
        similarityRegistration->AddObserver( itk::InitializeEvent(), similarityObserver );

        try
          {
          this->Logger() << std::endl << "*** Running similarity registration ***" << std::endl << std::endl;
          similarityObserver->Execute( similarityRegistration, itk::StartEvent() );
          similarityRegistration->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cout << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }

        // Add calculated transform to the composite transform or add it to the composite transform
        // which is incorporated into the fixed image header.
        if( this->m_ApplyLinearTransformsToFixedImageHeader && this->m_AllPreviousTransformsAreLinear )
          {
          this->m_CompositeLinearTransformForFixedImageHeader->AddTransform( const_cast<SimilarityTransformType *>(
                                                                               similarityRegistration->GetOutput()->Get() ) );
          }
        else
          {
          this->m_CompositeTransform->AddTransform( const_cast<SimilarityTransformType *>( similarityRegistration->
                                                                                           GetOutput()->Get() ) );
          }
        }
        break;
      case Translation:
        {
        typedef itk::TranslationTransform<RealType, VImageDimension> TranslationTransformType;

        typedef itk::ImageRegistrationMethodv4<ImageType, ImageType,
                                               TranslationTransformType> TranslationRegistrationType;
        typename TranslationRegistrationType::Pointer translationRegistration = TranslationRegistrationType::New();
        for( unsigned int n = 0; n < stageMetricList.size(); n++ )
          {
          translationRegistration->SetFixedImage( n, preprocessedFixedImagesPerStage[n] );
          translationRegistration->SetMovingImage( n, preprocessedMovingImagesPerStage[n] );
          }
        if( useMultiMetric )
          {
          translationRegistration->SetMetric( multiMetric );
          }
        else
          {
          translationRegistration->SetMetric( singleMetric );
          }

        translationRegistration->SetNumberOfLevels( numberOfLevels );
        for( unsigned int level = 0; level < numberOfLevels; ++level )
          {
          translationRegistration->SetShrinkFactorsPerDimension( level, shrinkFactorsPerDimensionForAllLevels[level] );
          }
        translationRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
        translationRegistration->SetSmoothingSigmasAreSpecifiedInPhysicalUnits(
          this->m_SmoothingSigmasAreInPhysicalUnits[currentStageNumber] );
        translationRegistration->SetMetricSamplingStrategy(
          static_cast<typename TranslationRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
        translationRegistration->SetMetricSamplingPercentage( samplingPercentage );
        translationRegistration->SetOptimizer( optimizer );
        if( this->m_CompositeTransform->GetNumberOfTransforms() > 0 )
          {
          translationRegistration->SetMovingInitialTransform( this->m_CompositeTransform );
          }
        if( this->m_FixedInitialTransform->GetNumberOfTransforms() > 0 )
          {
          translationRegistration->SetFixedInitialTransform( this->m_FixedInitialTransform );
          }

        typedef antsRegistrationCommandIterationUpdate<TranslationRegistrationType> TranslationCommandType;
        typename TranslationCommandType::Pointer translationObserver = TranslationCommandType::New();
        translationObserver->SetLogStream(*this->m_LogStream);
        translationObserver->SetNumberOfIterations( currentStageIterations );

        translationRegistration->AddObserver( itk::IterationEvent(), translationObserver );
        translationRegistration->AddObserver( itk::InitializeEvent(), translationObserver );

        try
          {
          this->Logger() << std::endl << "*** Running translation registration ***" << std::endl << std::endl;
          translationObserver->Execute( translationRegistration, itk::StartEvent() );
          translationRegistration->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cout << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }

        // Add calculated transform to the composite transform or add it to the composite transform
        // which is incorporated into the fixed image header.
        if( this->m_ApplyLinearTransformsToFixedImageHeader && this->m_AllPreviousTransformsAreLinear )
          {
          this->m_CompositeLinearTransformForFixedImageHeader->AddTransform( const_cast<TranslationTransformType *>(
                                                                               translationRegistration->GetOutput()->
                                                                               Get() ) );
          }
        else
          {
          this->m_CompositeTransform->AddTransform( const_cast<TranslationTransformType *>( translationRegistration->
                                                                                            GetOutput()->Get() ) );
          }
        }
        break;
      case GaussianDisplacementField:
        {
        typedef itk::Vector<RealType, VImageDimension> VectorType;
        VectorType zeroVector( 0.0 );
        //typedef itk::Image<VectorType, VImageDimension> DisplacementFieldType;

        // ORIENTATION ALERT: Original code set image size to
        // fixedImage buffered region, & if fixedImage BufferedRegion
        // != LargestPossibleRegion, this code would be wrong.

        typename DisplacementFieldType::Pointer displacementField = AllocImage<DisplacementFieldType>(
            preprocessedFixedImagesPerStage[0], zeroVector );

        typedef itk::GaussianSmoothingOnUpdateDisplacementFieldTransform<RealType, VImageDimension>
          GaussianDisplacementFieldTransformType;

        typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, GaussianDisplacementFieldTransformType>
          DisplacementFieldRegistrationType;
        typename DisplacementFieldRegistrationType::Pointer displacementFieldRegistration =
          DisplacementFieldRegistrationType::New();

        if( this->m_RestrictDeformationOptimizerWeights.size() == VImageDimension )
          {
          typename DisplacementFieldRegistrationType::OptimizerWeightsType optimizerWeights( VImageDimension );
          for( unsigned int d = 0; d < VImageDimension; d++ )
            {
            optimizerWeights[d] = this->m_RestrictDeformationOptimizerWeights[d];
            }
          displacementFieldRegistration->SetOptimizerWeights( optimizerWeights );
          }

        typename GaussianDisplacementFieldTransformType::Pointer outputDisplacementFieldTransform =
          const_cast<GaussianDisplacementFieldTransformType *>( displacementFieldRegistration->GetOutput()->Get() );

        // Create the transform adaptors

        typedef itk::GaussianSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor<
          GaussianDisplacementFieldTransformType> DisplacementFieldTransformAdaptorType;
        typename DisplacementFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;

        // Extract parameters

        RealType varianceForUpdateField =
          this->m_TransformMethods[currentStageNumber].m_UpdateFieldVarianceInVarianceSpace;
        RealType varianceForTotalField  =
          this->m_TransformMethods[currentStageNumber].m_TotalFieldVarianceInVarianceSpace;

        outputDisplacementFieldTransform->SetGaussianSmoothingVarianceForTheUpdateField( varianceForUpdateField );
        outputDisplacementFieldTransform->SetGaussianSmoothingVarianceForTheTotalField( varianceForTotalField );
        outputDisplacementFieldTransform->SetDisplacementField( displacementField );
        // Create the transform adaptors
        // For the gaussian displacement field, the specified variances are in image spacing terms
        // and, in normal practice, we typically don't change these values at each level.  However,
        // if the user wishes to add that option, they can use the class
        // GaussianSmoothingOnUpdateDisplacementFieldTransformAdaptor
        for( unsigned int level = 0; level < numberOfLevels; level++ )
          {
          // We use the shrink image filter to calculate the fixed parameters of the virtual
          // domain at each level.  To speed up calculation and avoid unnecessary memory
          // usage, we could calculate these fixed parameters directly.

          typedef itk::ShrinkImageFilter<DisplacementFieldType, DisplacementFieldType> ShrinkFilterType;
          typename ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
          shrinkFilter->SetShrinkFactors( shrinkFactorsPerDimensionForAllLevels[level] );
          shrinkFilter->SetInput( displacementField );
          shrinkFilter->Update();

          typename DisplacementFieldTransformAdaptorType::Pointer fieldTransformAdaptor =
            DisplacementFieldTransformAdaptorType::New();
          fieldTransformAdaptor->SetRequiredSpacing( shrinkFilter->GetOutput()->GetSpacing() );
          fieldTransformAdaptor->SetRequiredSize( shrinkFilter->GetOutput()->GetBufferedRegion().GetSize() );
          fieldTransformAdaptor->SetRequiredDirection( shrinkFilter->GetOutput()->GetDirection() );
          fieldTransformAdaptor->SetRequiredOrigin( shrinkFilter->GetOutput()->GetOrigin() );
          fieldTransformAdaptor->SetTransform( outputDisplacementFieldTransform );

          adaptors.push_back( fieldTransformAdaptor.GetPointer() );
          }
        for( unsigned int n = 0; n < stageMetricList.size(); n++ )
          {
          displacementFieldRegistration->SetFixedImage( n, preprocessedFixedImagesPerStage[n] );
          displacementFieldRegistration->SetMovingImage( n, preprocessedMovingImagesPerStage[n] );
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
        displacementFieldRegistration->SetSmoothingSigmasAreSpecifiedInPhysicalUnits(
          this->m_SmoothingSigmasAreInPhysicalUnits[currentStageNumber] );
        displacementFieldRegistration->SetMetricSamplingStrategy(
          static_cast<typename DisplacementFieldRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
        displacementFieldRegistration->SetMetricSamplingPercentage( samplingPercentage );
        displacementFieldRegistration->SetOptimizer( optimizer );
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
        displacementFieldRegistrationObserver->SetLogStream(*this->m_LogStream);
        displacementFieldRegistrationObserver->SetNumberOfIterations( currentStageIterations );

        displacementFieldRegistration->AddObserver( itk::IterationEvent(), displacementFieldRegistrationObserver );
        displacementFieldRegistration->AddObserver( itk::InitializeEvent(), displacementFieldRegistrationObserver );

        try
          {
          this->Logger() << std::endl
                         << "*** Running gaussian displacement field registration (varianceForUpdateField = "
                         << varianceForUpdateField << ", varianceForTotalField = " << varianceForTotalField << ") ***"
                         << std::endl << std::endl;
          displacementFieldRegistrationObserver->Execute( displacementFieldRegistration, itk::StartEvent() );
          displacementFieldRegistration->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cout << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }

        // Add calculated transform to the composite transform
        this->m_CompositeTransform->AddTransform( outputDisplacementFieldTransform );

        this->m_AllPreviousTransformsAreLinear = false;
        }
        break;
      case BSplineDisplacementField:
        {
        typedef itk::Vector<RealType, VImageDimension> VectorType;
        VectorType zeroVector( 0.0 );
        //typedef itk::Image<VectorType, VImageDimension> DisplacementFieldType;

        // ORIENTATION ALERT -- see comment above.

        typename DisplacementFieldType::Pointer displacementField = AllocImage<DisplacementFieldType>(
            preprocessedFixedImagesPerStage[0], zeroVector );

        typedef itk::BSplineSmoothingOnUpdateDisplacementFieldTransform<RealType, VImageDimension>
          BSplineDisplacementFieldTransformType;

        typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, BSplineDisplacementFieldTransformType>
          DisplacementFieldRegistrationType;
        typename DisplacementFieldRegistrationType::Pointer displacementFieldRegistration =
          DisplacementFieldRegistrationType::New();

        if( this->m_RestrictDeformationOptimizerWeights.size() == VImageDimension )
          {
          typename DisplacementFieldRegistrationType::OptimizerWeightsType optimizerWeights( VImageDimension );
          for( unsigned int d = 0; d < VImageDimension; d++ )
            {
            optimizerWeights[d] = this->m_RestrictDeformationOptimizerWeights[d];
            }
          displacementFieldRegistration->SetOptimizerWeights( optimizerWeights );
          }

        typename BSplineDisplacementFieldTransformType::Pointer outputDisplacementFieldTransform =
          const_cast<BSplineDisplacementFieldTransformType *>( displacementFieldRegistration->GetOutput()->Get() );
        outputDisplacementFieldTransform->SetDisplacementField( displacementField );

        // Create the transform adaptors

        typedef itk::BSplineSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor<
            BSplineDisplacementFieldTransformType> DisplacementFieldTransformAdaptorType;
        typename DisplacementFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;

        // Extract parameters

        const std::vector<unsigned int> & meshSizeForTheUpdateField =
          this->m_TransformMethods[currentStageNumber].m_UpdateFieldMeshSizeAtBaseLevel;
        std::vector<unsigned int> meshSizeForTheTotalField =
          this->m_TransformMethods[currentStageNumber].m_TotalFieldMeshSizeAtBaseLevel;

        outputDisplacementFieldTransform->SetSplineOrder( this->m_TransformMethods[currentStageNumber].m_SplineOrder );

        if( meshSizeForTheUpdateField.size() != VImageDimension || meshSizeForTheTotalField.size() != VImageDimension )
          {
          std::cout << "ERROR:  The mesh size(s) don't match the ImageDimension." << std::endl;
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
          // We use the shrink image filter to calculate the fixed parameters of the virtual
          // domain at each level.  To speed up calculation and avoid unnecessary memory
          // usage, we could calculate these fixed parameters directly.

          typedef itk::ShrinkImageFilter<DisplacementFieldType, DisplacementFieldType> ShrinkFilterType;
          typename ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
          shrinkFilter->SetShrinkFactors( shrinkFactorsPerDimensionForAllLevels[level] );
          shrinkFilter->SetInput( displacementField );
          shrinkFilter->Update();

          typedef itk::BSplineSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor<
              BSplineDisplacementFieldTransformType> BSplineDisplacementFieldTransformAdaptorType;
          typename BSplineDisplacementFieldTransformAdaptorType::Pointer bsplineFieldTransformAdaptor =
            BSplineDisplacementFieldTransformAdaptorType::New();
          bsplineFieldTransformAdaptor->SetRequiredSpacing( shrinkFilter->GetOutput()->GetSpacing() );
          bsplineFieldTransformAdaptor->SetRequiredSize( shrinkFilter->GetOutput()->GetBufferedRegion().GetSize() );
          bsplineFieldTransformAdaptor->SetRequiredDirection( shrinkFilter->GetOutput()->GetDirection() );
          bsplineFieldTransformAdaptor->SetRequiredOrigin( shrinkFilter->GetOutput()->GetOrigin() );
          bsplineFieldTransformAdaptor->SetTransform( outputDisplacementFieldTransform );

          // A good heuristic is to RealType the b-spline mesh resolution at each level
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
        for( unsigned int n = 0; n < stageMetricList.size(); n++ )
          {
          displacementFieldRegistration->SetFixedImage( n, preprocessedFixedImagesPerStage[n] );
          displacementFieldRegistration->SetMovingImage( n, preprocessedMovingImagesPerStage[n] );
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
        displacementFieldRegistration->SetSmoothingSigmasAreSpecifiedInPhysicalUnits(
          this->m_SmoothingSigmasAreInPhysicalUnits[currentStageNumber] );
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
        displacementFieldRegistration->SetOptimizer( optimizer );
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
                         << "*** Running bspline displacement field registration (updateMeshSizeAtBaseLevel = "
                         << updateMeshSize << ", totalMeshSizeAtBaseLevel = " << totalMeshSize << ") ***" << std::endl
                         << std::endl;
          displacementFieldRegistrationObserver->Execute( displacementFieldRegistration, itk::StartEvent() );
          displacementFieldRegistration->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cout << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }

        // Add calculated transform to the composite transform
        this->m_CompositeTransform->AddTransform( outputDisplacementFieldTransform );

        this->m_AllPreviousTransformsAreLinear = false;
        }
        break;
      case BSpline:
        {
        const unsigned int SplineOrder = 3;
        typedef itk::BSplineTransform<RealType, VImageDimension, SplineOrder> BSplineTransformType;

        typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, BSplineTransformType> BSplineRegistrationType;
        typename BSplineRegistrationType::Pointer bsplineRegistration = BSplineRegistrationType::New();

        if( this->m_RestrictDeformationOptimizerWeights.size() == VImageDimension )
          {
          typename BSplineRegistrationType::OptimizerWeightsType optimizerWeights( VImageDimension );
          for( unsigned int d = 0; d < VImageDimension; d++ )
            {
            optimizerWeights[d] = this->m_RestrictDeformationOptimizerWeights[d];
            }
          bsplineRegistration->SetOptimizerWeights( optimizerWeights );
          }

        typename BSplineTransformType::Pointer outputBSplineTransform =
          const_cast<BSplineTransformType *>( bsplineRegistration->GetOutput()->Get() );

        const std::vector<unsigned int> & size = this->m_TransformMethods[currentStageNumber].m_MeshSizeAtBaseLevel;

        typename BSplineTransformType::PhysicalDimensionsType physicalDimensions;
        typename BSplineTransformType::MeshSizeType meshSize;
        for( unsigned int d = 0; d < VImageDimension; d++ )
          {
          physicalDimensions[d] = preprocessedFixedImagesPerStage[0]->GetSpacing()[d]
            * static_cast<RealType>( preprocessedFixedImagesPerStage[0]->GetLargestPossibleRegion().GetSize()[d] - 1 );
          meshSize[d] = size[d];
          }

        // Create the transform adaptors

        typedef itk::BSplineTransformParametersAdaptor<BSplineTransformType> BSplineTransformAdaptorType;
        typename BSplineRegistrationType::TransformParametersAdaptorsContainerType adaptors;
        // Create the transform adaptors specific to B-splines
        for( unsigned int level = 0; level < numberOfLevels; level++ )
          {
          typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkFilterType;
          typename ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
          shrinkFilter->SetShrinkFactors( shrinkFactorsPerDimensionForAllLevels[level] );
          shrinkFilter->SetInput( preprocessedFixedImagesPerStage[0] );
          shrinkFilter->Update();

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
          bsplineAdaptor->SetRequiredTransformDomainOrigin( shrinkFilter->GetOutput()->GetOrigin() );
          bsplineAdaptor->SetRequiredTransformDomainDirection( shrinkFilter->GetOutput()->GetDirection() );
          bsplineAdaptor->SetRequiredTransformDomainPhysicalDimensions( physicalDimensions );

          adaptors.push_back( bsplineAdaptor.GetPointer() );
          }
        for( unsigned int n = 0; n < stageMetricList.size(); n++ )
          {
          bsplineRegistration->SetFixedImage( n, preprocessedFixedImagesPerStage[n] );
          bsplineRegistration->SetMovingImage( n, preprocessedMovingImagesPerStage[n] );
          }
        if( useMultiMetric )
          {
          bsplineRegistration->SetMetric( multiMetric );
          }
        else
          {
          bsplineRegistration->SetMetric( singleMetric );
          }

        bsplineRegistration->SetNumberOfLevels( numberOfLevels );
        for( unsigned int level = 0; level < numberOfLevels; ++level )
          {
          bsplineRegistration->SetShrinkFactorsPerDimension( level, shrinkFactorsPerDimensionForAllLevels[level] );
          }
        bsplineRegistration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
        bsplineRegistration->SetSmoothingSigmasAreSpecifiedInPhysicalUnits( this->m_SmoothingSigmasAreInPhysicalUnits[
                                                                              currentStageNumber] );
        bsplineRegistration->SetMetricSamplingStrategy(
          static_cast<typename BSplineRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
        bsplineRegistration->SetMetricSamplingPercentage( samplingPercentage );
        bsplineRegistration->SetOptimizer( optimizer );
        if( this->m_CompositeTransform->GetNumberOfTransforms() > 0 )
          {
          bsplineRegistration->SetMovingInitialTransform( this->m_CompositeTransform );
          }
        if( this->m_FixedInitialTransform->GetNumberOfTransforms() > 0 )
          {
          bsplineRegistration->SetFixedInitialTransform( this->m_FixedInitialTransform );
          }
        bsplineRegistration->SetTransformParametersAdaptorsPerLevel( adaptors );
        outputBSplineTransform->SetTransformDomainOrigin( preprocessedFixedImagesPerStage[0]->GetOrigin() );
        outputBSplineTransform->SetTransformDomainPhysicalDimensions( physicalDimensions );
        outputBSplineTransform->SetTransformDomainMeshSize( meshSize );
        outputBSplineTransform->SetTransformDomainDirection( preprocessedFixedImagesPerStage[0]->GetDirection() );
        outputBSplineTransform->SetIdentity();

        typedef antsRegistrationCommandIterationUpdate<BSplineRegistrationType> BSplineCommandType;
        typename BSplineCommandType::Pointer bsplineObserver = BSplineCommandType::New();
        bsplineObserver->SetLogStream(*this->m_LogStream);
        bsplineObserver->SetNumberOfIterations( currentStageIterations );

        bsplineRegistration->AddObserver( itk::IterationEvent(), bsplineObserver );
        bsplineRegistration->AddObserver( itk::InitializeEvent(), bsplineObserver );

        try
          {
          this->Logger() << std::endl << "*** Running bspline registration (meshSizeAtBaseLevel = " << meshSize
                         << ") ***"
                         << std::endl << std::endl;
          bsplineObserver->Execute( bsplineRegistration, itk::StartEvent() );
          bsplineRegistration->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cout << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }
        // Add calculated transform to the composite transform
        this->m_CompositeTransform->AddTransform( outputBSplineTransform );

        this->m_AllPreviousTransformsAreLinear = false;
        }
        break;
      case TimeVaryingVelocityField:
        {
        typedef itk::Vector<RealType, VImageDimension> VectorType;
        VectorType zeroVector( 0.0 );

        // Determine the parameters (size, spacing, etc) for the time-varying velocity field

        typedef itk::Image<VectorType, VImageDimension + 1> TimeVaryingVelocityFieldType;

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
                                                                       TimeVaryingVelocityFieldOutputTransformType>
                                                                                            VelocityFieldRegistrationType;
        typename VelocityFieldRegistrationType::Pointer velocityFieldRegistration =
          VelocityFieldRegistrationType::New();

        if( this->m_RestrictDeformationOptimizerWeights.size() == VImageDimension )
          {
          typename VelocityFieldRegistrationType::OptimizerWeightsType optimizerWeights( VImageDimension );
          for( unsigned int d = 0; d < VImageDimension; d++ )
            {
            optimizerWeights[d] = this->m_RestrictDeformationOptimizerWeights[d];
            }
          velocityFieldRegistration->SetOptimizerWeights( optimizerWeights );
          }

        typedef typename VelocityFieldRegistrationType::OutputTransformType OutputTransformType;
        typename OutputTransformType::Pointer outputTransform =
          const_cast<OutputTransformType *>( velocityFieldRegistration->GetOutput()->Get() );
        for( unsigned int n = 0; n < stageMetricList.size(); n++ )
          {
          velocityFieldRegistration->SetFixedImage( n, preprocessedFixedImagesPerStage[n] );
          velocityFieldRegistration->SetMovingImage( n, preprocessedMovingImagesPerStage[n] );
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
          typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkFilterType;
          typename ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
          shrinkFilter->SetShrinkFactors( shrinkFactorsPerDimensionForAllLevels[level] );
          shrinkFilter->SetInput( preprocessedFixedImagesPerStage[0] );
          shrinkFilter->Update();

          // Although we shrink the images for the given levels,
          // we keep the size in time the same

          velocityFieldSize.Fill( numberOfTimeIndices );
          velocityFieldOrigin.Fill( 0.0 );
          velocityFieldSpacing.Fill( 1.0 );
          velocityFieldDirection.SetIdentity();

          fixedImageSize = shrinkFilter->GetOutput()->GetBufferedRegion().GetSize();
          fixedImageOrigin = shrinkFilter->GetOutput()->GetOrigin();
          fixedImageSpacing = shrinkFilter->GetOutput()->GetSpacing();
          fixedImageDirection = shrinkFilter->GetOutput()->GetDirection();
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
          std::cout << "Exception caught: " << e << std::endl;
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
          std::cout << "The transform domain mesh size does not have the correct number of elements."
                           << "For image dimension = " << VImageDimension << ", you need " << VImageDimension + 1
                           << "elements. " << std::endl;
          return EXIT_FAILURE;
          }

        unsigned int numberOfTimePointSamples =
          this->m_TransformMethods[currentStageNumber].m_NumberOfTimePointSamples;
        unsigned int splineOrder = this->m_TransformMethods[currentStageNumber].m_SplineOrder;

        typedef itk::Image<VectorType, VImageDimension + 1> TimeVaryingVelocityFieldControlPointLatticeType;

        typename ImageType::SizeType fixedImageSize = preprocessedFixedImagesPerStage[0]->GetBufferedRegion().GetSize();
        typename ImageType::PointType fixedImageOrigin = preprocessedFixedImagesPerStage[0]->GetOrigin();
        typename ImageType::SpacingType fixedImageSpacing = preprocessedFixedImagesPerStage[0]->GetSpacing();
        typename ImageType::DirectionType fixedImageDirection = preprocessedFixedImagesPerStage[0]->GetDirection();

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

        typedef itk::TimeVaryingBSplineVelocityFieldImageRegistrationMethod<ImageType, ImageType,
          TimeVaryingBSplineVelocityFieldOutputTransformType>
          VelocityFieldRegistrationType;
        typename VelocityFieldRegistrationType::Pointer velocityFieldRegistration =
          VelocityFieldRegistrationType::New();

        if( this->m_RestrictDeformationOptimizerWeights.size() == VImageDimension )
          {
          typename VelocityFieldRegistrationType::OptimizerWeightsType optimizerWeights( VImageDimension );
          for( unsigned int d = 0; d < VImageDimension; d++ )
            {
            optimizerWeights[d] = this->m_RestrictDeformationOptimizerWeights[d];
            }
          velocityFieldRegistration->SetOptimizerWeights( optimizerWeights );
          }

        typedef typename VelocityFieldRegistrationType::OutputTransformType OutputTransformType;
        typename OutputTransformType::Pointer outputTransform =
          const_cast<OutputTransformType *>( velocityFieldRegistration->GetOutput()->Get() );
        for( unsigned int n = 0; n < stageMetricList.size(); n++ )
          {
          velocityFieldRegistration->SetFixedImage( n, preprocessedFixedImagesPerStage[n] );
          velocityFieldRegistration->SetMovingImage( n, preprocessedMovingImagesPerStage[n] );
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
        velocityFieldRegistration->SetNumberOfTimePointSamples( numberOfTimePointSamples );
        velocityFieldRegistration->SetMetricSamplingStrategy(
          static_cast<typename VelocityFieldRegistrationType::MetricSamplingStrategyType>( metricSamplingStrategy ) );
        velocityFieldRegistration->SetMetricSamplingPercentage( samplingPercentage );
        velocityFieldRegistration->SetLearningRate( learningRate );
        velocityFieldRegistration->SetConvergenceThreshold( convergenceThreshold );
        velocityFieldRegistration->SetConvergenceWindowSize( convergenceWindowSize );
        outputTransform->SetSplineOrder( splineOrder );
        outputTransform->SetLowerTimeBound( 0.0 );
        outputTransform->SetUpperTimeBound( 1.0 );

        typedef itk::TimeVaryingBSplineVelocityFieldTransformParametersAdaptor<OutputTransformType>
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

        typename OutputTransformType::VelocityFieldPointType        sampledVelocityFieldOrigin;
        typename OutputTransformType::VelocityFieldSpacingType      sampledVelocityFieldSpacing;
        typename OutputTransformType::VelocityFieldSizeType         sampledVelocityFieldSize;
        typename OutputTransformType::VelocityFieldDirectionType    sampledVelocityFieldDirection;

        sampledVelocityFieldOrigin.Fill( 0.0 );
        sampledVelocityFieldSpacing.Fill( 1.0 );
        sampledVelocityFieldSize.Fill( numberOfTimePointSamples );
        sampledVelocityFieldDirection.SetIdentity();
        for( unsigned int i = 0; i < VImageDimension; i++ )
          {
          sampledVelocityFieldOrigin[i] = preprocessedFixedImagesPerStage[0]->GetOrigin()[i];
          sampledVelocityFieldSpacing[i] = preprocessedFixedImagesPerStage[0]->GetSpacing()[i];
          sampledVelocityFieldSize[i] = preprocessedFixedImagesPerStage[0]->GetRequestedRegion().GetSize()[i];
          for( unsigned int j = 0; j < VImageDimension; j++ )
            {
            sampledVelocityFieldDirection[i][j] = preprocessedFixedImagesPerStage[0]->GetDirection()[i][j];
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
          std::cout << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }
        // Add calculated transform to the composite transform
        this->m_CompositeTransform->AddTransform( outputTransform );

        this->m_AllPreviousTransformsAreLinear = false;
        }
        break;
      case SyN:
        {
        typedef itk::Vector<RealType, VImageDimension> VectorType;
        VectorType zeroVector( 0.0 );
        //typedef itk::Image<VectorType, VImageDimension> DisplacementFieldType;

        typename DisplacementFieldType::Pointer displacementField = AllocImage<DisplacementFieldType>(
            preprocessedFixedImagesPerStage[0], zeroVector );

        typename DisplacementFieldType::Pointer inverseDisplacementField = AllocImage<DisplacementFieldType>(
            preprocessedFixedImagesPerStage[0], zeroVector );

        typedef itk::SyNImageRegistrationMethod<ImageType, ImageType,
                                                DisplacementFieldTransformType> DisplacementFieldRegistrationType;
        typename DisplacementFieldRegistrationType::Pointer displacementFieldRegistration =
          DisplacementFieldRegistrationType::New();

        if( this->m_RestrictDeformationOptimizerWeights.size() == VImageDimension )
          {
          typename DisplacementFieldRegistrationType::OptimizerWeightsType optimizerWeights( VImageDimension );
          for( unsigned int d = 0; d < VImageDimension; d++ )
            {
            optimizerWeights[d] = this->m_RestrictDeformationOptimizerWeights[d];
            }
          displacementFieldRegistration->SetOptimizerWeights( optimizerWeights );
          }

        typename DisplacementFieldTransformType::Pointer outputDisplacementFieldTransform =
          const_cast<DisplacementFieldTransformType *>( displacementFieldRegistration->GetOutput()->Get() );

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
          // TODO:
          // We use the shrink image filter to calculate the fixed parameters of the virtual
          // domain at each level.  To speed up calculation and avoid unnecessary memory
          // usage, we could calculate these fixed parameters directly.
          typedef itk::ShrinkImageFilter<DisplacementFieldType, DisplacementFieldType> ShrinkFilterType;
          typename ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
          shrinkFilter->SetShrinkFactors( shrinkFactorsPerDimensionForAllLevels[level] );
          shrinkFilter->SetInput( displacementField );
          shrinkFilter->Update();

          typename DisplacementFieldTransformAdaptorType::Pointer fieldTransformAdaptor =
            DisplacementFieldTransformAdaptorType::New();
          fieldTransformAdaptor->SetRequiredSpacing( shrinkFilter->GetOutput()->GetSpacing() );
          fieldTransformAdaptor->SetRequiredSize( shrinkFilter->GetOutput()->GetBufferedRegion().GetSize() );
          fieldTransformAdaptor->SetRequiredDirection( shrinkFilter->GetOutput()->GetDirection() );
          fieldTransformAdaptor->SetRequiredOrigin( shrinkFilter->GetOutput()->GetOrigin() );
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
          displacementFieldRegistration->SetFixedImage( n, preprocessedFixedImagesPerStage[n] );
          displacementFieldRegistration->SetMovingImage( n, preprocessedMovingImagesPerStage[n] );
          }
        if( useMultiMetric )
          {
          displacementFieldRegistration->SetMetric( multiMetric );
          }
        else
          {
          displacementFieldRegistration->SetMetric( singleMetric );
          }

        if( this->m_CompositeTransform->GetNumberOfTransforms() > 0 )
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
          displacementFieldRegistrationObserver2->SetWriteInterationsOutputsInIntervals( this->m_WriteIntervalVolumes );
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
          std::cout << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }

        // Add calculated transform to the composite transform
        this->m_CompositeTransform->AddTransform( outputDisplacementFieldTransform );
        this->m_AllPreviousTransformsAreLinear = false;
        }
        break;
      case BSplineSyN:
        {
        typedef itk::Vector<RealType, VImageDimension> VectorType;
        VectorType zeroVector( 0.0 );
        //typedef itk::Image<VectorType, VImageDimension> DisplacementFieldType;

        typename DisplacementFieldType::Pointer displacementField = AllocImage<DisplacementFieldType>(
            preprocessedFixedImagesPerStage[0], zeroVector );

        typename DisplacementFieldType::Pointer inverseDisplacementField = AllocImage<DisplacementFieldType>(
            preprocessedFixedImagesPerStage[0], zeroVector );

        typedef itk::BSplineSmoothingOnUpdateDisplacementFieldTransform<RealType, VImageDimension>
          BSplineDisplacementFieldTransformType;

        typedef itk::BSplineSyNImageRegistrationMethod<ImageType, ImageType, BSplineDisplacementFieldTransformType>
          DisplacementFieldRegistrationType;
        typename DisplacementFieldRegistrationType::Pointer displacementFieldRegistration =
          DisplacementFieldRegistrationType::New();

        if( this->m_RestrictDeformationOptimizerWeights.size() == VImageDimension )
          {
          typename DisplacementFieldRegistrationType::OptimizerWeightsType optimizerWeights( VImageDimension );
          for( unsigned int d = 0; d < VImageDimension; d++ )
            {
            optimizerWeights[d] = this->m_RestrictDeformationOptimizerWeights[d];
            }
          displacementFieldRegistration->SetOptimizerWeights( optimizerWeights );
          }

        typename BSplineDisplacementFieldTransformType::Pointer outputDisplacementFieldTransform =
          const_cast<BSplineDisplacementFieldTransformType *>( displacementFieldRegistration->GetOutput()->Get() );

        // Create the transform adaptors

        typedef itk::BSplineSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor<
            BSplineDisplacementFieldTransformType> DisplacementFieldTransformAdaptorType;
        typename DisplacementFieldRegistrationType::TransformParametersAdaptorsContainerType adaptors;

        // Extract parameters

        const std::vector<unsigned int> & meshSizeForTheUpdateField =
          this->m_TransformMethods[currentStageNumber].m_UpdateFieldMeshSizeAtBaseLevel;
        std::vector<unsigned int> meshSizeForTheTotalField =
          this->m_TransformMethods[currentStageNumber].m_TotalFieldMeshSizeAtBaseLevel;

        outputDisplacementFieldTransform->SetSplineOrder( this->m_TransformMethods[currentStageNumber].m_SplineOrder );

        if( meshSizeForTheUpdateField.size() != VImageDimension || meshSizeForTheTotalField.size() != VImageDimension )
          {
          std::cout << "ERROR:  The mesh size(s) don't match the ImageDimension." << std::endl;
          return EXIT_FAILURE;
          }

        typename BSplineDisplacementFieldTransformType::ArrayType updateMeshSize;
        typename BSplineDisplacementFieldTransformType::ArrayType totalMeshSize;
        for( unsigned int d = 0; d < VImageDimension; d++ )
          {
          updateMeshSize[d] = meshSizeForTheUpdateField[d];
          totalMeshSize[d] = meshSizeForTheTotalField[d];
          }
        // Create the transform adaptors
        for( unsigned int level = 0; level < numberOfLevels; level++ )
          {
          // We use the shrink image filter to calculate the fixed parameters of the virtual
          // domain at each level.  To speed up calculation and avoid unnecessary memory
          // usage, we could calculate these fixed parameters directly.

          typedef itk::ShrinkImageFilter<DisplacementFieldType, DisplacementFieldType> ShrinkFilterType;
          typename ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
          shrinkFilter->SetShrinkFactors( shrinkFactorsPerDimensionForAllLevels[level] );
          shrinkFilter->SetInput( displacementField );
          shrinkFilter->Update();

          typedef itk::BSplineSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor<
              BSplineDisplacementFieldTransformType>
            BSplineDisplacementFieldTransformAdaptorType;
          typename BSplineDisplacementFieldTransformAdaptorType::Pointer bsplineFieldTransformAdaptor =
            BSplineDisplacementFieldTransformAdaptorType::New();
          bsplineFieldTransformAdaptor->SetRequiredSpacing( shrinkFilter->GetOutput()->GetSpacing() );
          bsplineFieldTransformAdaptor->SetRequiredSize( shrinkFilter->GetOutput()->GetBufferedRegion().GetSize() );
          bsplineFieldTransformAdaptor->SetRequiredDirection( shrinkFilter->GetOutput()->GetDirection() );
          bsplineFieldTransformAdaptor->SetRequiredOrigin( shrinkFilter->GetOutput()->GetOrigin() );
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

        // Extract parameters
        typename DisplacementFieldRegistrationType::NumberOfIterationsArrayType numberOfIterationsPerLevel;
        numberOfIterationsPerLevel.SetSize( numberOfLevels );
        for( unsigned int d = 0; d < numberOfLevels; d++ )
          {
          numberOfIterationsPerLevel[d] = currentStageIterations[d];
          }
        for( unsigned int n = 0; n < stageMetricList.size(); n++ )
          {
          displacementFieldRegistration->SetFixedImage( n, preprocessedFixedImagesPerStage[n] );
          displacementFieldRegistration->SetMovingImage( n, preprocessedMovingImagesPerStage[n] );
          }
        if( useMultiMetric )
          {
          displacementFieldRegistration->SetMetric( multiMetric );
          }
        else
          {
          displacementFieldRegistration->SetMetric( singleMetric );
          }

        if( this->m_CompositeTransform->GetNumberOfTransforms() > 0 )
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
        outputDisplacementFieldTransform->SetDisplacementField( displacementField );
        outputDisplacementFieldTransform->SetInverseDisplacementField( inverseDisplacementField );

        typedef antsRegistrationCommandIterationUpdate<DisplacementFieldRegistrationType> DisplacementFieldCommandType;
        typename DisplacementFieldCommandType::Pointer displacementFieldRegistrationObserver =
          DisplacementFieldCommandType::New();
        displacementFieldRegistrationObserver->SetLogStream(*this->m_LogStream);
        displacementFieldRegistrationObserver->SetNumberOfIterations( currentStageIterations );

        displacementFieldRegistration->AddObserver( itk::IterationEvent(), displacementFieldRegistrationObserver );
        displacementFieldRegistration->AddObserver( itk::InitializeEvent(), displacementFieldRegistrationObserver );

        try
          {
          this->Logger() << std::endl << "*** Running B-spline SyN registration (updateMeshSizeAtBaseLevel = "
                         << updateMeshSize << ", totalMeshSizeAtBaseLevel = " << totalMeshSize << ") ***" << std::endl
                         << std::endl;
          displacementFieldRegistrationObserver->Execute( displacementFieldRegistration, itk::StartEvent() );
          displacementFieldRegistration->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cout << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }

        // Add calculated transform to the composite transform
        this->m_CompositeTransform->AddTransform( outputDisplacementFieldTransform );

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

        typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, GaussianDisplacementFieldTransformType>
          DisplacementFieldRegistrationType;
        typename DisplacementFieldRegistrationType::Pointer displacementFieldRegistration =
          DisplacementFieldRegistrationType::New();

        if( this->m_RestrictDeformationOptimizerWeights.size() == VImageDimension )
          {
          typename DisplacementFieldRegistrationType::OptimizerWeightsType optimizerWeights( VImageDimension );
          for( unsigned int d = 0; d < VImageDimension; d++ )
            {
            optimizerWeights[d] = this->m_RestrictDeformationOptimizerWeights[d];
            }
          displacementFieldRegistration->SetOptimizerWeights( optimizerWeights );
          }

        typename GaussianDisplacementFieldTransformType::Pointer outputDisplacementFieldTransform =
          const_cast<GaussianDisplacementFieldTransformType *>( displacementFieldRegistration->GetOutput()->Get() );

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
        // Create the transform adaptors
        // For the gaussian displacement field, the specified variances are in image spacing terms
        // and, in normal practice, we typically don't change these values at each level.  However,
        // if the user wishes to add that option, they can use the class
        // GaussianSmoothingOnUpdateDisplacementFieldTransformAdaptor
        for( unsigned int level = 0; level < numberOfLevels; level++ )
          {
          // We use the shrink image filter to calculate the fixed parameters of the virtual
          // domain at each level.  To speed up calculation and avoid unnecessary memory
          // usage, we could calculate these fixed parameters directly.

          typedef itk::ShrinkImageFilter<ConstantVelocityFieldType, ConstantVelocityFieldType> ShrinkFilterType;
          typename ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
          shrinkFilter->SetShrinkFactors( shrinkFactorsPerDimensionForAllLevels[level] );
          shrinkFilter->SetInput( constantVelocityField );
          shrinkFilter->Update();

          typename DisplacementFieldTransformAdaptorType::Pointer fieldTransformAdaptor =
            DisplacementFieldTransformAdaptorType::New();
          fieldTransformAdaptor->SetRequiredSpacing( shrinkFilter->GetOutput()->GetSpacing() );
          fieldTransformAdaptor->SetRequiredSize( shrinkFilter->GetOutput()->GetBufferedRegion().GetSize() );
          fieldTransformAdaptor->SetRequiredDirection( shrinkFilter->GetOutput()->GetDirection() );
          fieldTransformAdaptor->SetRequiredOrigin( shrinkFilter->GetOutput()->GetOrigin() );
          fieldTransformAdaptor->SetTransform( outputDisplacementFieldTransform );

          fieldTransformAdaptor->SetGaussianSmoothingVarianceForTheUpdateField( varianceForUpdateField );
          fieldTransformAdaptor->SetGaussianSmoothingVarianceForTheConstantVelocityField( varianceForVelocityField );

          adaptors.push_back( fieldTransformAdaptor.GetPointer() );
          }
        for( unsigned int n = 0; n < stageMetricList.size(); n++ )
          {
          displacementFieldRegistration->SetFixedImage( n, preprocessedFixedImagesPerStage[n] );
          displacementFieldRegistration->SetMovingImage( n, preprocessedMovingImagesPerStage[n] );
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
          std::cout << "Exception caught: " << e << std::endl;
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

        typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, BSplineDisplacementFieldTransformType>
          DisplacementFieldRegistrationType;
        typename DisplacementFieldRegistrationType::Pointer displacementFieldRegistration =
          DisplacementFieldRegistrationType::New();

        if( this->m_RestrictDeformationOptimizerWeights.size() == VImageDimension )
          {
          typename DisplacementFieldRegistrationType::OptimizerWeightsType optimizerWeights( VImageDimension );
          for( unsigned int d = 0; d < VImageDimension; d++ )
            {
            optimizerWeights[d] = this->m_RestrictDeformationOptimizerWeights[d];
            }
          displacementFieldRegistration->SetOptimizerWeights( optimizerWeights );
          }

        typename BSplineDisplacementFieldTransformType::Pointer outputDisplacementFieldTransform =
          const_cast<BSplineDisplacementFieldTransformType *>( displacementFieldRegistration->GetOutput()->Get() );

        // Create the transform adaptors

        typedef itk::BSplineSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor<
            BSplineDisplacementFieldTransformType>
          DisplacementFieldTransformAdaptorType;
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

        if( meshSizeForTheUpdateField.size() != VImageDimension || meshSizeForTheVelocityField.size() !=
            VImageDimension )
          {
          std::cout << "ERROR:  The mesh size(s) don't match the ImageDimension." << std::endl;
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
          // We use the shrink image filter to calculate the fixed parameters of the virtual
          // domain at each level.  To speed up calculation and avoid unnecessary memory
          // usage, we could calculate these fixed parameters directly.

          typedef itk::ShrinkImageFilter<ConstantVelocityFieldType, ConstantVelocityFieldType> ShrinkFilterType;
          typename ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
          shrinkFilter->SetShrinkFactors( shrinkFactorsPerDimensionForAllLevels[level] );
          shrinkFilter->SetInput( constantVelocityField );
          shrinkFilter->Update();

          typedef itk::BSplineExponentialDiffeomorphicTransformParametersAdaptor<BSplineDisplacementFieldTransformType>
            BSplineDisplacementFieldTransformAdaptorType;
          typename BSplineDisplacementFieldTransformAdaptorType::Pointer bsplineFieldTransformAdaptor =
            BSplineDisplacementFieldTransformAdaptorType::New();
          bsplineFieldTransformAdaptor->SetRequiredSpacing( shrinkFilter->GetOutput()->GetSpacing() );
          bsplineFieldTransformAdaptor->SetRequiredSize( shrinkFilter->GetOutput()->GetBufferedRegion().GetSize() );
          bsplineFieldTransformAdaptor->SetRequiredDirection( shrinkFilter->GetOutput()->GetDirection() );
          bsplineFieldTransformAdaptor->SetRequiredOrigin( shrinkFilter->GetOutput()->GetOrigin() );
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
          displacementFieldRegistration->SetFixedImage( n, preprocessedFixedImagesPerStage[n] );
          displacementFieldRegistration->SetMovingImage( n, preprocessedMovingImagesPerStage[n] );
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
          std::cout << "Exception caught: " << e << std::endl;
          return EXIT_FAILURE;
          }

        // Add calculated transform to the composite transform
        this->m_CompositeTransform->AddTransform( outputDisplacementFieldTransform );

        this->m_AllPreviousTransformsAreLinear = false;
        }
        break;
      default:
        std::cout << "ERROR:  Unrecognized transform option - " << whichTransform << std::endl;
        return EXIT_FAILURE;
      }
    timer.Stop();
    this->Logger() << "  Elapsed time (stage " << currentStageNumber << "): " << timer.GetMean() << std::endl
                   << std::endl;
    }

  if( this->m_ApplyLinearTransformsToFixedImageHeader &&
      this->m_CompositeLinearTransformForFixedImageHeader->GetNumberOfTransforms() > 0 )
    {
    this->m_CompositeTransform->PrependTransform( this->m_CompositeLinearTransformForFixedImageHeader );
    this->m_CompositeTransform->FlattenTransformQueue();
    }

  totalTimer.Stop();
  this->Logger() << std::endl << "Total elapsed time: " << totalTimer.GetMean() << std::endl;

  return EXIT_SUCCESS;
}

template <class TComputeType, unsigned VImageDimension>
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

template <class TComputeType, unsigned VImageDimension>
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

    if( this->m_ApplyLinearTransformsToFixedImageHeader && compXfrm->IsLinear() )
      {
      this->m_CompositeLinearTransformForFixedImageHeader = compToAdd;
      }
    else
      {
      this->m_FixedInitialTransform = compToAdd;
      this->m_AllPreviousTransformsAreLinear = false;
      }
    }
  else
    {
    compToAdd = CompositeTransformType::New();
    typename TransformType::Pointer xfrm = initialTransform->Clone();
    compToAdd->AddTransform( xfrm );
    if( this->m_ApplyLinearTransformsToFixedImageHeader && initialTransform->IsLinear() )
      {
      this->m_CompositeLinearTransformForFixedImageHeader = compToAdd;
      }
    else
      {
      this->m_FixedInitialTransform = compToAdd;
      this->m_AllPreviousTransformsAreLinear = false;
      }
    }
}

template <class TComputeType, unsigned VImageDimension>
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
      inputImage->GetLargestPossibleRegion().GetSize()[d] - 1 ) * inputImage->GetSpacing()[d];
    meshSize.push_back( static_cast<unsigned int>( vcl_ceil( domain / knotSpacing ) ) );
//     unsigned long extraPadding = static_cast<unsigned long>(
//       ( numberOfSpans * splineDistance - domain ) / inputImage->GetSpacing()[d] + 0.5 );
//     lowerBound[d] = static_cast<unsigned long>( 0.5 * extraPadding );
//     upperBound[d] = extraPadding - lowerBound[d];
//     numberOfControlPoints[d] = meshSize[d] + splineOrder;
    }

  return meshSize;
}

template <class TComputeType, unsigned VImageDimension>
typename RegistrationHelper<TComputeType, VImageDimension>::AffineTransformType::Pointer
RegistrationHelper<TComputeType, VImageDimension>
::CollapseLinearTransforms( const CompositeTransformType * compositeTransform )
{
  if( !compositeTransform->IsLinear() )
    {
    itkExceptionMacro( "The composite transform is not linear." );
    }

  typedef itk::TranslationTransform<RealType, VImageDimension> TranslationTransformType;

  typename AffineTransformType::Pointer totalTransform = AffineTransformType::New();
  for( unsigned int n = 0; n < compositeTransform->GetNumberOfTransforms(); n++ )
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
        dynamic_cast<MatrixOffsetTransformBaseType * const>( transform.GetPointer() );
      nthTransform->SetMatrix( matrixOffsetTransform->GetMatrix() );
      nthTransform->SetOffset( matrixOffsetTransform->GetOffset() );
      }
    totalTransform->Compose( nthTransform, true );
    }
  return totalTransform;
}

template <class TComputeType, unsigned VImageDimension>
typename RegistrationHelper<TComputeType, VImageDimension>::DisplacementFieldTransformPointer
RegistrationHelper<TComputeType, VImageDimension>
::CollapseDisplacementFieldTransforms( const CompositeTransformType * compositeTransform )
{
  if( compositeTransform->GetTransformCategory() != TransformType::DisplacementField  )
    {
    itkExceptionMacro( "The composite transform is not composed strictly of displacement fields." );
    }

  DisplacementFieldTransformPointer totalTransform = DisplacementFieldTransformType::New();

  if( compositeTransform->GetNumberOfTransforms() == 0 )
    {
    itkWarningMacro( "The composite transform is empty.  Returning empty displacement field transform." );

    return totalTransform;
    }

  bool hasInverse = true;
  for( unsigned int n = 0; n < compositeTransform->GetNumberOfTransforms(); n++ )
    {
    typename TransformType::Pointer transform = compositeTransform->GetNthTransform( n );

    typename DisplacementFieldTransformType::Pointer nthTransform =
      dynamic_cast<DisplacementFieldTransformType *>( transform.GetPointer() );

    if( n == 0 )
      {
      totalTransform->SetDisplacementField( nthTransform->GetModifiableDisplacementField() );
      if( nthTransform->GetInverseDisplacementField() )
        {
        totalTransform->SetInverseDisplacementField( nthTransform->GetModifiableInverseDisplacementField() );
        }
      else
        {
        hasInverse = false;
        }
      }
    else
      {
      typedef itk::ComposeDisplacementFieldsImageFilter<DisplacementFieldType> ComposerType;

      typename ComposerType::Pointer composer = ComposerType::New();
      composer->SetWarpingField( nthTransform->GetDisplacementField() );
      composer->SetDisplacementField( totalTransform->GetDisplacementField() );

      typename DisplacementFieldType::Pointer totalField = composer->GetOutput();
      totalField->Update();
      totalField->DisconnectPipeline();

      typename DisplacementFieldType::Pointer totalInverseField = NULL;

      if( hasInverse && nthTransform->GetInverseDisplacementField() )
        {
        typename ComposerType::Pointer inverseComposer = ComposerType::New();
        inverseComposer->SetWarpingField( totalTransform->GetInverseDisplacementField() );
        inverseComposer->SetDisplacementField( nthTransform->GetInverseDisplacementField() );

        totalInverseField = inverseComposer->GetOutput();
        totalInverseField->Update();
        totalInverseField->DisconnectPipeline();
        }
      else
        {
        hasInverse = false;
        }
      totalTransform->SetDisplacementField( totalField );
      totalTransform->SetInverseDisplacementField( totalInverseField );
      }
    }

  return totalTransform;
}

template <class TComputeType, unsigned VImageDimension>
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
  else if( compositeTransform->GetTransformCategory() == TransformType::DisplacementField )
    {
    collapsedCompositeTransform->AddTransform( this->CollapseDisplacementFieldTransforms( compositeTransform ) );
    return collapsedCompositeTransform;
    }

  // Find the first linear or displacement field transform
  typename TransformType::TransformCategoryType currentTransformCategory = TransformType::UnknownTransformCategory;
  unsigned int startIndex = 0;
  for( unsigned int n = 0; n < compositeTransform->GetNumberOfTransforms(); n++ )
    {
    typename TransformType::TransformCategoryType transformCategory =
      compositeTransform->GetNthTransform( n )->GetTransformCategory();
    if( transformCategory == TransformType::Linear || transformCategory == TransformType::DisplacementField )
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
  if( currentTransformCategory != TransformType::UnknownTransformCategory )
    {
    CompositeTransformPointer currentCompositeTransform = CompositeTransformType::New();
    currentCompositeTransform->AddTransform( compositeTransform->GetNthTransform( startIndex ) );
    for( unsigned int n = startIndex + 1; n < compositeTransform->GetNumberOfTransforms(); n++ )
      {
      typename TransformType::TransformCategoryType transformCategory =
        compositeTransform->GetNthTransform( n )->GetTransformCategory();
      if( transformCategory == currentTransformCategory )
        {
        currentCompositeTransform->AddTransform( compositeTransform->GetNthTransform( n ) );
        if( n == compositeTransform->GetNumberOfTransforms() - 1 )
          {
          if( currentTransformCategory == TransformType::Linear )
            {
            collapsedCompositeTransform->AddTransform( this->CollapseLinearTransforms( currentCompositeTransform ) );
            }
          else if( currentTransformCategory == TransformType::DisplacementField )
            {
            collapsedCompositeTransform->AddTransform( this->CollapseDisplacementFieldTransforms(
                                                         currentCompositeTransform ) );
            }
          }
        }
      else
        {
        if( currentTransformCategory == TransformType::Linear )
          {
          collapsedCompositeTransform->AddTransform( this->CollapseLinearTransforms( currentCompositeTransform ) );
          currentCompositeTransform->ClearTransformQueue();
          }
        else if( currentTransformCategory == TransformType::DisplacementField )
          {
          collapsedCompositeTransform->AddTransform( this->CollapseDisplacementFieldTransforms(
                                                       currentCompositeTransform ) );
          currentCompositeTransform->ClearTransformQueue();
          }
        currentTransformCategory = transformCategory;

        if( ( transformCategory == TransformType::Linear || transformCategory == TransformType::DisplacementField ) &&
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

  return collapsedCompositeTransform;
}

template <class TComputeType, unsigned VImageDimension>
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

template <class TComputeType, unsigned VImageDimension>
void
RegistrationHelper<TComputeType, VImageDimension>
::PrintState() const
{
  this->Logger() << "Dimension = " << Self::ImageDimension << std::endl
                 << "Number of stages = " << this->m_NumberOfStages << std::endl
                 << "Use Histogram Matching " << ( this->m_UseHistogramMatching ? "true" : "false" )
                 << std::endl
                 << "Winsorize Image Intensities "
                 << ( this->m_WinsorizeImageIntensities ? "true" : "false" ) << std::endl
                 << "Lower Quantile = " << this->m_LowerQuantile << std::endl
                 << "Upper Quantile = " << this->m_UpperQuantile << std::endl;

  for( unsigned i = 0; i < this->m_NumberOfStages; i++ )
    {
    this->Logger() << "Stage " << i + 1 << " State" << std::endl; // NOTE: + 1 for consistency.
    const Metric &          curMetric = this->m_Metrics[i];
    const TransformMethod & curTransform = this->m_TransformMethods[i];
    this->Logger() << "   Metric = " << curMetric.GetMetricAsString() << std::endl
                   << "     Fixed Image = " << curMetric.m_FixedImage << std::endl
                   << "     Moving Image = " << curMetric.m_MovingImage << std::endl
                   << "     Weighting = " << curMetric.m_Weighting << std::endl
                   << "     Sampling Strategy = "
                   << (curMetric.m_SamplingStrategy ==
        random ? "random" : (curMetric.m_SamplingStrategy == regular ) ? "regular" : (curMetric.m_SamplingStrategy ==
                                                                                      none ) ? "none" :
        "WARNING: UNKNOWN")
                   << std::endl
                   << "     NumberOfBins = " << curMetric.m_NumberOfBins << std::endl
                   << "     Radius = " << curMetric.m_Radius << std::endl
                   << "     Sampling percentage  = " << curMetric.m_SamplingPercentage << std::endl
                   << "   Transform = " << curTransform.XfrmMethodAsString() << std::endl
                   << "     Gradient Step = " << curTransform.m_GradientStep << std::endl
                   << "     Update Field Sigma (voxel space) = "
                   << curTransform.m_UpdateFieldVarianceInVarianceSpace << std::endl
                   << "     Total Field Sigma (voxel space) = "
                   << curTransform.m_TotalFieldVarianceInVarianceSpace << std::endl
                   << "     Update Field Time Sigma = " << curTransform.m_UpdateFieldTimeSigma << std::endl
                   << "     Total Field Time Sigma  = " << curTransform.m_TotalFieldTimeSigma << std::endl
                   << "     Number of Time Indices = " << curTransform.m_NumberOfTimeIndices << std::endl
                   << "     Number of Time Point Samples = " << curTransform.m_NumberOfTimeIndices << std::endl;
    }
}
} // namespace ants

#endif // __itkantsRegistrationHelper_hxx
