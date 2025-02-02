/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __antsRegistrationHelper_h
#define __antsRegistrationHelper_h

#include "ReadWriteData.h"
#include "itkANTSAffine3DTransform.h"
#include "itkANTSCenteredAffine2DTransform.h"
#include "itkANTSNeighborhoodCorrelationImageToImageMetricv4.h"
#include "itkAffineTransform.h"
#include "itkArray.h"
#include "itkBSplineExponentialDiffeomorphicTransform.h"
#include "itkBSplineExponentialDiffeomorphicTransformParametersAdaptor.h"
#include "itkBSplineSmoothingOnUpdateDisplacementFieldTransform.h"
#include "itkBSplineSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor.h"
#include "itkBSplineSyNImageRegistrationMethod.h"
#include "itkBSplineTransform.h"
#include "itkBSplineTransformParametersAdaptor.h"
#include "itkCastImageFilter.h"
#include "itkCenteredAffineTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkCommand.h"
#include "itkComposeDisplacementFieldsImageFilter.h"
#include "itkCompositeTransform.h"
#include "itkConjugateGradientLineSearchOptimizerv4.h"
#include "itkCorrelationImageToImageMetricv4.h"
#include "itkDemonsImageToImageMetricv4.h"
#include "itkDisplacementFieldTransform.h"
#include "itkEuclideanDistancePointSetToPointSetMetricv4.h"
#include "itkEuler2DTransform.h"
#include "itkEuler3DTransform.h"
#include "itkExpectationBasedPointSetToPointSetMetricv4.h"
#include "itkGaussianExponentialDiffeomorphicTransform.h"
#include "itkGaussianExponentialDiffeomorphicTransformParametersAdaptor.h"
#include "itkGaussianSmoothingOnUpdateDisplacementFieldTransform.h"
#include "itkGaussianSmoothingOnUpdateDisplacementFieldTransformParametersAdaptor.h"
#include "itkGradientDescentOptimizerv4.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImageRegistrationMethodv4.h"
#include "itkImageToHistogramFilter.h"
#include "itkImageToImageMetricv4.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkJensenHavrdaCharvatTsallisPointSetToPointSetMetricv4.h"
#include "itkJointHistogramMutualInformationImageToImageMetricv4.h"
#include "itkLabeledPointSetToPointSetMetricv4.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkMeanSquaresPointSetToPointSetIntensityMetricv4.h"
#include "itkMultiGradientOptimizerv4.h"
#include "itkObject.h"
#include "itkObjectToObjectMetric.h"
#include "itkObjectToObjectMultiMetricv4.h"
#include "itkQuaternionRigidTransform.h"
#include "itkRegistrationParameterScalesFromPhysicalShift.h"
#include "itkSimilarity2DTransform.h"
#include "itkSimilarity3DTransform.h"
#include "itkSyNImageRegistrationMethod.h"
#include "itkTimeProbe.h"
#include "itkTimeVaryingBSplineVelocityFieldImageRegistrationMethod.h"
#include "itkTimeVaryingBSplineVelocityFieldTransformParametersAdaptor.h"
#include "itkTimeVaryingVelocityFieldImageRegistrationMethodv4.h"
#include "itkTimeVaryingVelocityFieldTransform.h"
#include "itkTimeVaryingVelocityFieldTransformParametersAdaptor.h"
#include "itkTransform.h"
#include "itkTransformFactory.h"
#include "itkTranslationTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "itkWeakPointer.h"

#include "itkantsReadWriteTransform.h"

#include "antsAllocImage.h"
#include "antsCommandLineParser.h"

#include <string>
#include <iostream>
#include <sstream>
#include <deque>
#include <iomanip>

#include "antsRegistrationCommandIterationUpdate.h"
#include "antsRegistrationOptimizerCommandIterationUpdate.h"
#include "antsDisplacementAndVelocityFieldRegistrationCommandIterationUpdate.h"

namespace ants
{
typedef itk::ants::CommandLineParser ParserType;
typedef ParserType::OptionType       OptionType;

template <typename TComputeType, unsigned VImageDimension>
class RegistrationHelper final : public itk::Object
{
public:
  /** Standard class typedefs */
  typedef RegistrationHelper            Self;
  typedef itk::Object                   Superclass;
  typedef itk::SmartPointer<Self>       Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;
  typedef itk::WeakPointer<const Self>  ConstWeakPointer;

  typedef TComputeType                           RealType;
  typedef TComputeType                           PixelType;
  typedef itk::Image<PixelType, VImageDimension> ImageType;
  typedef typename ImageType::Pointer            ImagePointer;
  typedef itk::ImageBase<VImageDimension>        ImageBaseType;
  typedef typename ImageType::SpacingType        ImageSpacingType;

  typedef itk::Array<TComputeType>                                      IntensityAndGradientArrayType;
  typedef itk::PointSet<unsigned int, VImageDimension>                  LabeledPointSetType;
  typedef typename LabeledPointSetType::Pointer                         LabeledPointSetPointer;
  typedef itk::PointSet<IntensityAndGradientArrayType, VImageDimension> IntensityPointSetType;
  typedef typename IntensityPointSetType::Pointer                       IntensityPointSetPointer;

  typedef itk::Transform<TComputeType, VImageDimension, VImageDimension> TransformType;
  typedef itk::AffineTransform<RealType, VImageDimension>                AffineTransformType;
  typedef itk::ImageRegistrationMethodv4<ImageType, ImageType, AffineTransformType, ImageType, LabeledPointSetType>
    AffineRegistrationType;
  typedef
    typename AffineRegistrationType::ShrinkFactorsPerDimensionContainerType ShrinkFactorsPerDimensionContainerType;
  typedef typename AffineTransformType::Superclass                          MatrixOffsetTransformBaseType;
  typedef typename MatrixOffsetTransformBaseType::Pointer                   MatrixOffsetTransformBasePointer;
  typedef itk::CompositeTransform<RealType, VImageDimension>                CompositeTransformType;
  typedef typename CompositeTransformType::Pointer                          CompositeTransformPointer;
  typedef itk::DisplacementFieldTransform<RealType, VImageDimension>        DisplacementFieldTransformType;
  typedef typename DisplacementFieldTransformType::Pointer                  DisplacementFieldTransformPointer;
  typedef typename DisplacementFieldTransformType::DisplacementFieldType    DisplacementFieldType;
  typedef itk::TimeVaryingVelocityFieldTransform<RealType, VImageDimension> TimeVaryingVelocityFieldTransformType;
  typedef itk::ObjectToObjectMetric<VImageDimension, VImageDimension, ImageType, RealType>        ObjectMetricType;
  typedef itk::ObjectToObjectMultiMetricv4<VImageDimension, VImageDimension, ImageType, RealType> MultiMetricType;
  typedef itk::ImageToImageMetricv4<ImageType, ImageType, ImageType, RealType>                    ImageMetricType;
  typedef itk::ImageMaskSpatialObject<VImageDimension>   ImageMaskSpatialObjectType;
  typedef typename ImageMaskSpatialObjectType::ImageType MaskImageType;

  typedef itk::InterpolateImageFunction<ImageType, RealType> InterpolatorType;

  typedef itk::ConjugateGradientLineSearchOptimizerv4Template<TComputeType> ConjugateGradientDescentOptimizerType;
  typedef itk::GradientDescentOptimizerv4Template<TComputeType>             GradientDescentOptimizerType;


  enum MetricEnumeration
  {
    CC = 0,
    MI = 1,
    Mattes = 2,
    MeanSquares = 3,
    Demons = 4,
    GC = 5,
    ICP = 6,
    PSE = 7,
    JHCT = 8,
    IGDM = 9,
    IllegalMetric = 10
  };
  enum SamplingStrategy
  {
    none = 0,    // aka dense
    regular = 1, // regularly spaced sub-sampling
    random = 2,  // irregularly spaced sub-sampling
    invalid = 17
  };

  bool
  IsPointSetMetric(const MetricEnumeration metricType) const
  {
    if (metricType == ICP || metricType == PSE || metricType == JHCT || metricType == IGDM)
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  class Metric
  {
  public:
    Metric(MetricEnumeration       metricType,
           ImageType *             fixedImage,
           ImageType *             movingImage,
           LabeledPointSetType *   fixedLabeledPointSet,
           LabeledPointSetType *   movingLabeledPointSet,
           IntensityPointSetType * fixedIntensityPointSet,
           IntensityPointSetType * movingIntensityPointSet,
           unsigned int            stageID,
           RealType                weighting,
           SamplingStrategy        samplingStrategy,
           int                     numberOfBins,
           unsigned int            radius,
           bool                    useGradientFilter,
           bool                    useBoundaryPointsOnly,
           RealType                pointSetSigma,
           unsigned int            evaluationKNeighborhood,
           RealType                alpha,
           bool                    useAnisotropicCovariances,
           RealType                samplingPercentage,
           RealType                intensityDistanceSigma,
           RealType                euclideanDistanceSigma
           )
      : m_MetricType(metricType)
      , m_FixedImage(fixedImage)
      , m_MovingImage(movingImage)
      , m_StageID(stageID)
      , m_Weighting(weighting)
      , m_SamplingStrategy(samplingStrategy)
      , m_NumberOfBins(numberOfBins)
      , m_Radius(radius)
      , m_UseGradientFilter(useGradientFilter)
      , m_FixedLabeledPointSet(fixedLabeledPointSet)
      , m_MovingLabeledPointSet(movingLabeledPointSet)
      , m_FixedIntensityPointSet(fixedIntensityPointSet)
      , m_MovingIntensityPointSet(movingIntensityPointSet)
      , m_UseBoundaryPointsOnly(useBoundaryPointsOnly)
      , m_PointSetSigma(pointSetSigma)
      , m_EvaluationKNeighborhood(evaluationKNeighborhood)
      , m_Alpha(alpha)
      , m_UseAnisotropicCovariances(useAnisotropicCovariances)
      , m_SamplingPercentage(samplingPercentage)
      , m_IntensityDistanceSigma(intensityDistanceSigma)
      , m_EuclideanDistanceSigma(euclideanDistanceSigma)
    {}

    const std::string
    GetMetricAsString() const
    {
      switch (this->m_MetricType)
      {
        case CC:
        {
          return std::string("CC");
        }
        case MI:
        {
          return std::string("JointHistogramMI");
        }
        case Mattes:
        {
          return std::string("MattesMI");
        }
        case MeanSquares:
        {
          return std::string("MeanSquares");
        }
        case Demons:
        {
          return std::string("Demons");
        }
        case GC:
        {
          return std::string("GC");
        }
        case ICP:
        {
          return std::string("ICP");
        }
        case PSE:
        {
          return std::string("PSE");
        }
        case JHCT:
        {
          return std::string("JHCT");
        }
        case IGDM:
        {
          return std::string("IGDM");
        }
        default:
        {
        }
        break;
      }
      return std::string("Illegal");
    }

    MetricEnumeration m_MetricType;

    // Variables for image metrics

    ImagePointer     m_FixedImage;
    ImagePointer     m_MovingImage;
    unsigned int     m_StageID;
    RealType         m_Weighting;
    SamplingStrategy m_SamplingStrategy;
    int              m_NumberOfBins;
    unsigned int     m_Radius; // Only for CC metric
    bool             m_UseGradientFilter;

    // Variables for point-set metrics

    LabeledPointSetPointer   m_FixedLabeledPointSet;
    LabeledPointSetPointer   m_MovingLabeledPointSet;
    IntensityPointSetPointer m_FixedIntensityPointSet;
    IntensityPointSetPointer m_MovingIntensityPointSet;

    bool     m_UseBoundaryPointsOnly;
    RealType m_PointSetSigma;             // Only for PSE,JHCT metrics
    RealType m_EvaluationKNeighborhood;   // Only for PSE,JHCT metrics
    RealType m_Alpha;                     // Only for JHCT metric
    bool     m_UseAnisotropicCovariances; // Only for JHCT metric

    RealType m_SamplingPercentage; // Only for PSE,JHCT metrics

    RealType m_IntensityDistanceSigma; // Only for IGDM metric
    RealType m_EuclideanDistanceSigma; // Only for IGDM metric
  };

  typedef std::deque<Metric> MetricListType;

  enum XfrmMethod
  {
    Rigid = 0,
    Affine = 1,
    CompositeAffine = 2,
    Similarity = 3,
    Translation = 4,
    BSpline = 5,
    GaussianDisplacementField = 6,
    BSplineDisplacementField = 7,
    TimeVaryingVelocityField = 8,
    TimeVaryingBSplineVelocityField = 9,
    SyN = 10,
    BSplineSyN = 11,
    Exponential = 12,
    BSplineExponential = 13,
    UnknownXfrm = 14
  };

  class TransformMethod
  {
  public:
    TransformMethod()
      : m_XfrmMethod(Rigid)
      , m_GradientStep(0)
      , m_UpdateFieldVarianceInVarianceSpace(0.0)
      , m_TotalFieldVarianceInVarianceSpace(0.0)
      , m_SplineOrder(3)
      , m_UpdateFieldTimeSigma(0.0)
      , m_TotalFieldTimeSigma(0.0)
      , m_NumberOfTimeIndices(0)
      , m_NumberOfTimePointSamples(4)
      , m_VelocityFieldVarianceInVarianceSpace(0.0)
    {}

    std::string
    XfrmMethodAsString() const
    {
      switch (this->m_XfrmMethod)
      {
        case Rigid:
        {
          return std::string("Rigid");
        }
        case Affine:
        {
          return std::string("Affine");
        }
        case CompositeAffine:
        {
          return std::string("CompositeAffine");
        }
        case Similarity:
        {
          return std::string("Similarity");
        }
        case Translation:
        {
          return std::string("Translation");
        }
        case BSpline:
        {
          return std::string("BSpline");
        }
        case GaussianDisplacementField:
        {
          return std::string("GaussianDisplacementField");
        }
        case BSplineDisplacementField:
        {
          return std::string("BSplineDisplacementField");
        }
        case TimeVaryingVelocityField:
        {
          return std::string("TimeVaryingVelocityField");
        }
        case TimeVaryingBSplineVelocityField:
        {
          return std::string("TimeVaryingBSplineVelocityField");
        }
        case SyN:
        {
          return std::string("SyN");
        }
        case BSplineSyN:
        {
          return std::string("BSplineSyN");
        }
        case Exponential:
        {
          return std::string("Exponential");
        }
        case BSplineExponential:
        {
          return std::string("BSplineExponential");
        }
        case UnknownXfrm:
          return std::string("UnknownXfrm");
      }
      return std::string("Impossible");
    }

    XfrmMethod m_XfrmMethod;
    // all transforms
    RealType m_GradientStep;
    // BSpline
    std::vector<unsigned int> m_MeshSizeAtBaseLevel;
    // GaussianDisplacementField
    RealType m_UpdateFieldVarianceInVarianceSpace;
    RealType m_TotalFieldVarianceInVarianceSpace;
    // BSplineDisplacementField
    std::vector<unsigned int> m_TotalFieldMeshSizeAtBaseLevel;
    std::vector<unsigned int> m_UpdateFieldMeshSizeAtBaseLevel;
    unsigned int              m_SplineOrder; // also anything B-spline
    // TimeVaryingVelocityField
    RealType     m_UpdateFieldTimeSigma;
    RealType     m_TotalFieldTimeSigma;
    unsigned int m_NumberOfTimeIndices;
    // TimeVaryingBSplineVelocityField
    std::vector<unsigned int> m_VelocityFieldMeshSize;
    unsigned int              m_NumberOfTimePointSamples;
    // Exponential
    RealType m_VelocityFieldVarianceInVarianceSpace;
    // BSplineExponential
    std::vector<unsigned int> m_VelocityFieldMeshSizeAtBaseLevel;
  };

  typedef std::deque<TransformMethod> TransformMethodListType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(RegistrationHelper);

  /** Dimension of the image.  This constant is used by functions that are
   * templated over image type (as opposed to being templated over pixel type
   * and dimension) when they need compile time access to the dimension of
   * the image. */
  static constexpr unsigned int ImageDimension = VImageDimension;

  /**
   * add a metric, corresponding to the registration stage
   */
  void
  AddMetric(MetricEnumeration       metricType,
            ImageType *             fixedImage,
            ImageType *             movingImage,
            LabeledPointSetType *   fixedLabeledPointSet,
            LabeledPointSetType *   movingLabeledPointSet,
            IntensityPointSetType * fixedIntensityPointSet,
            IntensityPointSetType * movingIntensityPointSet,
            unsigned int            stageID,
            RealType                weighting,
            SamplingStrategy        samplingStrategy,
            int                     numberOfBins,
            unsigned int            radius,
            bool                    useGradientFilter,
            bool                    useBoundaryPointsOnly,
            RealType                pointSetSigma,
            unsigned int            evaluationKNeighborhood,
            RealType                alpha,
            bool                    useAnisotropicCovariances,
            RealType                samplingPercentage,
            RealType                intensityDistanceSigma,
            RealType                euclideanDistanceSigma);

  /** For backwards compatibility */
  inline void
  AddMetric(MetricEnumeration metricType,
            ImageType *       fixedImage,
            ImageType *       movingImage,
            unsigned int      stageID,
            RealType          weighting,
            SamplingStrategy  samplingStrategy,
            int               numberOfBins,
            unsigned int      radius,
            bool              useGradientFilter,
            RealType          samplingPercentage)
  {
    this->AddMetric(metricType,
                    fixedImage,
                    movingImage,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    stageID,
                    weighting,
                    samplingStrategy,
                    numberOfBins,
                    radius,
                    useGradientFilter,
                    false,
                    1.0,
                    50,
                    1.1,
                    false,
                    samplingPercentage,
                    std::sqrt(5),
                    std::sqrt(5));
  }


  /**
   * Get set of metrics per stage.  If we have more than one, we have
   * to use the MultiMetricType
   */
  MetricListType
  GetMetricListPerStage(unsigned int stageID);

  /**
   * return the enumerated Metric type based on a string representation
   */
  MetricEnumeration
  StringToMetricType(const std::string & str) const;

  /**
   * return the enumerated transform method specified by str
   */
  XfrmMethod
  StringToXfrmMethod(const std::string & str) const;

  /**
   * set the fixed initial transform.
   */
  void
  SetFixedInitialTransform(const TransformType * initialTransform);

  /**
   * set the moving initial transform.
   */
  void
  SetMovingInitialTransform(const TransformType * initialTransform);

  /**
   * restore the state transform and set that as an initial transform.
   */
  void
  SetRestoreStateTransform(const TransformType * initialTransform);

  /**
   * add a rigid transform
   */
  void
  AddRigidTransform(RealType GradientStep);

  /**
   * add an affine transform
   */
  void
  AddAffineTransform(RealType GradientStep);

  /**
   * add a composite affine transform
   */
  void
  AddCompositeAffineTransform(RealType GradientStep);

  /**
   * add a similarity transform
   */
  void
  AddSimilarityTransform(RealType GradientStep);

  /**
   * add a translation transform
   */
  void
  AddTranslationTransform(RealType GradientStep);

  /**
   * add a spline transform
   */
  void
  AddBSplineTransform(RealType GradientStep, std::vector<unsigned int> & MeshSizeAtBaseLevel);

  /**
   * add gaussian displacement transform
   */
  void
  AddGaussianDisplacementFieldTransform(RealType GradientStep,
                                        RealType UpdateFieldVarianceInVarianceSpace,
                                        RealType TotalFieldVarianceInVarianceSpace);

  /**
   * add bspline displacement transform
   */
  void
  AddBSplineDisplacementFieldTransform(RealType                    GradientStep,
                                       std::vector<unsigned int> & UpdateFieldMeshSizeAtBaseLevel,
                                       std::vector<unsigned int> & TotalFieldMeshSizeAtBaseLevel,
                                       unsigned int                SplineOrder);

  /**
   * add a time varying velocity field transform
   */
  void
  AddTimeVaryingVelocityFieldTransform(RealType     GradientStep,
                                       unsigned int NumberOfTimeIndices,
                                       RealType     UpdateFieldVarianceInVarianceSpace,
                                       RealType     UpdateFieldTimeSigma,
                                       RealType     TotalFieldVarianceInVarianceSpace,
                                       RealType     TotalFieldTimeSigma);

  /**
   * add a time varying b spline velocity field transform
   */
  void
  AddTimeVaryingBSplineVelocityFieldTransform(RealType                  GradientStep,
                                              std::vector<unsigned int> VelocityFieldMeshSize,
                                              unsigned int              NumberOfTimePointSamples,
                                              unsigned int              SplineOrder);

  /**
   * add a SyN transform
   */
  void
  AddSyNTransform(RealType GradientStep,
                  RealType UpdateFieldVarianceInVarianceSpace,
                  RealType TotalFieldVarianceInVarianceSpace);

  /**
   * add a B-spline SyN transform
   */
  void
  AddBSplineSyNTransform(RealType                    GradientStep,
                         std::vector<unsigned int> & UpdateFieldMeshSizeAtBaseLevel,
                         std::vector<unsigned int> & TotalFieldMeshSizeAtBaseLevel,
                         unsigned int                SplineOrder);

  /**
   * add an exponential transform
   */
  void
  AddExponentialTransform(RealType     GradientStep,
                          RealType     UpdateFieldVarianceInVarianceSpace,
                          RealType     VelocityFieldVarianceInVarianceSpace,
                          unsigned int NumberOfIntegrationSteps);

  /**
   * add a B-spline exponential transform
   */
  void
  AddBSplineExponentialTransform(RealType                    GradientStep,
                                 std::vector<unsigned int> & UpdateFieldMeshSizeAtBaseLevel,
                                 std::vector<unsigned int> & VelocityFieldMeshSizeAtBaseLevel,
                                 unsigned int                NumberOfIntegrationSteps,
                                 unsigned int                SplineOrder);
  /**
   * Add the collected output prefix
   */
  void
  SetOutputPrefix(const std::string & outputPrefix);

  /**
   * Add the collected iterations list
   */
  void
  SetIterations(const std::vector<std::vector<unsigned int>> & Iterations);

  /**
   * Add the collected convergence thresholds
   */
  void
  SetConvergenceThresholds(const std::vector<RealType> & thresholds);

  /**
   * Add the collected convergence window sizes
   */
  void
  SetConvergenceWindowSizes(const std::vector<unsigned int> & windowSizes);

  /**
   * Add the collected smoothing sigmas list
   */
  void
  SetSmoothingSigmas(const std::vector<std::vector<float>> & SmoothingSigmas);

  /**
   * Add the restrict deformation optimizer weights
   */
  void
  SetRestrictDeformationOptimizerWeights(const std::vector<std::vector<RealType>> & restrictDeformationWeights);

  /**
   * Add the collected bool smoothing sigmas in voxel units list
   */
  void
  SetSmoothingSigmasAreInPhysicalUnits(const std::vector<bool> & SmoothingSigmasAreInPhysicalUnits);

  /**
   * Add the collected shrink factors list
   */
  void
  SetShrinkFactors(const std::vector<std::vector<unsigned int>> & ShrinkFactors);

  /**
   * Given an image and a specified shrink factor, calculate the shrink factor for
   * each dimension.  The heuristic we use is to apply the specified shrink factor
   * to the image dimension with the highest resolution (i.e. smallest spacing).
   * Then, for each other dimension, we find the factor ( < specified factor ) which
   * makes the subsampled image as close to isotropic (in terms of spacing) as possible.
   */
  ShrinkFactorsPerDimensionContainerType
  CalculateShrinkFactorsPerDimension(unsigned int, ImageSpacingType);

  /**
   * turn on histogram matching of the input images
   */
  itkSetMacro(UseHistogramMatching, bool);
  itkGetConstMacro(UseHistogramMatching, bool);
  itkBooleanMacro(UseHistogramMatching);

  /**
   * turn on the option that lets you estimate the learning rate step size only at the beginning of each level.
   * useful as a second stage of fine-scale registration.
   */
  itkSetMacro(DoEstimateLearningRateAtEachIteration, bool);
  itkGetConstMacro(DoEstimateLearningRateAtEachIteration, bool);
  itkBooleanMacro(DoEstimateLearningRateAtEachIteration);

  /**
   * turn on the option that prints the CC similarity measure
   * between the full-size fixed and moving input images at each iteraton.
   */
  itkSetMacro(PrintSimilarityMeasureInterval, unsigned int);
  itkGetConstMacro(PrintSimilarityMeasureInterval, unsigned int);

  /**
   * turn on the option that writes the output volume in intervals of iterations.
   */
  itkSetMacro(WriteIntervalVolumes, unsigned int);
  itkGetConstMacro(WriteIntervalVolumes, unsigned int);

  /**
   * turn on the option that cause the direct initialization of the linear transforms at each stage.
   */
  itkSetMacro(InitializeTransformsPerStage, bool);
  itkGetConstMacro(InitializeTransformsPerStage, bool);
  itkBooleanMacro(InitializeTransformsPerStage);

  /**
   * Set a constant random seed with an int != 0
   */
  itkSetMacro(RegistrationRandomSeed, int);
  itkGetConstMacro(RegistrationRandomSeed, int);

  /**
   * turn on winsorize image intensity normalization
   */
  void
  SetWinsorizeImageIntensities(bool Winsorize, float LowerQuantile = 0.0, float UpperQuantile = 1.0);

  itkGetModifiableObjectMacro(CompositeTransform, CompositeTransformType);
  itkGetModifiableObjectMacro(RegistrationState, CompositeTransformType);
  /**
   * Set/Get the interpolator.  Linear is default.
   */
  itkSetObjectMacro(Interpolator, InterpolatorType);
  itkGetModifiableObjectMacro(Interpolator, InterpolatorType);

  /**
   *  Get the Warped Image & Inverse Warped Images
   */
  typename ImageType::Pointer
  GetWarpedImage() const;

  typename ImageType::Pointer
  GetInverseWarpedImage() const;

  /**
   * Set the fixed/moving image masks with a spatial object
   */
  void
  AddFixedImageMask(typename ImageMaskSpatialObjectType::Pointer & fixedImageMask)
  {
    this->m_FixedImageMasks.push_back(fixedImageMask);
  }

  void
  AddMovingImageMask(typename ImageMaskSpatialObjectType::Pointer & movingImageMask)
  {
    this->m_MovingImageMasks.push_back(movingImageMask);
  }

  /**
   * Set the fixed/moving mask image. this will be used to instantiate
   * ImageMaskSpatialObject masks.
   */
  void
  AddFixedImageMask(typename MaskImageType::Pointer & fixedImageMask);
  void
  AddMovingImageMask(typename MaskImageType::Pointer & movingImageMask);

  /**
   * Collapse a composite transform by adjacent linear or displacement field transforms.
   */
  CompositeTransformPointer
  CollapseCompositeTransform(const CompositeTransformType *);

  void
  ApplyCompositeLinearTransformToImageHeader(const CompositeTransformType *,
                                             ImageBaseType * const,
                                             const bool applyInverse);

  /**
   * Collapse a composite transform composed of displacement field transforms.
   * We return a composite transform since we don't combine mixed displacement field
   * transforms (i.e., transforms that do and do not have inverses).
   */
  typename CompositeTransformType::Pointer
  CollapseDisplacementFieldTransforms(const CompositeTransformType *);

  /**
   * Collapse a composite linear transform to a generic affine transform.
   */
  typename AffineTransformType::Pointer
  CollapseLinearTransforms(const CompositeTransformType *);

  /**
   *
   */
  template <typename TTransformType>
  bool
  InitializeWithPreviousLinearTransform(const CompositeTransformType *,
                                        const std::string,
                                        typename TTransformType::Pointer &);

  /**
   * Compute approximate mesh size for a specified isotropic knot spacing
   */
  std::vector<unsigned int>
  CalculateMeshSizeForSpecifiedKnotSpacing(ImageBaseType * const, const RealType, const unsigned int);

  /**
   * Do the registration. Will return EXIT_FAILURE if there is any
   * problem completing the registration.
   */
  int
  DoRegistration();

  /**
   * print out the internal registration helper state
   */
  void
  PrintState() const;

  void
  SetLogStream(std::ostream & logStream)
  {
    this->m_LogStream = &logStream;
  }

protected:
  RegistrationHelper();
  ~RegistrationHelper() override;

private:
  typename itk::ImageBase<VImageDimension>::Pointer
  GetShrinkImageOutputInformation(
    const itk::ImageBase<VImageDimension> * inputImageInformation,
    const typename RegistrationHelper<TComputeType, VImageDimension>::ShrinkFactorsPerDimensionContainerType &
      shrinkFactorsPerDimensionForCurrentLevel) const;

  int
  ValidateParameters();

  std::ostream &
  Logger() const
  {
    return *m_LogStream;
  }

  template <typename RegistrationMethodType>
  typename RegistrationMethodType::Pointer
  PrepareRegistrationMethod(
    CompositeTransformType *                                               compositeTransform,
    const unsigned int                                                     currentStageNumber,
    const unsigned int                                                     parametersDimensionSize,
    std::vector<typename RegistrationMethodType::FixedImageType::Pointer>  preprocessedFixedImagesPerStage,
    std::vector<typename RegistrationMethodType::MovingImageType::Pointer> preprocessedMovingImagesPerStage,
    std::vector<typename RegistrationMethodType::PointSetType::Pointer>    fixedPointSetsPerStage,
    std::vector<typename RegistrationMethodType::PointSetType::Pointer>    movingPointSetsPerStage,
    const MetricListType                                                   stageMetricList,
    ObjectMetricType *                                                     singleMetric,
    MultiMetricType *                                                      multiMetric,
    ConjugateGradientDescentOptimizerType *                                optimizer,
    const unsigned int                                                     numberOfLevels,
    const std::vector<ShrinkFactorsPerDimensionContainerType>              shrinkFactorsPerDimensionForAllLevels,
    const typename RegistrationMethodType::SmoothingSigmasArrayType        smoothingSigmasPerLevel,
    typename AffineRegistrationType::MetricSamplingStrategyEnum            metricSamplingStrategy,
    const float                                                            samplingPercentage)
  {
    typename RegistrationMethodType::Pointer                     registrationMethod = RegistrationMethodType::New();
    typedef typename RegistrationMethodType::OutputTransformType RegistrationMethodTransformType;

    if (this->m_RegistrationRandomSeed != 0)
    {
      registrationMethod->MetricSamplingReinitializeSeed(this->m_RegistrationRandomSeed);
    }

    for (unsigned int n = 0; n < stageMetricList.size(); n++)
    {
      if (!this->IsPointSetMetric(stageMetricList[n].m_MetricType))
      {
        registrationMethod->SetFixedImage(n, preprocessedFixedImagesPerStage[n]);
        registrationMethod->SetMovingImage(n, preprocessedMovingImagesPerStage[n]);
      }
      else
      {
        registrationMethod->SetFixedPointSet(n, fixedPointSetsPerStage[n]);
        registrationMethod->SetMovingPointSet(n, movingPointSetsPerStage[n]);
      }
    }

    if (multiMetric)
    {
      registrationMethod->SetMetric(multiMetric);
    }
    else
    {
      registrationMethod->SetMetric(singleMetric);
    }

    registrationMethod->SetNumberOfLevels(numberOfLevels);

    for (unsigned int level = 0; level < numberOfLevels; ++level)
    {
      registrationMethod->SetShrinkFactorsPerDimension(level, shrinkFactorsPerDimensionForAllLevels[level]);
    }
    registrationMethod->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);
    registrationMethod->SetSmoothingSigmasAreSpecifiedInPhysicalUnits(
      this->m_SmoothingSigmasAreInPhysicalUnits[currentStageNumber]);
    registrationMethod->SetMetricSamplingStrategy(
      static_cast<typename RegistrationMethodType::MetricSamplingStrategyEnum>(metricSamplingStrategy));
    registrationMethod->SetMetricSamplingPercentage(samplingPercentage);

    if (this->m_RestrictDeformationOptimizerWeights.size() > currentStageNumber)
    {
      if (this->m_RestrictDeformationOptimizerWeights[currentStageNumber].size() == parametersDimensionSize)
      {
        typename RegistrationMethodType::OptimizerWeightsType optimizerWeights(parametersDimensionSize);
        for (unsigned int d = 0; d < parametersDimensionSize; d++)
        {
          optimizerWeights[d] = this->m_RestrictDeformationOptimizerWeights[currentStageNumber][d];
        }
        registrationMethod->SetOptimizerWeights(optimizerWeights);
      }
    }

    registrationMethod->SetOptimizer(optimizer);

    typename RegistrationMethodTransformType::Pointer currentTransform = RegistrationMethodTransformType::New();

    std::string            t = currentTransform->GetNameOfClass();
    std::string            s = "Transform";
    std::string::size_type index = t.find(s);

    if (index != std::string::npos)
    {
      t.erase(index, s.length());
    }

    if (compositeTransform->GetNumberOfTransforms() > 0)
    {
      if (this->m_InitializeTransformsPerStage)
      {
        const unsigned int numOfTransforms = compositeTransform->GetNumberOfTransforms();
        this->Logger() << "Current number of transforms in the composite transform: " << numOfTransforms << std::endl;
        for (unsigned int i = 0; i < numOfTransforms; i++)
        {
          this->Logger() << i + 1 << ") " << compositeTransform->GetNthTransform(i)->GetNameOfClass() << std::endl;
        }

        if (this->InitializeWithPreviousLinearTransform<RegistrationMethodTransformType>(
              compositeTransform, t.c_str(), currentTransform))
        {
          this->Logger() << "Registration process is run using direct initialization!" << std::endl;
          compositeTransform->RemoveTransform(); // Remove previous initial transform,
                                                 // since it is included in current results.
          registrationMethod->SetInitialTransform(currentTransform);
        }
      }
    }

    if (compositeTransform->GetNumberOfTransforms() > 0)
    {
      registrationMethod->SetMovingInitialTransform(compositeTransform);
    }
    if (this->m_FixedInitialTransform->GetNumberOfTransforms() > 0)
    {
      registrationMethod->SetFixedInitialTransform(this->m_FixedInitialTransform);
    }

    return registrationMethod;
  }

  template <typename RegistrationMethodType>
  int
  AddLinearTransformToCompositeTransform(
    CompositeTransformType *                                               compositeTransform,
    const unsigned int                                                     currentStageNumber,
    const unsigned int                                                     parametersDimensionSize,
    std::vector<typename RegistrationMethodType::FixedImageType::Pointer>  preprocessedFixedImagesPerStage,
    std::vector<typename RegistrationMethodType::MovingImageType::Pointer> preprocessedMovingImagesPerStage,
    std::vector<typename RegistrationMethodType::PointSetType::Pointer>    fixedPointSetsPerStage,
    std::vector<typename RegistrationMethodType::PointSetType::Pointer>    movingPointSetsPerStage,
    const MetricListType                                                   stageMetricList,
    ObjectMetricType *                                                     singleMetric,
    MultiMetricType *                                                      multiMetric,
    ConjugateGradientDescentOptimizerType *                                optimizer,
    const unsigned int                                                     numberOfLevels,
    const std::vector<ShrinkFactorsPerDimensionContainerType>              shrinkFactorsPerDimensionForAllLevels,
    const typename RegistrationMethodType::SmoothingSigmasArrayType        smoothingSigmasPerLevel,
    typename AffineRegistrationType::MetricSamplingStrategyEnum            metricSamplingStrategy,
    const float                                                            samplingPercentage)
  {
    typename RegistrationMethodType::Pointer registrationMethod =
      this->PrepareRegistrationMethod<RegistrationMethodType>(compositeTransform,
                                                              currentStageNumber,
                                                              parametersDimensionSize,
                                                              preprocessedFixedImagesPerStage,
                                                              preprocessedMovingImagesPerStage,
                                                              fixedPointSetsPerStage,
                                                              movingPointSetsPerStage,
                                                              stageMetricList,
                                                              singleMetric,
                                                              multiMetric,
                                                              optimizer,
                                                              numberOfLevels,
                                                              shrinkFactorsPerDimensionForAllLevels,
                                                              smoothingSigmasPerLevel,
                                                              metricSamplingStrategy,
                                                              samplingPercentage);

    typedef antsRegistrationCommandIterationUpdate<RegistrationMethodType> TransformCommandType;
    typename TransformCommandType::Pointer transformObserver = TransformCommandType::New();
    transformObserver->SetLogStream(*this->m_LogStream);
    transformObserver->SetNumberOfIterations(this->m_Iterations[currentStageNumber]);
    registrationMethod->AddObserver(itk::IterationEvent(), transformObserver);
    registrationMethod->AddObserver(itk::InitializeEvent(), transformObserver);

    try
    {
      typedef typename RegistrationMethodType::OutputTransformType RegistrationMethodTransformType;
      typename RegistrationMethodTransformType::Pointer currentTransform = RegistrationMethodTransformType::New();

      this->Logger() << std::endl
                     << "*** Running " << currentTransform->GetNameOfClass() << " registration ***" << std::endl
                     << std::endl;
      transformObserver->Execute(registrationMethod, itk::StartEvent());
      registrationMethod->Update();
    }
    catch (const itk::ExceptionObject & e)
    {
      this->Logger() << "Exception caught: " << e << std::endl;
      return EXIT_FAILURE;
    }


    compositeTransform->AddTransform(registrationMethod->GetModifiableTransform());

    return EXIT_SUCCESS;
  }


  typename CompositeTransformType::Pointer m_CompositeTransform;
  typename CompositeTransformType::Pointer m_RegistrationState;
  typename CompositeTransformType::Pointer m_FixedInitialTransform;

  std::vector<typename ImageMaskSpatialObjectType::Pointer> m_FixedImageMasks;
  std::vector<typename ImageMaskSpatialObjectType::Pointer> m_MovingImageMasks;

  typename InterpolatorType::Pointer m_Interpolator;

  unsigned int                           m_NumberOfStages;
  MetricListType                         m_Metrics;
  TransformMethodListType                m_TransformMethods;
  std::string                            m_OutputPrefix;
  std::vector<std::vector<unsigned int>> m_Iterations;
  std::vector<RealType>                  m_ConvergenceThresholds;
  std::vector<unsigned int>              m_ConvergenceWindowSizes;
  std::vector<std::vector<float>>        m_SmoothingSigmas;
  std::vector<bool>                      m_SmoothingSigmasAreInPhysicalUnits;
  std::vector<std::vector<RealType>>     m_RestrictDeformationOptimizerWeights;
  std::vector<std::vector<unsigned int>> m_ShrinkFactors;
  bool                                   m_UseHistogramMatching;
  bool                                   m_WinsorizeImageIntensities;
  bool                                   m_DoEstimateLearningRateAtEachIteration;
  RealType                               m_LowerQuantile;
  RealType                               m_UpperQuantile;
  std::ostream *                         m_LogStream;

  int m_RegistrationRandomSeed;

  unsigned int m_PrintSimilarityMeasureInterval;
  unsigned int m_WriteIntervalVolumes;
  bool         m_InitializeTransformsPerStage;
  bool         m_AllPreviousTransformsAreLinear;
};

// ##########################################################################
// ##########################################################################

/**
 * Transform traits to generalize the rigid transform
 */
template <typename TComputeType, unsigned int ImageDimension>
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

template <typename TComputeType, unsigned int ImageDimension>
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

// ##########################################################################
// ##########################################################################

// Provide common way of reading transforms.
template <typename TComputeType, unsigned VImageDimension>
typename ants::RegistrationHelper<TComputeType, VImageDimension>::CompositeTransformType::Pointer
GetCompositeTransformFromParserOption(typename ParserType::Pointer &           parser,
                                      typename ParserType::OptionType::Pointer initialTransformOption,
                                      std::vector<bool> &                      derivedTransforms,
                                      bool                                     useStaticCastForR = false)
{
  typedef typename ants::RegistrationHelper<TComputeType, VImageDimension> RegistrationHelperType;
  typedef typename RegistrationHelperType::CompositeTransformType          CompositeTransformType;
  typename CompositeTransformType::Pointer compositeTransform = CompositeTransformType::New();

  typedef typename RegistrationHelperType::ImageType ImageType;
  typename ImageType::Pointer                        fixedImage = nullptr;
  typename ImageType::Pointer                        movingImage = nullptr;

  bool                                                       verbose = false;
  typename itk::ants::CommandLineParser::OptionType::Pointer verboseOption = parser->GetOption("verbose");
  if (verboseOption && verboseOption->GetNumberOfFunctions())
  {
    verbose = parser->Convert<bool>(verboseOption->GetFunction(0)->GetName());
  }

  std::deque<std::string> transformNames;
  std::deque<std::string> transformTypes;
  derivedTransforms.resize(initialTransformOption->GetNumberOfFunctions());
  for (unsigned int n = 0; n < initialTransformOption->GetNumberOfFunctions(); n++)
  {
    std::string initialTransformName;
    bool        useInverse = false;
    bool        calculatedTransformFromImages = false;

    derivedTransforms[n] = false;

    if (initialTransformOption->GetFunction(n)->GetNumberOfParameters() == 0)
    {
      initialTransformName = initialTransformOption->GetFunction(n)->GetName();
    }
    else if (initialTransformOption->GetFunction(n)->GetNumberOfParameters() > 2)
    {
      std::string fixedImageFileName = initialTransformOption->GetFunction(n)->GetParameter(0);
      ReadImage<ImageType>(fixedImage, fixedImageFileName.c_str());
      if (fixedImage)
      {
        std::string movingImageFileName = initialTransformOption->GetFunction(n)->GetParameter(1);
        ReadImage<ImageType>(movingImage, movingImageFileName.c_str());
      }

      // initialization feature types:
      //  0: image centers
      //  1: center of (intensity) mass
      //  2: image origins

      unsigned short initializationFeature =
        parser->Convert<unsigned short>(initialTransformOption->GetFunction(n)->GetParameter(2));

      typedef typename RigidTransformTraits<TComputeType, VImageDimension>::TransformType TransformType;

      typename TransformType::Pointer transform = TransformType::New();

      if (initializationFeature == 0 || initializationFeature == 1)
      {
        typedef itk::CenteredTransformInitializer<TransformType, ImageType, ImageType> TransformInitializerType;
        typename TransformInitializerType::Pointer initializer = TransformInitializerType::New();

        initializer->SetTransform(transform);

        initializer->SetFixedImage(fixedImage);
        initializer->SetMovingImage(movingImage);

        if (initializationFeature == 0)
        {
          initializer->GeometryOn();
          initialTransformName = std::string("Image center alignment using ");
        }
        else
        {
          initializer->MomentsOn();
          initialTransformName = std::string("Center of mass alignment using ");
        }
        initializer->InitializeTransform();
      }
      else
      {
        initialTransformName = std::string("Image origin alignment using ");

        typename TransformType::InputPointType   rotationCenter;
        typename TransformType::OutputVectorType translationVector;

        typename ImageType::PointType fixedOrigin = fixedImage->GetOrigin();
        typename ImageType::PointType movingOrigin = movingImage->GetOrigin();

        for (unsigned int d = 0; d < VImageDimension; d++)
        {
          rotationCenter[d] = fixedOrigin[d];
          translationVector[d] = movingOrigin[d] - fixedOrigin[d];
        }
        transform->SetCenter(rotationCenter);
        transform->SetTranslation(translationVector);
      }

      initialTransformName += std::string("fixed image: ") + initialTransformOption->GetFunction(n)->GetParameter(0) +
                              std::string(" and moving image: ") +
                              initialTransformOption->GetFunction(n)->GetParameter(1);

      compositeTransform->AddTransform(transform);

      calculatedTransformFromImages = true;
      derivedTransforms[n] = true;

      transformNames.push_back(initialTransformName);
      transformTypes.push_back(transform->GetNameOfClass());
    }

    if (!calculatedTransformFromImages)
    {
      if (initialTransformOption->GetFunction(n)->GetNumberOfParameters() == 0)
      {
        initialTransformName = initialTransformOption->GetFunction(n)->GetName();
        useInverse = false;
      }
      else
      {
        initialTransformName = initialTransformOption->GetFunction(n)->GetParameter(0);
        if (initialTransformOption->GetFunction(n)->GetNumberOfParameters() > 1)
        {
          useInverse = parser->Convert<bool>(initialTransformOption->GetFunction(n)->GetParameter(1));
        }
      }
    }

    static bool MatOffRegistered(false); // Only register once for each template dimension.
    if (!MatOffRegistered)
    {
      MatOffRegistered = true;
      // Register the matrix offset transform base class to the
      // transform factory for compatibility with the current ANTs.
      typedef itk::MatrixOffsetTransformBase<TComputeType, VImageDimension, VImageDimension> MatrixOffsetTransformType;
      itk::TransformFactory<MatrixOffsetTransformType>::RegisterTransform();
    }

    if (!calculatedTransformFromImages)
    {
      typedef typename RegistrationHelperType::TransformType TransformType;
      typename TransformType::Pointer                        initialTransform;
      if (std::strcmp(initialTransformName.c_str(), "identity") == 0 ||
          std::strcmp(initialTransformName.c_str(), "Identity") == 0)
      {
        typedef itk::MatrixOffsetTransformBase<TComputeType, VImageDimension, VImageDimension>
                                                    MatrixOffsetTransformType;
        typename MatrixOffsetTransformType::Pointer identityTransform = MatrixOffsetTransformType::New();
        identityTransform->SetIdentity();

        initialTransform = identityTransform;
      }
      else
      {
        initialTransform =
          itk::ants::ReadTransform<TComputeType, VImageDimension>(initialTransformName, useStaticCastForR);
      }
      if (initialTransform.IsNull())
      {
        if (verbose)
        {
          std::cout << "Can't read initial transform " << initialTransformName << std::endl;
        }
        return nullptr;
      }
      if (useInverse)
      {
        initialTransform = dynamic_cast<TransformType *>(initialTransform->GetInverseTransform().GetPointer());
        if (initialTransform.IsNull())
        {
          if (verbose)
          {
            std::cout << "Inverse does not exist for " << initialTransformName << std::endl;
          }
          return nullptr;
        }
        initialTransformName = std::string("inverse of ") + initialTransformName;
      }
      static const std::string CompositeTransformID("CompositeTransform");
      if (initialTransform->GetNameOfClass() == CompositeTransformID)
      {
        const typename itk::CompositeTransform<TComputeType, VImageDimension>::ConstPointer tempComp =
          dynamic_cast<const itk::CompositeTransform<TComputeType, VImageDimension> *>(initialTransform.GetPointer());
        for (unsigned int i = 0; i < tempComp->GetNumberOfTransforms(); ++i)
        {
          std::stringstream tempstream;
          tempstream << initialTransformName << "[" << i << "]";

          compositeTransform->AddTransform(tempComp->GetNthTransform(i));
          transformNames.push_back(tempstream.str());
          transformTypes.push_back(tempComp->GetNthTransform(i)->GetNameOfClass());
        }
      }
      else
      {
        compositeTransform->AddTransform(initialTransform);
        transformNames.push_back(initialTransformName);
        transformTypes.push_back(initialTransform->GetNameOfClass());
      }
    }
  }
  if (verbose)
  {
    std::cout << "=============================================================================" << std::endl;
    std::cout << "The composite transform comprises the following transforms (in order): " << std::endl;
    for (unsigned int n = 0; n < transformNames.size(); n++)
    {
      std::cout << "  " << n + 1 << ". " << transformNames[n] << " (type = " << transformTypes[n] << ")" << std::endl;
    }
    std::cout << "=============================================================================" << std::endl;
  }
  return compositeTransform;
}

// ####################
} // namespace ants

#include "itkantsRegistrationHelper.hxx"

#endif
