/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkExpectationBasedPointSetRegistrationFunction_h_
#define _itkExpectationBasedPointSetRegistrationFunction_h_

#include "itkPointSetFunction.h"
#include "itkGaussianMembershipFunction.h"
#include "itkKdTreeGenerator.h"
#include "itkListSample.h"
#include "itkMatrix.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkMeshSource.h"
#include "itkPointSet.h"
#include "itkVector.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkPDEDeformableRegistrationFunction.h"
#include "itkPoint.h"
#include "itkCovariantVector.h"
#include "itkInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkCentralDifferenceImageFunction.h"

namespace itk
{
/**
 * \class ExpectationBasedPointSetRegistrationFunction
 *
 * This class encapsulate the PDE which drives the demons registration
 * algorithm. It is used by ExpectationBasedPointSetRegistrationFilter to compute the
 * output deformation field which will map a moving image onto a
 * a fixed image.
 *
 * Non-integer moving image values are obtained by using
 * interpolation. The default interpolator is of type
 * LinearInterpolateImageFunction. The user may set other
 * interpolators via method SetMovingImageInterpolator. Note that the input
 * interpolator must derive from baseclass InterpolateImageFunction.
 *
 * This class is templated over the fixed image type, moving image type,
 * and the deformation field type.
 *
 * \warning This filter assumes that the fixed image type, moving image type
 * and deformation field type all have the same number of dimensions.
 *
 * \sa ExpectationBasedPointSetRegistrationFilter
 * \ingroup FiniteDifferenceFunctions
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField, typename TPointSet>
class ExpectationBasedPointSetRegistrationFunction final
  : public AvantsPDEDeformableRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>
{
public:
  /** Standard class typedefs. */
  using Self = ExpectationBasedPointSetRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField, TPointSet>;
  using Superclass = PDEDeformableRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(ExpectationBasedPointSetRegistrationFunction);

  /** MovingImage image type. */
  using MovingImageType = typename Superclass::MovingImageType;
  using MovingImagePointer = typename Superclass::MovingImagePointer;

  /** FixedImage image type. */
  using FixedImageType = typename Superclass::FixedImageType;
  using FixedImagePointer = typename Superclass::FixedImagePointer;
  using ImagePointType = typename FixedImageType::PointType;
  using IndexType = typename FixedImageType::IndexType;
  using SizeType = typename FixedImageType::SizeType;
  using SpacingType = typename FixedImageType::SpacingType;

  /** Deformation field type. */
  using DisplacementFieldType = typename Superclass::DisplacementFieldType;
  using DisplacementFieldTypePointer = typename Superclass::DisplacementFieldTypePointer;
  using VectorType = typename DisplacementFieldType::PixelType;

  /** Inherit some enums from the superclass. */
  static constexpr unsigned int ImageDimension = Superclass::ImageDimension;
  static constexpr unsigned int MeasurementDimension = Superclass::ImageDimension;

  /** Inherit some enums from the superclass. */
  using PixelType = typename Superclass::PixelType;
  using RadiusType = typename Superclass::RadiusType;
  using NeighborhoodType = typename Superclass::NeighborhoodType;
  using FloatOffsetType = typename Superclass::FloatOffsetType;
  using TimeStepType = typename Superclass::TimeStepType;

  /** Covariant vector type. */
  using CovariantVectorType = CovariantVector<double, (Self::ImageDimension)>;

  /**  PointSet Types */
  using PointSetType = TPointSet;
  using PointSetPointer = typename PointSetType::Pointer;
  using PointType = typename PointSetType::PointType;
  using PointDataType = typename PointSetType::PixelType;
  using LabelSetType = std::vector<PointDataType>;
  //  typedef long PointDataType;
  using MeasurementVectorType = Vector<typename PointSetType::CoordRepType, MeasurementDimension>;
  using SampleType = typename Statistics::ListSample<MeasurementVectorType>;
  using TreeGeneratorType = typename Statistics::WeightedCentroidKdTreeGenerator<SampleType>;
  using NeighborhoodIdentifierType = typename TreeGeneratorType::KdTreeType::InstanceIdentifierVectorType;

  /** Bspline stuff */
  using BSplinePointSetType = PointSet<VectorType, (Self::ImageDimension)>;
  using BSplineFilterType = BSplineScatteredDataPointSetToImageFilter<BSplinePointSetType, DisplacementFieldType>;
  using BSplineWeightsType = typename BSplineFilterType::WeightsContainerType;
  using ControlPointLatticeType = typename BSplineFilterType::PointDataImageType;
  using ArrayType = typename BSplineFilterType::ArrayType;

  /** Other typedef */
  using RealType = float;
  using OutputType = float;
  using RandomizerType = Statistics::MersenneTwisterRandomVariateGenerator;
  using GaussianType = typename Statistics::GaussianMembershipFunction<VectorType>;

  /** Fixed image gradient calculator type. */
  using GradientCalculatorType = CentralDifferenceImageFunction<FixedImageType>;
  using GradientCalculatorPointer = typename GradientCalculatorType::Pointer;

  /** Moving image gradient calculator type. */

  using MovingImageGradientCalculatorType = CentralDifferenceImageFunction<MovingImageType>;
  using MovingImageGradientCalculatorPointer = typename MovingImageGradientCalculatorType::Pointer;

  /** This class uses a constant timestep of 1. */
  TimeStepType
  ComputeGlobalTimeStep(void * itkNotUsed(GlobalData)) const override
  {
    return m_TimeStep;
  }

  /** Return a pointer to a global data structure that is passed to
   * this object from the solver at each calculation.  */
  void *
  GetGlobalDataPointer() const override
  {
    auto * global = new GlobalDataStruct();

    global->m_SumOfSquaredDifference = 0.0;
    global->m_NumberOfPixelsProcessed = 0L;
    global->m_SumOfSquaredChange = 0;
    return global;
  }

  /** Release memory for global data structure. */
  void
  ReleaseGlobalDataPointer(void * gd) const override;

  void
  ExpectationLandmarkField(float weight, bool whichdirection);

  void
  FastExpectationLandmarkField(float weight, bool whichdirection, long whichlabel, bool dobspline);

  /** Set the object's state before each iteration. */
  void
  InitializeIteration() override;

  /** This method is called by a finite difference solver image filter at
   * each pixel that does not lie on a data set boundary */
  PixelType
  ComputeUpdate(const NeighborhoodType & it,
                void *                   globalData,
                const FloatOffsetType &  offset = FloatOffsetType(0.0)) override;

  PixelType
  ComputeUpdateInv(const NeighborhoodType & it,
                   void *                   globalData,
                   const FloatOffsetType &  offset = FloatOffsetType(0.0)) override;

  /** Get the metric value. The metric value is the mean square difference
   * in intensity between the fixed image and transforming moving image
   * computed over the the overlapping region between the two images. */
  virtual double
  GetMetric() const
  {
    return m_Metric;
  }

  /** Get the rms change in deformation field. */
  virtual double
  GetRMSChange() const
  {
    return m_RMSChange;
  }

  /** Set/Get the threshold below which the absolute difference of
   * intensity yields a match. When the intensities match between a
   * moving and fixed image pixel, the update vector (for that
   * iteration) will be the zero vector. Default is 0.001. */
  virtual void
  SetEuclideanDistanceThreshold(double);

  virtual double
  GetEuclideanDistanceThreshold() const;

  void
  SetFixedPointSetSigma(float f)
  {
    this->m_FixedPointSetSigma = f;
  }

  float
  GetFixedPointSetSigma()
  {
    return this->m_FixedPointSetSigma;
  }

  void
  SetMovingPointSetSigma(float f)
  {
    this->m_MovingPointSetSigma = f;
  }

  float
  GetMovingPointSetSigma()
  {
    return this->m_MovingPointSetSigma;
  }

  void
  SetKNeighborhood(unsigned int n)
  {
    this->m_KNeighborhood = n;
  }

  unsigned int
  GetKNeighborhood()
  {
    return this->m_KNeighborhood;
  }

  void
  SetUseSymmetricMatching(unsigned int b)
  {
    this->m_UseSymmetricMatching = b;
  }

protected:
  ExpectationBasedPointSetRegistrationFunction();
  ~ExpectationBasedPointSetRegistrationFunction() override = default;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /** FixedImage image neighborhood iterator type. */
  using FixedImageNeighborhoodIteratorType = ConstNeighborhoodIterator<FixedImageType>;

  /** A global data type for this class of equation. Used to store
   * information for computing the metric. */
  struct GlobalDataStruct
  {
    double        m_SumOfSquaredDifference;
    unsigned long m_NumberOfPixelsProcessed;
    double        m_SumOfSquaredChange;
  };

  void
  SetUpKDTrees(long whichlabel);

private:
  ExpectationBasedPointSetRegistrationFunction(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  /** Cache fixed image information. */
  SpacingType    m_FixedImageSpacing;
  ImagePointType m_FixedImageOrigin;
  double         m_Normalizer;

  /** Function to compute derivatives of the fixed image. */
  GradientCalculatorPointer m_FixedImageGradientCalculator;

  /** Function to compute derivatives of the moving image. */
  MovingImageGradientCalculatorPointer m_MovingImageGradientCalculator;
  bool                                 m_UseMovingImageGradient;

  /** The global timestep. */
  TimeStepType m_TimeStep;

  /** Threshold below which the denominator term is considered zero. */
  double m_DenominatorThreshold;

  /** Threshold below which two intensity value are assumed to match. */
  double m_EuclideanDistanceThreshold;

  /** The metric value is the mean square difference in intensity between
   * the fixed image and transforming moving image computed over the
   * the overlapping region between the two images. */
  mutable double        m_Metric;
  mutable double        m_SumOfSquaredDifference;
  mutable unsigned long m_NumberOfPixelsProcessed;
  mutable double        m_RMSChange;
  mutable double        m_SumOfSquaredChange;

  /** Mutex lock to protect modification to metric. */
  //  mutable SimpleFastMutexLock     m_MetricCalculationLock;

  DisplacementFieldTypePointer m_DerivativeFixedField;
  DisplacementFieldTypePointer m_DerivativeMovingField;

  float m_FixedPointSetSigma;
  float m_MovingPointSetSigma;

  float m_LandmarkEnergy;

  unsigned int m_KNeighborhood;
  unsigned int m_BucketSize;
  RealType     m_Sigma;

  typename TreeGeneratorType::Pointer m_FixedKdTreeGenerator;
  typename SampleType::Pointer        m_FixedSamplePoints;
  typename TreeGeneratorType::Pointer m_MovingKdTreeGenerator;
  typename SampleType::Pointer        m_MovingSamplePoints;
  typename RandomizerType::Pointer    m_Randomizer;
  bool                                m_Normalize;
  LabelSetType                        m_LabelSet;
  unsigned int                        m_UseSymmetricMatching;

  typename BSplinePointSetType::Pointer m_bpoints;
  typename BSplineWeightsType::Pointer  m_bweights;
  unsigned int                          m_bcount;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkExpectationBasedPointSetRegistrationFunction.hxx"
#endif

#endif
