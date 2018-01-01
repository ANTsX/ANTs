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
template <class TFixedImage, class TMovingImage, class TDisplacementField, class TPointSet>
class ExpectationBasedPointSetRegistrationFunction :
  public         AvantsPDEDeformableRegistrationFunction<TFixedImage,
                                                         TMovingImage,
                                                         TDisplacementField>
{
public:
  /** Standard class typedefs. */
  typedef ExpectationBasedPointSetRegistrationFunction Self;
  typedef PDEDeformableRegistrationFunction<TFixedImage,
                                            TMovingImage, TDisplacementField
                                            >    Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( ExpectationBasedPointSetRegistrationFunction,
                PDEDeformableRegistrationFunction );

  /** MovingImage image type. */
  typedef typename Superclass::MovingImageType    MovingImageType;
  typedef typename Superclass::MovingImagePointer MovingImagePointer;

  /** FixedImage image type. */
  typedef typename Superclass::FixedImageType    FixedImageType;
  typedef typename Superclass::FixedImagePointer FixedImagePointer;
  typedef typename FixedImageType::PointType     ImagePointType;
  typedef typename FixedImageType::IndexType     IndexType;
  typedef typename FixedImageType::SizeType      SizeType;
  typedef typename FixedImageType::SpacingType   SpacingType;

  /** Deformation field type. */
  typedef typename Superclass::DisplacementFieldType DisplacementFieldType;
  typedef typename Superclass::DisplacementFieldTypePointer
    DisplacementFieldTypePointer;
  typedef typename DisplacementFieldType::PixelType VectorType;

  /** Inherit some enums from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned  int, Superclass::ImageDimension);
  itkStaticConstMacro(MeasurementDimension, unsigned  int, Superclass::ImageDimension);

  /** Inherit some enums from the superclass. */
  typedef typename Superclass::PixelType        PixelType;
  typedef typename Superclass::RadiusType       RadiusType;
  typedef typename Superclass::NeighborhoodType NeighborhoodType;
  typedef typename Superclass::FloatOffsetType  FloatOffsetType;
  typedef typename Superclass::TimeStepType     TimeStepType;

  /** Covariant vector type. */
  typedef CovariantVector<double, itkGetStaticConstMacro(ImageDimension)> CovariantVectorType;

  /**  PointSet Types */
  typedef  TPointSet                       PointSetType;
  typedef  typename PointSetType::Pointer  PointSetPointer;
  typedef typename PointSetType::PointType PointType;
  typedef typename PointSetType::PixelType PointDataType;
  typedef std::vector<PointDataType>       LabelSetType;
//  typedef long PointDataType;
  typedef Vector<typename PointSetType::CoordRepType, MeasurementDimension> MeasurementVectorType;
  typedef typename Statistics::ListSample<MeasurementVectorType>            SampleType;
  typedef typename
    Statistics::WeightedCentroidKdTreeGenerator<SampleType>   TreeGeneratorType;
  typedef typename TreeGeneratorType::KdTreeType::
    InstanceIdentifierVectorType                              NeighborhoodIdentifierType;

  /** Bspline stuff */
  typedef PointSet<VectorType,
                   itkGetStaticConstMacro( ImageDimension )>              BSplinePointSetType;
  typedef BSplineScatteredDataPointSetToImageFilter
    <BSplinePointSetType, DisplacementFieldType>            BSplineFilterType;
  typedef typename BSplineFilterType::WeightsContainerType BSplineWeightsType;
  typedef typename BSplineFilterType::PointDataImageType   ControlPointLatticeType;
  typedef typename BSplineFilterType::ArrayType            ArrayType;

  /** Other typedef */
  typedef  float RealType;
  typedef  float OutputType;
  typedef typename Statistics
    ::MersenneTwisterRandomVariateGenerator                  RandomizerType;
  typedef typename Statistics
    ::GaussianMembershipFunction<VectorType>          GaussianType;

  /** Fixed image gradient calculator type. */
  typedef CentralDifferenceImageFunction<FixedImageType> GradientCalculatorType;
  typedef typename GradientCalculatorType::Pointer       GradientCalculatorPointer;

  /** Moving image gradient calculator type. */

  typedef CentralDifferenceImageFunction<MovingImageType> MovingImageGradientCalculatorType;
  typedef typename MovingImageGradientCalculatorType::Pointer
    MovingImageGradientCalculatorPointer;

  /** This class uses a constant timestep of 1. */
  TimeStepType ComputeGlobalTimeStep(void * itkNotUsed(GlobalData) ) const ITK_OVERRIDE
  {
    return m_TimeStep;
  }

  /** Return a pointer to a global data structure that is passed to
   * this object from the solver at each calculation.  */
  void * GetGlobalDataPointer() const ITK_OVERRIDE
  {
    GlobalDataStruct *global = new GlobalDataStruct();

    global->m_SumOfSquaredDifference  = 0.0;
    global->m_NumberOfPixelsProcessed = 0L;
    global->m_SumOfSquaredChange      = 0;
    return global;
  }

  /** Release memory for global data structure. */
  void ReleaseGlobalDataPointer( void *GlobalData ) const ITK_OVERRIDE;

  void ExpectationLandmarkField(float weight, bool whichdirection);

  void FastExpectationLandmarkField(float weight, bool whichdirection, long whichlabel, bool dobsp);

  /** Set the object's state before each iteration. */
  void InitializeIteration() ITK_OVERRIDE;

  /** This method is called by a finite difference solver image filter at
   * each pixel that does not lie on a data set boundary */
  PixelType  ComputeUpdate(const NeighborhoodType & neighborhood, void *globalData, const FloatOffsetType & offset = FloatOffsetType(
                                       0.0) ) ITK_OVERRIDE;

  PixelType  ComputeUpdateInv(const NeighborhoodType & neighborhood, void *globalData, const FloatOffsetType & offset = FloatOffsetType(
                                          0.0) ) ITK_OVERRIDE;

  /** Get the metric value. The metric value is the mean square difference
   * in intensity between the fixed image and transforming moving image
   * computed over the the overlapping region between the two images. */
  virtual double GetMetric() const
  {
    return m_Metric;
  }

  /** Get the rms change in deformation field. */
  virtual double GetRMSChange() const
  {
    return m_RMSChange;
  }

  /** Set/Get the threshold below which the absolute difference of
   * intensity yields a match. When the intensities match between a
   * moving and fixed image pixel, the update vector (for that
   * iteration) will be the zero vector. Default is 0.001. */
  virtual void SetEuclideanDistanceThreshold(double);

  virtual double GetEuclideanDistanceThreshold() const;

  void SetFixedPointSetSigma( float f )
  {
    this->m_FixedPointSetSigma = f;
  }

  float GetFixedPointSetSigma()
  {
    return this->m_FixedPointSetSigma;
  }

  void SetMovingPointSetSigma( float f )
  {
    this->m_MovingPointSetSigma = f;
  }

  float GetMovingPointSetSigma()
  {
    return this->m_MovingPointSetSigma;
  }

  void SetKNeighborhood( unsigned int n )
  {
    this->m_KNeighborhood = n;
  }

  unsigned int GetKNeighborhood()
  {
    return this->m_KNeighborhood;
  }

  void SetUseSymmetricMatching(unsigned int b)
  {
    this->m_UseSymmetricMatching = b;
  }

protected:
  ExpectationBasedPointSetRegistrationFunction();
  virtual ~ExpectationBasedPointSetRegistrationFunction() ITK_OVERRIDE
  {
  }

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;

  /** FixedImage image neighborhood iterator type. */
  typedef ConstNeighborhoodIterator<FixedImageType> FixedImageNeighborhoodIteratorType;

  /** A global data type for this class of equation. Used to store
   * information for computing the metric. */
  struct GlobalDataStruct
    {
    double m_SumOfSquaredDifference;
    unsigned long m_NumberOfPixelsProcessed;
    double m_SumOfSquaredChange;
    };

  void SetUpKDTrees(long whichlabel);

private:
  ExpectationBasedPointSetRegistrationFunction(const Self &); // purposely not implemented
  void operator=(const Self &);                               // purposely not implemented

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

  typename TreeGeneratorType::Pointer                        m_FixedKdTreeGenerator;
  typename SampleType::Pointer                               m_FixedSamplePoints;
  typename TreeGeneratorType::Pointer                        m_MovingKdTreeGenerator;
  typename SampleType::Pointer                               m_MovingSamplePoints;
  typename RandomizerType::Pointer                           m_Randomizer;
  bool         m_Normalize;
  LabelSetType m_LabelSet;
  unsigned int m_UseSymmetricMatching;

  typename BSplinePointSetType::Pointer m_bpoints;
  typename BSplineWeightsType::Pointer m_bweights;
  unsigned int m_bcount;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkExpectationBasedPointSetRegistrationFunction.hxx"
#endif

#endif
