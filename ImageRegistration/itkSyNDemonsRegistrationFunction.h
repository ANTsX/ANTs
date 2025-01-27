/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkSyNDemonsRegistrationFunction_h_
#define _itkSyNDemonsRegistrationFunction_h_

#include "itkPDEDeformableRegistrationFunction.h"
#include "itkPoint.h"
#include "itkCovariantVector.h"
#include "itkInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkCentralDifferenceImageFunction.h"

namespace itk
{
/**
 * \class SyNDemonsRegistrationFunction
 *
 * This class encapsulate the PDE which drives the demons registration
 * algorithm. It is used by SyNDemonsRegistrationFilter to compute the
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
 * \sa SyNDemonsRegistrationFilter
 * \ingroup FiniteDifferenceFunctions
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
class SyNDemonsRegistrationFunction final
  : public AvantsPDEDeformableRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>
{
public:
  /** Standard class typedefs. */
  typedef SyNDemonsRegistrationFunction                                                    Self;
  typedef PDEDeformableRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField> Superclass;
  typedef SmartPointer<Self>                                                               Pointer;
  typedef SmartPointer<const Self>                                                         ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(SyNDemonsRegistrationFunction);

  /** MovingImage image type. */
  typedef typename Superclass::MovingImageType    MovingImageType;
  typedef typename Superclass::MovingImagePointer MovingImagePointer;

  /** FixedImage image type. */
  typedef typename Superclass::FixedImageType    FixedImageType;
  typedef typename Superclass::FixedImagePointer FixedImagePointer;
  typedef typename FixedImageType::IndexType     IndexType;
  typedef typename FixedImageType::SizeType      SizeType;
  typedef typename FixedImageType::SpacingType   SpacingType;

  /** Deformation field type. */
  typedef typename Superclass::DisplacementFieldType        DisplacementFieldType;
  typedef typename Superclass::DisplacementFieldTypePointer DisplacementFieldTypePointer;

  /** Inherit some enums from the superclass. */
  static constexpr unsigned int ImageDimension = Superclass::ImageDimension;

  /** Inherit some enums from the superclass. */
  typedef typename Superclass::PixelType        PixelType;
  typedef typename Superclass::RadiusType       RadiusType;
  typedef typename Superclass::NeighborhoodType NeighborhoodType;
  typedef typename Superclass::FloatOffsetType  FloatOffsetType;
  typedef typename Superclass::TimeStepType     TimeStepType;

  /** Interpolator type. */
  typedef double                                                        CoordRepType;
  typedef InterpolateImageFunction<MovingImageType, CoordRepType>       InterpolatorType;
  typedef typename InterpolatorType::Pointer                            InterpolatorPointer;
  typedef typename InterpolatorType::PointType                          PointType;
  typedef LinearInterpolateImageFunction<MovingImageType, CoordRepType> DefaultInterpolatorType;

  /** Covariant vector type. */
  typedef CovariantVector<double, Self::ImageDimension> CovariantVectorType;

  /** Fixed image gradient calculator type. */
  typedef CentralDifferenceImageFunction<FixedImageType> GradientCalculatorType;
  typedef typename GradientCalculatorType::Pointer       GradientCalculatorPointer;

  /** Moving image gradient calculator type. */
  typedef CentralDifferenceImageFunction<MovingImageType, CoordRepType> MovingImageGradientCalculatorType;
  typedef typename MovingImageGradientCalculatorType::Pointer           MovingImageGradientCalculatorPointer;

  /** Set the moving image interpolator. */
  void
  SetMovingImageInterpolator(InterpolatorType * ptr)
  {
    m_MovingImageInterpolator = ptr;
  }

  /** Get the moving image interpolator. */
  InterpolatorType *
  GetMovingImageInterpolator(void)
  {
    return m_MovingImageInterpolator;
  }

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
    GlobalDataStruct * global = new GlobalDataStruct();

    global->m_SumOfSquaredDifference = 0.0;
    global->m_NumberOfPixelsProcessed = 0L;
    global->m_SumOfSquaredChange = 0;
    return global;
  }

  /** Release memory for global data structure. */
  void
  ReleaseGlobalDataPointer(void * GlobalData) const override;

  /** Set the object's state before each iteration. */
  void
  InitializeIteration() override;

  /** This method is called by a finite difference solver image filter at
   * each pixel that does not lie on a data set boundary */
  PixelType
  ComputeUpdate(const NeighborhoodType & neighborhood,
                void *                   globalData,
                const FloatOffsetType &  offset = FloatOffsetType(0.0)) override;

  PixelType
  ComputeUpdateInv(const NeighborhoodType & neighborhood,
                   void *                   globalData,
                   const FloatOffsetType &  offset = FloatOffsetType(0.0)) override;

  void
  SetUseSSD(bool b)
  {
    this->m_UseSSD = b;
  }

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

  /** Select if the fixed image or moving image gradient is used for
   * the computating the demon forces. The fixed image gradient is used
   * by default. */
  virtual void
  SetUseMovingImageGradient(bool flag)
  {
    m_UseMovingImageGradient = flag;
  }

  virtual bool
  GetUseMovingImageGradient() const
  {
    return m_UseMovingImageGradient;
  }

  /** Set/Get the threshold below which the absolute difference of
   * intensity yields a match. When the intensities match between a
   * moving and fixed image pixel, the update vector (for that
   * iteration) will be the zero vector. Default is 0.001. */
  virtual void
  SetIntensityDifferenceThreshold(double);

  virtual double
  GetIntensityDifferenceThreshold() const;

protected:
  SyNDemonsRegistrationFunction();
  virtual ~SyNDemonsRegistrationFunction() override {}

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /** FixedImage image neighborhood iterator type. */
  typedef ConstNeighborhoodIterator<FixedImageType> FixedImageNeighborhoodIteratorType;

  /** A global data type for this class of equation. Used to store
   * information for computing the metric. */
  struct GlobalDataStruct
  {
    double        m_SumOfSquaredDifference;
    unsigned long m_NumberOfPixelsProcessed;
    double        m_SumOfSquaredChange;
  };

private:
  SyNDemonsRegistrationFunction(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  /** Cache fixed image information. */
  SpacingType m_FixedImageSpacing;
  PointType   m_FixedImageOrigin;
  double      m_Normalizer;

  /** Function to compute derivatives of the fixed image. */
  GradientCalculatorPointer m_FixedImageGradientCalculator;

  /** Function to compute derivatives of the moving image. */
  MovingImageGradientCalculatorPointer m_MovingImageGradientCalculator;
  bool                                 m_UseMovingImageGradient;

  /** Function to interpolate the moving image. */
  InterpolatorPointer m_MovingImageInterpolator;

  /** The global timestep. */
  TimeStepType m_TimeStep;

  /** Threshold below which the denominator term is considered zero. */
  double m_DenominatorThreshold;

  /** Threshold below which two intensity value are assumed to match. */
  double m_IntensityDifferenceThreshold;

  /** The metric value is the mean square difference in intensity between
   * the fixed image and transforming moving image computed over the
   * the overlapping region between the two images. */
  mutable double        m_Metric;
  mutable double        m_SumOfSquaredDifference;
  mutable unsigned long m_NumberOfPixelsProcessed;
  mutable double        m_RMSChange;
  mutable double        m_SumOfSquaredChange;

  bool m_UseSSD;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkSyNDemonsRegistrationFunction.cxx"
#endif

#endif
