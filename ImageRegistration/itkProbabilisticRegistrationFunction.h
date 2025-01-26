/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkProbabilisticRegistrationFunction_h_
#define _itkProbabilisticRegistrationFunction_h_

#include "antsAllocImage.h"
#include "itkAvantsPDEDeformableRegistrationFunction.h"
#include "itkPoint.h"
#include "itkCovariantVector.h"
#include "itkInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkCentralDifferenceImageFunction.h"
#include "itkGradientRecursiveGaussianImageFilter.h"

#include "itkAvantsMutualInformationRegistrationFunction.h"

namespace itk
{
/**
 * \class ProbabilisticRegistrationFunction
 *
 * This class encapsulate the PDE which drives the demons registration
 * algorithm. It is used by ProbabilisticRegistrationFilter to compute the
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
 * \sa ProbabilisticRegistrationFilter
 * \ingroup FiniteDifferenceFunctions
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
class ProbabilisticRegistrationFunction final
  : public AvantsPDEDeformableRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>
{
public:
  /** Standard class typedefs. */
  typedef ProbabilisticRegistrationFunction                                                      Self;
  typedef AvantsPDEDeformableRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField> Superclass;
  typedef SmartPointer<Self>                                                                     Pointer;
  typedef SmartPointer<const Self>                                                               ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(ProbabilisticRegistrationFunction);

  /** MovingImage image type. */
  typedef typename Superclass::MovingImageType    MovingImageType;
  typedef typename Superclass::MovingImagePointer MovingImagePointer;

  /** FixedImage image type. */
  typedef typename Superclass::MetricImageType          MetricImageType;
  typedef typename Superclass::MetricImageType::Pointer MetricImagePointer;
  typedef typename Superclass::FixedImageType           FixedImageType;
  typedef typename Superclass::FixedImagePointer        FixedImagePointer;
  typedef typename FixedImageType::IndexType            IndexType;
  typedef typename FixedImageType::SizeType             SizeType;

  /** Deformation field type. */
  typedef typename Superclass::DisplacementFieldType        DisplacementFieldType;
  typedef typename Superclass::DisplacementFieldTypePointer DisplacementFieldTypePointer;
  typedef typename TDisplacementField::PixelType            VectorType;

  typedef CovariantVector<float, Self::ImageDimension>           GradientPixelType;
  typedef Image<GradientPixelType, Self::ImageDimension>         GradientImageType;
  typedef SmartPointer<GradientImageType>                                          GradientImagePointer;
  typedef GradientRecursiveGaussianImageFilter<MetricImageType, GradientImageType> GradientImageFilterType;
  typedef typename GradientImageFilterType::Pointer                                GradientImageFilterPointer;
  typedef Image<float, Self::ImageDimension>                     BinaryImageType;
  typedef typename BinaryImageType::Pointer                                        BinaryImagePointer;

  /** Inherit some enums from the superclass. */
  static constexpr unsigned int ImageDimension = Superclass::ImageDimension;

  /** Inherit some enums from the superclass. */
  typedef typename Superclass::PixelType        PixelType;
  typedef typename Superclass::RadiusType       RadiusType;
  typedef typename Superclass::NeighborhoodType NeighborhoodType;
  //  typedef typename Superclass::NeighborhoodType    BoundaryNeighborhoodType;
  typedef typename Superclass::FloatOffsetType FloatOffsetType;
  typedef typename Superclass::TimeStepType    TimeStepType;

  /** Interpolator type. */
  typedef double                                                        CoordRepType;
  typedef InterpolateImageFunction<MovingImageType, CoordRepType>       InterpolatorType;
  typedef typename InterpolatorType::Pointer                            InterpolatorPointer;
  typedef typename InterpolatorType::PointType                          PointType;
  typedef LinearInterpolateImageFunction<MovingImageType, CoordRepType> DefaultInterpolatorType;

  /** Covariant vector type. */
  typedef CovariantVector<double, Self::ImageDimension> CovariantVectorType;

  /** Gradient calculator type. */
  typedef CentralDifferenceImageFunction<FixedImageType> GradientCalculatorType;
  typedef typename GradientCalculatorType::Pointer       GradientCalculatorPointer;

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

  typename TDisplacementField::PixelType
  ComputeMetricAtPairB(IndexType fixedindex, typename TDisplacementField::PixelType vec);
  typename TDisplacementField::PixelType
  ComputeMetricAtPairC(IndexType fixedindex, typename TDisplacementField::PixelType vec);

  /** This class uses a constant timestep of 1. */
  TimeStepType
  ComputeGlobalTimeStep(void * /* GlobalData */) const override
  {
    return m_TimeStep;
  }

  /** Return a pointer to a global data structure that is passed to
   * this object from the solver at each calculation.  */
  void *
  GetGlobalDataPointer() const override
  {
    GlobalDataStruct * global = new GlobalDataStruct();

    return global;
  }

  /** Release memory for global data structure. */
  void
  ReleaseGlobalDataPointer(void * GlobalData) const override
  {
    // HACK: The signature of this function should be reconsidered
    //      a delete of a re-interpret cast is a dangerous
    //      proposition.
    delete reinterpret_cast<GlobalDataStruct *>(GlobalData);
  }

  /** Set the object's state before each iteration. */
  void
  InitializeIteration() override;

  /** Set the object's state before each iteration. */
  void
  InitializeIterationOld();

  double
  ComputeCrossCorrelation()
  {
    if (finitediffimages[0])
    {
      double                                                totalcc = 0;
      unsigned long                                         ct = 0;
      typedef ImageRegionIteratorWithIndex<MetricImageType> ittype;
      ittype it(this->finitediffimages[0], this->finitediffimages[0]->GetLargestPossibleRegion().GetSize());
      for (it.GoToBegin(); !it.IsAtEnd(); ++it)
      {
        IndexType oindex = it.GetIndex();
        double    sfm = finitediffimages[2]->GetPixel(oindex);
        double    sff = finitediffimages[3]->GetPixel(oindex);
        double    smm = finitediffimages[4]->GetPixel(oindex);
        double    cc = 0;
        if (fabs(sff * smm) > 0)
        {
          cc += sfm * sfm / (sff * smm);
          ct++;
        }
        totalcc += cc;
      }
      this->m_Energy = totalcc / static_cast<double>(ct) * (-1.0);
      return this->m_Energy;
    }
    else
    {
      return 0;
    }
  }

  virtual VectorType
  OpticalFlowUpdate(const NeighborhoodType & neighborhood)
  {
    // Get fixed image related information
    IndexType index = neighborhood.GetIndex();

    typename TDisplacementField::PixelType vec = Superclass::m_DisplacementField->GetPixel(index);
    VectorType                             update;
    update.Fill(0.0);
    double              fixedValue;
    CovariantVectorType fixedGradient;
    double              fixedGradientSquaredMagnitude = 0;
    fixedValue = (double)Superclass::Superclass::m_FixedImage->GetPixel(index);
    fixedGradient = m_FixedImageGradientCalculator->EvaluateAtIndex(index);
    unsigned int j = 0;
    for (j = 0; j < ImageDimension; j++)
    {
      fixedGradientSquaredMagnitude += itk::Math::sqr(fixedGradient[j]);
    }
    double movingValue;

    PointType mappedPoint;
    for (j = 0; j < ImageDimension; j++)
    {
      mappedPoint[j] = double(index[j]) * m_FixedImageSpacing[j] + m_FixedImageOrigin[j];
      mappedPoint[j] += static_cast<double>(vec[j]);
    }
    if (m_MovingImageInterpolator->IsInsideBuffer(mappedPoint))
    {
      movingValue = m_MovingImageInterpolator->Evaluate(mappedPoint);
    }
    else
    {
      for (j = 0; j < ImageDimension; j++)
      {
        update[j] = 0.0;
      }
      return update;
    }
    double speedValue = fixedValue - movingValue;
    if (std::fabs(speedValue) < static_cast<double>(this->m_RobustnessParameter))
    {
      speedValue = 0;
    }
    double denominator = itk::Math::sqr(speedValue) / static_cast<double>(m_Normalizer) + fixedGradientSquaredMagnitude;
    double DenominatorThreshold = 1e-9;
    double IntensityDifferenceThreshold = 0.001;
    if (itk::Math::abs(speedValue) < IntensityDifferenceThreshold || denominator < DenominatorThreshold)
    {
      for (j = 0; j < ImageDimension; j++)
      {
        update[j] = 0.0;
      }
      return update;
    }
    for (j = 0; j < ImageDimension; j++)
    {
      update[j] = speedValue * fixedGradient[j] / denominator;
    }
    return update;
  }

  void
  SetFixedImageMask(MetricImageType * img)
  {
    m_FixedImageMask = img;
  }

  VectorType
  ComputeUpdate(const NeighborhoodType & neighborhood,
                void * /* globalData */,
                const FloatOffsetType & /* offset */ = FloatOffsetType(0.0)) override
  {
    VectorType update;

    update.Fill(0.0);
    IndexType        oindex = neighborhood.GetIndex();
    FixedImageType * img = const_cast<FixedImageType *>(Superclass::m_FixedImage.GetPointer());
    if (!img)
    {
      return update;
    }
    update = this->ComputeMetricAtPairB(oindex, update);

    return update;
  }

  VectorType
  ComputeUpdateInv(const NeighborhoodType & neighborhood,
                   void * /* globalData */,
                   const FloatOffsetType & /* offset */ = FloatOffsetType(0.0)) override
  {
    VectorType update;

    update.Fill(0.0);
    IndexType        oindex = neighborhood.GetIndex();
    FixedImageType * img = const_cast<FixedImageType *>(Superclass::m_FixedImage.GetPointer());
    if (!img)
    {
      return update;
    }
    update = this->ComputeMetricAtPairC(oindex, update);

    return update;
  }

  void
  SetFullyRobust(bool b)
  {
    m_FullyRobust = b;
  }

  void
  GetProbabilities();

  double localProbabilistic;
  float  m_TEMP;

protected:
  ProbabilisticRegistrationFunction();
  ~ProbabilisticRegistrationFunction() override = default;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /** FixedImage image neighborhood iterator type. */
  typedef ConstNeighborhoodIterator<FixedImageType> FixedImageNeighborhoodIteratorType;

  /** A global data type for this class of equation. Used to store
   * iterators for the fixed image. */
  struct GlobalDataStruct
  {
    FixedImageNeighborhoodIteratorType m_FixedImageIterator;
  };

  MetricImagePointer
  MakeImage()
  {
    FixedImageType * img = const_cast<FixedImageType *>(Superclass::m_FixedImage.GetPointer());

    this->m_MetricImage = AllocImage<MetricImageType>(img->GetLargestPossibleRegion(), 0);
    this->m_MetricImage->SetSpacing(img->GetSpacing());
    this->m_MetricImage->SetOrigin(img->GetOrigin());

    bool makebinimg = false;
    if (makebinimg)
    {
      m_Iteration = 0;
      binaryimage = AllocImage<BinaryImageType>(img->GetLargestPossibleRegion(), 1);
      binaryimage->SetSpacing(img->GetSpacing());
      binaryimage->SetOrigin(img->GetOrigin());
    }
    return this->m_MetricImage;
  }

private:
  ProbabilisticRegistrationFunction(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  /** Cache fixed image information. */
  typename TFixedImage::SpacingType m_FixedImageSpacing;
  typename TFixedImage::PointType   m_FixedImageOrigin;

  /** Function to compute derivatives of the fixed image. */
  GradientCalculatorPointer m_FixedImageGradientCalculator;

  GradientCalculatorPointer m_MovingImageGradientCalculator;

  /** Function to interpolate the moving image. */
  InterpolatorPointer m_MovingImageInterpolator;

  /** The global timestep. */
  TimeStepType m_TimeStep;

  /** Threshold below which the denominator term is considered zero. */
  double m_DenominatorThreshold;

  /** Threshold below which two intensity value are assumed to match. */
  double m_IntensityDifferenceThreshold;

  mutable double m_MetricTotal;
  mutable float  m_MinMag;
  mutable float  m_MaxMag;
  mutable float  m_AvgMag;
  mutable float  m_Thresh;

  GradientImagePointer m_MetricGradientImage;

  MetricImagePointer finitediffimages[5];
  BinaryImagePointer binaryimage;

  MetricImagePointer m_FixedImageMask;
  MetricImagePointer m_MovingImageMask;

  typedef itk::AvantsMutualInformationRegistrationFunction<FixedImageType, MovingImageType, DisplacementFieldType>
                                MetricType2;
  typename MetricType2::Pointer m_MIFunct;
  unsigned int                  m_NumberOfHistogramBins;
  bool                          m_FullyRobust;
  unsigned int                  m_Iteration;
  float                         m_Normalizer;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkProbabilisticRegistrationFunction.cxx"
#endif

#endif
