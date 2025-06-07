/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkJensenTsallisBSplineRegistrationFunction_h_
#define _itkJensenTsallisBSplineRegistrationFunction_h_

#include "itkAvantsPDEDeformableRegistrationFunction.h"
#include "itkJensenHavrdaCharvatTsallisLabeledPointSetMetric.h"

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkPointSet.h"

namespace itk
{
/**
 * \class JensenTsallisBSplineRegistrationFunction
 * \ingroup FiniteDifferenceFunctions
 */
template <typename TFixedImage,
          typename TFixedPointSet,
          typename TMovingImage,
          typename TMovingPointSet,
          typename TDisplacementField>
class JensenTsallisBSplineRegistrationFunction
  : public AvantsPDEDeformableRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>
{
public:
  /** Standard class typedefs. */
  typedef JensenTsallisBSplineRegistrationFunction                                               Self;
  typedef AvantsPDEDeformableRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField> Superclass;
  typedef SmartPointer<Self>                                                                     Pointer;
  typedef SmartPointer<const Self>                                                               ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkOverrideGetNameOfClassMacro(JensenTsallisBSplineRegistrationFunction);

  /**
   * Inherit some enums from the superclass.
   */
  static constexpr unsigned int ImageDimension = Superclass::ImageDimension;
  static constexpr unsigned int PointDimension = TFixedPointSet::PointDimension;

  typedef typename Superclass::NeighborhoodType NeighborhoodType;
  typedef typename Superclass::FloatOffsetType  FloatOffsetType;
  typedef typename Superclass::TimeStepType     TimeStepType;

  typedef JensenHavrdaCharvatTsallisLabeledPointSetMetric<TMovingPointSet> PointSetMetricType;
  typedef typename PointSetMetricType::RealType                            RealType;

  /**
   * Moving typedefs
   */
  typedef TMovingImage                         MovingImageType;
  typedef TMovingPointSet                      MovingPointSetType;
  typedef typename MovingImageType::Pointer    MovingImagePointer;
  typedef typename MovingPointSetType::Pointer MovingPointSetPointer;

  /**
   * Fixed typedefs
   */
  typedef TFixedImage                         FixedImageType;
  typedef TFixedPointSet                      FixedPointSetType;
  typedef typename FixedImageType::Pointer    FixedImagePointer;
  typedef typename FixedPointSetType::Pointer FixedPointSetPointer;

  /**
   * Deformation field typedefs
   */
  typedef TDisplacementField                        DisplacementFieldType;
  typedef typename DisplacementFieldType::Pointer   DisplacementFieldPointer;
  typedef typename DisplacementFieldType::PixelType VectorType;

  /**
   * BSpline typedefs
   */
  /** Typedefs for B-spline filter */
  typedef PointSet<VectorType, Self::ImageDimension>                          BSplinePointSetType;
  typedef BSplineScatteredDataPointSetToImageFilter<BSplinePointSetType, DisplacementFieldType> BSplineFilterType;
  typedef typename BSplineFilterType::WeightsContainerType                                      BSplineWeightsType;
  typedef typename BSplineFilterType::PointDataImageType                                        ControlPointLatticeType;
  typedef typename BSplineFilterType::ArrayType                                                 ArrayType;

  void
  SetAlpha(RealType r)
  {
    this->m_Alpha = r;
  }

  RealType
  GetAlpha()
  {
    return this->m_Alpha;
  }

  void
  SetUseRegularizationTerm(bool b)
  {
    this->m_UseRegularizationTerm = b;
  }

  RealType
  GetUseRegularizationTerm()
  {
    return this->m_UseRegularizationTerm;
  }

  void
  SetUseInputAsSamples(bool b)
  {
    this->m_UseInputAsSamples = b;
  }

  RealType
  GetUseInputAsSamples()
  {
    return this->m_UseInputAsSamples;
  }

  void
  SetUseAnisotropicCovariances(bool b)
  {
    this->m_UseAnisotropicCovariances = b;
  }

  RealType
  GetUseAnisotropicCovariances()
  {
    return this->m_UseAnisotropicCovariances;
  }

  void
  SetFixedPointSetSigma(RealType r)
  {
    this->m_FixedPointSetSigma = r;
  }

  RealType
  GetFixedPointSetSigma()
  {
    return this->m_FixedPointSetSigma;
  }

  void
  SetFixedKernelSigma(RealType r)
  {
    this->m_FixedKernelSigma = r;
  }

  RealType
  GetFixedKernelSigma()
  {
    return this->m_FixedKernelSigma;
  }

  void
  SetFixedEvaluationKNeighborhood(unsigned int r)
  {
    this->m_FixedEvaluationKNeighborhood = r;
  }

  unsigned int
  GetFixedEvaluationKNeighborhood()
  {
    return this->m_FixedEvaluationKNeighborhood;
  }

  void
  SetFixedCovarianceKNeighborhood(unsigned int r)
  {
    this->m_FixedCovarianceKNeighborhood = r;
  }

  unsigned int
  GetFixedCovarianceKNeighborhood()
  {
    return this->m_FixedCovarianceKNeighborhood;
  }

  void
  SetNumberOfFixedSamples(unsigned long r)
  {
    this->m_NumberOfFixedSamples = r;
  }

  unsigned long
  GetNumberOfFixedSamples()
  {
    return this->m_NumberOfFixedSamples;
  }

  void
  SetMovingPointSetSigma(RealType r)
  {
    this->m_MovingPointSetSigma = r;
  }

  RealType
  GetMovingPointSetSigma()
  {
    return this->m_MovingPointSetSigma;
  }

  void
  SetMovingKernelSigma(RealType r)
  {
    this->m_MovingKernelSigma = r;
  }

  RealType
  GetMovingKernelSigma()
  {
    return this->m_MovingKernelSigma;
  }

  void
  SetMovingEvaluationKNeighborhood(unsigned int r)
  {
    this->m_MovingEvaluationKNeighborhood = r;
  }

  unsigned int
  GetMovingEvaluationKNeighborhood()
  {
    return this->m_MovingEvaluationKNeighborhood;
  }

  void
  SetMovingCovarianceKNeighborhood(unsigned int r)
  {
    this->m_MovingCovarianceKNeighborhood = r;
  }

  unsigned int
  GetMovingCovarianceKNeighborhood()
  {
    return this->m_MovingCovarianceKNeighborhood;
  }

  void
  SetNumberOfMovingSamples(unsigned long r)
  {
    this->m_NumberOfMovingSamples = r;
  }

  unsigned long
  GetNumberOfMovingSamples()
  {
    return this->m_NumberOfMovingSamples;
  }

  void
  SetMeshResolution(ArrayType a)
  {
    this->m_MeshResolution = a;
  }

  ArrayType
  GetMeshResolution()
  {
    return this->m_MeshResolution;
  }

  void
  SetNumberOfLevels(unsigned int r)
  {
    this->m_NumberOfLevels = r;
  }

  unsigned int
  GetNumberOfLevels()
  {
    return this->m_NumberOfLevels;
  }

  void
  SetSplineOrder(unsigned int r)
  {
    this->m_SplineOrder = r;
  }

  unsigned int
  GetSplineOrder()
  {
    return this->m_SplineOrder;
  }

  /** This class uses a constant timestep of 1. */
  virtual TimeStepType
  ComputeGlobalTimeStep(void * itkNotUsed(GlobalData)) const
  {
    return m_TimeStep;
  }

  /** Return a pointer to a global data structure that is passed to
   * this object from the solver at each calculation.  */
  virtual void *
  GetGlobalDataPointer() const
  {
    GlobalDataStruct * global = new GlobalDataStruct();

    return global;
  }

  /** Release memory for global data structure. */
  virtual void
  ReleaseGlobalDataPointer(void * GlobalData) const
  {
    delete (GlobalDataStruct *)GlobalData;
  }

  /** Set the object's state before each iteration. */
  virtual void
  InitializeIteration();

  virtual VectorType
  ComputeUpdate(const NeighborhoodType & neighborhood,
                void *                   globalData,
                const FloatOffsetType &  offset = FloatOffsetType(0.0));

  virtual VectorType
  ComputeUpdateInv(const NeighborhoodType & neighborhood,
                   void *                   globalData,
                   const FloatOffsetType &  offset = FloatOffsetType(0.0));

protected:
  JensenTsallisBSplineRegistrationFunction();
  ~JensenTsallisBSplineRegistrationFunction() {}

  /** FixedImage image neighborhood iterator type. */
  typedef ConstNeighborhoodIterator<FixedImageType> FixedImageNeighborhoodIteratorType;

  /** A global data type for this class of equation. Used to store
   * iterators for the fixed image. */
  struct GlobalDataStruct
  {
    FixedImageNeighborhoodIteratorType m_FixedImageIterator;
  };

  void
  PrintSelf(std::ostream & os, Indent indent) const;

private:
  //  typename ControlPointLatticeType::Pointer                m_FixedControlPointLattice;
  //  typename ControlPointLatticeType::Pointer                m_MovingControlPointLattice;

  /** The global timestep. */
  TimeStepType m_TimeStep;

  /**
   * point-set related variables
   */
  bool m_UseRegularizationTerm;
  bool m_UseInputAsSamples;
  bool m_UseAnisotropicCovariances;

  typename DisplacementFieldType::Pointer m_DerivativeFixedField;
  typename DisplacementFieldType::Pointer m_DerivativeMovingField;

  RealType      m_FixedPointSetSigma;
  RealType      m_FixedKernelSigma;
  unsigned int  m_FixedEvaluationKNeighborhood;
  unsigned int  m_FixedCovarianceKNeighborhood;
  unsigned long m_NumberOfFixedSamples;

  RealType      m_MovingPointSetSigma;
  RealType      m_MovingKernelSigma;
  unsigned int  m_MovingEvaluationKNeighborhood;
  unsigned int  m_MovingCovarianceKNeighborhood;
  unsigned long m_NumberOfMovingSamples;

  RealType m_Alpha;

  /**
   * Bspline related variables
   */
  unsigned int m_SplineOrder;
  unsigned int m_NumberOfLevels;
  ArrayType    m_MeshResolution;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkJensenTsallisBSplineRegistrationFunction.hxx"
#endif

#endif
