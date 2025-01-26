/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkJensenHavrdaCharvatTsallisLabeledPointSetMetric_h
#define _itkJensenHavrdaCharvatTsallisLabeledPointSetMetric_h

#include "itkPointSetToPointSetMetric.h"

#include "itkIdentityTransform.h"

namespace itk
{
/** \class JensenHavrdaCharvatTsallisLabeledPointSetMetric
 *
 *
 * This class is templated over the fixed point-set type, moving
 * point-set type.
 *
 */
template <typename TPointSet>
class JensenHavrdaCharvatTsallisLabeledPointSetMetric : public PointSetToPointSetMetric<TPointSet, TPointSet>
{
public:
  /** Standard class typedefs. */
  typedef JensenHavrdaCharvatTsallisLabeledPointSetMetric Self;
  typedef PointSetToPointSetMetric<TPointSet, TPointSet>  Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkOverrideGetNameOfClassMacro(JensenHavrdaCharvatTsallisLabeledPointSetMetric);

  static constexpr unsigned int PointDimension = TPointSet::PointDimension;

  /** Types transferred from the base class */
  typedef typename Superclass::TransformType           TransformType;
  typedef typename Superclass::TransformPointer        TransformPointer;
  typedef typename Superclass::TransformParametersType TransformParametersType;
  typedef typename Superclass::TransformJacobianType   TransformJacobianType;

  typedef typename Superclass::MeasureType    MeasureType;
  typedef typename Superclass::DerivativeType DerivativeType;

  typedef TPointSet                        PointSetType;
  typedef typename PointSetType::PointType PointType;
  typedef typename PointSetType::PixelType PixelType;

  /**
   * Other typedefs
   */
  typedef double                                      RealType;
  typedef IdentityTransform<RealType, PointDimension> DefaultTransformType;

  typedef std::vector<PixelType> LabelSetType;

  /**
   * Public function definitions
   */

  itkSetConstObjectMacro(FixedPointSet, PointSetType);
  itkGetConstObjectMacro(FixedPointSet, PointSetType);

  itkSetConstObjectMacro(MovingPointSet, PointSetType);
  itkGetConstObjectMacro(MovingPointSet, PointSetType);

  itkSetObjectMacro(Transform, TransformType);
  itkGetModifiableObjectMacro(Transform, TransformType);

  /**
   * The only transform type used is Identity
   */
  unsigned int
  GetNumberOfParameters(void) const
  {
    return m_Transform->GetNumberOfParameters();
  }

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void
  Initialize(void);

  /** Get the number of values */
  unsigned int
  GetNumberOfValues() const;

  /** Get the derivatives of the match measure. */
  void
  GetDerivative(const TransformParametersType & parameters, DerivativeType & Derivative) const;

  /**  Get the value for single valued optimizers. */
  MeasureType
  GetValue(const TransformParametersType & parameters) const;

  /**  Get value and derivatives for multiple valued optimizers. */
  void
  GetValueAndDerivative(const TransformParametersType & parameters,
                        MeasureType &                   Value,
                        DerivativeType &                Derivative) const;

  itkSetClampMacro(Alpha, RealType, 1.0, 2.0);
  itkGetConstMacro(Alpha, RealType);

  itkSetMacro(UseRegularizationTerm, bool);
  itkGetConstMacro(UseRegularizationTerm, bool);
  itkBooleanMacro(UseRegularizationTerm);

  itkSetMacro(UseWithRespectToTheMovingPointSet, bool);
  itkGetConstMacro(UseWithRespectToTheMovingPointSet, bool);
  itkBooleanMacro(UseWithRespectToTheMovingPointSet);

  itkSetMacro(MovingPointSetSigma, RealType);
  itkGetConstMacro(MovingPointSetSigma, RealType);

  itkSetMacro(MovingEvaluationKNeighborhood, unsigned int);
  itkGetConstMacro(MovingEvaluationKNeighborhood, unsigned int);

  itkSetMacro(FixedPointSetSigma, RealType);
  itkGetConstMacro(FixedPointSetSigma, RealType);

  itkSetMacro(FixedEvaluationKNeighborhood, unsigned int);
  itkGetConstMacro(FixedEvaluationKNeighborhood, unsigned int);

  itkSetMacro(UseInputAsSamples, bool);
  itkGetConstMacro(UseInputAsSamples, bool);
  itkBooleanMacro(UseInputAsSamples);

  /**
   * If this->m_UseInputAsSamples = true, the following
   * two variables are not used.
   */

  itkSetMacro(NumberOfMovingSamples, unsigned long);
  itkGetConstMacro(NumberOfMovingSamples, unsigned long);

  itkSetMacro(NumberOfFixedSamples, unsigned long);
  itkGetConstMacro(NumberOfFixedSamples, unsigned long);

  itkSetMacro(UseAnisotropicCovariances, bool);
  itkGetConstMacro(UseAnisotropicCovariances, bool);
  itkBooleanMacro(UseAnisotropicCovariances);

  /**
   * If this->m_UseAnisotropicCovariances = false, the
   * following four variables are not used.
   */

  itkSetMacro(FixedCovarianceKNeighborhood, unsigned int);
  itkGetConstMacro(FixedCovarianceKNeighborhood, unsigned int);

  itkSetMacro(FixedKernelSigma, RealType);
  itkGetConstMacro(FixedKernelSigma, RealType);

  itkSetMacro(MovingCovarianceKNeighborhood, unsigned int);
  itkGetConstMacro(MovingCovarianceKNeighborhood, unsigned int);

  itkSetMacro(MovingKernelSigma, RealType);
  itkGetConstMacro(MovingKernelSigma, RealType);

  void
  SetFixedLabelSet(LabelSetType labels)
  {
    typename LabelSetType::const_iterator iter;
    for (iter = labels.begin(); iter != labels.end(); ++iter)
    {
      this->m_FixedLabelSet.push_back(*iter);
    }
    this->Modified();
  }

  LabelSetType *
  GetFixedLabelSet()
  {
    return &this->m_FixedLabelSet;
  }

  void
  SetMovingLabelSet(LabelSetType labels)
  {
    typename LabelSetType::const_iterator iter;
    for (iter = labels.begin(); iter != labels.end(); ++iter)
    {
      this->m_MovingLabelSet.push_back(*iter);
    }
    this->Modified();
  }

  LabelSetType *
  GetMovingLabelSet()
  {
    return &this->m_MovingLabelSet;
  }

protected:
  JensenHavrdaCharvatTsallisLabeledPointSetMetric();
  ~JensenHavrdaCharvatTsallisLabeledPointSetMetric() {}

  void
  PrintSelf(std::ostream & os, Indent indent) const;

private:
  JensenHavrdaCharvatTsallisLabeledPointSetMetric(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  bool m_UseRegularizationTerm;
  bool m_UseInputAsSamples;
  bool m_UseAnisotropicCovariances;
  bool m_UseWithRespectToTheMovingPointSet;

  RealType      m_MovingPointSetSigma;
  RealType      m_MovingKernelSigma;
  unsigned int  m_MovingCovarianceKNeighborhood;
  unsigned int  m_MovingEvaluationKNeighborhood;
  unsigned long m_NumberOfMovingSamples;

  RealType      m_FixedPointSetSigma;
  RealType      m_FixedKernelSigma;
  unsigned int  m_FixedCovarianceKNeighborhood;
  unsigned int  m_FixedEvaluationKNeighborhood;
  unsigned long m_NumberOfFixedSamples;

  RealType m_Alpha;

  TransformPointer m_Transform;

  LabelSetType m_FixedLabelSet;
  LabelSetType m_MovingLabelSet;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkJensenHavrdaCharvatTsallisLabeledPointSetMetric.hxx"
#endif

#endif
