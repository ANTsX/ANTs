/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
  https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt
  for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __antsJointHistogramParzenShapeAndOrientationListSampleFunction_h
#define __antsJointHistogramParzenShapeAndOrientationListSampleFunction_h

#include "antsListSampleFunction.h"

#include "itkImage.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
namespace itk
{
namespace ants
{
namespace Statistics
{
/** \class JointHistogramParzenShapeAndOrientationListSampleFunction.h
 * \brief point set filter.
 */

template <typename TListSample, typename TOutput = double, typename TCoordRep = double>
class JointHistogramParzenShapeAndOrientationListSampleFunction final
  : public ListSampleFunction<TListSample, TOutput, TCoordRep>
{
public:
  typedef JointHistogramParzenShapeAndOrientationListSampleFunction Self;
  typedef ListSampleFunction<TListSample, TOutput, TCoordRep>       Superclass;
  typedef SmartPointer<Self>                                        Pointer;
  typedef SmartPointer<const Self>                                  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(JointHistogramParzenShapeAndOrientationListSampleFunction);

  typedef typename Superclass::InputListSampleType        InputListSampleType;
  typedef typename Superclass::InputMeasurementVectorType InputMeasurementVectorType;
  typedef typename Superclass::InputMeasurementType       InputMeasurementType;
  /** List sample typedef support. */
  typedef TListSample ListSampleType;

  /** Other typedef */
  typedef TOutput RealType;
  typedef TOutput OutputType;

  typedef Image<RealType, 2>                          JointHistogramImageType;
  typedef typename JointHistogramImageType::PointType JointHistogramImagePointType;
  typedef typename JointHistogramImageType::Pointer   JointHistogramImagePointer;
  typedef Vector<RealType, 2>                         ThetaPsiType;

  typedef typename JointHistogramImageType::IndexType              IndexType;
  typedef typename IndexType::IndexValueType                       IndexValueType;
  typedef BSplineInterpolateImageFunction<JointHistogramImageType> InterpolatorType;
  typedef LinearInterpolateImageFunction<JointHistogramImageType>  LInterpolatorType;
  typedef typename InterpolatorType::Pointer                       InterpolatorPointer;

  /** Helper functions */

  void
  SetInputListSample(const InputListSampleType * ptr) override;

  itkSetMacro(ShapeSigma, RealType);
  itkGetConstMacro(ShapeSigma, RealType);
  itkSetMacro(OrientationSigma, RealType);
  itkGetConstMacro(OrientationSigma, RealType);

  itkSetMacro(NumberOfShapeJointHistogramBins, unsigned int);
  itkGetConstMacro(NumberOfShapeJointHistogramBins, unsigned int);
  itkSetMacro(NumberOfOrientationJointHistogramBins, unsigned int);
  itkGetConstMacro(NumberOfOrientationJointHistogramBins, unsigned int);

  TOutput
  Evaluate(const InputMeasurementVectorType &) const override;

protected:
  JointHistogramParzenShapeAndOrientationListSampleFunction();
  ~JointHistogramParzenShapeAndOrientationListSampleFunction() override;
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  GenerateData();

  void IncrementJointHistogramForShape(RealType, RealType);
  void
  IncrementJointHistogramForOrientation(RealType, RealType, RealType, unsigned int);

private:
  // purposely not implemented
  JointHistogramParzenShapeAndOrientationListSampleFunction(const Self &);
  void
  operator=(const Self &);

  unsigned int               m_NumberOfShapeJointHistogramBins;
  unsigned int               m_NumberOfOrientationJointHistogramBins;
  RealType                   m_ShapeSigma;
  RealType                   m_OrientationSigma;
  RealType                   m_MaximumEigenvalue1;
  RealType                   m_MaximumEigenvalue2;
  RealType                   m_MinimumEigenvalue1;
  RealType                   m_MinimumEigenvalue2;
  JointHistogramImagePointer m_JointHistogramImages[3];
  InterpolatorPointer        m_Interpolator;
  bool                       m_UseNearestNeighborIncrements;
};
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "antsJointHistogramParzenShapeAndOrientationListSampleFunction.hxx"
#endif

#endif
