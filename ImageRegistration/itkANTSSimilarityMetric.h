/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/ANTsX/ANTs/blob/main/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkANTSSimilarityMetric_h
#define __itkANTSSimilarityMetric_h

#include "ITKGetterSetterMacroShims.h"
#include "itkObject.h"
#include "itkObjectFactory.h"

#include "itkAvantsPDEDeformableRegistrationFunction.h"
#include "itkMacro.h"
#include "itkVector.h"
#include "itkANTSLabeledPointSet.h"

namespace itk
{
template <unsigned int TDimension = 3, typename TReal = float>
class ANTSSimilarityMetric final : public Object
{
public:
  /** Standard class typedefs. */
  typedef ANTSSimilarityMetric     Self;
  typedef Object                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(ANTSSimilarityMetric);
  static constexpr unsigned int Dimension = TDimension;

  typedef TReal                                                RealType;
  typedef Image<RealType, Self::Dimension>   ImageType;
  typedef typename ImageType::Pointer                          ImagePointer;
  typedef Vector<RealType, Self::Dimension>  VectorType;
  typedef Image<VectorType, Self::Dimension> DisplacementFieldType;

  /** Point Types  for landmarks and labeled point-sets */
  typedef itk::ANTSLabeledPointSet<Dimension>        LabeledPointSetType;
  typedef typename LabeledPointSetType::Pointer      LabeledPointSetPointer;
  typedef typename LabeledPointSetType::PointSetType PointSetType;
  typedef typename PointSetType::Pointer             PointSetPointer;
  typedef typename PointSetType::PointType           PointType;
  typedef typename PointSetType::PixelType           PointDataType;
  typedef typename ImageType::PointType              ImagePointType;

  /** Typedefs for similarity metrics */
  typedef AvantsPDEDeformableRegistrationFunction<ImageType, ImageType, DisplacementFieldType> MetricType;
  typedef typename MetricType::Pointer                                                         MetricPointer;
  typedef typename MetricType::RadiusType                                                      RadiusType;

  itkNonVirtualSetMacro(FixedImage, ImagePointer);
  itkNonVirtualGetConstMacro(FixedImage, ImagePointer);
  itkNonVirtualSetMacro(MovingImage, ImagePointer);
  itkNonVirtualGetConstMacro(MovingImage, ImagePointer);
  itkNonVirtualSetMacro(FixedPointSet, PointSetPointer);
  itkNonVirtualGetConstMacro(FixedPointSet, PointSetPointer);
  itkNonVirtualSetMacro(MovingPointSet, PointSetPointer);
  itkNonVirtualGetConstMacro(MovingPointSet, PointSetPointer);
  itkNonVirtualSetMacro(WeightImage, ImagePointer);
  itkNonVirtualGetConstMacro(WeightImage, ImagePointer);
  itkNonVirtualSetClampMacro(WeightScalar, RealType, 0.0, NumericTraits<RealType>::max());
  itkNonVirtualGetConstMacro(WeightScalar, RealType);

  itkNonVirtualSetObjectMacro(Metric, MetricType);
  itkNonVirtualGetModifiableObjectMacro(Metric, MetricType);

  itkNonVirtualSetMacro(MaximizeMetric, bool);
  itkNonVirtualGetConstMacro(MaximizeMetric, bool);
  itkNonVirtualBooleanMacro(MaximizeMetric);

private:
  MetricPointer m_Metric;
  bool          m_MaximizeMetric;

  ImagePointer m_FixedImage;
  ImagePointer m_MovingImage;

  PointSetPointer m_FixedPointSet;
  PointSetPointer m_MovingPointSet;

  ImagePointer m_WeightImage;
  RealType     m_WeightScalar;
};
} // end namespace itk

#endif
