/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkANTSSimilarityMetric.h,v $
  Language:  C++
  Date:      $Date: 2008/11/15 23:46:06 $
  Version:   $Revision: 1.18 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkANTSSimilarityMetric_h
#define __itkANTSSimilarityMetric_h

#include "itkObject.h"
#include "itkObjectFactory.h"

#include "itkAvantsPDEDeformableRegistrationFunction.h"
#include "itkMacro.h"
#include "itkVector.h"
#include "itkANTSLabeledPointSet.h"

namespace itk
{
template <unsigned int TDimension = 3, class TReal = float>
class ANTSSimilarityMetric
  : public       Object
{
public:
  /** Standard class typedefs. */
  typedef ANTSSimilarityMetric     Self;
  typedef Object                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ANTSSimilarityMetric, Object );
  itkStaticConstMacro( Dimension, unsigned int, TDimension );

  typedef TReal RealType;
  typedef Image<RealType,
                itkGetStaticConstMacro( Dimension )>                   ImageType;
  typedef typename ImageType::Pointer ImagePointer;
  typedef Vector<RealType,
                 itkGetStaticConstMacro( Dimension )>                   VectorType;
  typedef Image<VectorType,
                itkGetStaticConstMacro( Dimension )>                   DisplacementFieldType;

  /** Point Types  for landmarks and labeled point-sets */
  typedef itk::ANTSLabeledPointSet<Dimension>        LabeledPointSetType;
  typedef typename LabeledPointSetType::Pointer      LabeledPointSetPointer;
  typedef typename LabeledPointSetType::PointSetType PointSetType;
  typedef typename PointSetType::Pointer             PointSetPointer;
  typedef typename PointSetType::PointType           PointType;
  typedef typename PointSetType::PixelType           PointDataType;
  typedef typename ImageType::PointType              ImagePointType;

  /** Typedefs for similarity metrics */
  typedef AvantsPDEDeformableRegistrationFunction
    <ImageType, ImageType, DisplacementFieldType>           MetricType;
  typedef typename MetricType::Pointer    MetricPointer;
  typedef typename MetricType::RadiusType RadiusType;

  itkSetMacro( FixedImage, ImagePointer );
  itkGetConstMacro( FixedImage, ImagePointer );
  itkSetMacro( MovingImage, ImagePointer );
  itkGetConstMacro( MovingImage, ImagePointer );
  itkSetMacro( FixedPointSet, PointSetPointer );
  itkGetConstMacro( FixedPointSet, PointSetPointer );
  itkSetMacro( MovingPointSet, PointSetPointer );
  itkGetConstMacro( MovingPointSet, PointSetPointer );
  itkSetMacro( WeightImage, ImagePointer );
  itkGetConstMacro( WeightImage, ImagePointer );
  itkSetClampMacro( WeightScalar, RealType, 0.0, NumericTraits<RealType>::max() );
  itkGetConstMacro( WeightScalar, RealType );

  itkSetObjectMacro( Metric, MetricType );
  itkGetModifiableObjectMacro( Metric, MetricType );

  itkSetMacro( MaximizeMetric, bool );
  itkGetConstMacro( MaximizeMetric, bool );
  itkBooleanMacro( MaximizeMetric );
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
