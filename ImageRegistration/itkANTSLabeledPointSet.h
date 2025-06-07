/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkANTSLabeledPointSet_h
#define __itkANTSLabeledPointSet_h

#include "itkObject.h"
#include "itkObjectFactory.h"

#include "itkMacro.h"
#include "itkVector.h"
#include "itkPointSet.h"

namespace itk
{
template <unsigned int TDimension = 3>
class ANTSLabeledPointSet : public Object
{
public:
  /** Standard class typedefs. */
  typedef ANTSLabeledPointSet      Self;
  typedef Object                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(ANTSLabeledPointSet);
  static constexpr unsigned int Dimension = TDimension;

  typedef float                                                RealType;
  typedef Image<RealType, Self::Dimension>   ImageType;
  typedef typename ImageType::Pointer                          ImagePointer;
  typedef Vector<RealType, Self::Dimension>  VectorType;
  typedef Image<VectorType, Self::Dimension> DisplacementFieldType;

  /** Point Types  for landmarks and labeled point-sets */
  typedef long                                          PointDataVectorType;
  typedef itk::PointSet<PointDataVectorType, Dimension> PointSetType;
  typedef typename PointSetType::Pointer                PointSetPointer;
  typedef typename PointSetType::PointType              PointType;
  typedef typename PointSetType::PixelType              PointDataType;
  typedef typename ImageType::PointType                 ImagePointType;

  itkSetMacro(PointSet, PointSetPointer);
  itkGetConstMacro(PointSet, PointSetPointer);

  PointType
  GetPoint(unsigned long ii)
  {
    PointType point;

    this->m_PointSet->GetPoint(ii, &point);
    return point;
  }

  PointDataType
  GetPointData(unsigned long ii)
  {
    PointDataType data;

    this->m_PointSet->GetPointData(ii, &data);
    return data;
  }

  void
  SetPoint(unsigned long ii, PointType point)
  {
    this->m_PointSet->SetPoint(ii, point);
  }

  void
  SetPointData(unsigned long ii, PointDataType label)
  {
    this->m_PointSet->SetPointData(ii, label);
  }

  void
  SetPointAndData(unsigned long ii, PointType point, PointDataType label)
  {
    this->m_PointSet->SetPoint(ii, point);
    this->m_PointSet->SetPointData(ii, label);
  }

private:
  PointSetPointer m_PointSet;
};
} // end namespace itk

#endif
