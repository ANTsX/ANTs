/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkFEMElement3DC0LinearTriangular_h
#define __itkFEMElement3DC0LinearTriangular_h

#include "itkFEMElementStd.h"

namespace itk
{
namespace fem
{
/**
 * \class Element3DC0LinearTriangular
 * \brief 3-noded, linear, C0 continuous finite element in 2D space.
 */
class Element3DC0LinearTriangular : public ElementStd<3, 3>
{
  typedef ElementStd<3, 3> TemplatedParentClass;
  FEM_ABSTRACT_CLASS(Element3DC0LinearTriangular, TemplatedParentClass)
public:
  // ////////////////////////////////////////////////////////////////////////
  /*
   * Methods related to numeric integration
   */

  enum
  {
    DefaultIntegrationOrder = 1
  };

  virtual void
  GetIntegrationPointAndWeight(unsigned int i, VectorType & pt, Float & w, unsigned int order) const;

  virtual unsigned int
  GetNumberOfIntegrationPoints(unsigned int order) const;

  // ////////////////////////////////////////////////////////////////////////
  /*
   * Methods related to the geometry of an element
   */

  virtual VectorType
  ShapeFunctions(const VectorType & pt) const;

  virtual void
  ShapeFunctionDerivatives(const VectorType & pt, MatrixType & shapeD) const;

  // FIXME: Write a proper implementation
  virtual bool
  GetLocalFromGlobalCoordinates(const VectorType & globalPt, VectorType & localPt) const;

  virtual Float
  JacobianDeterminant(const VectorType & pt, const MatrixType * pJ = 0) const;

  virtual void
  JacobianInverse(const VectorType & pt, MatrixType & invJ, const MatrixType * pJ = 0) const;

  /**
   * Draw the element on the specified device context
   */
#ifdef FEM_BUILD_VISUALIZATION
  void
  Draw(CDC * pDC, Solution::ConstPointer sol) const;

#endif

  /**
   * Constants for integration rules.
   */
  static const Float trigGaussRuleInfo[6][7][4];

  /**
   * Array that holds number of integration point for each order
   * of numerical integration.
   */
  static const unsigned int Nip[6];
};
} // namespace fem
} // namespace itk

#endif // #ifndef __itkFEMElement3DC0LinearTriangular_h
