/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkFEMElement3DC0LinearTriangularMembrane_h
#define __itkFEMElement3DC0LinearTriangularMembrane_h

#include "itkFEMElement3DC0LinearTriangular.h"
#include "itkFEMElement3DMembrane.h"
#include "itkFEMElement3DMembrane1DOF.h"

namespace itk
{
namespace fem
{
/**
 * \class Element3DC0LinearTriangularMembrane
 * \brief 3-noded finite element class in 3D space for surface membrane problem.
 *
 * This element is combined from Element3DC0LinearTriangular and Element3DMembrane.
 */
class Element3DC0LinearTriangularMembrane : public Element3DMembrane<Element3DC0LinearTriangular>
{
  FEM_CLASS(Element3DC0LinearTriangularMembrane, Element3DMembrane<Element3DC0LinearTriangular>)
public:
  HANDLE_ELEMENT_LOADS();

  /**
   * Default constructor only clears the internal storage
   */
  Element3DC0LinearTriangularMembrane();

  /**
   * Construct an element by specifying pointers to
   * 3 points and a material.
   */
  Element3DC0LinearTriangularMembrane(NodeIDType n1_, NodeIDType n2_, NodeIDType n3_, Material::ConstPointer p_);

  //  virtual unsigned int GetNumberOfDegreesOfFreedomPerNode( void ) const
  //  { return 1; }

  //  virtual void GetStiffnessMatrix( MatrixType& Ke ) const;

  //  void Read( std::istream&, void* info );
}; // class Element3DC0LinearTriangularMembrane

FEM_CLASS_INIT(Element3DC0LinearTriangularMembrane)
} // namespace fem
} // namespace itk

#endif // #ifndef __itkFEMElement3DC0LinearTriangularMembrane_h
