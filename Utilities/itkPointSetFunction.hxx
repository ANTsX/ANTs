/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPointSetFunction_hxx
#define __itkPointSetFunction_hxx

#include "itkPointSetFunction.h"

namespace itk
{
/**
 * Constructor
 */
template <class TInputPointSet, class TOutput, class TCoordRep>
PointSetFunction<TInputPointSet, TOutput, TCoordRep>
::PointSetFunction()
{
  m_PointSet = ITK_NULLPTR;
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputPointSet, class TOutput, class TCoordRep>
void
PointSetFunction<TInputPointSet, TOutput, TCoordRep>
::PrintSelf(
  std::ostream& os,
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "InputPointSet: " << m_PointSet.GetPointer() << std::endl;
}

/**
 * Initialize by setting the input point set
 */
template <class TInputPointSet, class TOutput, class TCoordRep>
void
PointSetFunction<TInputPointSet, TOutput, TCoordRep>
::SetInputPointSet(
  const InputPointSetType * ptr )
{
  // set the input image
  m_PointSet = ptr;
}
} // end namespace itk

#endif
