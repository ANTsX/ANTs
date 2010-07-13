/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: antsListSampleFunction.txx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:20:04 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
  http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt
  for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __antsListSampleFunction_txx
#define __antsListSampleFunction_txx

#include "antsListSampleFunction.h"

namespace itk
{
namespace ants
{
namespace Statistics
{
/**
 * Constructor
 */
template <class TInputListSample, class TOutput, class TCoordRep>
ListSampleFunction<TInputListSample, TOutput, TCoordRep>
::ListSampleFunction()
{
  this->m_ListSample = NULL;
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputListSample, class TOutput, class TCoordRep>
void
ListSampleFunction<TInputListSample, TOutput, TCoordRep>
::PrintSelf(
  std::ostream& os,
  Indent indent) const
{
  os << indent << "InputListSample: " << m_ListSample.GetPointer() << std::endl;
}

template <class TInputListSample, class TOutput, class TCoordRep>
void
ListSampleFunction<TInputListSample, TOutput, TCoordRep>
::SetWeights( WeightArrayType* array )
{
  this->m_Weights = *array;
  this->Modified();
}

template <class TInputListSample, class TOutput, class TCoordRep>
typename ListSampleFunction<TInputListSample, TOutput, TCoordRep>::WeightArrayType
* ListSampleFunction<TInputListSample, TOutput, TCoordRep>
::GetWeights()
  {
  return &this->m_Weights;
  }

/**
 * Initialize by setting the input point set
 */
template <class TInputListSample, class TOutput, class TCoordRep>
void
ListSampleFunction<TInputListSample, TOutput, TCoordRep>
::SetInputListSample(
  const InputListSampleType * ptr )
{
  // set the input image
  m_ListSample = ptr;
}
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#endif
