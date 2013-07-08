/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: antsListSampleFunction.hxx,v $
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
#ifndef __antsListSampleFunction_hxx
#define __antsListSampleFunction_hxx

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
  this->m_ListSamples.clear();
  this->m_ListSampleWeights.clear();
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputListSample, class TOutput, class TCoordRep>
void
ListSampleFunction<TInputListSample, TOutput, TCoordRep>
::PrintSelf( std::ostream& os, Indent indent) const
{
  for( unsigned int d = 0; d < this->m_ListSamples.size(); d++ )
    {
    os << indent << "InputListSample: " << this->m_ListSamples[d] << std::endl;
    }
}

template <class TInputListSample, class TOutput, class TCoordRep>
void
ListSampleFunction<TInputListSample, TOutput, TCoordRep>
::SetListSampleWeights( const unsigned int idx, ListSampleWeightArrayType* array )
{
  if( idx >= this->m_ListSampleWeights.size() )
    {
    this->m_ListSampleWeights.resize( idx + 1 );
    this->m_ListSampleWeights[idx] = array;
    this->Modified();
    }
  if( this->m_ListSampleWeights[idx] != array )
    {
    this->m_ListSampleWeights[idx] = array;
    this->Modified();
    }
}

template <class TInputListSample, class TOutput, class TCoordRep>
typename ListSampleFunction<TInputListSample, TOutput, TCoordRep>::ListSampleWeightArrayType
* ListSampleFunction<TInputListSample, TOutput, TCoordRep>
::GetListSampleWeights( const unsigned int idx )
  {
  if( idx < this->m_ListSampleWeights.size() )
    {
    return this->m_ListSampleWeights[idx];
    }
  else
    {
    return NULL;
    }
  }

/**
 * Initialize by setting the input point set
 */
template <class TInputListSample, class TOutput, class TCoordRep>
void
ListSampleFunction<TInputListSample, TOutput, TCoordRep>
::SetIndexedInputListSample( const unsigned int idx, const InputListSampleType * ptr )
{
  if( idx >= this->m_ListSamples.size() )
    {
    this->m_ListSamples.resize( idx + 1 );
    this->m_ListSamples[idx] = ptr;
    this->Modified();
    }
  if( this->m_ListSamples[idx] != ptr )
    {
    this->m_ListSamples[idx] = ptr;
    this->Modified();
    }
}

template <class TInputListSample, class TOutput, class TCoordRep>
const
typename ListSampleFunction<TInputListSample, TOutput, TCoordRep>::InputListSampleType
* ListSampleFunction<TInputListSample, TOutput, TCoordRep>
::GetInputListSample( const unsigned int idx ) const
  {
  if( idx < this->m_ListSamples.size() )
    {
    return this->m_ListSamples[idx];
    }
  else
    {
    return NULL;
    }
  }
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#endif
