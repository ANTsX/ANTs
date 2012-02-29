/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: antsListSampleToListSampleFilter.hxx,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
  http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt
  for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkListSampleToListSampleFilter_hxx
#define _itkListSampleToListSampleFilter_hxx

#include "antsListSampleToListSampleFilter.h"

namespace itk
{
namespace ants
{
namespace Statistics
{
/**
 *
 */
template <class TInputListSample, class TOutputListSample>
ListSampleToListSampleFilter<TInputListSample, TOutputListSample>
::ListSampleToListSampleFilter()
{
  // Modify superclass default values, can be overridden by subclasses
  this->SetNumberOfRequiredInputs( 1 );
  this->SetNumberOfRequiredOutputs( 1 );
}

template <class TInputListSample, class TOutputListSample>
void
ListSampleToListSampleFilter<TInputListSample, TOutputListSample>
::SetInputListSample( const TInputListSample *input )
{
//   this->m_InputListSample = const_cast<InputListSampleType *>( input );
  this->ProcessObject::SetNthInput( 0,
                                    reinterpret_cast<DataObject *>(
                                      const_cast<InputListSampleType *>( input ) ) );
}

template <class TInputListSample, class TOutputListSample>
void
ListSampleToListSampleFilter<TInputListSample, TOutputListSample>
::AllocateOutput()
{
  typename DataObject::Pointer obj =
    reinterpret_cast<DataObject *>(TOutputListSample::New().GetPointer() );

//   typename TOutputListSample::Pointer output
//     = reinterpret_cast<TOutputListSample*>(obj.GetPointer());
  this->ProcessObject::SetNumberOfRequiredOutputs( 1 );
  this->ProcessObject::SetNthOutput( 0, obj.GetPointer() );

//   this->m_OutputListSample = OutputListSampleType::New();
}

/**
 *
 */
template <class TInputListSample, class TOutputListSample>
typename ListSampleToListSampleFilter<TInputListSample, TOutputListSample>::InputListSampleType
* ListSampleToListSampleFilter<TInputListSample, TOutputListSample>
::GetInput()
  {
  return reinterpret_cast<TInputListSample *>(
    this->ProcessObject::GetInput( 0 ) );
  }

template <class TInputListSample, class TOutputListSample>
typename ListSampleToListSampleFilter<TInputListSample, TOutputListSample>::OutputListSampleType
* ListSampleToListSampleFilter<TInputListSample, TOutputListSample>
::GetOutput()
  {
  if( this->GetNumberOfOutputs() < 1 )
    {
    return 0;
    }

  // we assume that the first output is of the templated type
  return reinterpret_cast<TOutputListSample *>( this->ProcessObject::GetOutput( 0 ) );
//   return this->m_OutputListSample.GetPointer();
  }
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#endif
