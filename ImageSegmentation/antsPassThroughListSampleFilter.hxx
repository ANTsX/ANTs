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
#ifndef __antsPassThroughListSampleFilter_hxx
#define __antsPassThroughListSampleFilter_hxx

#include "antsPassThroughListSampleFilter.h"

namespace itk
{
namespace ants
{
namespace Statistics
{
template <class TListSample>
PassThroughListSampleFilter<TListSample>
::PassThroughListSampleFilter()
{
  this->AllocateOutput();
  this->GetOutput()->SetMeasurementVectorSize(
    this->GetInput()->GetMeasurementVectorSize() );
}

template <class TListSample>
PassThroughListSampleFilter<TListSample>
::~PassThroughListSampleFilter()
{
}

template <class TListSample>
void
PassThroughListSampleFilter<TListSample>
::GenerateData()
{
  /**
   * Simply pass the input to the output.
   */
  typename ListSampleType::ConstIterator It = this->GetInput()->Begin();
  while( It != this->GetInput()->End() )
    {
    this->GetOutput()->PushBack( It.GetMeasurementVector() );
    ++It;
    }
}

template <class TListSample>
void
PassThroughListSampleFilter<TListSample>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );
}
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#endif
