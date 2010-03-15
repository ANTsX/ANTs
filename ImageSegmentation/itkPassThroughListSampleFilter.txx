/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPassThroughListSampleFilter.txx,v $
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPassThroughListSampleFilter_txx
#define __itkPassThroughListSampleFilter_txx

#include "itkPassThroughListSampleFilter.h"

namespace itk
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
  const unsigned int measurementVectorSize =
    this->GetOutput()->GetMeasurementVectorSize();

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
} // end of namespace itk

#endif
