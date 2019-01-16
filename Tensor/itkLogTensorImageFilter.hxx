/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkLogTensorImageFilter_hxx
#define _itkLogTensorImageFilter_hxx

#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkOffset.h"
#include "itkProgressReporter.h"
#include "itkObjectFactory.h"
#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "itkLogTensorImageFilter.h"
#include "TensorFunctions.h"

namespace itk
{
template <typename TInputImage, typename TOutputImage>
LogTensorImageFilter<TInputImage, TOutputImage>
::LogTensorImageFilter()
= default;

template <typename TInputImage, typename TOutputImage>
void
LogTensorImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  InputImagePointer  input = this->GetInput();
  OutputImagePointer output = this->GetOutput();

  ImageRegionConstIterator<InputImageType> inputIt( input, input->GetLargestPossibleRegion() );

  output->SetRegions( input->GetLargestPossibleRegion() );
  output->Allocate();

  ImageRegionIteratorWithIndex<OutputImageType> outputIt( output, output->GetLargestPossibleRegion() );
  for( inputIt.GoToBegin(), outputIt.GoToBegin();
       !inputIt.IsAtEnd() && !outputIt.IsAtEnd();
       ++inputIt, ++outputIt )
    {
    InputPixelType result = TensorLog<InputPixelType>(inputIt.Value() );
    outputIt.Set( result );
    }
}

/**
 * Standard "PrintSelf" method
 */
template <typename TInputImage, typename TOutput>
void
LogTensorImageFilter<TInputImage, TOutput>
::PrintSelf(
  std::ostream& os,
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}
} // end namespace itk

#endif
