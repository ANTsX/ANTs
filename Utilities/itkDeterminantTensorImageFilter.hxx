/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkDeterminantTensorImageFilter_hxx
#define _itkDeterminantTensorImageFilter_hxx


#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkProgressReporter.h"

#include "vnl/vnl_det.h"

namespace itk
{

template <typename TInputImage, typename TRealType, typename TOutputImage>
DeterminantTensorImageFilter<TInputImage, TRealType, TOutputImage>::DeterminantTensorImageFilter()
{
  this->DynamicMultiThreadingOff();
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
void
DeterminantTensorImageFilter<TInputImage, TRealType, TOutputImage>::ThreadedGenerateData(
  const OutputImageRegionType & outputRegionForThread,
  ThreadIdType                  threadId)
{
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  ImageRegionIterator<OutputImageType>     ItD(this->GetOutput(), outputRegionForThread);
  ImageRegionConstIterator<InputImageType> ItM(this->GetInput(), outputRegionForThread);

  ItD.GoToBegin();
  ItM.GoToBegin();
  while (!ItD.IsAtEnd() && !ItM.IsAtEnd())
  {
    ItD.Set(static_cast<OutputPixelType>(vnl_det((ItM.Get()).GetVnlMatrix())));
    progress.CompletedPixel();
    ++ItD;
    ++ItM;
  }
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
void
DeterminantTensorImageFilter<TInputImage, TRealType, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

} // end namespace itk

#endif
