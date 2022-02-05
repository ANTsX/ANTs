/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkSplitAlternatingTimeSeriesImageFilter_hxx
#define __itkSplitAlternatingTimeSeriesImageFilter_hxx

#include "itkProgressReporter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{
template <typename TInputImage, typename TOutputImage>
SplitAlternatingTimeSeriesImageFilter<TInputImage, TOutputImage>::SplitAlternatingTimeSeriesImageFilter()
{
  this->SetNumberOfRequiredOutputs(2);
  this->SetNthOutput(0, this->MakeOutput(0));
  this->SetNthOutput(1, this->MakeOutput(1));
}

template <typename TInputImage, typename TOutputImage>
void
SplitAlternatingTimeSeriesImageFilter<TInputImage, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

/**
 * \sa UnaryFunctorImageFilter::GenerateOutputInformation()
 */
template <typename TInputImage, typename TOutputImage>
void
SplitAlternatingTimeSeriesImageFilter<TInputImage, TOutputImage>::GenerateOutputInformation()
{
  // do not call the superclass' implementation of this method since
  // this filter allows the input the output to be of different dimensions

  // get pointers to the input and output
  typename Superclass::InputImageConstPointer inputPtr = this->GetInput();
  typename Superclass::OutputImagePointer     outputPtr0 = this->GetOutput(0);
  typename Superclass::OutputImagePointer     outputPtr1 = this->GetOutput(1);

  if (!outputPtr0 || !outputPtr1 || !inputPtr)
  {
    return;
  }

  // Set the output image largest possible region.  Use a RegionCopier
  // so that the input and output images can be different dimensions.
  OutputImageRegionType outputLargestPossibleRegion;
  this->CallCopyInputRegionToOutputRegion(outputLargestPossibleRegion, inputPtr->GetLargestPossibleRegion());

  // for the last dimension, assumed to be time
  unsigned int halfTime = inputPtr->GetLargestPossibleRegion().GetSize()[InputImageDimension - 1] / 2;
  outputLargestPossibleRegion.SetSize(InputImageDimension - 1, halfTime);

  outputPtr0->SetLargestPossibleRegion(outputLargestPossibleRegion);
  outputPtr1->SetLargestPossibleRegion(outputLargestPossibleRegion);
  outputPtr0->SetSpacing(inputPtr->GetSpacing());
  outputPtr1->SetSpacing(inputPtr->GetSpacing());

  outputPtr0->SetDirection(inputPtr->GetDirection());
  outputPtr1->SetDirection(inputPtr->GetDirection());

  typename InputImageType::PointType origin = inputPtr->GetOrigin();
  outputPtr0->SetOrigin(origin);
  origin[InputImageDimension - 1] += inputPtr->GetSpacing()[InputImageDimension - 1];
  outputPtr1->SetOrigin(origin);

  // Support VectorImages by setting number of components on output.
  const unsigned int numComponents = inputPtr->GetNumberOfComponentsPerPixel();
  if (numComponents != outputPtr0->GetNumberOfComponentsPerPixel())
  {
    outputPtr0->SetNumberOfComponentsPerPixel(numComponents);
  }
  if (numComponents != outputPtr1->GetNumberOfComponentsPerPixel())
  {
    outputPtr1->SetNumberOfComponentsPerPixel(numComponents);
  }
}

template <typename TInputImage, typename TOutputImage>
void
SplitAlternatingTimeSeriesImageFilter<TInputImage, TOutputImage>::ThreadedGenerateData(
  const OutputImageRegionType & outputRegionForThread,
  ThreadIdType                  threadId)
{
  itkDebugMacro(<< "Actually executing");

  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  OutputImageRegionType outputRegion = outputRegionForThread;

  ImageRegionIterator<OutputImageType> outIt0(this->GetOutput(0), outputRegion);
  ImageRegionIterator<OutputImageType> outIt1(this->GetOutput(1), outputRegion);

  while (!outIt0.IsAtEnd())
  {
    typename InputImageType::IndexType idx = outIt0.GetIndex();
    idx[InputImageDimension - 1] *= 2;
    outIt0.Set(this->GetInput()->GetPixel(idx));
    idx[InputImageDimension - 1] += 1;
    outIt1.Set(this->GetInput()->GetPixel(idx));

    ++outIt0;
    ++outIt1;
  }
}
} // end namespace itk

#endif
