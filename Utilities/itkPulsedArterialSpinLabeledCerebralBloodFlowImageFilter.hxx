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
#ifndef __itkPulsedArterialSpinLabeledCerebralBloodFlowImageFilter_hxx
#define __itkPulsedArterialSpinLabeledCerebralBloodFlowImageFilter_hxx

#include "itkProgressReporter.h"
#include "itkImageRegionIterator.h"

namespace itk
{
template <typename TInputImage, typename TReferenceImage, typename TOutputImage>
PulsedArterialSpinLabeledCerebralBloodFlowImageFilter<TInputImage, TReferenceImage, TOutputImage>::
  PulsedArterialSpinLabeledCerebralBloodFlowImageFilter()
{
  this->m_TI1 = 700;
  this->m_TI2 = 1700;
  this->m_T1blood = 1664;
  this->m_Lambda = 0.9;
  this->m_Alpha = 0.95;
  this->m_SliceDelay = 45;

  this->SetNumberOfRequiredInputs(2);
}

template <typename TInputImage, typename TReferenceImage, typename TOutputImage>
void
PulsedArterialSpinLabeledCerebralBloodFlowImageFilter<TInputImage, TReferenceImage, TOutputImage>::PrintSelf(
  std::ostream & os,
  Indent         indent) const
{
  Superclass::PrintSelf(os, indent);
}

template <typename TInputImage, typename TReferenceImage, typename TOutputImage>
void
PulsedArterialSpinLabeledCerebralBloodFlowImageFilter<TInputImage, TReferenceImage, TOutputImage>::SetDifferenceImage(
  const TInputImage * img)
{
  this->SetNthInput(0, const_cast<TInputImage *>(img));
}

template <typename TInputImage, typename TReferenceImage, typename TOutputImage>
typename TInputImage::ConstPointer
PulsedArterialSpinLabeledCerebralBloodFlowImageFilter<TInputImage, TReferenceImage, TOutputImage>::GetDifferenceImage()
{
  return static_cast<const TInputImage *>(this->ProcessObject::GetInput(0));
}

template <typename TInputImage, typename TReferenceImage, typename TOutputImage>
void
PulsedArterialSpinLabeledCerebralBloodFlowImageFilter<TInputImage, TReferenceImage, TOutputImage>::SetReferenceImage(
  const TReferenceImage * ref)
{
  this->SetNthInput(1, const_cast<TReferenceImage *>(ref));
}

template <typename TInputImage, typename TReferenceImage, typename TOutputImage>
typename TReferenceImage::ConstPointer
PulsedArterialSpinLabeledCerebralBloodFlowImageFilter<TInputImage, TReferenceImage, TOutputImage>::GetReferenceImage()
{
  return static_cast<const TReferenceImage *>(this->ProcessObject::GetInput(1));
}

template <typename TInputImage, typename TReferenceImage, typename TOutputImage>
void
PulsedArterialSpinLabeledCerebralBloodFlowImageFilter<TInputImage, TReferenceImage, TOutputImage>::ThreadedGenerateData(
  const OutputImageRegionType & outputRegionForThread,
  ThreadIdType                  threadId)
{
  itkDebugMacro(<< "Actually executing");

  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  OutputImageRegionType outputRegion = outputRegionForThread;

  InputImageRegionType inputRegion;
  this->CallCopyOutputRegionToInputRegion(inputRegion, outputRegionForThread);

  ImageRegionIterator<OutputImageType> outIt(this->GetOutput(), outputRegion);

  ImageRegionConstIterator<InputImageType> inIt(this->GetInput(), inputRegion);

  while (!outIt.IsAtEnd())
  {
    // account for delay between acquisition of slices
    float TI = (this->m_TI2 - this->m_TI1) + this->m_SliceDelay * (outIt.GetIndex()[2] - 1);

    typename ReferenceImageType::IndexType idx;
    for (unsigned int i = 0; i < ReferenceImageDimension; i++)
    {
      idx[i] = inIt.GetIndex()[i];
    }

    float ratio = inIt.Value() / this->GetReferenceImage()->GetPixel(idx);

    // 540,000 is a unit conversion to give ml/100g/min
    float cbf = ratio * 5400000.0f * this->m_Lambda /
                (2.0f * this->m_Alpha * this->m_TI1 * static_cast<float>(exp(-TI / this->m_T1blood)));

    outIt.Set(cbf);
    ++outIt;
    ++inIt;
  }
}
} // end namespace itk

#endif
