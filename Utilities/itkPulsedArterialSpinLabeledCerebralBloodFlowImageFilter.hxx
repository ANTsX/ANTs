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

#include "itkPulsedArterialSpinLabeledCerebralBloodFlowImageFilter.h"
#include "itkProgressReporter.h"
#include "itkImageRegionIterator.h"

namespace itk
{
template <class TInputImage, class TOutputImage>
PulsedArterialSpinLabeledCerebralBloodFlowImageFilter<TInputImage, TOutputImage>
::PulsedArterialSpinLabeledCerebralBloodFlowImageFilter()
{
  this->m_TI1 = 700;
  this->m_TI2 = 1700;
  this->m_T1blood = 1664;
  this->m_Lambda = 0.9;
  this->m_Alpha = 0.95;
  this->m_SliceDelay = 45;
}

template <class TInputImage, class TOutputImage>
void
PulsedArterialSpinLabeledCerebralBloodFlowImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

template <class TInputImage, class TOutputImage>
void
PulsedArterialSpinLabeledCerebralBloodFlowImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                       ThreadIdType threadId)
{
  itkDebugMacro(<< "Actually executing");

  ProgressReporter progress( this, threadId, outputRegionForThread.GetNumberOfPixels() );

  OutputImageRegionType outputRegion = outputRegionForThread;

  InputImageRegionType inputRegion;
  this->CallCopyOutputRegionToInputRegion(inputRegion, outputRegionForThread);

  ImageRegionIterator<OutputImageType> outIt(
    this->GetOutput(), outputRegion);

  ImageRegionConstIterator<InputImageType> inIt(
    this->GetInput(), inputRegion );

  while( !outIt.IsAtEnd() )
    {
    // account for delay between acquisition of slices
    float TI = ( this->m_TI2 - this->m_TI1)
      + this->m_SliceDelay * (outIt.GetIndex()[2] - 1);

    // 540,000 is a unit conversion to give ml/100g/min
    float cbf = inIt.Value() * 5400000.0 * this->m_Lambda
      / ( 2.0 * this->m_Alpha * this->m_TI1 * exp( -TI / this->m_T1blood ) );

    outIt.Set( cbf );
    ++outIt;
    ++inIt;
    }
}
} // end namespace itk

#endif
