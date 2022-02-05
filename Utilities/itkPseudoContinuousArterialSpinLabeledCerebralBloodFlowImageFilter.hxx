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
#ifndef __itkPseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter_hxx
#define __itkPseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter_hxx

#include "itkProgressReporter.h"
#include "itkImageRegionIterator.h"

namespace itk
{
template <typename TInputImage, typename TReferenceImage, typename TOutputImage>
PseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter<TInputImage, TReferenceImage, TOutputImage>::
  PseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter()
{
  this->m_TI1 = 700;
  this->m_TI2 = 1700;
  this->m_T1blood = 1664;
  this->m_Lambda = 0.9;
  this->m_Alpha = 0.85; // 0.95;
  this->m_SliceDelay = 45;

  this->SetNumberOfRequiredInputs(2);
}

template <typename TInputImage, typename TReferenceImage, typename TOutputImage>
void
PseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter<TInputImage, TReferenceImage, TOutputImage>::PrintSelf(
  std::ostream & os,
  Indent         indent) const
{
  Superclass::PrintSelf(os, indent);
}

template <typename TInputImage, typename TReferenceImage, typename TOutputImage>
void
PseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter<TInputImage, TReferenceImage, TOutputImage>::
  SetDifferenceImage(const TInputImage * img)
{
  SetNthInput(0, const_cast<TInputImage *>(img));
}

template <typename TInputImage, typename TReferenceImage, typename TOutputImage>
typename TInputImage::ConstPointer
PseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter<TInputImage, TReferenceImage, TOutputImage>::
  GetDifferenceImage()
{
  return static_cast<const TInputImage *>(this->ProcessObject::GetInput(0));
}

template <typename TInputImage, typename TReferenceImage, typename TOutputImage>
void
PseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter<TInputImage, TReferenceImage, TOutputImage>::
  SetReferenceImage(const TReferenceImage * ref)
{
  SetNthInput(1, const_cast<TReferenceImage *>(ref));
}

template <typename TInputImage, typename TReferenceImage, typename TOutputImage>
typename TReferenceImage::ConstPointer
PseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter<TInputImage, TReferenceImage, TOutputImage>::
  GetReferenceImage()
{
  return static_cast<const TReferenceImage *>(this->ProcessObject::GetInput(1));
}

template <typename TInputImage, typename TReferenceImage, typename TOutputImage>
void
PseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter<TInputImage, TReferenceImage, TOutputImage>::
  ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId)
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
    // const float TI = ( this->m_TI2 - this->m_TI1) + this->m_SliceDelay * (outIt.GetIndex()[2] - 1);

    typename ReferenceImageType::IndexType idx;
    for (unsigned int i = 0; i < ReferenceImageDimension; i++)
    {
      idx[i] = inIt.GetIndex()[i];
    }

    // float ratio = inIt.Value() / this->GetReferenceImage()->GetPixel( idx );


    // 1.06 * M0W *  exp( 1 / 40.0 - 1 / 80.0) * TE;


    // RealType deltaM = itIt.Value();                  //  control - tagged if  control is odd

    // RealType T_1a = 1650; // 3T  from ASL_whyNhow.pdf
    // RealType T_1t = 1300; // 3T

    // from "Impact of equilibrium magnetization of blood on ASL quantification"
    // RealType tau  = 2100; // FIXME milliseconds from PMC3049525
    // RealType w    = 700;  // FIXME milliseconds from PMC3049525

    // Label width: Not in dicom, but sequence-specific -- magic parameter. Reference values: pCASL 1.5, CASL 1.6,
    // PASL 0.7.
    // from PMC3049525
    // const RealType scaling = 4 * this->Alpha * M_0 * T_1t
    //  * ( exp( -1.0 * ( tau + w ) / T_1a ) - exp( -1.0 * w / T_1t )  );

    // float cbf = this->Lambda * inIt.Value() * ( -1.0 )  / scaling;
    float cbf = 0;

    outIt.Set(cbf);
    ++outIt;
    ++inIt;
  }
}
} // end namespace itk

#endif
