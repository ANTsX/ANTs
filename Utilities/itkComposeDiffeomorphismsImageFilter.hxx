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
#ifndef __itkComposeDiffeomorphismsImageFilter_hxx
#define __itkComposeDiffeomorphismsImageFilter_hxx


#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVectorLinearInterpolateImageFunction.h"

namespace itk
{
/*
 * ComposeDiffeomorphismsImageFilter class definitions
 */
template <typename InputImage, typename TOutputImage>
ComposeDiffeomorphismsImageFilter<InputImage, TOutputImage>::ComposeDiffeomorphismsImageFilter()
{
  this->SetNumberOfRequiredInputs(2);

  typedef VectorLinearInterpolateImageFunction<InputFieldType, RealType> DefaultInterpolatorType;
  typename DefaultInterpolatorType::Pointer                              interpolator = DefaultInterpolatorType::New();
  this->m_Interpolator = interpolator;
}

template <typename InputImage, typename TOutputImage>
ComposeDiffeomorphismsImageFilter<InputImage, TOutputImage>::~ComposeDiffeomorphismsImageFilter()
{}

template <typename InputImage, typename TOutputImage>
void
ComposeDiffeomorphismsImageFilter<InputImage, TOutputImage>::SetInterpolator(InterpolatorType * interpolator)
{
  itkDebugMacro("setting Interpolator to " << interpolator);
  if (this->m_Interpolator != interpolator)
  {
    this->m_Interpolator = interpolator;
    this->Modified();
    if (!this->GetDisplacementField())
    {
      this->m_Interpolator->SetInputImage(this->GetInput(0));
    }
  }
}

template <typename InputImage, typename TOutputImage>
void
ComposeDiffeomorphismsImageFilter<InputImage, TOutputImage>::BeforeThreadedGenerateData()
{
  VectorType zeroVector(0.0);

  this->AllocateOutputs();
  this->GetOutput()->FillBuffer(zeroVector);
}

template <typename InputImage, typename TOutputImage>
void
ComposeDiffeomorphismsImageFilter<InputImage, TOutputImage>::ThreadedGenerateData(const RegionType & region,
                                                                                  ThreadIdType itkNotUsed(threadId))
{
  VectorType zeroVector(0.0);

  typename OutputFieldType::Pointer output = this->GetOutput();

  ImageRegionConstIterator<InputFieldType>      ItW(this->GetWarpingField(), region);
  ImageRegionIteratorWithIndex<OutputFieldType> ItF(output, region);
  for (ItW.GoToBegin(), ItF.GoToBegin(); !ItW.IsAtEnd(); ++ItW, ++ItF)
  {
    PointType point1;
    output->TransformIndexToPhysicalPoint(ItF.GetIndex(), point1);

    typename InputFieldType::PixelType tmpWarp = ItW.Get();
    PointType                          point2 = point1;
    for (unsigned int d = 0; d < point2.PointDimension; d++)
    {
      point2[d] += tmpWarp[d];
    }

    typename InterpolatorType::OutputType displacement;
    if (this->m_Interpolator->IsInsideBuffer(point2))
    {
      displacement = this->m_Interpolator->Evaluate(point2);
      typename OutputFieldType::PixelType tmpOut;
      tmpOut = (point2 + displacement) - point1;
      ItF.Set(tmpOut);
    }
  }
}

template <typename InputImage, typename TOutputImage>
void
ComposeDiffeomorphismsImageFilter<InputImage, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}
} // end namespace itk

#endif
