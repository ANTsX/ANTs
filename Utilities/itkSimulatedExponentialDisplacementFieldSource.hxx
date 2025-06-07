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
#ifndef itkSimulatedExponentialDisplacementFieldSource_hxx
#define itkSimulatedExponentialDisplacementFieldSource_hxx


#include "itkAddImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkExponentialDisplacementFieldImageFilter.h"
#include "itkGaussianOperator.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkMultiplyImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVectorNeighborhoodOperatorImageFilter.h"


namespace itk
{
template <typename TOutputImage>
SimulatedExponentialDisplacementFieldSource<TOutputImage>::SimulatedExponentialDisplacementFieldSource()
  : m_NumberOfIntegrationSteps(0)
  , m_GradientStep(0.1)
  , m_DilationRadius(5)
  , m_SmoothingStandardDeviation(1.0)
{
  this->m_DisplacementNoiseStandardDeviation.Fill(1.0);
}

template <typename TOutputImage>
void
SimulatedExponentialDisplacementFieldSource<TOutputImage>::GenerateData()
{
  VectorType zeroVector(0.0);

  typename OutputImageType::Pointer randomField = OutputImageType::New();
  randomField->SetRegions(this->GetOutputSize());
  randomField->SetDirection(this->GetOutputDirection());
  randomField->SetOrigin(this->GetOutputOrigin());
  randomField->SetSpacing(this->GetOutputSpacing());
  randomField->AllocateInitialized();

  SizeType outputSize = this->GetOutputSize();

  for (SizeValueType i = 0; i < ImageDimension; i++)
  {

    // The grayscale morphology filter ignores the negative values so we
    // have to handle the positive and negative values separately and put
    // them together at the end.

    typename RealImageType::Pointer componentImagePositive = RealImageType::New();
    componentImagePositive->SetRegions(this->GetOutputSize());
    componentImagePositive->CopyInformation(randomField);
    componentImagePositive->AllocateInitialized();

    typename RealImageType::Pointer componentImageNegative = RealImageType::New();
    componentImageNegative->SetRegions(this->GetOutputSize());
    componentImageNegative->CopyInformation(randomField);
    componentImageNegative->AllocateInitialized();

    for (SizeValueType n = 0; n < this->GetNumberOfRandomPoints(); n++)
    {
      typename RealImageType::IndexType randomIndex;
      for (SizeValueType d = 0; d < ImageDimension; d++)
      {
        randomIndex[d] = this->GetRandomizer()->GetIntegerVariate(outputSize[d] - 1);
      }
      RealType voxelValue = static_cast<RealType>(this->GetRandomizer()->GetNormalVariate(
        NumericTraits<double>::ZeroValue(), std::pow(this->m_DisplacementNoiseStandardDeviation[i], 2)));

      if (voxelValue > 0)
      {
        componentImagePositive->SetPixel(randomIndex, voxelValue);
      }
      else
      {
        componentImageNegative->SetPixel(randomIndex, -voxelValue);
      }
    }
    typedef BinaryBallStructuringElement<RealType, ImageDimension> StructuringElementType;
    StructuringElementType                                         structuringElement;
    structuringElement.SetRadius(this->m_DilationRadius);
    structuringElement.CreateStructuringElement();

    using GrayscaleDilateImageFilterType =
      GrayscaleDilateImageFilter<RealImageType, RealImageType, StructuringElementType>;

    typename GrayscaleDilateImageFilterType::Pointer dilateFilterPositive = GrayscaleDilateImageFilterType::New();
    dilateFilterPositive->SetInput(componentImagePositive);
    dilateFilterPositive->SetKernel(structuringElement);

    typename GrayscaleDilateImageFilterType::Pointer dilateFilterNegative = GrayscaleDilateImageFilterType::New();
    dilateFilterNegative->SetInput(componentImageNegative);
    dilateFilterNegative->SetKernel(structuringElement);

    using MultiplierType = MultiplyImageFilter<RealImageType, RealImageType, RealImageType>;
    typename MultiplierType::Pointer multiplier = MultiplierType::New();
    multiplier->SetInput(dilateFilterNegative->GetOutput());
    multiplier->SetConstant(-1.0);

    using AdderType = AddImageFilter<RealImageType, RealImageType, RealImageType>;
    typename AdderType::Pointer adder = AdderType::New();
    adder->SetInput1(dilateFilterPositive->GetOutput());
    adder->SetInput2(multiplier->GetOutput());
    adder->Update();

    ImageRegionIterator<RealImageType>   ItC(adder->GetOutput(), adder->GetOutput()->GetRequestedRegion());
    ImageRegionIterator<OutputImageType> ItO(randomField, randomField->GetRequestedRegion());
    for (ItC.GoToBegin(), ItO.GoToBegin(); !ItC.IsAtEnd(); ++ItC, ++ItO)
    {
      VectorType vector = ItO.Get();
      vector[i] = ItC.Get();
      ItO.Set(vector);
    }
  }

  using DuplicatorType = ImageDuplicator<OutputImageType>;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(randomField);
  duplicator->Update();

  OutputImagePointer randomFieldSmooth = duplicator->GetOutput();

  RealType variance = std::pow(this->m_SmoothingStandardDeviation, 2);
  if (variance > 0.0)
  {
    using GaussianSmoothingOperatorType = GaussianOperator<OutputPixelComponentType, ImageDimension>;
    GaussianSmoothingOperatorType gaussianSmoothingOperator;

    using GaussianSmoothingSmootherType = VectorNeighborhoodOperatorImageFilter<OutputImageType, OutputImageType>;
    typename GaussianSmoothingSmootherType::Pointer smoother = GaussianSmoothingSmootherType::New();

    for (SizeValueType d = 0; d < ImageDimension; d++)
    {
      gaussianSmoothingOperator.SetDirection(d);
      gaussianSmoothingOperator.SetVariance(variance);
      gaussianSmoothingOperator.SetMaximumError(0.001);
      gaussianSmoothingOperator.SetMaximumKernelWidth(randomFieldSmooth->GetRequestedRegion().GetSize()[d]);
      gaussianSmoothingOperator.CreateDirectional();

      // todo: make sure we only smooth within the buffered region
      smoother->SetOperator(gaussianSmoothingOperator);
      smoother->SetInput(randomFieldSmooth);
      try
      {
        smoother->Update();
      }
      catch (const ExceptionObject & exc)
      {
        std::string msg("Caught exception: ");
        msg += exc.what();
        itkExceptionMacro(<< msg);
      }

      randomFieldSmooth = smoother->GetOutput();
      randomFieldSmooth->Update();
      randomFieldSmooth->DisconnectPipeline();
    }
  }

  using ExponentiatorType = ExponentialDisplacementFieldImageFilter<OutputImageType, OutputImageType>;
  typename ExponentiatorType::Pointer exponentiator = ExponentiatorType::New();
  exponentiator->SetInput(randomFieldSmooth);

  if (this->m_NumberOfIntegrationSteps == 0)
  {
    exponentiator->SetAutomaticNumberOfIterations(true);
  }
  else
  {
    exponentiator->SetAutomaticNumberOfIterations(false);
    exponentiator->SetMaximumNumberOfIterations(this->m_NumberOfIntegrationSteps);
  }
  exponentiator->SetComputeInverse(false);

  OutputImagePointer exponentiatedField = exponentiator->GetOutput();
  exponentiatedField->Update();
  exponentiatedField->DisconnectPipeline();

  typename DuplicatorType::Pointer duplicator2 = DuplicatorType::New();
  duplicator2->SetInputImage(exponentiatedField);
  duplicator2->Update();

  OutputImagePointer exponentiatedFieldSmooth = duplicator2->GetOutput();

  if (this->GetEnforceStationaryBoundary())
  {
    RealType weight1 = 1.0;
    if (variance < 0.5)
    {
      weight1 = 1.0 - 1.0 * (variance / 0.5);
    }
    RealType weight2 = 1.0 - weight1;

    const RegionType                          region = exponentiatedField->GetLargestPossibleRegion();
    const SizeType                            size = region.GetSize();
    const typename OutputImageType::IndexType startIndex = region.GetIndex();

    ImageRegionConstIteratorWithIndex<OutputImageType> ItF(exponentiatedField,
                                                           exponentiatedField->GetLargestPossibleRegion());
    ImageRegionIteratorWithIndex<OutputImageType>      ItS(exponentiatedFieldSmooth,
                                                      exponentiatedFieldSmooth->GetLargestPossibleRegion());
    for (ItF.GoToBegin(), ItS.GoToBegin(); !ItF.IsAtEnd(); ++ItF, ++ItS)
    {
      typename OutputImageType::IndexType index = ItF.GetIndex();
      bool                                isOnBoundary = false;
      for (unsigned int d = 0; d < ImageDimension; d++)
      {
        if (index[d] == startIndex[d] || index[d] == static_cast<IndexValueType>(size[d]) - startIndex[d] - 1)
        {
          isOnBoundary = true;
          break;
        }
      }
      if (isOnBoundary)
      {
        ItS.Set(zeroVector);
      }
      else
      {
        ItS.Set(ItS.Get() * weight1 + ItF.Get() * weight2);
      }
    }
  }

  this->SetNthOutput(0, exponentiatedFieldSmooth);
}

template <typename TOutputImage>
void
SimulatedExponentialDisplacementFieldSource<TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Number of compositions: " << this->m_NumberOfIntegrationSteps << std::endl;
  os << indent << "Gradient step: " << this->m_GradientStep << std::endl;
  os << indent << "Dilation radius: " << this->m_DilationRadius << std::endl;
  os << indent << "Displacement noise standard deviation: " << this->m_DisplacementNoiseStandardDeviation << std::endl;
  os << indent << "Smoothing standard deviation: " << this->m_SmoothingStandardDeviation << std::endl;
}


} // end namespace itk

#endif
