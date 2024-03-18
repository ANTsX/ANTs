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
#ifndef itkNonLocalSuperresolutionImageFilter_hxx
#define itkNonLocalSuperresolutionImageFilter_hxx


#include "itkAbsoluteValueDifferenceImageFilter.h"
#include "itkArray.h"
#include "itkBoxMeanImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkCastImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkIterationReporter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkMath.h"
#include "itkNeighborhoodIterator.h"
#include "itkResampleImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkSubtractImageFilter.h"

#include "itkProgressReporter.h"

#include <algorithm>
#include <numeric>

namespace itk
{

template <typename TInputImage, typename TOutputImage>
NonLocalSuperresolutionImageFilter<TInputImage, TOutputImage>::NonLocalSuperresolutionImageFilter()
  : m_EpsilonThreshold(0.1)
  , m_PatchSimilaritySigma(1.0)
  , m_IntensityDifferenceSigma(1.0)
  , m_PerformInitialMeanCorrection(false)
  , m_CurrentIteration(0)
{
  this->SetNumberOfRequiredInputs(2);

  this->m_WeightSumImage = nullptr;

  this->m_InterpolatedLowResolutionInputImage = nullptr;

  // Interpolator --- default to linear
  typedef LinearInterpolateImageFunction<InputImageType, RealType> LinearInterpolatorType;
  this->m_Interpolator = LinearInterpolatorType::New();

  this->SetSimilarityMetric(NonLocalPatchBasedImageFilterEnums::SimilarityMetric::MEAN_SQUARES);
}

template <typename TInputImage, typename TOutputImage>
void
NonLocalSuperresolutionImageFilter<TInputImage, TOutputImage>::VerifyInputInformation() const
{}

template <typename TInputImage, typename TOutputImage>
void
NonLocalSuperresolutionImageFilter<TInputImage, TOutputImage>::AllocateOutputs()
{
  typename OutputImageType::Pointer outputImage = this->GetOutput();
  outputImage->CopyInformation(this->GetHighResolutionReferenceImage());
  outputImage->SetRegions(this->GetHighResolutionReferenceImage()->GetBufferedRegion());
  outputImage->AllocateInitialized();
}

template <typename TInputImage, typename TOutputImage>
void
NonLocalSuperresolutionImageFilter<TInputImage, TOutputImage>::GenerateOutputInformation()
{
  // Call the superclass' implementation of this method
  Superclass::GenerateOutputInformation();

  // Get pointers to the input and output
  OutputImageType * outputImage = this->GetOutput();
  if (!outputImage)
  {
    return;
  }

  const InputImageType * referenceImage = this->GetHighResolutionReferenceImage();

  // Set the size of the output region
  if (referenceImage)
  {
    outputImage->SetLargestPossibleRegion(referenceImage->GetLargestPossibleRegion());
    outputImage->SetSpacing(referenceImage->GetSpacing());
    outputImage->SetOrigin(referenceImage->GetOrigin());
    outputImage->SetDirection(referenceImage->GetDirection());
  }
}

template <typename TInputImage, typename TOutputImage>
void
NonLocalSuperresolutionImageFilter<TInputImage, TOutputImage>::GenerateInputRequestedRegion()
{
  // Call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  if (!this->GetInput())
  {
    return;
  }

  // Get pointers to the input
  InputImagePointer inputPtr = const_cast<TInputImage *>(this->GetInput());

  // Determining the actual input region is non-trivial, especially
  // when we cannot assume anything about the transform being used.
  // So we do the easy thing and request the entire input image.
  //
  inputPtr->SetRequestedRegionToLargestPossibleRegion();
}

template <typename TInputImage, typename TOutputImage>
void
NonLocalSuperresolutionImageFilter<TInputImage, TOutputImage>::GenerateData()
{
  if (this->m_ScaleLevels.size() == 0)
  {
    itkExceptionMacro("There are no scale levels.");
  }

  IterationReporter reporter(this, 0, 1);

  bool isConverged = false;
  this->m_CurrentIteration = 0;
  this->m_CurrentEpsilon = NumericTraits<RealType>::max();

  while (this->m_CurrentIteration < this->m_ScaleLevels.size() && isConverged == false)
  {
    reporter.CompletedStep();

    this->BeforeThreadedGenerateData();

    typename ImageSource<TOutputImage>::ThreadStruct str1;
    str1.Filter = this;

    this->GetMultiThreader()->SetSingleMethod(this->ThreaderCallback, &str1);

    this->GetMultiThreader()->SingleMethodExecute();

    this->AfterThreadedGenerateData();

    OutputImageType * outputImage = this->GetOutput();

    typedef DivideImageFilter<OutputImageType, RealImageType, InputImageType> DividerType;
    typename DividerType::Pointer                                             divider = DividerType::New();
    divider->SetInput1(outputImage);
    divider->SetInput2(this->m_WeightSumImage);
    divider->Update();

    InputImagePointer meanCorrectedImage = this->PerformMeanCorrection(divider->GetOutput());

    typedef AbsoluteValueDifferenceImageFilter<InputImageType, InputImageType, InputImageType> AbsoluterType;
    typename AbsoluterType::Pointer absoluter = AbsoluterType::New();
    absoluter->SetInput1(outputImage);
    absoluter->SetInput2(meanCorrectedImage);

    typedef StatisticsImageFilter<InputImageType> StatsFilterType;
    typename StatsFilterType::Pointer             stats = StatsFilterType::New();
    stats->SetInput(absoluter->GetOutput());
    stats->Update();

    this->m_CurrentEpsilon = stats->GetMean();

    if (this->m_CurrentEpsilon < this->m_EpsilonThreshold)
    {
      isConverged = true;
    }

    typedef CastImageFilter<InputImageType, OutputImageType> CasterType;
    typename CasterType::Pointer                             caster = CasterType::New();
    caster->SetInput(meanCorrectedImage);
    caster->Update();

    this->SetNthOutput(0, caster->GetOutput());

    this->m_CurrentIteration++;
  }
}

template <typename TInputImage, typename TOutputImage>
void
NonLocalSuperresolutionImageFilter<TInputImage, TOutputImage>::BeforeThreadedGenerateData()
{
  if (this->m_CurrentIteration == 0)
  {
    Superclass::BeforeThreadedGenerateData();

    this->m_Interpolator->SetInputImage(this->GetLowResolutionInputImage());

    typedef IdentityTransform<RealType, ImageDimension> IdentityTransformType;
    typename IdentityTransformType::Pointer             identityTransform = IdentityTransformType::New();
    identityTransform->SetIdentity();

    typedef ResampleImageFilter<InputImageType, InputImageType, RealType> ResamplerType;
    typename ResamplerType::Pointer                                       resampler = ResamplerType::New();
    resampler->SetInterpolator(this->m_Interpolator);
    resampler->SetInput(this->GetLowResolutionInputImage());
    resampler->SetTransform(identityTransform);
    resampler->SetOutputParametersFromImage(this->GetHighResolutionReferenceImage());
    resampler->Update();

    this->m_InterpolatedLowResolutionInputImage = resampler->GetOutput();

    // Initialize the weight sum image
    this->m_WeightSumImage = RealImageType::New();
    this->m_WeightSumImage->CopyInformation(this->GetHighResolutionReferenceImage());
    this->m_WeightSumImage->SetRegions(this->GetHighResolutionReferenceImage()->GetBufferedRegion());
    this->m_WeightSumImage->SetLargestPossibleRegion(
      this->GetHighResolutionReferenceImage()->GetLargestPossibleRegion());
    this->m_WeightSumImage->Allocate();
    this->m_WeightSumImage->FillBuffer(itk::NumericTraits<typename RealImageType::PixelType>::OneValue());

    Superclass::SetTargetImageRegion(this->GetHighResolutionReferenceImage()->GetBufferedRegion());

    if (this->m_PerformInitialMeanCorrection)
    {
      InputImageType *  referenceImage = const_cast<InputImageType *>(this->GetHighResolutionReferenceImage());
      InputImagePointer meanCorrectedImage = this->PerformMeanCorrection(referenceImage);
      this->SetHighResolutionReferenceImage(meanCorrectedImage);
    }

    this->AllocateOutputs();
  }
  else
  {
    typedef CastImageFilter<OutputImageType, InputImageType> CasterType;
    typename CasterType::Pointer                             caster = CasterType::New();
    caster->SetInput(this->GetOutput());
    caster->Update();

    this->m_InterpolatedLowResolutionInputImage = caster->GetOutput();

    this->m_WeightSumImage->FillBuffer(1.0);
  }
}

template <typename TInputImage, typename TOutputImage>
void
NonLocalSuperresolutionImageFilter<TInputImage, TOutputImage>::ThreadedGenerateData(const RegionType & region,
                                                                                    ThreadIdType       threadId)
{
  ProgressReporter progress(this, threadId, region.GetNumberOfPixels(), 100);

  const InputImageType * highResolutionInputImage = this->GetInput(1);

  OutputImageType * outputImage = this->GetOutput();

  NeighborhoodOffsetListType searchNeighborhoodOffsetList = this->GetNeighborhoodSearchOffsetList();
  SizeValueType              searchNeighborhoodSize = searchNeighborhoodOffsetList.size();

  // This is used for future extensions to include multiple high reference images

  InputImageList highResolutionInputImageList;
  highResolutionInputImageList.push_back(const_cast<InputImageType *>(highResolutionInputImage));

  ConstNeighborhoodIteratorType It(this->GetNeighborhoodPatchRadius(), highResolutionInputImage, region);

  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    progress.CompletedPixel();

    IndexType currentCenterIndex = It.GetIndex();

    InputImagePixelVectorType highResolutionPatch =
      this->VectorizeImageListPatch(highResolutionInputImageList, currentCenterIndex, true);

    for (SizeValueType i = 0; i < searchNeighborhoodSize; i++)
    {
      IndexType searchIndex = currentCenterIndex + searchNeighborhoodOffsetList[i];

      if (searchIndex == currentCenterIndex)
      {
        continue;
      }

      if (!outputImage->GetBufferedRegion().IsInside(searchIndex))
      {
        continue;
      }

      RealType intensityDifference = It.GetCenterPixel() - highResolutionInputImage->GetPixel(searchIndex);

      if (std::fabs(intensityDifference) >
          static_cast<RealType>(3.0) * this->m_IntensityDifferenceSigma *
            static_cast<RealType>(itk::Math::sqr(this->m_ScaleLevels[this->m_CurrentIteration])))
      {
        continue;
      }

      RealType patchSimilarity =
        this->ComputeNeighborhoodPatchSimilarity(highResolutionInputImageList, searchIndex, highResolutionPatch, true);

      RealType intensityWeight = itk::Math::sqr(
        intensityDifference / (this->m_IntensityDifferenceSigma * this->m_ScaleLevels[this->m_CurrentIteration]));

      RealType patchWeight = itk::Math::sqr(
        patchSimilarity / (this->m_PatchSimilaritySigma * this->m_ScaleLevels[this->m_CurrentIteration]));

      RealType weight = std::exp(-(intensityWeight + patchWeight));

      outputImage->SetPixel(currentCenterIndex,
                            outputImage->GetPixel(currentCenterIndex) +
                              weight * this->m_InterpolatedLowResolutionInputImage->GetPixel(searchIndex));

      this->m_WeightSumImage->SetPixel(currentCenterIndex,
                                       this->m_WeightSumImage->GetPixel(currentCenterIndex) + weight);
    }
  }
}

template <typename TInputImage, typename TOutputImage>
void
NonLocalSuperresolutionImageFilter<TInputImage, TOutputImage>::AfterThreadedGenerateData()
{}

template <typename TInputImage, typename TOutputImage>
typename NonLocalSuperresolutionImageFilter<TInputImage, TOutputImage>::InputImagePointer
NonLocalSuperresolutionImageFilter<TInputImage, TOutputImage>::PerformMeanCorrection(InputImageType * image)
{
  typedef BoxMeanImageFilter<TInputImage, TInputImage> BoxMeanFilterType;
  typename BoxMeanFilterType::Pointer                  boxMeanFilter = BoxMeanFilterType::New();
  boxMeanFilter->SetInput(image);

  typename InputImageType::SpacingType lowResolutionSpacing = this->GetLowResolutionInputImage()->GetSpacing();
  typename InputImageType::SpacingType highResolutionSpacing = this->GetHighResolutionReferenceImage()->GetSpacing();

  typename BoxMeanFilterType::RadiusType boxRadius;
  for (SizeValueType d = 0; d < ImageDimension; d++)
  {
    boxRadius[d] = static_cast<SizeValueType>(std::ceil(highResolutionSpacing[d] / lowResolutionSpacing[d])) - 1;
  }
  boxMeanFilter->SetRadius(boxRadius);
  boxMeanFilter->Update();

  typedef NearestNeighborInterpolateImageFunction<InputImageType, RealType> NearestNeighborInterpolatorType;
  typename NearestNeighborInterpolatorType::Pointer                         nearestNeighborInterpolator =
    NearestNeighborInterpolatorType::New();

  nearestNeighborInterpolator->SetInputImage(boxMeanFilter->GetOutput());

  typedef IdentityTransform<RealType, ImageDimension> IdentityTransformType;
  typename IdentityTransformType::Pointer             identityTransform = IdentityTransformType::New();
  identityTransform->SetIdentity();

  typedef ResampleImageFilter<InputImageType, InputImageType, RealType> ResamplerType;
  typename ResamplerType::Pointer                                       resampler = ResamplerType::New();
  resampler->SetInterpolator(nearestNeighborInterpolator);
  resampler->SetInput(boxMeanFilter->GetOutput());
  resampler->SetTransform(identityTransform);
  resampler->SetOutputParametersFromImage(this->GetLowResolutionInputImage());

  typedef SubtractImageFilter<InputImageType> SubtracterType;
  typename SubtracterType::Pointer            subtracter = SubtracterType::New();
  subtracter->SetInput1(resampler->GetOutput());
  subtracter->SetInput2(this->GetLowResolutionInputImage());
  subtracter->Update();

  nearestNeighborInterpolator->SetInputImage(subtracter->GetOutput());

  typedef ResampleImageFilter<InputImageType, InputImageType, RealType> ResamplerType2;
  typename ResamplerType2::Pointer                                      resampler2 = ResamplerType2::New();
  resampler2->SetInterpolator(nearestNeighborInterpolator);
  resampler2->SetInput(subtracter->GetOutput());
  resampler2->SetTransform(identityTransform);
  resampler2->SetOutputParametersFromImage(this->GetHighResolutionReferenceImage());

  typedef SubtractImageFilter<InputImageType> SubtracterType2;
  typename SubtracterType2::Pointer           subtracter2 = SubtracterType2::New();
  subtracter2->SetInput1(image);
  subtracter2->SetInput2(resampler2->GetOutput());
  subtracter2->Update();

  ImageRegionIteratorWithIndex<InputImageType> It(subtracter2->GetOutput(),
                                                  subtracter2->GetOutput()->GetRequestedRegion());
  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    if (It.Get() < NumericTraits<OutputPixelType>::ZeroValue())
    {
      It.Set(this->m_InterpolatedLowResolutionInputImage->GetPixel(It.GetIndex()));
    }
  }

  InputImagePointer meanCorrectedImage = subtracter2->GetOutput();
  meanCorrectedImage->DisconnectPipeline();

  return meanCorrectedImage;
}

template <typename TInputImage, typename TOutputImage>
void
NonLocalSuperresolutionImageFilter<TInputImage, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << "Interpolator: " << std::endl;
  this->m_Interpolator->Print(os, indent);

  os << indent << "Intensity difference sigma = " << this->m_IntensityDifferenceSigma << std::endl;
  os << indent << "Patch similarity sigma = " << this->m_PatchSimilaritySigma << std::endl;
}

} // end namespace itk

#endif
