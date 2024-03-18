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
#ifndef itkWeightedVotingFusionImageFilter_hxx
#define itkWeightedVotingFusionImageFilter_hxx


#include "itkImageRegionIteratorWithIndex.h"
#include "itkProgressReporter.h"

#include <algorithm>
#include <numeric>

#include <vnl/algo/vnl_cholesky.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_inverse.h>

namespace itk
{

template <typename TInputImage, typename TOutputImage>
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>::WeightedVotingFusionImageFilter()
  : m_IsWeightedAveragingComplete(false)
  , m_NumberOfAtlases(0)
  , m_NumberOfAtlasSegmentations(0)
  , m_NumberOfAtlasModalities(0)
  , m_Alpha(0.1)
  , m_Beta(2.0)
  , m_RetainLabelPosteriorProbabilityImages(false)
  , m_RetainAtlasVotingWeightImages(false)
  , m_ConstrainSolutionToNonnegativeWeights(false)
{
  this->m_MaskImage = nullptr;

  this->m_CountImage = nullptr;

  this->m_NeighborhoodSearchRadiusImage = nullptr;

  this->SetSimilarityMetric(itk::NonLocalPatchBasedImageFilterEnums::SimilarityMetric::PEARSON_CORRELATION);
}

template <typename TInputImage, typename TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>::UpdateInputs()
{
  // Set all the inputs

  this->SetNumberOfIndexedInputs(this->m_NumberOfAtlases * this->m_NumberOfAtlasModalities +
                                 this->m_NumberOfAtlasSegmentations + this->m_TargetImage.size() +
                                 this->m_LabelExclusionImages.size());

  SizeValueType nthInput = 0;

  for (SizeValueType i = 0; i < this->m_TargetImage.size(); i++)
  {
    this->SetNthInput(nthInput++, this->m_TargetImage[i]);
  }

  for (SizeValueType i = 0; i < this->m_NumberOfAtlases; i++)
  {
    for (SizeValueType j = 0; j < this->m_NumberOfAtlasModalities; j++)
    {
      this->SetNthInput(nthInput++, this->m_AtlasImages[i][j]);
    }
  }

  for (SizeValueType i = 0; i < this->m_NumberOfAtlasSegmentations; i++)
  {
    this->SetNthInput(nthInput++, this->m_AtlasSegmentations[i]);
  }

  typename LabelExclusionMap::const_iterator it;
  for (it = m_LabelExclusionImages.begin(); it != m_LabelExclusionImages.end(); ++it)
  {
    this->SetNthInput(nthInput++, it->second);
  }

  if (this->m_MaskImage.IsNotNull())
  {
    this->SetNthInput(nthInput++, this->m_MaskImage);
  }

  this->Modified();
}

template <typename TInputImage, typename TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

  // Get the output requested region
  RegionType outRegion = this->GetOutput()->GetRequestedRegion();

  // Pad this region by the search window and patch size
  if (this->m_NeighborhoodSearchRadiusImage.IsNull())
  {
    outRegion.PadByRadius(this->GetNeighborhoodSearchRadius());
  }
  else
  {
    NeighborhoodRadiusType maxNeighborhoodSearchRadius;
    maxNeighborhoodSearchRadius.Fill(0);

    ImageRegionConstIterator<RadiusImageType> ItR(this->m_NeighborhoodSearchRadiusImage,
                                                  this->m_NeighborhoodSearchRadiusImage->GetRequestedRegion());
    for (ItR.GoToBegin(); !ItR.IsAtEnd(); ++ItR)
    {
      RadiusValueType localSearchRadius = ItR.Get();
      if (localSearchRadius > maxNeighborhoodSearchRadius[0])
      {
        maxNeighborhoodSearchRadius.Fill(localSearchRadius);
      }
    }
    outRegion.PadByRadius(maxNeighborhoodSearchRadius);
  }
  outRegion.PadByRadius(this->GetNeighborhoodPatchRadius());

  // Iterate over all the inputs to this filter

  for (SizeValueType i = 0; i < this->m_TargetImage.size(); i++)
  {
    InputImageType * input = this->m_TargetImage[i];
    if (i == 0)
    {
      this->SetTargetImageRegion(input->GetRequestedRegion());
    }
    RegionType region = outRegion;
    region.Crop(input->GetLargestPossibleRegion());
    input->SetRequestedRegion(region);
  }

  for (SizeValueType i = 0; i < this->m_NumberOfAtlases; i++)
  {
    for (SizeValueType j = 0; j < this->m_NumberOfAtlasModalities; j++)
    {
      InputImageType * input = this->m_AtlasImages[i][j];
      RegionType       region = outRegion;
      region.Crop(input->GetLargestPossibleRegion());
      input->SetRequestedRegion(region);
    }
  }

  for (SizeValueType i = 0; i < this->m_NumberOfAtlasSegmentations; i++)
  {
    LabelImageType * input = this->m_AtlasSegmentations[i];
    RegionType       region = outRegion;
    region.Crop(input->GetLargestPossibleRegion());
    input->SetRequestedRegion(region);
  }

  typename LabelExclusionMap::const_iterator it;
  for (it = m_LabelExclusionImages.begin(); it != m_LabelExclusionImages.end(); ++it)
  {
    LabelImageType * input = it->second;
    RegionType       region = outRegion;
    region.Crop(input->GetLargestPossibleRegion());
    input->SetRequestedRegion(region);
  }

  if (this->m_MaskImage.IsNotNull())
  {
    MaskImageType * input = this->m_MaskImage;
    RegionType      region = outRegion;
    region.Crop(input->GetLargestPossibleRegion());
    input->SetRequestedRegion(region);
  }
}

template <typename TInputImage, typename TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>::GenerateData()
{
  this->BeforeThreadedGenerateData();

  /**
   * Multithread processing for the weighted averaging
   */
  typename ImageSource<TOutputImage>::ThreadStruct str1;
  str1.Filter = this;

  //  this->GetMultiThreader()->SetGlobalDefaultNumberOfThreads( this->GetNumberOfThreads() );
  this->GetMultiThreader()->SetSingleMethod(this->ThreaderCallback, &str1);

  this->GetMultiThreader()->SingleMethodExecute();

  this->m_IsWeightedAveragingComplete = true;

  /**
   * Multithread processing for the image(s) reconstruction
   */

  typename ImageSource<TOutputImage>::ThreadStruct str2;
  str2.Filter = this;

  //  this->GetMultiThreader()->SetGlobalDefaultNumberOfThreads( this->GetNumberOfThreads() );
  this->GetMultiThreader()->SetSingleMethod(this->ThreaderCallback, &str2);

  this->GetMultiThreader()->SingleMethodExecute();

  this->AfterThreadedGenerateData();
}

template <typename TInputImage, typename TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>::BeforeThreadedGenerateData()
{
  Superclass::BeforeThreadedGenerateData();

  if (this->m_NumberOfAtlasSegmentations != this->m_NumberOfAtlases)
  {
    // Set the number of atlas segmentations to 0 since we're just going to
    // doing joint intensity fusion
    this->m_NumberOfAtlasSegmentations = 0;
  }

  // Check to see if the number of target images is equal to 1 or equal to the number
  // of atlas modalities
  if (this->m_TargetImage.size() != 1 && this->m_TargetImage.size() != this->m_NumberOfAtlasModalities)
  {
    itkExceptionMacro("The number of target images must be 1 or must be the number of atlas modalities.");
  }

  // Find all the unique labels in the atlas segmentations
  this->m_LabelSet.clear();
  for (unsigned int i = 0; i < this->m_NumberOfAtlasSegmentations; i++)
  {
    ImageRegionConstIteratorWithIndex<LabelImageType> It(this->m_AtlasSegmentations[i],
                                                         this->m_AtlasSegmentations[i]->GetRequestedRegion());
    for (It.GoToBegin(); !It.IsAtEnd(); ++It)
    {
      if (!this->m_MaskImage || this->m_MaskImage->GetPixel(It.GetIndex()) != NumericTraits<LabelType>::ZeroValue())
      {
        this->m_LabelSet.insert(It.Get());
      }
    }
  }

  // Initialize the posterior maps
  this->m_LabelPosteriorProbabilityImages.clear();

  typename LabelSetType::const_iterator labelIt;
  for (labelIt = this->m_LabelSet.begin(); labelIt != this->m_LabelSet.end(); ++labelIt)
  {
    typename ProbabilityImageType::Pointer labelProbabilityImage = ProbabilityImageType::New();
    labelProbabilityImage->CopyInformation(this->m_TargetImage[0]);
    labelProbabilityImage->SetRegions(this->m_TargetImage[0]->GetRequestedRegion());
    labelProbabilityImage->SetLargestPossibleRegion(this->m_TargetImage[0]->GetLargestPossibleRegion());
    labelProbabilityImage->AllocateInitialized();

    this->m_LabelPosteriorProbabilityImages.insert(
      std::pair<LabelType, ProbabilityImagePointer>(*labelIt, labelProbabilityImage));
  }

  // Initialize the atlas voting weight images
  if (this->m_RetainAtlasVotingWeightImages)
  {
    this->m_AtlasVotingWeightImages.clear();
    this->m_AtlasVotingWeightImages.resize(this->m_NumberOfAtlases);

    for (SizeValueType i = 0; i < this->m_NumberOfAtlases; i++)
    {
      this->m_AtlasVotingWeightImages[i] = ProbabilityImageType::New();
      this->m_AtlasVotingWeightImages[i]->CopyInformation(this->m_TargetImage[0]);
      this->m_AtlasVotingWeightImages[i]->SetRegions(this->m_TargetImage[0]->GetRequestedRegion());
      this->m_AtlasVotingWeightImages[i]->SetLargestPossibleRegion(this->m_TargetImage[0]->GetLargestPossibleRegion());
      this->m_AtlasVotingWeightImages[i]->AllocateInitialized();
    }
  }

  // Do the joint intensity fusion
  this->m_JointIntensityFusionImage.clear();
  this->m_JointIntensityFusionImage.resize(this->m_NumberOfAtlasModalities);

  for (SizeValueType i = 0; i < this->m_NumberOfAtlasModalities; i++)
  {
    this->m_JointIntensityFusionImage[i] = InputImageType::New();
    this->m_JointIntensityFusionImage[i]->CopyInformation(this->m_TargetImage[0]);
    this->m_JointIntensityFusionImage[i]->SetRegions(this->m_TargetImage[0]->GetRequestedRegion());
    this->m_JointIntensityFusionImage[i]->SetLargestPossibleRegion(this->m_TargetImage[0]->GetLargestPossibleRegion());
    this->m_JointIntensityFusionImage[i]->AllocateInitialized();
  }

  // Initialize the weight sum image
  this->m_WeightSumImage = ProbabilityImageType::New();
  this->m_WeightSumImage->CopyInformation(this->m_TargetImage[0]);
  this->m_WeightSumImage->SetRegions(this->m_TargetImage[0]->GetRequestedRegion());
  this->m_WeightSumImage->SetLargestPossibleRegion(this->m_TargetImage[0]->GetLargestPossibleRegion());
  this->m_WeightSumImage->AllocateInitialized();

  // Initialize the count image
  this->m_CountImage = CountImageType::New();
  this->m_CountImage->CopyInformation(this->m_TargetImage[0]);
  this->m_CountImage->SetRegions(this->m_TargetImage[0]->GetRequestedRegion());
  this->m_CountImage->SetLargestPossibleRegion(this->m_TargetImage[0]->GetLargestPossibleRegion());
  this->m_CountImage->AllocateInitialized();

  // Determine the ordered search offset list (or map if an search radius image is specified)

  typename InputImageType::SpacingType spacing = this->m_TargetImage[0]->GetSpacing();

  NeighborhoodOffsetListType orderedNeighborhoodSearchOffsetList;
  orderedNeighborhoodSearchOffsetList.clear();

  this->m_NeighborhoodSearchOffsetSetsMap.clear();

  if (this->m_NeighborhoodSearchRadiusImage.IsNull())
  {
    ConstNeighborhoodIterator<InputImageType> It(
      this->GetNeighborhoodSearchRadius(), this->GetInput(), this->GetInput()->GetRequestedRegion());

    DistanceIndexVectorType squaredDistances;
    squaredDistances.resize(this->GetNeighborhoodSearchSize());

    for (unsigned int n = 0; n < this->GetNeighborhoodSearchSize(); n++)
    {
      NeighborhoodOffsetType offset = (It.GetNeighborhood()).GetOffset(n);

      squaredDistances[n].first = n;
      squaredDistances[n].second = 0.0;
      for (unsigned int d = 0; d < ImageDimension; d++)
      {
        squaredDistances[n].second += itk::Math::sqr(offset[d] * spacing[d]);
      }
    }
    std::sort(squaredDistances.begin(), squaredDistances.end(), DistanceIndexComparator());

    for (unsigned int n = 0; n < this->GetNeighborhoodSearchSize(); n++)
    {
      orderedNeighborhoodSearchOffsetList.push_back((It.GetNeighborhood()).GetOffset(squaredDistances[n].first));
    }
    this->SetNeighborhoodSearchOffsetList(orderedNeighborhoodSearchOffsetList);
  }
  else
  {
    ImageRegionConstIterator<RadiusImageType> ItR(this->m_NeighborhoodSearchRadiusImage,
                                                  this->m_NeighborhoodSearchRadiusImage->GetRequestedRegion());

    for (ItR.GoToBegin(); !ItR.IsAtEnd(); ++ItR)
    {
      RadiusValueType localSearchRadius = ItR.Get();
      if (localSearchRadius > 0 && this->m_NeighborhoodSearchOffsetSetsMap.find(localSearchRadius) ==
                                     this->m_NeighborhoodSearchOffsetSetsMap.end())
      {
        NeighborhoodRadiusType localNeighborhoodSearchRadius;
        localNeighborhoodSearchRadius.Fill(localSearchRadius);

        std::vector<NeighborhoodOffsetType> localNeighborhoodSearchOffsetList;

        ConstNeighborhoodIterator<InputImageType> It(
          localNeighborhoodSearchRadius, this->GetInput(), this->GetInput()->GetRequestedRegion());

        RadiusValueType localNeighborhoodSearchSize = (It.GetNeighborhood()).Size();

        DistanceIndexVectorType squaredDistances;
        squaredDistances.resize(localNeighborhoodSearchSize);

        for (unsigned int n = 0; n < localNeighborhoodSearchSize; n++)
        {
          NeighborhoodOffsetType offset = (It.GetNeighborhood()).GetOffset(n);

          squaredDistances[n].first = n;
          squaredDistances[n].second = 0.0;
          for (unsigned int d = 0; d < ImageDimension; d++)
          {
            squaredDistances[n].second += itk::Math::sqr(offset[d] * spacing[d]);
          }
        }
        std::sort(squaredDistances.begin(), squaredDistances.end(), DistanceIndexComparator());

        for (unsigned int n = 0; n < localNeighborhoodSearchSize; n++)
        {
          localNeighborhoodSearchOffsetList.push_back((It.GetNeighborhood()).GetOffset(squaredDistances[n].first));
        }
        this->m_NeighborhoodSearchOffsetSetsMap[localSearchRadius] = localNeighborhoodSearchOffsetList;
      }
    }
  }

  this->AllocateOutputs();
}

template <typename TInputImage, typename TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>::ThreadedGenerateData(const RegionType & region,
                                                                                 ThreadIdType       threadId)
{
  if (!this->m_IsWeightedAveragingComplete)
  {
    this->ThreadedGenerateDataForWeightedAveraging(region, threadId);
  }
  else
  {
    this->ThreadedGenerateDataForReconstruction(region, threadId);
  }
}

template <typename TInputImage, typename TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>::ThreadedGenerateDataForWeightedAveraging(
  const RegionType & region,
  ThreadIdType       threadId)
{
  ProgressReporter progress(this, threadId, region.GetNumberOfPixels(), 100);

  typename OutputImageType::Pointer output = this->GetOutput();

  SizeValueType numberOfTargetModalities = this->m_TargetImage.size();

  MatrixType absoluteAtlasPatchDifferences(this->m_NumberOfAtlases,
                                           this->GetNeighborhoodPatchSize() * numberOfTargetModalities);

  MatrixType originalAtlasPatchIntensities(this->m_NumberOfAtlases,
                                           this->GetNeighborhoodPatchSize() * this->m_NumberOfAtlasModalities);

  std::vector<SizeValueType> minimumAtlasOffsetIndices(this->m_NumberOfAtlases);

  bool useOnlyFirstAtlasImage = true;
  if (numberOfTargetModalities == this->m_NumberOfAtlasModalities)
  {
    useOnlyFirstAtlasImage = false;
  }

  std::vector<NeighborhoodOffsetType> searchNeighborhoodOffsetList = this->GetNeighborhoodSearchOffsetList();

  // Iterate over the input region
  ConstNeighborhoodIteratorType ItN(this->GetNeighborhoodPatchRadius(), this->m_TargetImage[0], region);
  for (ItN.GoToBegin(); !ItN.IsAtEnd(); ++ItN)
  {
    progress.CompletedPixel();

    IndexType currentCenterIndex = ItN.GetIndex();

    if (this->m_MaskImage && this->m_MaskImage->GetPixel(currentCenterIndex) == NumericTraits<LabelType>::ZeroValue())
    {
      continue;
    }

    // Do not do the following check from Paul's original code.  Since we're incorporating
    // joint intensity fusion, we want to calculate at every voxel (except outside of a
    // possible mask) even if there are no segmentation labels at that voxel.

    //     // Check to see if there are any non-zero labels
    //     if( this->m_NumberOfAtlasSegmentations > 0 )
    //       {
    //       bool nonBackgroundLabelExistAtThisVoxel = false;
    //       for( SizeValueType i = 0; i < this->m_NumberOfAtlasSegmentations; i++ )
    //         {
    //         if( this->m_AtlasSegmentations[i]->GetPixel( currentCenterIndex ) > 0 )
    //           {
    //           nonBackgroundLabelExistAtThisVoxel = true;
    //           break;
    //           }
    //         }
    //       if( ! nonBackgroundLabelExistAtThisVoxel )
    //         {
    //         continue;
    //         }
    //       }

    // Determine the search neighborhood offset list for the current center voxel
    if (this->m_NeighborhoodSearchRadiusImage.IsNotNull())
    {
      RadiusValueType localSearchRadius = this->m_NeighborhoodSearchRadiusImage->GetPixel(currentCenterIndex);
      if (localSearchRadius <= 0)
      {
        continue;
      }
      searchNeighborhoodOffsetList = this->m_NeighborhoodSearchOffsetSetsMap[localSearchRadius];
    }
    SizeValueType searchNeighborhoodSize = searchNeighborhoodOffsetList.size();

    // if Metric is MSQ, create both target patch and normalized target patch for use.
    // target patch is for metric calculation, and normalized patch is for weight calculation.
    // If metric is PC, only need to use normalizedTargetPatch

    const InputImagePixelVectorType normalizedTargetPatch =
      this->VectorizeImageListPatch(this->m_TargetImage, currentCenterIndex, true);

    InputImagePixelVectorType targetPatch;
    if (this->m_SimilarityMetric == NonLocalPatchBasedImageFilterEnums::SimilarityMetric::MEAN_SQUARES)
    {
      targetPatch = this->VectorizeImageListPatch(this->m_TargetImage, currentCenterIndex, false);
    }

    absoluteAtlasPatchDifferences.fill(0.0);
    originalAtlasPatchIntensities.fill(0.0);

    // In each atlas, search for a patch that matches the target patch
    for (SizeValueType i = 0; i < this->m_NumberOfAtlases; i++)
    {

      RealType      minimumPatchSimilarity = NumericTraits<RealType>::max();
      SizeValueType minimumPatchOffsetIndex = 0;

      for (SizeValueType j = 0; j < searchNeighborhoodSize; j++)
      {
        IndexType searchIndex = currentCenterIndex + searchNeighborhoodOffsetList[j];

        if (!output->GetRequestedRegion().IsInside(searchIndex))
        {
          continue;
        }

        const RealType patchSimilarity = [&]() -> RealType {
          // use non-normalized vector for MSQ
          switch (this->m_SimilarityMetric)
          {
            case NonLocalPatchBasedImageFilterEnums::SimilarityMetric::MEAN_SQUARES:
              // use non-normalized vectors for MSE
              return this->ComputeNeighborhoodPatchSimilarity(
                this->m_AtlasImages[i], searchIndex, targetPatch, useOnlyFirstAtlasImage);
              break;
            case NonLocalPatchBasedImageFilterEnums::SimilarityMetric::PEARSON_CORRELATION:
              // use normalized vector for PC
              return this->ComputeNeighborhoodPatchSimilarity(
                this->m_AtlasImages[i], searchIndex, normalizedTargetPatch, useOnlyFirstAtlasImage);
              break;
            default:
              itkGenericExceptionMacro("Invalid SimilarityMetric Chosen.");
          }
          return 0.0;
        }();

        if (patchSimilarity < minimumPatchSimilarity)
        {
          minimumPatchSimilarity = patchSimilarity;
          minimumPatchOffsetIndex = j;
        }
      }

      // Once the patch has been found, normalize it and then compute the absolute
      // difference with target patch

      IndexType minimumIndex = currentCenterIndex + searchNeighborhoodOffsetList[minimumPatchOffsetIndex];
      InputImagePixelVectorType normalizedMinimumAtlasPatch;
      if (numberOfTargetModalities == this->m_NumberOfAtlasModalities)
      {
        normalizedMinimumAtlasPatch = this->VectorizeImageListPatch(this->m_AtlasImages[i], minimumIndex, true);
      }
      else
      {
        normalizedMinimumAtlasPatch = this->VectorizeImagePatch(this->m_AtlasImages[i][0], minimumIndex, true);
      }

      typename InputImagePixelVectorType::const_iterator itA = normalizedMinimumAtlasPatch.begin();
      typename InputImagePixelVectorType::const_iterator itT = normalizedTargetPatch.begin();
      while (itA != normalizedMinimumAtlasPatch.end())
      {
        RealType value = std::fabs(*itA - *itT);
        absoluteAtlasPatchDifferences(i, itA - normalizedMinimumAtlasPatch.begin()) = value;

        ++itA;
        ++itT;
      }

      InputImagePixelVectorType originalMinimumAtlasPatch =
        this->VectorizeImageListPatch(this->m_AtlasImages[i], minimumIndex, false);

      typename InputImagePixelVectorType::const_iterator itO = originalMinimumAtlasPatch.begin();
      while (itO != originalMinimumAtlasPatch.end())
      {
        originalAtlasPatchIntensities(i, itO - originalMinimumAtlasPatch.begin()) = *itO;
        ++itO;
      }

      minimumAtlasOffsetIndices[i] = minimumPatchOffsetIndex;
    }

    // Allocate Mx
    MatrixType Mx(this->m_NumberOfAtlases, this->m_NumberOfAtlases);

    // Compute Mx values
    for (SizeValueType i = 0; i < this->m_NumberOfAtlases; i++)
    {
      for (SizeValueType j = 0; j <= i; j++)
      {
        RealType mxValue = 0.0;

        for (unsigned int k = 0; k < this->GetNeighborhoodPatchSize() * numberOfTargetModalities; k++)
        {
          mxValue += absoluteAtlasPatchDifferences[i][k] * absoluteAtlasPatchDifferences[j][k];
        }
        mxValue /= static_cast<RealType>(this->GetNeighborhoodPatchSize() - 1);

        if (!itk::Math::FloatAlmostEqual(this->m_Beta, NumericTraits<RealType>::OneValue()))
        {
          if (itk::Math::FloatAlmostEqual(this->m_Beta, static_cast<RealType>(2.0)))
          {
            mxValue *= mxValue;
          }
          else
          {
            mxValue = std::pow(mxValue, this->m_Beta);
          }
        }

        if (!std::isfinite(mxValue))
        {
          mxValue = 0.0;
        }

        Mx(i, j) = Mx(j, i) = mxValue;
      }
    }

    // Compute the weights by solving for the inverse of Mx
    MatrixType MxBar(this->m_NumberOfAtlases, this->m_NumberOfAtlases, 0.0);
    MxBar.fill_diagonal(this->m_Alpha);
    MxBar += Mx;

    // Define a vector of all ones
    VectorType ones(this->m_NumberOfAtlases, 1.0);
    VectorType W(this->m_NumberOfAtlases, 1.0);

    if (this->m_ConstrainSolutionToNonnegativeWeights)
    {
      W = this->NonNegativeLeastSquares(MxBar, ones, 1e-6);
    }
    else
    {
      vnl_cholesky cholesky(MxBar, vnl_cholesky::estimate_condition);
      if (cholesky.rcond() > itk::Math::sqrteps)
      {
        // well-conditioned matrix
        W = cholesky.solve(ones);
      }
      else
      {
        // ill-conditioned matrix
        W = vnl_svd<RealType>(MxBar).solve(ones);
      }

      for (double & i : W)
      {
        if (i < 0.0)
        {
          i = 0.0;
        }
      }
    }

    // Normalize the weights
    W *= 1.0 / dot_product(W, ones);

    // Do joint intensity fusion
    VectorType estimatedNeighborhoodIntensities = W;

    estimatedNeighborhoodIntensities.post_multiply(originalAtlasPatchIntensities);

    for (SizeValueType i = 0; i < this->m_NumberOfAtlasModalities; i++)
    {
      for (SizeValueType j = 0; j < this->GetNeighborhoodPatchSize(); j++)
      {
        IndexType neighborhoodIndex = ItN.GetIndex(j);

        if (!output->GetRequestedRegion().IsInside(neighborhoodIndex))
        {
          continue;
        }

        if (this->m_MaskImage &&
            this->m_MaskImage->GetPixel(neighborhoodIndex) == NumericTraits<LabelType>::ZeroValue())
        {
          continue;
        }

        RealType estimatedValue =
          (static_cast<RealType>(estimatedNeighborhoodIntensities[i * this->GetNeighborhoodPatchSize() + j]) +
           static_cast<RealType>(this->m_JointIntensityFusionImage[i]->GetPixel(neighborhoodIndex)));

        if (!std::isfinite(estimatedValue))
        {
          estimatedValue = 0.0;
        }

        this->m_JointIntensityFusionImage[i]->SetPixel(neighborhoodIndex,
                                                       static_cast<InputImagePixelType>(estimatedValue));
        if (i == 0)
        {
          this->m_CountImage->SetPixel(neighborhoodIndex, this->m_CountImage->GetPixel(neighborhoodIndex) + 1);
        }
      }
    }

    if (this->m_NumberOfAtlasSegmentations > 0)
    {
      // Perform voting using Hongzhi's averaging scheme. Iterate over all segmentation patches
      for (SizeValueType n = 0; n < this->GetNeighborhoodPatchSize(); n++)
      {
        IndexType neighborhoodIndex = ItN.GetIndex(n);
        if (!output->GetRequestedRegion().IsInside(neighborhoodIndex))
        {
          continue;
        }

        for (SizeValueType i = 0; i < this->m_NumberOfAtlasSegmentations; i++)
        {
          // The segmentation at the corresponding patch location in atlas i
          IndexType minimumIndex = neighborhoodIndex + searchNeighborhoodOffsetList[minimumAtlasOffsetIndices[i]];

          if (!output->GetRequestedRegion().IsInside(minimumIndex))
          {
            continue;
          }

          LabelType label = this->m_AtlasSegmentations[i]->GetPixel(minimumIndex);

          if (this->m_LabelSet.find(label) == this->m_LabelSet.end())
          {
            continue;
          }

          // Add that weight the posterior map for voxel at idx
          this->m_LabelPosteriorProbabilityImages[label]->SetPixel(
            neighborhoodIndex,
            this->m_LabelPosteriorProbabilityImages[label]->GetPixel(neighborhoodIndex) + static_cast<float>(W[i]));
          this->m_WeightSumImage->SetPixel(
            neighborhoodIndex, this->m_WeightSumImage->GetPixel(neighborhoodIndex) + static_cast<float>(W[i]));

          if (this->m_RetainAtlasVotingWeightImages)
          {
            this->m_AtlasVotingWeightImages[i]->SetPixel(
              neighborhoodIndex,
              this->m_AtlasVotingWeightImages[i]->GetPixel(neighborhoodIndex) + static_cast<float>(W[i]));
          }
        }
      }
    }
  }
}

template <typename TInputImage, typename TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>::ThreadedGenerateDataForReconstruction(
  const RegionType & region,
  ThreadIdType       threadId)
{
  ProgressReporter progress(this, threadId, 2 * region.GetNumberOfPixels(), 100);

  typename OutputImageType::Pointer output = this->GetOutput();

  // Perform voting at each voxel
  ImageRegionIteratorWithIndex<OutputImageType> It(output, region);

  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    progress.CompletedPixel();

    IndexType index = It.GetIndex();

    if (this->m_MaskImage && this->m_MaskImage->GetPixel(It.GetIndex()) == NumericTraits<LabelType>::ZeroValue())
    {
      continue;
    }

    RealType  maxPosteriorProbability = 0.0;
    LabelType winningLabel = NumericTraits<LabelType>::ZeroValue();

    typename LabelSetType::const_iterator labelIt;
    for (labelIt = this->m_LabelSet.begin(); labelIt != this->m_LabelSet.end(); ++labelIt)
    {
      // check if the label is excluded
      typename LabelExclusionMap::const_iterator xIt = this->m_LabelExclusionImages.find(*labelIt);
      bool isLabelExcluded = (xIt != m_LabelExclusionImages.end() && xIt->second->GetPixel(index) != 0);

      if (!isLabelExcluded)
      {
        typename ProbabilityImageType::PixelType posteriorProbability =
          this->m_LabelPosteriorProbabilityImages[*labelIt]->GetPixel(index);

        // Vote!
        if (maxPosteriorProbability < static_cast<RealType>(posteriorProbability))
        {
          maxPosteriorProbability = static_cast<RealType>(posteriorProbability);
          winningLabel = *labelIt;
        }
      }
    }
    It.Set(winningLabel);
  }

  if (this->m_RetainLabelPosteriorProbabilityImages || this->m_RetainAtlasVotingWeightImages)
  {
    ImageRegionIteratorWithIndex<ProbabilityImageType> ItW(this->m_WeightSumImage, region);

    for (ItW.GoToBegin(); !ItW.IsAtEnd(); ++ItW)
    {
      progress.CompletedPixel();

      typename ProbabilityImageType::PixelType weightSum = ItW.Get();

      IndexType index = ItW.GetIndex();

      if (weightSum < static_cast<typename ProbabilityImageType::PixelType>(0.1))
      {
        continue;
      }

      if (this->m_RetainLabelPosteriorProbabilityImages)
      {
        typename LabelSetType::const_iterator labelIt;
        for (labelIt = this->m_LabelSet.begin(); labelIt != this->m_LabelSet.end(); ++labelIt)
        {
          typename ProbabilityImageType::PixelType labelProbability =
            this->m_LabelPosteriorProbabilityImages[*labelIt]->GetPixel(index);
          this->m_LabelPosteriorProbabilityImages[*labelIt]->SetPixel(index, labelProbability / weightSum);
        }
      }

      if (this->m_RetainAtlasVotingWeightImages)
      {
        for (SizeValueType i = 0; i < this->m_NumberOfAtlases; i++)
        {
          typename ProbabilityImageType::PixelType votingWeight = this->m_AtlasVotingWeightImages[i]->GetPixel(index);
          this->m_AtlasVotingWeightImages[i]->SetPixel(index, votingWeight / weightSum);
        }
      }
    }
  }
}

template <typename TInputImage, typename TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>::AfterThreadedGenerateData()
{
  // Clear posterior maps if not kept
  if (!this->m_RetainLabelPosteriorProbabilityImages)
  {
    this->m_LabelPosteriorProbabilityImages.clear();
  }

  // Normalize the joint intensity fusion images.
  for (SizeValueType i = 0; i < this->m_NumberOfAtlasModalities; i++)
  {
    ImageRegionIterator<InputImageType> ItJ(this->m_JointIntensityFusionImage[i],
                                            this->m_JointIntensityFusionImage[i]->GetRequestedRegion());
    ImageRegionIterator<CountImageType> ItC(this->m_CountImage, this->m_CountImage->GetRequestedRegion());

    for (ItJ.GoToBegin(), ItC.GoToBegin(); !ItJ.IsAtEnd(); ++ItJ, ++ItC)
    {
      typename CountImageType::PixelType count = ItC.Get();

      if (count > 0)
      {
        ItJ.Set(ItJ.Get() / static_cast<InputImagePixelType>(count));
      }
    }
  }
}

template <typename TInputImage, typename TOutputImage>
typename WeightedVotingFusionImageFilter<TInputImage, TOutputImage>::VectorType
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>::NonNegativeLeastSquares(const MatrixType & A,
                                                                                    const VectorType & y,
                                                                                    const RealType     tolerance)
{
  // Algorithm based on
  // Lawson, Charles L.; Hanson, Richard J. (1995). Solving Least Squares Problems. SIAM.
  // cf https://en.wikipedia.org/wiki/Non-negative_least_squares

  SizeValueType m = A.rows();
  SizeValueType n = A.cols();

  // This fortran implementation sets a maximum iteration number of 3 times the
  // number of columns:
  //    http://www.netlib.org/lawson-hanson/all

  const SizeValueType maximumNumberOfIterations = 3 * n;

  // Initialization

  VectorType P(n, 0);
  VectorType R(n, 1);
  VectorType x(n, 0);
  VectorType s(n, 0);
  VectorType w = A.transpose() * (y - A * x);

  RealType      wMaxValue = w.max_value();
  SizeValueType maxIndex = NumericTraits<SizeValueType>::max();
  wMaxValue = NumericTraits<RealType>::NonpositiveMin();
  for (SizeValueType i = 0; i < n; i++)
  {
    if (itk::Math::FloatAlmostEqual(R[i], NumericTraits<RealType>::OneValue()) && wMaxValue < w[i])
    {
      maxIndex = i;
      wMaxValue = w[i];
    }
  }

  // Outer loop

  SizeValueType numberOfIterations = 0;
  while (R.sum() > 0 && wMaxValue > tolerance && numberOfIterations++ < maximumNumberOfIterations)
  {
    P[maxIndex] = 1;
    R[maxIndex] = 0;

    SizeValueType sizeP = P.sum();

    MatrixType    AP(m, sizeP, 0);
    SizeValueType jIndex = 0;
    for (SizeValueType j = 0; j < n; j++)
    {
      if (itk::Math::FloatAlmostEqual(P[j], NumericTraits<RealType>::OneValue()))
      {
        AP.set_column(jIndex++, A.get_column(j));
      }
    }

    VectorType sP = vnl_svd<RealType>(AP).pinverse() * y;

    SizeValueType iIndex = 0;

    for (SizeValueType i = 0; i < n; i++)
    {
      if (!itk::Math::FloatAlmostEqual(R[i], NumericTraits<RealType>::ZeroValue()))
      {
        s[i] = 0;
      }
      else
      {
        s[i] = sP[iIndex++];
      }
    }

    // Inner loop
    while (sP.min_value() <= tolerance && sizeP > 0)
    {
      RealType alpha = NumericTraits<RealType>::max();

      for (SizeValueType i = 0; i < n; i++)
      {
        if (itk::Math::FloatAlmostEqual(P[i], NumericTraits<RealType>::OneValue()) && s[i] <= tolerance)
        {
          RealType value = x[i] / (x[i] - s[i]);
          if (value < alpha)
          {
            alpha = value;
          }
        }
      }

      x += alpha * (s - x);

      for (SizeValueType i = 0; i < n; i++)
      {
        if (itk::Math::FloatAlmostEqual(P[i], NumericTraits<RealType>::OneValue()) && std::fabs(x[i]) < tolerance)
        {
          P[i] = 0;
          R[i] = 1;
        }
      }

      sizeP = P.sum();
      if (sizeP == 0)
      {
        break;
      }

      AP.set_size(m, sizeP);
      jIndex = 0;
      for (SizeValueType j = 0; j < n; j++)
      {
        if (itk::Math::FloatAlmostEqual(P[j], NumericTraits<RealType>::OneValue()))
        {
          AP.set_column(jIndex++, A.get_column(j));
        }
      }

      sP = vnl_svd<RealType>(AP).pinverse() * y;

      iIndex = 0;
      for (SizeValueType i = 0; i < n; i++)
      {
        if (!itk::Math::FloatAlmostEqual(R[i], NumericTraits<RealType>::ZeroValue()))
        {
          s[i] = 0;
        }
        else
        {
          s[i] = sP[iIndex++];
        }
      }
    }

    x = s;
    w = A.transpose() * (y - A * x);

    maxIndex = NumericTraits<SizeValueType>::max();
    wMaxValue = NumericTraits<RealType>::NonpositiveMin();
    for (SizeValueType i = 0; i < n; i++)
    {
      if (itk::Math::FloatAlmostEqual(R[i], NumericTraits<RealType>::OneValue()) && wMaxValue < w[i])
      {
        maxIndex = i;
        wMaxValue = w[i];
      }
    }
  }

  return x;
}

template <typename TInputImage, typename TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << "Number of atlases = " << this->m_NumberOfAtlases << std::endl;
  os << "Number of atlas segmentations = " << this->m_NumberOfAtlasSegmentations << std::endl;
  os << "Number of atlas modalities = " << this->m_NumberOfAtlasModalities << std::endl;
  os << "Alpha = " << this->m_Alpha << std::endl;
  os << "Beta = " << this->m_Beta << std::endl;
  if (this->m_ConstrainSolutionToNonnegativeWeights)
  {
    os << "Constrain solution to positive weights using NNLS." << std::endl;
  }

  os << "Label set: ";
  typename LabelSetType::const_iterator labelIt;
  for (labelIt = this->m_LabelSet.begin(); labelIt != this->m_LabelSet.end(); ++labelIt)
  {
    os << *labelIt << " ";
  }
  os << std::endl;
}

} // end namespace itk

#endif
