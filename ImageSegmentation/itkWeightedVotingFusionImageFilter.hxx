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

#include "itkWeightedVotingFusionImageFilter.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkProgressReporter.h"

#include <algorithm>
#include <numeric>

#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_inverse.h>

namespace itk {

template<typename TInputImage, typename TOutputImage>
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::WeightedVotingFusionImageFilter() :
  m_IsWeightedAveragingComplete( false ),
  m_NumberOfAtlases( 0 ),
  m_NumberOfAtlasSegmentations( 0 ),
  m_NumberOfAtlasModalities( 0 ),
  m_PatchNeighborhoodSize( 0 ),
  m_Alpha( 0.1 ),
  m_Beta( 2.0 ),
  m_RetainLabelPosteriorProbabilityImages( false ),
  m_RetainAtlasVotingWeightImages( false ),
  m_ConstrainSolutionToNonnegativeWeights( false )
{
  this->m_MaskImage = ITK_NULLPTR;

  this->m_CountImage = ITK_NULLPTR;

  this->m_SearchNeighborhoodRadiusImage = ITK_NULLPTR;
  this->m_SearchNeighborhoodRadius.Fill( 3 );

  this->m_PatchNeighborhoodRadius.Fill( 2 );

  this->m_SimilarityMetric = PEARSON_CORRELATION;
}

template <class TInputImage, class TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::UpdateInputs()
{
  // Set all the inputs

  this->SetNumberOfIndexedInputs( this->m_NumberOfAtlases * this->m_NumberOfAtlasModalities +
    this->m_NumberOfAtlasSegmentations + this->m_TargetImage.size() + this->m_LabelExclusionImages.size() );

  SizeValueType nthInput = 0;

  for( SizeValueType i = 0; i < this->m_TargetImage.size(); i++ )
    {
    this->SetNthInput( nthInput++, this->m_TargetImage[i] );
    }

  for( SizeValueType i = 0; i < this->m_NumberOfAtlases; i++ )
    {
    for( SizeValueType j = 0; j < this->m_NumberOfAtlasModalities; j++ )
      {
      this->SetNthInput( nthInput++, this->m_AtlasImages[i][j] );
      }
    }

  for( SizeValueType i = 0; i < this->m_NumberOfAtlasSegmentations; i++ )
    {
    this->SetNthInput( nthInput++, this->m_AtlasSegmentations[i] );
    }

  typename LabelExclusionMap::const_iterator it;
  for( it = m_LabelExclusionImages.begin(); it != m_LabelExclusionImages.end(); ++it )
    {
    this->SetNthInput( nthInput++, it->second );
    }

  if( this->m_MaskImage.IsNotNull() )
    {
    this->SetNthInput( nthInput++, this->m_MaskImage );
    }

  this->Modified();
}

template <class TInputImage, class TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

  // Get the output requested region
  RegionType outRegion = this->GetOutput()->GetRequestedRegion();

  // Pad this region by the search window and patch size
  if( this->m_SearchNeighborhoodRadiusImage.IsNull() )
    {
    outRegion.PadByRadius( this->m_SearchNeighborhoodRadius );
    }
  else
    {
    NeighborhoodRadiusType maxSearchNeighborhoodRadius;
    maxSearchNeighborhoodRadius.Fill( 0 );

    ImageRegionConstIterator<RadiusImageType> ItR( this->m_SearchNeighborhoodRadiusImage,
      this->m_SearchNeighborhoodRadiusImage->GetRequestedRegion() );
    for( ItR.GoToBegin(); !ItR.IsAtEnd(); ++ItR )
      {
      RadiusValueType localSearchRadius = ItR.Get();
      if( localSearchRadius > maxSearchNeighborhoodRadius[0] )
        {
        maxSearchNeighborhoodRadius.Fill( localSearchRadius );
        }
      }
    outRegion.PadByRadius( maxSearchNeighborhoodRadius );
    }
  outRegion.PadByRadius( this->m_PatchNeighborhoodRadius );

  // Iterate over all the inputs to this filter

  for( SizeValueType i = 0; i < this->m_TargetImage.size(); i++ )
    {
    InputImageType *input = this->m_TargetImage[i];
    if( i == 0 )
      {
      this->m_TargetImageRequestedRegion = input->GetRequestedRegion();
      }
    RegionType region = outRegion;
    region.Crop( input->GetLargestPossibleRegion() );
    input->SetRequestedRegion( region );
    }

  for( SizeValueType i = 0; i < this->m_NumberOfAtlases; i++ )
    {
    for( SizeValueType j = 0; j < this->m_NumberOfAtlasModalities; j++ )
      {
      InputImageType *input = this->m_AtlasImages[i][j];
      RegionType region = outRegion;
      region.Crop( input->GetLargestPossibleRegion() );
      input->SetRequestedRegion( region );
      }
    }

  for( SizeValueType i = 0; i < this->m_NumberOfAtlasSegmentations; i++ )
    {
    LabelImageType *input = this->m_AtlasSegmentations[i];
    RegionType region = outRegion;
    region.Crop( input->GetLargestPossibleRegion() );
    input->SetRequestedRegion( region );
    }

  typename LabelExclusionMap::const_iterator it;
  for( it = m_LabelExclusionImages.begin(); it != m_LabelExclusionImages.end(); ++it )
    {
    LabelImageType *input = it->second;
    RegionType region = outRegion;
    region.Crop( input->GetLargestPossibleRegion() );
    input->SetRequestedRegion( region );
    }

  if( this->m_MaskImage.IsNotNull() )
    {
    MaskImageType *input = this->m_MaskImage;
    RegionType region = outRegion;
    region.Crop( input->GetLargestPossibleRegion() );
    input->SetRequestedRegion( region );
    }
}

template <class TInputImage, class TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  this->BeforeThreadedGenerateData();

  /**
   * Multithread processing for the weighted averaging
   */
  typename ImageSource<TOutputImage>::ThreadStruct str1;
  str1.Filter = this;

  this->GetMultiThreader()->SetNumberOfThreads( this->GetNumberOfThreads() );
  this->GetMultiThreader()->SetSingleMethod( this->ThreaderCallback, &str1 );

  this->GetMultiThreader()->SingleMethodExecute();

  this->m_IsWeightedAveragingComplete = true;

  /**
   * Multithread processing for the image(s) reconstruction
   */

  typename ImageSource<TOutputImage>::ThreadStruct str2;
  str2.Filter = this;

  this->GetMultiThreader()->SetNumberOfThreads( this->GetNumberOfThreads() );
  this->GetMultiThreader()->SetSingleMethod( this->ThreaderCallback, &str2 );

  this->GetMultiThreader()->SingleMethodExecute();

  this->AfterThreadedGenerateData();
}

template <class TInputImage, class TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
  if( this->m_NumberOfAtlasSegmentations != this->m_NumberOfAtlases )
    {
    // Set the number of atlas segmentations to 0 since we're just going to
    // doing joint intensity fusion
    this->m_NumberOfAtlasSegmentations = 0;
    }

  // Check to see if the number of target images is equal to 1 or equal to the number
  // of atlas modalities
  if( this->m_TargetImage.size() != 1 && this->m_TargetImage.size() != this->m_NumberOfAtlasModalities )
    {
    itkExceptionMacro( "The number of target images must be 1 or must be the number of atlas modalities." );
    }

  // Find all the unique labels in the atlas segmentations
  this->m_LabelSet.clear();
  for( unsigned int i = 0; i < this->m_NumberOfAtlasSegmentations; i++ )
    {
    ImageRegionConstIteratorWithIndex<LabelImageType> It(
      this->m_AtlasSegmentations[i], this->m_AtlasSegmentations[i]->GetRequestedRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      if( !this->m_MaskImage ||
          this->m_MaskImage->GetPixel( It.GetIndex() ) != NumericTraits<LabelType>::ZeroValue() )
        {
        this->m_LabelSet.insert( It.Get() );
        }
      }
    }

  // Initialize the posterior maps
  this->m_LabelPosteriorProbabilityImages.clear();

  typename LabelSetType::const_iterator labelIt;
  for( labelIt = this->m_LabelSet.begin(); labelIt != this->m_LabelSet.end(); ++labelIt )
    {
    typename ProbabilityImageType::Pointer labelProbabilityImage = ProbabilityImageType::New();
    labelProbabilityImage->CopyInformation( this->m_TargetImage[0] );
    labelProbabilityImage->SetRegions( this->m_TargetImage[0]->GetRequestedRegion() );
    labelProbabilityImage->SetLargestPossibleRegion( this->m_TargetImage[0]->GetLargestPossibleRegion() );
    labelProbabilityImage->Allocate();
    labelProbabilityImage->FillBuffer( 0.0 );

    this->m_LabelPosteriorProbabilityImages.insert(
      std::pair<LabelType, ProbabilityImagePointer>( *labelIt, labelProbabilityImage ) );
    }

  // Initialize the atlas voting weight images
  if( this->m_RetainAtlasVotingWeightImages )
    {
    this->m_AtlasVotingWeightImages.clear();
    this->m_AtlasVotingWeightImages.resize( this->m_NumberOfAtlases );

    for( SizeValueType i = 0; i < this->m_NumberOfAtlases; i++ )
      {
      this->m_AtlasVotingWeightImages[i] = ProbabilityImageType::New();
      this->m_AtlasVotingWeightImages[i]->CopyInformation( this->m_TargetImage[0] );
      this->m_AtlasVotingWeightImages[i]->SetRegions( this->m_TargetImage[0]->GetRequestedRegion() );
      this->m_AtlasVotingWeightImages[i]->SetLargestPossibleRegion( this->m_TargetImage[0]->GetLargestPossibleRegion() );
      this->m_AtlasVotingWeightImages[i]->Allocate();
      this->m_AtlasVotingWeightImages[i]->FillBuffer( 0.0 );
      }
    }

  // Do the joint intensity fusion
  this->m_JointIntensityFusionImage.clear();
  this->m_JointIntensityFusionImage.resize( this->m_NumberOfAtlasModalities );

  for( SizeValueType i = 0; i < this->m_NumberOfAtlasModalities; i++ )
    {
    this->m_JointIntensityFusionImage[i] = InputImageType::New();
    this->m_JointIntensityFusionImage[i]->CopyInformation( this->m_TargetImage[0] );
    this->m_JointIntensityFusionImage[i]->SetRegions( this->m_TargetImage[0]->GetRequestedRegion() );
    this->m_JointIntensityFusionImage[i]->SetLargestPossibleRegion( this->m_TargetImage[0]->GetLargestPossibleRegion() );
    this->m_JointIntensityFusionImage[i]->Allocate();
    this->m_JointIntensityFusionImage[i]->FillBuffer( 0.0 );
    }

  // Initialize the weight sum image
  this->m_WeightSumImage = ProbabilityImageType::New();
  this->m_WeightSumImage->CopyInformation( this->m_TargetImage[0] );
  this->m_WeightSumImage->SetRegions( this->m_TargetImage[0]->GetRequestedRegion() );
  this->m_WeightSumImage->SetLargestPossibleRegion( this->m_TargetImage[0]->GetLargestPossibleRegion() );
  this->m_WeightSumImage->Allocate();
  this->m_WeightSumImage->FillBuffer( 0.0 );

  // Initialize the count image
  this->m_CountImage = CountImageType::New();
  this->m_CountImage->CopyInformation( this->m_TargetImage[0] );
  this->m_CountImage->SetRegions( this->m_TargetImage[0]->GetRequestedRegion() );
  this->m_CountImage->SetLargestPossibleRegion( this->m_TargetImage[0]->GetLargestPossibleRegion() );
  this->m_CountImage->Allocate();
  this->m_CountImage->FillBuffer( 0 );

  // Determine the search offset list (or map if an search radius image is specified)

  typename InputImageType::SpacingType spacing = this->m_TargetImage[0]->GetSpacing();

  this->m_SearchNeighborhoodOffsetList.clear();
  this->m_SearchNeighborhoodOffsetSetsMap.clear();

  if( this->m_SearchNeighborhoodRadiusImage.IsNull() )
    {
    ConstNeighborhoodIterator<InputImageType> It( this->m_SearchNeighborhoodRadius,
      this->GetInput(), this->GetInput()->GetRequestedRegion() );

    SizeValueType searchNeighborhoodSize = ( It.GetNeighborhood() ).Size();

    DistanceIndexVectorType squaredDistances;
    squaredDistances.resize( searchNeighborhoodSize );

    for( unsigned int n = 0; n < searchNeighborhoodSize; n++ )
      {
      NeighborhoodOffsetType offset = ( It.GetNeighborhood() ).GetOffset( n );

      squaredDistances[n].first = n;
      squaredDistances[n].second = 0.0;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        squaredDistances[n].second += vnl_math_sqr( offset[d] * spacing[d] );
        }
      }
    std::sort( squaredDistances.begin(), squaredDistances.end(), DistanceIndexComparator() );

    for( unsigned int n = 0; n < searchNeighborhoodSize; n++ )
      {
      this->m_SearchNeighborhoodOffsetList.push_back( ( It.GetNeighborhood() ).GetOffset( squaredDistances[n].first ) );
      }
    }
  else
    {
    ImageRegionConstIterator<RadiusImageType> ItR( this->m_SearchNeighborhoodRadiusImage,
      this->m_SearchNeighborhoodRadiusImage->GetRequestedRegion() );

    for( ItR.GoToBegin(); !ItR.IsAtEnd(); ++ItR )
      {
      RadiusValueType localSearchRadius = ItR.Get();
      if( localSearchRadius > 0 &&
        this->m_SearchNeighborhoodOffsetSetsMap.find( localSearchRadius ) ==
          this->m_SearchNeighborhoodOffsetSetsMap.end() )
        {
        NeighborhoodRadiusType localSearchNeighborhoodRadius;
        localSearchNeighborhoodRadius.Fill( localSearchRadius );

        std::vector<NeighborhoodOffsetType> localSearchNeighborhoodOffsetList;

        ConstNeighborhoodIterator<InputImageType> It( localSearchNeighborhoodRadius,
          this->GetInput(), this->GetInput()->GetRequestedRegion() );

        RadiusValueType localSearchNeighborhoodSize = ( It.GetNeighborhood() ).Size();

        DistanceIndexVectorType squaredDistances;
        squaredDistances.resize( localSearchNeighborhoodSize );

        for( unsigned int n = 0; n < localSearchNeighborhoodSize; n++ )
          {
          NeighborhoodOffsetType offset = ( It.GetNeighborhood() ).GetOffset( n );

          squaredDistances[n].first = n;
          squaredDistances[n].second = 0.0;
          for( unsigned int d = 0; d < ImageDimension; d++ )
            {
            squaredDistances[n].second += vnl_math_sqr( offset[d] * spacing[d] );
            }
          }
        std::sort( squaredDistances.begin(), squaredDistances.end(), DistanceIndexComparator() );

        for( unsigned int n = 0; n < localSearchNeighborhoodSize; n++ )
          {
          localSearchNeighborhoodOffsetList.push_back( ( It.GetNeighborhood() ).GetOffset( squaredDistances[n].first ) );
          }
        this->m_SearchNeighborhoodOffsetSetsMap[localSearchRadius] = localSearchNeighborhoodOffsetList;
        }
      }
    }

  // Determine the patch offset list

  ConstNeighborhoodIterator<InputImageType> It2( this->m_PatchNeighborhoodRadius,
    this->GetInput(), this->GetInput()->GetRequestedRegion() );

  this->m_PatchNeighborhoodOffsetList.clear();

  this->m_PatchNeighborhoodSize = ( It2.GetNeighborhood() ).Size();
  for( unsigned int n = 0; n < this->m_PatchNeighborhoodSize; n++ )
    {
    this->m_PatchNeighborhoodOffsetList.push_back( ( It2.GetNeighborhood() ).GetOffset( n ) );
    }

  this->AllocateOutputs();
}

template <class TInputImage, class TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateData( const RegionType &region, ThreadIdType threadId )
{
  if( !this->m_IsWeightedAveragingComplete )
    {
    this->ThreadedGenerateDataForWeightedAveraging( region, threadId );
    }
  else
    {
    this->ThreadedGenerateDataForReconstruction( region, threadId );
    }
}

template <class TInputImage, class TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateDataForWeightedAveraging( const RegionType & region, ThreadIdType threadId )
{
  ProgressReporter progress( this, threadId, region.GetNumberOfPixels(), 100 );

  typename OutputImageType::Pointer output = this->GetOutput();

  SizeValueType numberOfTargetModalities = this->m_TargetImage.size();

  MatrixType absoluteAtlasPatchDifferences( this->m_NumberOfAtlases,
    this->m_PatchNeighborhoodSize * numberOfTargetModalities );

  MatrixType originalAtlasPatchIntensities( this->m_NumberOfAtlases,
    this->m_PatchNeighborhoodSize * this->m_NumberOfAtlasModalities );

  std::vector<SizeValueType> minimumAtlasOffsetIndices( this->m_NumberOfAtlases );

  bool useOnlyFirstAtlasImage = true;
  if( numberOfTargetModalities == this->m_NumberOfAtlasModalities )
    {
    useOnlyFirstAtlasImage = false;
    }

  // Iterate over the input region
  ConstNeighborhoodIteratorType ItN( this->m_PatchNeighborhoodRadius, this->m_TargetImage[0], region );
  for( ItN.GoToBegin(); !ItN.IsAtEnd(); ++ItN )
    {
    progress.CompletedPixel();

    IndexType currentCenterIndex = ItN.GetIndex();

    if( this->m_MaskImage &&
        this->m_MaskImage->GetPixel( currentCenterIndex ) == NumericTraits<LabelType>::ZeroValue() )
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
    std::vector<NeighborhoodOffsetType> searchNeighborhoodOffsetList;
    if( this->m_SearchNeighborhoodRadiusImage.IsNull() )
      {
      searchNeighborhoodOffsetList = this->m_SearchNeighborhoodOffsetList;
      }
    else
      {
      RadiusValueType localSearchRadius =
        this->m_SearchNeighborhoodRadiusImage->GetPixel( currentCenterIndex );
      if( localSearchRadius <= 0 )
        {
        continue;
        }
      searchNeighborhoodOffsetList = this->m_SearchNeighborhoodOffsetSetsMap[localSearchRadius];
      }
    SizeValueType searchNeighborhoodSize = searchNeighborhoodOffsetList.size();

    InputImagePixelVectorType normalizedTargetPatch =
      this->VectorizeImageListPatch( this->m_TargetImage, currentCenterIndex, true );

    absoluteAtlasPatchDifferences.fill( 0.0 );
    originalAtlasPatchIntensities.fill( 0.0 );

    // In each atlas, search for a patch that matches the target patch
    for( SizeValueType i = 0; i < this->m_NumberOfAtlases; i++ )
      {

      RealType minimumPatchSimilarity = NumericTraits<RealType>::max();
      SizeValueType minimumPatchOffsetIndex = 0;

      for( SizeValueType j = 0; j < searchNeighborhoodSize; j++ )
        {
        IndexType searchIndex = currentCenterIndex + searchNeighborhoodOffsetList[j];

        if( !output->GetRequestedRegion().IsInside( searchIndex ) )
          {
          continue;
          }

        RealType patchSimilarity = this->ComputeNeighborhoodPatchSimilarity(
          this->m_AtlasImages[i], searchIndex, normalizedTargetPatch, useOnlyFirstAtlasImage );

        if( patchSimilarity < minimumPatchSimilarity )
          {
          minimumPatchSimilarity = patchSimilarity;
          minimumPatchOffsetIndex = j;
          }
        }

      // Once the patch has been found, normalize it and then compute the absolute
      // difference with target patch

      IndexType minimumIndex = currentCenterIndex +
        searchNeighborhoodOffsetList[minimumPatchOffsetIndex];
      InputImagePixelVectorType normalizedMinimumAtlasPatch;
      if( numberOfTargetModalities == this->m_NumberOfAtlasModalities )
        {
        normalizedMinimumAtlasPatch =
          this->VectorizeImageListPatch( this->m_AtlasImages[i], minimumIndex, true );
        }
      else
        {
        normalizedMinimumAtlasPatch =
          this->VectorizeImagePatch( this->m_AtlasImages[i][0], minimumIndex, true );
        }

      typename InputImagePixelVectorType::const_iterator itA = normalizedMinimumAtlasPatch.begin();
      typename InputImagePixelVectorType::const_iterator itT = normalizedTargetPatch.begin();
      while( itA != normalizedMinimumAtlasPatch.end() )
        {
        RealType value = std::fabs( *itA - *itT );
        absoluteAtlasPatchDifferences(i, itA - normalizedMinimumAtlasPatch.begin()) = value;

        ++itA;
        ++itT;
        }

      InputImagePixelVectorType originalMinimumAtlasPatch =
        this->VectorizeImageListPatch( this->m_AtlasImages[i], minimumIndex, false );

      typename InputImagePixelVectorType::const_iterator itO = originalMinimumAtlasPatch.begin();
      while( itO != originalMinimumAtlasPatch.end() )
        {
        originalAtlasPatchIntensities(i, itO - originalMinimumAtlasPatch.begin()) = *itO;
        ++itO;
        }

      minimumAtlasOffsetIndices[i] = minimumPatchOffsetIndex;
      }

    // Allocate Mx
    MatrixType Mx( this->m_NumberOfAtlases, this->m_NumberOfAtlases );

    // Compute Mx values
    for( SizeValueType i = 0; i < this->m_NumberOfAtlases; i++ )
      {
      for( SizeValueType j = 0; j <= i; j++ )
        {
        RealType mxValue = 0.0;

        for( unsigned int k = 0; k < this->m_PatchNeighborhoodSize * numberOfTargetModalities; k++ )
          {
          mxValue += absoluteAtlasPatchDifferences[i][k] * absoluteAtlasPatchDifferences[j][k];
          }
        mxValue /= static_cast<RealType>( this->m_PatchNeighborhoodSize - 1 );

        if( this->m_Beta == 2.0 )
          {
          mxValue *= mxValue;
          }
        else
          {
          mxValue = std::pow( mxValue, this->m_Beta );
          }

        if( !std::isfinite( mxValue ) )
          {
          mxValue = 0.0;
          }

        Mx(i, j) = Mx(j, i) = mxValue;
        }
      }

    // Compute the weights by solving for the inverse of Mx
    MatrixType MxBar( this->m_NumberOfAtlases, this->m_NumberOfAtlases, 0.0 );
    MxBar.fill_diagonal( this->m_Alpha );
    MxBar += Mx;

    // Define a vector of all ones
    VectorType ones( this->m_NumberOfAtlases, 1.0 );
    VectorType W( this->m_NumberOfAtlases, 1.0 );

    if( this->m_ConstrainSolutionToNonnegativeWeights )
      {
      W = this->NonNegativeLeastSquares( MxBar, ones, 1e-6 );
      }
    else
      {
      W = vnl_svd<RealType>( MxBar ).solve( ones );

      for( SizeValueType i = 0; i < W.size(); i++ )
        {
        if( W[i] < 0.0 )
          {
          W[i] = 0.0;
          }
        }
      }

    // Normalize the weights
    W *= 1.0 / dot_product( W, ones );

    // Do joint intensity fusion
    VectorType estimatedNeighborhoodIntensities = W;

    estimatedNeighborhoodIntensities.post_multiply( originalAtlasPatchIntensities );

    for( SizeValueType i = 0; i < this->m_NumberOfAtlasModalities; i++ )
      {
      for( SizeValueType j = 0; j < this->m_PatchNeighborhoodSize; j++ )
        {
        IndexType neighborhoodIndex = ItN.GetIndex( j );

        if( !output->GetRequestedRegion().IsInside( neighborhoodIndex ) )
          {
          continue;
          }

        if( this->m_MaskImage &&
            this->m_MaskImage->GetPixel( neighborhoodIndex ) == NumericTraits<LabelType>::ZeroValue() )
          {
          continue;
          }

        RealType estimatedValue = (
          estimatedNeighborhoodIntensities[i * this->m_PatchNeighborhoodSize + j] +
          this->m_JointIntensityFusionImage[i]->GetPixel( neighborhoodIndex ) );

        if( !std::isfinite( estimatedValue ) )
          {
          estimatedValue = 0.0;
          }

        this->m_JointIntensityFusionImage[i]->SetPixel( neighborhoodIndex,
           static_cast<InputImagePixelType>( estimatedValue ) );
        if( i == 0 )
          {
          this->m_CountImage->SetPixel( neighborhoodIndex,
            this->m_CountImage->GetPixel( neighborhoodIndex ) + 1 );
          }
        }
      }

    if( this->m_NumberOfAtlasSegmentations > 0 )
      {
      // Perform voting using Hongzhi's averaging scheme. Iterate over all segmentation patches
      for( SizeValueType n = 0; n < this->m_PatchNeighborhoodSize; n++ )
        {
        IndexType neighborhoodIndex = ItN.GetIndex( n );
        if( !output->GetRequestedRegion().IsInside( neighborhoodIndex ) )
          {
          continue;
          }

        for( SizeValueType i = 0; i < this->m_NumberOfAtlasSegmentations; i++ )
          {
          // The segmentation at the corresponding patch location in atlas i
          IndexType minimumIndex = neighborhoodIndex +
            searchNeighborhoodOffsetList[minimumAtlasOffsetIndices[i]];

          if( !output->GetRequestedRegion().IsInside( minimumIndex ) )
            {
            continue;
            }

          LabelType label = this->m_AtlasSegmentations[i]->GetPixel( minimumIndex );

          if( this->m_LabelSet.find( label ) == this->m_LabelSet.end() )
            {
            continue;
            }

          // Add that weight the posterior map for voxel at idx
          this->m_LabelPosteriorProbabilityImages[label]->SetPixel( neighborhoodIndex,
            this->m_LabelPosteriorProbabilityImages[label]->GetPixel( neighborhoodIndex ) + W[i] );
          this->m_WeightSumImage->SetPixel( neighborhoodIndex,
            this->m_WeightSumImage->GetPixel( neighborhoodIndex ) + W[i] );

          if( this->m_RetainAtlasVotingWeightImages )
            {
            this->m_AtlasVotingWeightImages[i]->SetPixel( neighborhoodIndex,
              this->m_AtlasVotingWeightImages[i]->GetPixel( neighborhoodIndex ) + W[i] );
            }
          }
        }
      }
    }
}

template <class TInputImage, class TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateDataForReconstruction( const RegionType &region, ThreadIdType
  itkNotUsed( threadId ) )
{
  typename OutputImageType::Pointer output = this->GetOutput();

  // Perform voting at each voxel
  ImageRegionIteratorWithIndex<OutputImageType> It( output, region );

  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    IndexType index = It.GetIndex();

    if( this->m_MaskImage &&
        this->m_MaskImage->GetPixel( It.GetIndex() ) == NumericTraits<LabelType>::ZeroValue() )
      {
      continue;
      }

    RealType maxPosteriorProbability = 0.0;
    LabelType winningLabel = NumericTraits<LabelType>::ZeroValue();

    typename LabelSetType::const_iterator labelIt;
    for( labelIt = this->m_LabelSet.begin(); labelIt != this->m_LabelSet.end(); ++labelIt )
      {
      // check if the label is excluded
      typename LabelExclusionMap::const_iterator xIt = this->m_LabelExclusionImages.find( *labelIt );
      bool isLabelExcluded = ( xIt != m_LabelExclusionImages.end() && xIt->second->GetPixel( index ) != 0 );

      if( !isLabelExcluded )
        {
        typename ProbabilityImageType::PixelType posteriorProbability =
          this->m_LabelPosteriorProbabilityImages[*labelIt]->GetPixel( index );

        // Vote!
        if( maxPosteriorProbability < posteriorProbability )
          {
          maxPosteriorProbability = posteriorProbability;
          winningLabel = *labelIt;
          }
        }
      }
    It.Set( winningLabel );
    }

  ImageRegionIteratorWithIndex<ProbabilityImageType> ItW( this->m_WeightSumImage,
    region );

  for( ItW.GoToBegin(); !ItW.IsAtEnd(); ++ItW )
    {
    typename ProbabilityImageType::PixelType weightSum = ItW.Get();

    IndexType index = ItW.GetIndex();

    if( weightSum < 0.1 )
      {
      continue;
      }

    if( this->m_RetainLabelPosteriorProbabilityImages )
      {
      typename LabelSetType::const_iterator labelIt;
      for( labelIt = this->m_LabelSet.begin(); labelIt != this->m_LabelSet.end(); ++labelIt )
        {
        typename ProbabilityImageType::PixelType labelProbability =
          this->m_LabelPosteriorProbabilityImages[*labelIt]->GetPixel( index );
        this->m_LabelPosteriorProbabilityImages[*labelIt]->SetPixel( index, labelProbability / weightSum );
        }
      }

    if( this->m_RetainAtlasVotingWeightImages )
      {
      for( SizeValueType i = 0; i < this->m_NumberOfAtlases; i++ )
        {
        typename ProbabilityImageType::PixelType votingWeight =
          this->m_AtlasVotingWeightImages[i]->GetPixel( index );
        this->m_AtlasVotingWeightImages[i]->SetPixel( index, votingWeight / weightSum );
        }
      }
    }
}

template <class TInputImage, class TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::AfterThreadedGenerateData()
{
  // Clear posterior maps if not kept
  if( !this->m_RetainLabelPosteriorProbabilityImages )
    {
    this->m_LabelPosteriorProbabilityImages.clear();
    }

  // Normalize the joint intensity fusion images.
  for( SizeValueType i = 0; i < this->m_NumberOfAtlasModalities; i++ )
    {
    ImageRegionIterator<InputImageType> ItJ( this->m_JointIntensityFusionImage[i],
      this->m_JointIntensityFusionImage[i]->GetRequestedRegion() );
    ImageRegionIterator<CountImageType> ItC( this->m_CountImage,
      this->m_CountImage->GetRequestedRegion() );

    for( ItJ.GoToBegin(), ItC.GoToBegin(); !ItJ.IsAtEnd(); ++ItJ, ++ItC )
      {
      typename CountImageType::PixelType count = ItC.Get();

      if( count > 0 )
        {
        ItJ.Set( ItJ.Get() / static_cast<RealType>( count ) );
        }
      }
    }
}

template <class TInputImage, class TOutputImage>
typename WeightedVotingFusionImageFilter<TInputImage, TOutputImage>::InputImagePixelVectorType
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::VectorizeImageListPatch( const InputImageList &imageList, const IndexType index, const bool normalize )
{
  InputImagePixelVectorType patchVector( this->m_PatchNeighborhoodSize * imageList.size() );
  for( unsigned int i = 0; i < imageList.size(); i++ )
    {
    InputImagePixelVectorType patchVectorPerModality = this->VectorizeImagePatch( imageList[i], index, normalize );
    for( unsigned int j = 0; j < this->m_PatchNeighborhoodSize; j++ )
      {
      patchVector[i * this->m_PatchNeighborhoodSize + j] = patchVectorPerModality[j];
      }
    }
  return patchVector;
}

template <class TInputImage, class TOutputImage>
typename WeightedVotingFusionImageFilter<TInputImage, TOutputImage>::InputImagePixelVectorType
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::VectorizeImagePatch( const InputImagePointer image, const IndexType index, const bool normalize )
{
  InputImagePixelVectorType patchVector( this->m_PatchNeighborhoodSize );
  for( SizeValueType i = 0; i < this->m_PatchNeighborhoodSize; i++ )
    {
    IndexType neighborhoodIndex = index + this->m_PatchNeighborhoodOffsetList[i];

    bool isInBounds = this->m_TargetImageRequestedRegion.IsInside( neighborhoodIndex );
    if( isInBounds )
      {
      InputImagePixelType pixel = image->GetPixel( neighborhoodIndex );
      patchVector[i] = pixel;
      }
    else
      {
      patchVector[i] = std::numeric_limits<RealType>::quiet_NaN();
      }
    }

  if( normalize )
    {
    RealType mean = 0.0;
    RealType standardDeviation = 0.0;
    this->GetMeanAndStandardDeviationOfVectorizedImagePatch( patchVector, mean, standardDeviation );

    standardDeviation = std::max( standardDeviation, 1.0 );

    typename InputImagePixelVectorType::iterator it;
    for( it = patchVector.begin(); it != patchVector.end(); ++it )
      {
      *it = ( *it - mean ) / standardDeviation;
      }
    }
  return patchVector;
}

template <class TInputImage, class TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::GetMeanAndStandardDeviationOfVectorizedImagePatch(
  const InputImagePixelVectorType &patchVector, RealType &mean, RealType &standardDeviation )
{
  RealType sum = 0.0;
  RealType sumOfSquares = 0.0;
  RealType count = 0.0;

  typename InputImagePixelVectorType::const_iterator it;
  for( it = patchVector.begin(); it != patchVector.end(); ++it )
    {
    if( std::isfinite( *it ) )
      {
      sum += *it;
      sumOfSquares += vnl_math_sqr( *it );
      count += 1.0;
      }
    }

  mean = sum / count;
  standardDeviation = std::sqrt( ( sumOfSquares - count * vnl_math_sqr( mean ) ) / ( count - 1.0 ) );
}

template <class TInputImage, class TOutputImage>
typename WeightedVotingFusionImageFilter<TInputImage, TOutputImage>::RealType
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::ComputeNeighborhoodPatchSimilarity( const InputImageList &imageList, const IndexType index,
  const InputImagePixelVectorType &normalizedPatchVectorY, const bool useOnlyFirstImage )
{
  unsigned int numberOfImagesToUse = imageList.size();
  if( useOnlyFirstImage )
    {
    numberOfImagesToUse = 1;
    }

  RealType sumX = 0.0;
  RealType sumOfSquaresX = 0.0;
  RealType sumOfSquaredDifferencesXY = 0.0;
  RealType sumXY = 0.0;
  RealType N = 0.0;

  SizeValueType count = 0;
  for( SizeValueType i = 0; i < numberOfImagesToUse; i++ )
    {
    for( SizeValueType j = 0; j < this->m_PatchNeighborhoodSize; j++ )
      {
      IndexType neighborhoodIndex = index + this->m_PatchNeighborhoodOffsetList[j];

      bool isInBounds = this->m_TargetImageRequestedRegion.IsInside( neighborhoodIndex );
      if( isInBounds && std::isfinite( normalizedPatchVectorY[count] ) )
        {
        RealType x = static_cast<RealType>( imageList[i]->GetPixel( neighborhoodIndex ) );
        RealType y = static_cast<RealType>( normalizedPatchVectorY[count] );

        sumX += x;
        sumOfSquaresX += vnl_math_sqr( x );
        sumXY += ( x * y );

        sumOfSquaredDifferencesXY += vnl_math_sqr( y - x );
        N += 1.0;
        }
      ++count;
      }
    }

  // If we are on the boundary, a neighborhood patch might not overlap
  // with the image.  If we have 2 voxels or less for a neighborhood patch
  // we don't consider it to be a suitable match.
  if( N < 3.0 )
    {
    return NumericTraits<RealType>::max();
    }

  if( this->m_SimilarityMetric == PEARSON_CORRELATION )
    {
    RealType varianceX = sumOfSquaresX - vnl_math_sqr( sumX ) / N;
    varianceX = std::max( varianceX, 1.0e-6 );

    RealType measure = vnl_math_sqr( sumXY ) / varianceX;
    return ( sumXY > 0 ? -measure : measure );
    }
  else if( this->m_SimilarityMetric == MEAN_SQUARES )
    {
    return ( sumOfSquaredDifferencesXY / N );
    }
  else
    {
    itkExceptionMacro( "Unrecognized similarity metric." );
    }
}

template <class TInputImage, class TOutputImage>
typename WeightedVotingFusionImageFilter<TInputImage, TOutputImage>::VectorType
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::NonNegativeLeastSquares( const MatrixType A, const VectorType y, const RealType tolerance )
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

  VectorType P( n, 0 );
  VectorType R( n, 1 );
  VectorType x( n, 0 );
  VectorType s( n, 0 );
  VectorType w = A.transpose() * ( y - A * x );

  RealType wMaxValue = w.max_value();
  SizeValueType maxIndex = NumericTraits<SizeValueType>::max();
  wMaxValue = NumericTraits<RealType>::NonpositiveMin();
  for( SizeValueType i = 0; i < n; i++ )
    {
    if( R[i] == 1 && wMaxValue < w[i] )
      {
      maxIndex = i;
      wMaxValue = w[i];
      }
    }

  // Outer loop

  SizeValueType numberOfIterations = 0;
  while( R.sum() > 0 && wMaxValue > tolerance &&
    numberOfIterations++ < maximumNumberOfIterations )
    {
    P[maxIndex] = 1;
    R[maxIndex] = 0;

    SizeValueType sizeP = P.sum();

    MatrixType AP( m, sizeP, 0 );
    SizeValueType jIndex = 0;
    for( SizeValueType j = 0; j < n; j++ )
      {
      if( P[j] == 1 )
        {
        AP.set_column( jIndex++, A.get_column( j ) );
        }
      }

    VectorType sP = vnl_svd<RealType>( AP ).pinverse() * y;

    SizeValueType iIndex = 0;

    for( SizeValueType i = 0; i < n; i++ )
      {
      if( R[i] != 0 )
        {
        s[i] = 0;
        }
      else
        {
        s[i] = sP[iIndex++];
        }
      }

    // Inner loop
    while( sP.min_value() <= tolerance && sizeP > 0 )
      {
      RealType alpha = NumericTraits<RealType>::max();

      for( SizeValueType i = 0; i < n; i++ )
        {
        if( P[i] == 1 && s[i] <= tolerance )
          {
          RealType value = x[i] / ( x[i] - s[i] );
          if( value < alpha )
            {
            alpha = value;
            }
          }
        }

      x += alpha * ( s - x );

      for( SizeValueType i = 0; i < n; i++ )
        {
        if( P[i] == 1 && std::fabs( x[i] ) < tolerance )
          {
          P[i] = 0;
          R[i] = 1;
          }
        }

      sizeP = P.sum();
      if( sizeP == 0 )
        {
        break;
        }

      AP.set_size( m, sizeP );
      jIndex = 0;
      for( SizeValueType j = 0; j < n; j++ )
        {
        if( P[j] == 1 )
          {
          AP.set_column( jIndex++, A.get_column( j ) );
          }
        }

      sP = vnl_svd<RealType>( AP ).pinverse() * y;

      iIndex = 0;
      for( SizeValueType i = 0; i < n; i++ )
        {
        if( R[i] != 0 )
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
    w = A.transpose() * ( y - A * x );

    maxIndex = NumericTraits<SizeValueType>::max();
    wMaxValue = NumericTraits<RealType>::NonpositiveMin();
    for( SizeValueType i = 0; i < n; i++ )
      {
      if( R[i] == 1 && wMaxValue < w[i] )
        {
        maxIndex = i;
        wMaxValue = w[i];
        }
      }
    }

  return x;
}

template <class TInputImage, class TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::PrintSelf( std::ostream &os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << "Number of atlases = " << this->m_NumberOfAtlases << std::endl;
  os << "Number of atlas segmentations = " << this->m_NumberOfAtlasSegmentations << std::endl;
  os << "Number of atlas modalities = " << this->m_NumberOfAtlasModalities << std::endl;
  os << "Alpha = " << this->m_Alpha << std::endl;
  os << "Beta = " << this->m_Beta << std::endl;
  os << "Search neighborhood radius = " << this->m_SearchNeighborhoodRadius << std::endl;
  os << "Patch neighborhood radius = " << this->m_PatchNeighborhoodRadius << std::endl;
  if( this->m_SimilarityMetric == PEARSON_CORRELATION )
    {
    os << "Using Pearson correlation to measure the patch similarity." << std::endl;
    }
  else if( this->m_SimilarityMetric == MEAN_SQUARES )
    {
    os << "Using mean squares to measure the patch similarity." << std::endl;
    }
  if( this->m_ConstrainSolutionToNonnegativeWeights )
    {
    os << "Constrain solution to positive weights using NNLS." << std::endl;
    }

  os << "Label set: ";
  typename LabelSetType::const_iterator labelIt;
  for( labelIt = this->m_LabelSet.begin(); labelIt != this->m_LabelSet.end(); ++labelIt )
    {
    os << *labelIt << " ";
    }
  os << std::endl;
}

} // end namespace itk

#endif
