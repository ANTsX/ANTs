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

#include <vnl/algo/vnl_svd.h>

namespace itk {

template<typename TInputImage, typename TOutputImage>
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::WeightedVotingFusionImageFilter() :
  m_NumberOfAtlases( 0 ),
  m_NumberOfAtlasSegmentations( 0 ),
  m_NumberOfModalities( 0 ),
  m_PatchNeighborhoodSize( 0 ),
  m_Alpha( 0.1 ),
  m_Beta( 2.0 ),
  m_RetainLabelPosteriorProbabilityImages( false ),
  m_RetainAtlasVotingWeightImages( false )
{
  this->m_MaskImage = ITK_NULLPTR;
  this->m_MaskLabel = NumericTraits<LabelType>::OneValue();
}

template <class TInputImage, class TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::UpdateInputs()
{
  // Set all the inputs

  this->SetNumberOfIndexedInputs( ( this->m_NumberOfAtlases  + 1 ) * this->m_NumberOfModalities +
    this->m_NumberOfAtlasSegmentations + this->m_LabelExclusionImages.size() );

  SizeValueType nthInput = 0;

  for( SizeValueType i = 0; i < this->m_NumberOfModalities; i++ )
    {
    this->SetNthInput( nthInput++, this->m_TargetImage[i] );
    }

  for( SizeValueType i = 0; i < this->m_NumberOfAtlases; i++ )
    {
    for( SizeValueType j = 0; j < this->m_NumberOfModalities; j++ )
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
  outRegion.PadByRadius( this->m_SearchNeighborhoodRadius );
  outRegion.PadByRadius( this->m_PatchNeighborhoodRadius );

  // Iterate over all the inputs to this filter

  for( SizeValueType i = 0; i < this->m_NumberOfModalities; i++ )
    {
    InputImageType *input = this->m_TargetImage[i];
    RegionType region = outRegion;
    region.Crop( input->GetLargestPossibleRegion() );
    input->SetRequestedRegion( region );
    }

  for( SizeValueType i = 0; i < this->m_NumberOfAtlases; i++ )
    {
    for( SizeValueType j = 0; j < this->m_NumberOfModalities; j++ )
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
::BeforeThreadedGenerateData()
{
  if( this->m_NumberOfAtlasSegmentations != this->m_NumberOfAtlases )
    {
    // Set the number of atlas segmentations to 0 since we're just going to
    // doing joint intensity fusion
    this->m_NumberOfAtlasSegmentations = 0;
    }

  // Find all the unique labels in the atlas segmentations
  this->m_LabelSet.clear();
  for( unsigned int i = 0; i < this->m_NumberOfAtlasSegmentations; i++ )
    {
    ImageRegionConstIteratorWithIndex<LabelImageType> It(
      this->m_AtlasSegmentations[i], this->m_AtlasSegmentations[i]->GetRequestedRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      if( !this->m_MaskImage || this->m_MaskImage->GetPixel( It.GetIndex() ) == this->m_MaskLabel )
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
  this->m_JointIntensityFusionImage.resize( this->m_NumberOfModalities );

  for( SizeValueType i = 0; i < this->m_NumberOfModalities; i++ )
    {
    this->m_JointIntensityFusionImage[i] = InputImageType::New();
    this->m_JointIntensityFusionImage[i]->CopyInformation( this->m_TargetImage[i] );
    this->m_JointIntensityFusionImage[i]->SetRegions( this->m_TargetImage[i]->GetRequestedRegion() );
    this->m_JointIntensityFusionImage[i]->SetLargestPossibleRegion( this->m_TargetImage[i]->GetLargestPossibleRegion() );
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

  // Determine the offset list

  ConstNeighborhoodIterator<InputImageType> It( this->m_SearchNeighborhoodRadius,
    this->GetInput(), this->GetInput()->GetRequestedRegion() );

  this->m_SearchNeighborhoodOffsetList.clear();

  this->m_SearchNeighborhoodSize = ( It.GetNeighborhood() ).Size();
  for( unsigned int n = 0; n < this->m_SearchNeighborhoodSize; n++ )
    {
    this->m_SearchNeighborhoodOffsetList.push_back( ( It.GetNeighborhood() ).GetOffset( n ) );
    }

  // Determine neighborhood sizes

  this->m_PatchNeighborhoodSize = 1;
  this->m_SearchNeighborhoodSize = 1;
  for( SizeValueType d = 0; d < ImageDimension; d++ )
    {
    this->m_PatchNeighborhoodSize *= ( 2 * this->m_PatchNeighborhoodRadius[d] + 1 );
    this->m_SearchNeighborhoodSize *= ( 2 * this->m_SearchNeighborhoodRadius[d] + 1 );
    }

  this->AllocateOutputs();
}

template <class TInputImage, class TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateData( const RegionType & region, ThreadIdType threadId )
{
  ProgressReporter progress( this, threadId, region.GetNumberOfPixels(), 100 );

  typename OutputImageType::Pointer output = this->GetOutput();

  MatrixType absoluteAtlasPatchDifferences( this->m_NumberOfAtlases,
    this->m_PatchNeighborhoodSize * this->m_NumberOfModalities );

  MatrixType originalAtlasPatchIntensities( this->m_NumberOfAtlases,
    this->m_PatchNeighborhoodSize * this->m_NumberOfModalities );

  std::vector<SizeValueType> minimumAtlasOffsetIndices( this->m_NumberOfAtlases );

  // Iterate over the input region
  ConstNeighborhoodIteratorType ItN( this->m_PatchNeighborhoodRadius, this->m_TargetImage[0], region );
  for( ItN.GoToBegin(); !ItN.IsAtEnd(); ++ItN )
    {
    if( this->m_MaskImage && this->m_MaskImage->GetPixel( ItN.GetIndex() ) != this->m_MaskLabel )
      {
      progress.CompletedPixel();
      continue;
      }

    IndexType currentCenterIndex = ItN.GetIndex();

    InputImagePixelVectorType normalizedTargetPatch =
      this->VectorizeImageListPatch( this->m_TargetImage, currentCenterIndex, true );

    absoluteAtlasPatchDifferences.fill( 0.0 );
    originalAtlasPatchIntensities.fill( 0.0 );

    // In each atlas, search for a patch that matches the target patch
    for( SizeValueType i = 0; i < this->m_NumberOfAtlases; i++ )
      {

      RealType minimumPatchSimilarity = NumericTraits<RealType>::max();
      SizeValueType minimumPatchOffsetIndex = 0;

      for( SizeValueType j = 0; j < this->m_SearchNeighborhoodSize; j++ )
        {
        IndexType searchIndex = currentCenterIndex + this->m_SearchNeighborhoodOffsetList[j];

        if( !output->GetRequestedRegion().IsInside( searchIndex ) )
          {
          continue;
          }

        InputImagePixelVectorType individualAtlasPatch =
          this->VectorizeImageListPatch( this->m_AtlasImages[i], searchIndex, false );

        RealType patchSimilarity = this->ComputePatchSimilarity( individualAtlasPatch, normalizedTargetPatch  );

        if( patchSimilarity < minimumPatchSimilarity )
          {
          minimumPatchSimilarity = patchSimilarity;
          minimumPatchOffsetIndex = j;
          }
        }

      // Once the patch has been found, normalize it and then compute the absolute
      // difference with target patch
      IndexType minimumIndex = currentCenterIndex +
        this->m_SearchNeighborhoodOffsetList[minimumPatchOffsetIndex];
      InputImagePixelVectorType normalizedMinimumAtlasPatch =
        this->VectorizeImageListPatch( this->m_AtlasImages[i], minimumIndex, true );
      InputImagePixelVectorType originalMinimumAtlasPatch =
        this->VectorizeImageListPatch( this->m_AtlasImages[i], minimumIndex, false );

      typename InputImagePixelVectorType::const_iterator itA = normalizedMinimumAtlasPatch.begin();
      typename InputImagePixelVectorType::const_iterator itT = normalizedTargetPatch.begin();
      typename InputImagePixelVectorType::const_iterator itO = originalMinimumAtlasPatch.begin();
      while( itA != normalizedMinimumAtlasPatch.end() )
        {
        RealType value = std::fabs( *itA - *itT );
        absoluteAtlasPatchDifferences(i, itA - normalizedMinimumAtlasPatch.begin()) = value;
        originalAtlasPatchIntensities(i, itO - originalMinimumAtlasPatch.begin()) = *itO;

        ++itA;
        ++itT;
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

        for( unsigned int k = 0; k < this->m_PatchNeighborhoodSize * this->m_NumberOfModalities; k++ )
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

        if( vnl_math_isnan( mxValue ) || vnl_math_isinf( mxValue ) )
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
    VectorType W = vnl_svd<RealType>( MxBar ).solve( ones );

    // Normalize the weights
    W *= 1.0 / dot_product( W, ones );

    // Do joint intensity fusion
    VectorType estimatedNeighborhoodIntensities = W;

    for( SizeValueType i = 0; i < estimatedNeighborhoodIntensities.size(); i++ )
      {
      if( estimatedNeighborhoodIntensities[i] < 0.0 )
        {
        estimatedNeighborhoodIntensities[i] = 0.0;
        }
      }

    if( estimatedNeighborhoodIntensities.one_norm() > 0 )
      {
      estimatedNeighborhoodIntensities /= estimatedNeighborhoodIntensities.one_norm();
      }

    estimatedNeighborhoodIntensities.post_multiply( originalAtlasPatchIntensities );

    for( SizeValueType i = 0; i < this->m_NumberOfModalities; i++ )
      {
      for( SizeValueType j = 0; j < this->m_PatchNeighborhoodSize; j++ )
        {
        IndexType neighborhoodIndex = ItN.GetIndex( j );

        if( !output->GetRequestedRegion().IsInside( neighborhoodIndex ) )
          {
          continue;
          }

        if( this->m_MaskImage && this->m_MaskImage->GetPixel( neighborhoodIndex ) != this->m_MaskLabel )
          {
          continue;
          }

        RealType estimatedValue = (
          estimatedNeighborhoodIntensities[i * this->m_PatchNeighborhoodSize + j] +
          this->m_JointIntensityFusionImage[i]->GetPixel( neighborhoodIndex ) );

        if( vnl_math_isnan( estimatedValue ) || vnl_math_isinf( estimatedValue ) )
          {
          estimatedValue = 0.0;
          }

        this->m_JointIntensityFusionImage[i]->SetPixel( neighborhoodIndex,
           static_cast<InputImagePixelType>( estimatedValue ) );
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
            this->m_SearchNeighborhoodOffsetList[minimumAtlasOffsetIndices[i]];

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
    progress.CompletedPixel();
    }
}

template <class TInputImage, class TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::AfterThreadedGenerateData()
{
  typename OutputImageType::Pointer output = this->GetOutput();

  // Perform voting at each voxel
  ImageRegionIteratorWithIndex<OutputImageType> It( output, output->GetBufferedRegion() );

  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    IndexType index = It.GetIndex();

    if( this->m_MaskImage && this->m_MaskImage->GetPixel( It.GetIndex() ) != this->m_MaskLabel )
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

  // Clear posterior maps if not kept
  if( !this->m_RetainLabelPosteriorProbabilityImages )
    {
    this->m_LabelPosteriorProbabilityImages.clear();
    }

  ImageRegionIteratorWithIndex<ProbabilityImageType> ItW( this->m_WeightSumImage,
    this->m_WeightSumImage->GetBufferedRegion() );

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
typename WeightedVotingFusionImageFilter<TInputImage, TOutputImage>::InputImagePixelVectorType
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::VectorizeImageListPatch( const InputImageList &imageList, const IndexType index, const bool normalize )
{
  if( imageList.size() != this->m_NumberOfModalities )
    {
    itkExceptionMacro( "Input image list size does not equal the number of modalities." );
    }

  InputImagePixelVectorType patchVector( this->m_PatchNeighborhoodSize * this->m_NumberOfModalities );
  for( unsigned int i = 0; i < this->m_NumberOfModalities; i++ )
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
  ConstNeighborhoodIteratorType It( this->m_PatchNeighborhoodRadius, image, image->GetRequestedRegion() );
  It.SetLocation( index );

  InputImagePixelVectorType patchVector( this->m_PatchNeighborhoodSize );
  for( SizeValueType i = 0; i < this->m_PatchNeighborhoodSize; i++ )
    {
    bool isInBounds;
    InputImagePixelType pixel = It.GetPixel( i, isInBounds );
    if( isInBounds )
      {
      patchVector[i] = pixel;
      }
    }

  if( normalize )
    {
    RealType mean = 0.0;
    RealType standardDeviation = 0.0;
    this->GetMeanAndStandardDeviationOfVectorizedImagePatch( patchVector, mean, standardDeviation );

    standardDeviation = vnl_math_max( standardDeviation, 1.0 );

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
    sum += *it;
    sumOfSquares += vnl_math_sqr( *it );
    count += 1.0;
    }

  mean = sum / count;
  standardDeviation = std::sqrt( ( sumOfSquares - count * vnl_math_sqr( mean ) ) / ( count - 1.0 ) );
}

template <class TInputImage, class TOutputImage>
typename WeightedVotingFusionImageFilter<TInputImage, TOutputImage>::RealType
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::ComputePatchSimilarity( const InputImagePixelVectorType &patchVectorX,
  const InputImagePixelVectorType &normalizedPatchVectorY )
{
  RealType sumX = 0.0;
  RealType sumOfSquaresX = 0.0;
  RealType sumXY = 0.0;
  for( SizeValueType i = 0; i < patchVectorX.size(); i++ )
    {
    RealType x = static_cast<RealType>( patchVectorX[i] );
    RealType y = static_cast<RealType>( normalizedPatchVectorY[i] );
    sumX += x;
    sumOfSquaresX += vnl_math_sqr( x );
    sumXY += ( x * y );
    }

  RealType varianceX =
    sumOfSquaresX - vnl_math_sqr( sumX ) / static_cast<RealType>( patchVectorX.size() );
  varianceX = vnl_math_max( varianceX, 1.0 );

  return ( sumXY > 0 ? -vnl_math_sqr( sumXY ) / varianceX : vnl_math_sqr( sumXY ) / varianceX );
}

template <class TInputImage, class TOutputImage>
void
WeightedVotingFusionImageFilter<TInputImage, TOutputImage>
::PrintSelf( std::ostream &os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << "Number of atlases = " << this->m_NumberOfAtlases << std::endl;
  os << "Number of atlas segmentations = " << this->m_NumberOfAtlasSegmentations << std::endl;
  os << "Number of modalities = " << this->m_NumberOfModalities << std::endl;
  os << "Alpha = " << this->m_Alpha << std::endl;
  os << "Beta = " << this->m_Beta << std::endl;
  os << "Search neighborhood radius = " << this->m_SearchNeighborhoodRadius << std::endl;
  os << "Patch neighborhood radius = " << this->m_PatchNeighborhoodRadius << std::endl;

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
