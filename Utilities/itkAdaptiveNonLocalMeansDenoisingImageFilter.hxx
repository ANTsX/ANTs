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
#ifndef itkAdaptiveNonLocalMeansDenoisingImageFilter_hxx
#define itkAdaptiveNonLocalMeansDenoisingImageFilter_hxx

#include "itkAdaptiveNonLocalMeansDenoisingImageFilter.h"

#include "itkArray.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMeanImageFilter.h"
#include "itkNeighborhoodIterator.h"
#include "itkProgressReporter.h"
#include "itkStatisticsImageFilter.h"
#include "itkVarianceImageFilter.h"

namespace itk {

template <typename TInputImage, typename TOutputImage>
AdaptiveNonLocalMeansDenoisingImageFilter<TInputImage, TOutputImage>
::AdaptiveNonLocalMeansDenoisingImageFilter() :
  m_UseRicianNoiseModel( true ),
  m_Epsilon( 0.00001 ),
  m_MeanThreshold( 0.95 ),
  m_VarianceThreshold( 0.5 ),
  m_SmoothingFactor( 1.0 ),
  m_SmoothingVariance( 2.0 ),
  m_MaximumInputPixelIntensity( NumericTraits<RealType>::NonpositiveMin() ),
  m_MinimumInputPixelIntensity( NumericTraits<RealType>::max() )
{
  this->SetNumberOfRequiredInputs( 1 );

  this->m_MeanImage = ITK_NULLPTR;
  this->m_VarianceImage = ITK_NULLPTR;
  this->m_IntensitySquaredDistanceImage = ITK_NULLPTR;
  this->m_ThreadContributionCountImage = ITK_NULLPTR;

  this->m_RicianBiasImage = ITK_NULLPTR;

  this->m_BlockNeighborhoodRadius.Fill( 1 );
  this->m_NeighborhoodRadiusForLocalMeanAndVariance.Fill( 1 );
  this->m_NeighborhoodRadius2.Fill( 3 );
}

template<typename TInputImage, typename TOutputImage>
void
AdaptiveNonLocalMeansDenoisingImageFilter<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
  const InputImageType *inputImage = this->GetInput();

  typedef MeanImageFilter<InputImageType, RealImageType> MeanImageFilterType;
  typename MeanImageFilterType::Pointer meanImageFilter = MeanImageFilterType::New();
  meanImageFilter->SetInput( inputImage );
  meanImageFilter->SetRadius( this->GetNeighborhoodRadiusForLocalMeanAndVariance() );

  this->m_MeanImage = meanImageFilter->GetOutput();
  this->m_MeanImage->Update();
  this->m_MeanImage->DisconnectPipeline();

  typedef VarianceImageFilter<InputImageType, RealImageType> VarianceImageFilterType;
  typename VarianceImageFilterType::Pointer varianceImageFilter = VarianceImageFilterType::New();
  varianceImageFilter->SetInput( inputImage );
  varianceImageFilter->SetRadius( this->GetNeighborhoodRadiusForLocalMeanAndVariance() );

  this->m_VarianceImage = varianceImageFilter->GetOutput();
  this->m_VarianceImage->Update();
  this->m_VarianceImage->DisconnectPipeline();

  typedef StatisticsImageFilter<InputImageType> StatsFilterType;
  typename StatsFilterType::Pointer statsFilter = StatsFilterType::New();
  statsFilter->SetInput( inputImage );
  statsFilter->Update();

  this->m_MaximumInputPixelIntensity = static_cast<RealType>( statsFilter->GetMaximum() );
  this->m_MinimumInputPixelIntensity = static_cast<RealType>( statsFilter->GetMinimum() );

  this->m_ThreadContributionCountImage = RealImageType::New();
  this->m_ThreadContributionCountImage->CopyInformation( inputImage );
  this->m_ThreadContributionCountImage->SetRegions( inputImage->GetRequestedRegion() );
  this->m_ThreadContributionCountImage->Allocate();
  this->m_ThreadContributionCountImage->FillBuffer( 0.0 );

  if( this->m_UseRicianNoiseModel )
    {
    this->m_RicianBiasImage = RealImageType::New();
    this->m_RicianBiasImage->CopyInformation( inputImage );
    this->m_RicianBiasImage->SetRegions( inputImage->GetRequestedRegion() );
    this->m_RicianBiasImage->Allocate();
    this->m_RicianBiasImage->FillBuffer( 0.0 );
    }

  ConstNeighborhoodIterator<InputImageType> ItBI( this->m_BlockNeighborhoodRadius,
    this->GetInput(), this->GetInput()->GetRequestedRegion() );

  this->m_NeighborhoodOffsetList.clear();

  unsigned int neighborhoodBlockSize = ( ItBI.GetNeighborhood() ).Size();
  for( unsigned int n = 0; n < neighborhoodBlockSize; n++ )
    {
    this->m_NeighborhoodOffsetList.push_back( ( ItBI.GetNeighborhood() ).GetOffset( n ) );
    }

  this->AllocateOutputs();
}

template<typename TInputImage, typename TOutputImage>
void
AdaptiveNonLocalMeansDenoisingImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateData( const RegionType &region, ThreadIdType threadId )
{
  ProgressReporter progress( this, threadId, region.GetNumberOfPixels(), 100 );

  const InputImageType *inputImage = this->GetInput();

  OutputImageType *outputImage = this->GetOutput();

  ConstNeighborhoodIterator<RealImageType> ItV( this->m_NeighborhoodRadius2, this->m_VarianceImage, region );
  ConstNeighborhoodIterator<RealImageType> ItM( this->m_NeighborhoodRadius2, this->m_MeanImage, region );

  ConstNeighborhoodIterator<InputImageType> ItBI( this->m_BlockNeighborhoodRadius, inputImage, region );
  ConstNeighborhoodIterator<RealImageType> ItBM( this->m_BlockNeighborhoodRadius, this->m_MeanImage, region );

  NeighborhoodIterator<InputImageType> ItBO( this->m_BlockNeighborhoodRadius, outputImage, region );
  NeighborhoodIterator<RealImageType> ItBL( this->m_BlockNeighborhoodRadius, this->m_ThreadContributionCountImage, region );

  unsigned int neighborhoodSize = ( ItM.GetNeighborhood() ).Size();
  unsigned int neighborhoodBlockSize = ( ItBM.GetNeighborhood() ).Size();

  Array<RealType> weightedAverageIntensities( neighborhoodBlockSize );

  ItM.GoToBegin();
  ItV.GoToBegin();
  ItBI.GoToBegin();
  ItBM.GoToBegin();
  ItBO.GoToBegin();
  ItBL.GoToBegin();

  while( !ItM.IsAtEnd() )
    {
    progress.CompletedPixel();

    RealType meanCenterPixel = ItM.GetCenterPixel();
    RealType varianceCenterPixel = ItV.GetCenterPixel();

    RealType maxWeight = NumericTraits<RealType>::ZeroValue();
    RealType sumOfWeights = NumericTraits<RealType>::ZeroValue();

    weightedAverageIntensities.Fill( NumericTraits<RealType>::ZeroValue() );

    RealType meanNeighborhoodPixel = NumericTraits<RealType>::ZeroValue();
    RealType varianceNeighborhoodPixel = NumericTraits<RealType>::ZeroValue();

    if( meanCenterPixel > this->m_Epsilon && varianceCenterPixel > this->m_Epsilon )
      {
      RealType minimumDistance = NumericTraits<RealType>::max();
      for( unsigned int m = 0; m < neighborhoodSize; m++ )
        {
        if( ! ItM.IndexInBounds( m ) || m == static_cast<unsigned int>( 0.5 * neighborhoodSize ) )
          {
          continue;
          }

        meanNeighborhoodPixel = ItM.GetPixel( m );
        varianceNeighborhoodPixel = ItV.GetPixel( m );

        if( meanNeighborhoodPixel < this->m_Epsilon || varianceNeighborhoodPixel < this->m_Epsilon )
          {
          continue;
          }

        RealType meanRatio = meanCenterPixel / meanNeighborhoodPixel;
        RealType meanRatioInverse = ( this->m_MaximumInputPixelIntensity - meanCenterPixel ) /
          ( this->m_MaximumInputPixelIntensity - meanNeighborhoodPixel );

        RealType varianceRatio = varianceCenterPixel / varianceNeighborhoodPixel;

        if( ( ( meanRatio > this->m_MeanThreshold && meanRatio < 1.0 / this->m_MeanThreshold ) ||
            ( meanRatioInverse > this->m_MeanThreshold && meanRatioInverse < 1.0 / this->m_MeanThreshold ) ) &&
            varianceRatio > this->m_VarianceThreshold && varianceRatio < 1.0 / this->m_VarianceThreshold )
          {
          IndexType neighborhoodIndex = ItM.GetIndex( m );

          RealType averageDistance = 0.0;
          RealType count = 0.0;

          for( unsigned int n = 0; n < neighborhoodBlockSize; n++ )
            {
            IndexType neighborhoodBlockIndex = neighborhoodIndex + this->m_NeighborhoodOffsetList[n];

            if( ! inputImage->GetRequestedRegion().IsInside( neighborhoodBlockIndex ) )
              {
              continue;
              }
            averageDistance += vnl_math_sqr( ( ItBI.GetPixel( n ) - ItBM.GetPixel( n ) ) -
              ( inputImage->GetPixel( neighborhoodBlockIndex ) - this->m_MeanImage->GetPixel( neighborhoodBlockIndex ) ) );

            count += 1.0;
            }
          averageDistance /= count;
          minimumDistance = vnl_math_min( averageDistance, minimumDistance );
          }
        }

      if( this->m_UseRicianNoiseModel )
        {
        for( unsigned int n = 0; n < neighborhoodBlockSize; n++ )
          {
          if( ! ItBM.IndexInBounds( n ) )
            {
            continue;
            }
          if(  minimumDistance == NumericTraits<RealType>::max() )
            {
            this->m_RicianBiasImage->SetPixel( ItBM.GetIndex(), 0.0 );
            }
          else
            {
            this->m_RicianBiasImage->SetPixel( ItBM.GetIndex(), minimumDistance );
            }
          }
        }

      minimumDistance *= this->m_SmoothingFactor;
      if( minimumDistance == 0.0 )
        {
        minimumDistance = 0.00000000000001;
        }

      for( unsigned int m = 0; m < neighborhoodSize; m++ )
        {
        if( ! ItM.IndexInBounds( m ) || m == static_cast<unsigned int>( 0.5 * neighborhoodSize ) )
          {
          continue;
          }

        meanNeighborhoodPixel = ItM.GetPixel( m );
        varianceNeighborhoodPixel = ItV.GetPixel( m );

        if( meanNeighborhoodPixel < this->m_Epsilon || varianceNeighborhoodPixel < this->m_Epsilon )
          {
          continue;
          }

        RealType meanRatio = meanCenterPixel / meanNeighborhoodPixel;
        RealType meanRatioInverse = ( this->m_MaximumInputPixelIntensity - meanCenterPixel ) /
          ( this->m_MaximumInputPixelIntensity - meanNeighborhoodPixel );

        RealType varianceRatio = varianceCenterPixel / varianceNeighborhoodPixel;

        if( ( ( meanRatio > this->m_MeanThreshold && meanRatio < 1.0 / this->m_MeanThreshold ) ||
            ( meanRatioInverse > this->m_MeanThreshold && meanRatioInverse < 1.0 / this->m_MeanThreshold ) ) &&
            varianceRatio > this->m_VarianceThreshold && varianceRatio < 1.0 / this->m_VarianceThreshold )
          {
          IndexType neighborhoodIndex = ItM.GetIndex( m );

          RealType averageDistance = 0.0;
          RealType count = 0.0;
          for( unsigned int n = 0; n < neighborhoodBlockSize; n++ )
            {
            IndexType neighborhoodBlockIndex = neighborhoodIndex + this->m_NeighborhoodOffsetList[n];
            if( ! inputImage->GetRequestedRegion().IsInside( neighborhoodBlockIndex ) )
              {
              continue;
              }
            averageDistance += vnl_math_sqr( inputImage->GetPixel( neighborhoodBlockIndex ) - ItBI.GetPixel( n ) );
            count += 1.0;
            }
          averageDistance /= count;

          RealType weight = 0.0;
          if( averageDistance <= 3.0 * minimumDistance )
            {
            weight = std::exp( -averageDistance / minimumDistance );
            }
          if( weight > maxWeight )
            {
            maxWeight = weight;
            }

          if( weight > 0.0 && weight <= 1.0 )
            {
            for( unsigned int n = 0; n < neighborhoodBlockSize; n++ )
              {
              IndexType neighborhoodBlockIndex = neighborhoodIndex + this->m_NeighborhoodOffsetList[n];
              if( ! inputImage->GetRequestedRegion().IsInside( neighborhoodBlockIndex ) )
                {
                continue;
                }
              if( this->m_UseRicianNoiseModel )
                {
                weightedAverageIntensities[n] += weight * vnl_math_sqr( inputImage->GetPixel( neighborhoodBlockIndex ) );
                }
              else
                {
                weightedAverageIntensities[n] += weight * inputImage->GetPixel( neighborhoodBlockIndex );
                }
              }
            sumOfWeights += weight;
            }
          }
        }

      if( maxWeight == 0.0 )
        {
        maxWeight = 1.0;
        }

      for( unsigned int n = 0; n < neighborhoodBlockSize; n++ )
        {
        if( ! ItBM.IndexInBounds( n ) )
          {
          continue;
          }
        if( this->m_UseRicianNoiseModel )
          {
          weightedAverageIntensities[n] += maxWeight * vnl_math_sqr( ItBI.GetPixel( n ) );
          }
        else
          {
          weightedAverageIntensities[n] += maxWeight * ItBI.GetPixel( n );
          }
        }
      sumOfWeights += maxWeight;

      if( sumOfWeights > 0.0 )
        {
        for( unsigned int n = 0; n < neighborhoodBlockSize; n++ )
          {
          if( ! ItBO.IndexInBounds( n ) )
            {
            continue;
            }
          typename OutputImageType::PixelType estimate = ItBO.GetPixel( n );
          estimate += ( weightedAverageIntensities[n] / sumOfWeights );

          ItBO.SetPixel( n, estimate );
          }
        }
      }
    else
      {
      maxWeight = 1.0;

      for( unsigned int n = 0; n < neighborhoodBlockSize; n++ )
        {
        if( ! ItBM.IndexInBounds( n ) )
          {
          continue;
          }
        if( this->m_UseRicianNoiseModel )
          {
          weightedAverageIntensities[n] += maxWeight * vnl_math_sqr( ItBI.GetPixel( n ) );
          }
        else
          {
          weightedAverageIntensities[n] += maxWeight * ItBI.GetPixel( n );
          }
        }
      sumOfWeights += maxWeight;

      for( unsigned int n = 0; n < neighborhoodBlockSize; n++ )
        {
        if( ! ItBO.IndexInBounds( n ) )
          {
          continue;
          }
        typename OutputImageType::PixelType estimate = ItBO.GetPixel( n );
        estimate += ( weightedAverageIntensities[n] / sumOfWeights );

        ItBO.SetPixel( n, estimate );
        ItBL.SetPixel( n, ItBL.GetPixel( n ) + 1.0 );
        }
      }

    ++ItM;
    ++ItV;
    ++ItBI;
    ++ItBL;
    ++ItBM;
    ++ItBO;
    }
}

template<typename TInputImage, typename TOutputImage>
void
AdaptiveNonLocalMeansDenoisingImageFilter<TInputImage, TOutputImage>
::AfterThreadedGenerateData()
{
  if( this->m_UseRicianNoiseModel )
    {
    typedef DiscreteGaussianImageFilter<RealImageType, RealImageType> SmootherType;
    typename SmootherType::Pointer smoother = SmootherType::New();
    smoother->SetInput( this->m_RicianBiasImage );
    smoother->SetUseImageSpacingOn();
    smoother->Update();

    ImageRegionConstIterator<RealImageType> ItS( smoother->GetOutput(), smoother->GetOutput()->GetRequestedRegion() );
    ImageRegionConstIterator<RealImageType> ItM( this->m_MeanImage, this->m_MeanImage->GetRequestedRegion() );
    ImageRegionIterator<RealImageType> ItB( this->m_RicianBiasImage, this->m_RicianBiasImage->GetRequestedRegion() );

    for( ItB.GoToBegin(), ItS.GoToBegin(), ItM.GoToBegin(); !ItS.IsAtEnd(); ++ItS, ++ItB, ++ItM )
      {
      RealType snr = ItM.Get() / std::sqrt( ItS.Get() );

      RealType bias = 2.0 * ItS.Get() / this->CalculateCorrectionFactor( snr );

      if( vnl_math_isnan( bias ) || vnl_math_isinf( bias ) )
        {
        bias = 0.0;
        }
      ItB.Set( bias );
      }
    }

  ImageRegionIteratorWithIndex<OutputImageType> ItO( this->GetOutput(),
    this->GetOutput()->GetRequestedRegion() );
  ImageRegionConstIterator<RealImageType> ItL( this->m_ThreadContributionCountImage,
    this->m_ThreadContributionCountImage->GetRequestedRegion() );

  for( ItO.GoToBegin(), ItL.GoToBegin(); !ItO.IsAtEnd(); ++ItO, ++ItL )
    {
    RealType estimate = ItO.Get();

    if( ItL.Get() == 0.0 )
      {
      continue;
      }

//     estimate /= ItL.Get();

    if( this->m_UseRicianNoiseModel )
      {
      RealType bias = this->m_RicianBiasImage->GetPixel( ItO.GetIndex() );

      estimate -= bias;
      if( estimate < 0.0 )
        {
        estimate = 0.0;
        }
      estimate = std::sqrt( estimate );
      }

    ItO.Set( estimate );
    }
}

template<typename TInputImage, typename TOutputImage>
typename AdaptiveNonLocalMeansDenoisingImageFilter<TInputImage, TOutputImage>::RealType
AdaptiveNonLocalMeansDenoisingImageFilter<TInputImage, TOutputImage>
::CalculateCorrectionFactor( RealType snr )
{
   RealType snrSquared = vnl_math_sqr( snr );

   RealType value = 2.0 + snrSquared - 0.125 * vnl_math::pi * std::exp( -0.5 * snrSquared ) *
     vnl_math_sqr( ( 2.0 + snrSquared ) * this->m_ModifiedBesselCalculator.ModifiedBesselI0( 0.25 * snrSquared ) +
     snrSquared * this->m_ModifiedBesselCalculator.ModifiedBesselI1( 0.25 * snrSquared ) );

   if( value < 0.001 || value > 10.0 )
     {
     value = 1.0;
     }
   return value;
}

template<typename TInputImage, typename TOutputImage>
void
AdaptiveNonLocalMeansDenoisingImageFilter<TInputImage, TOutputImage>
::PrintSelf( std::ostream &os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( this->m_UseRicianNoiseModel )
    {
    os << indent << "Using Rician noise model." << std::endl;
    }
  else
    {
    os << indent << "Using Gaussian noise model." << std::endl;
    }

  os << indent << "Epsilon = " << this->m_Epsilon << std::endl;
  os << indent << "Mean threshold = " << this->m_MeanThreshold << std::endl;
  os << indent << "Variance threshold = " << this->m_VarianceThreshold << std::endl;
  os << indent << "Smoothing factor = " << this->m_SmoothingFactor << std::endl;
  os << indent << "Smoothing variance = " << this->m_SmoothingVariance << std::endl;

  os << indent << "Neighborhood radius for local mean and variance = " << this->m_NeighborhoodRadiusForLocalMeanAndVariance << std::endl;
  os << indent << "Neighborhood radius 2 = " << this->m_NeighborhoodRadius2 << std::endl;
  os << indent << "Block Neighborhood radius = " << this->m_BlockNeighborhoodRadius << std::endl;
}

} // end namespace itk

#endif
