/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: antsAtroposSegmentationImageFilter.txx,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
  http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt
  for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __antsAtroposSegmentationImageFilter_txx
#define __antsAtroposSegmentationImageFilter_txx

#include "antsAtroposSegmentationImageFilter.h"
#include "antsGaussianListSampleFunction.h"

#include "itkAddImageFilter.h"
#include "itkAddConstantToImageFilter.h"
#include "itkBinaryContourImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkDistanceToCentroidMembershipFunction.h"
#include "itkFastMarchingImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkIterationReporter.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkMinimumDecisionRule2.h"
#include "itkOtsuMultipleThresholdsCalculator.h"
#include "itkSampleClassifierFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkTimeProbe.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkWeightedCentroidKdTreeGenerator.h"

#include "vnl/vnl_vector.h"

namespace itk
{
namespace ants
{
template <class TInputImage, class TMaskImage, class TClassifiedImage>
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::AtroposSegmentationImageFilter()
{
  this->ProcessObject::SetNumberOfRequiredInputs( 1 );
  this->m_NumberOfIntensityImages = 1;

  this->m_NumberOfClasses = 3;
  this->m_MaximumNumberOfIterations = 5;
  this->m_ElapsedIterations = 0;
  this->m_ConvergenceThreshold = 0.001;

  this->m_MaskLabel = NumericTraits<LabelType>::One;

  this->m_InitializationStrategy = KMeans;

  this->m_PriorProbabilityWeight = 1.0;
  this->m_AdaptiveSmoothingWeights.clear();
  this->m_PriorLabelParameterMap.clear();
  this->m_PriorProbabilityThreshold = 0.0;
  this->m_PriorProbabilityImages.clear();
  this->m_PriorProbabilitySparseImages.clear();

  this->m_MRFSmoothingFactor = 0.3;
  this->m_MRFRadius.Fill( 1 );

  this->m_SplineOrder = 3;
  this->m_NumberOfLevels.Fill( 6 );
  this->m_NumberOfControlPoints.Fill( this->m_SplineOrder + 1 );

  this->m_MinimizeMemoryUsage = false;

  this->m_UseEuclideanDistanceForPriorLabels = false;
  this->m_PosteriorProbabilityImages.clear();
  this->m_DistancePriorProbabilityImages.clear();

  this->m_OutlierHandlingFilter = NULL;
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::~AtroposSegmentationImageFilter()
{
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::SetMaskImage( const MaskImageType * mask )
{
  this->SetNthInput( 1, const_cast<MaskImageType *>( mask ) );
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
const typename AtroposSegmentationImageFilter
<TInputImage, TMaskImage, TClassifiedImage>::MaskImageType
* AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GetMaskImage() const
  {
  const MaskImageType * maskImage =
    dynamic_cast<const MaskImageType *>( this->ProcessObject::GetInput( 1 ) );

  return maskImage;
  }

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::SetPriorLabelImage( const ClassifiedImageType * prior )
{
  this->m_InitializationStrategy = PriorLabelImage;
  this->SetNthInput( 2, const_cast<ClassifiedImageType *>( prior ) );
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
const typename AtroposSegmentationImageFilter
<TInputImage, TMaskImage, TClassifiedImage>::ClassifiedImageType
* AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GetPriorLabelImage() const
  {
  const ClassifiedImageType * prior =
    dynamic_cast<const ClassifiedImageType *>(
      this->ProcessObject::GetInput( 2 ) );

  return prior;
  }

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::SetPriorProbabilityImage(
  unsigned int whichClass, RealImageType * priorImage )
{
  if( whichClass < 1 || whichClass > this->m_NumberOfClasses )
    {
    itkExceptionMacro( "The requested prior probability image = "
                       << whichClass << " should be in the range [1, "
                       << this->m_NumberOfClasses << "]" );
    }
  if( this->m_MinimizeMemoryUsage )
    {
    // To make matters simpler, we force the index to be zero
    //   for each priorImage image.

    typename RealImageType::IndexType startIndex =
      priorImage->GetRequestedRegion().GetIndex();

    typename SparseImageType::Pointer sparsePriorImage =
      SparseImageType::New();
    sparsePriorImage->Initialize();

    unsigned long count = 0;

    ImageRegionConstIteratorWithIndex<RealImageType> It( priorImage,
                                                         priorImage->GetRequestedRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      if( It.Get() > this->m_PriorProbabilityThreshold )
        {
        typename RealImageType::IndexType index = It.GetIndex();
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          index[d] -= startIndex[d];
          }

        unsigned long number = this->IndexToNumber(
            index, priorImage->GetRequestedRegion().GetSize() );

        typename SparseImageType::PointType imageNumberIndex;
        imageNumberIndex[0] = number;

        sparsePriorImage->SetPoint( count, imageNumberIndex );
        sparsePriorImage->SetPointData( count, It.Get() );
        count++;
        }
      }
    if( this->m_PriorProbabilitySparseImages.size() < whichClass )
      {
      this->m_PriorProbabilitySparseImages.resize( whichClass );
      }
    this->m_PriorProbabilitySparseImages[whichClass - 1] = sparsePriorImage;
    }
  else
    {
    if( this->m_PriorProbabilityImages.size() < whichClass )
      {
      this->m_PriorProbabilityImages.resize( whichClass );
      }
    this->m_PriorProbabilityImages[whichClass - 1] = priorImage;
    }
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::RealImageType::Pointer
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GetPriorProbabilityImage( unsigned int whichClass ) const
{
  if( whichClass < 1 || whichClass > this->m_NumberOfClasses )
    {
    itkExceptionMacro( "The requested prior probability image = "
                       << whichClass << " should be in the range [1, "
                       << this->m_NumberOfClasses << "]" );
    }
  if( this->m_InitializationStrategy != PriorProbabilityImages )
    {
    return NULL;
    }

  if( ( this->m_MinimizeMemoryUsage &&
        this->m_PriorProbabilitySparseImages.size() != this->m_NumberOfClasses ) ||
      ( !this->m_MinimizeMemoryUsage &&
        this->m_PriorProbabilityImages.size() != this->m_NumberOfClasses ) )
    {
    itkExceptionMacro( "The number of prior probability images does not "
                       << "equal the number of classes." );
    }

  if( this->m_MinimizeMemoryUsage )
    {
    // To make matters simpler, we forced the prior probability images
    //   to all have indices of zeros.
    typename RealImageType::RegionType region;
    region.SetSize( this->GetIntensityImage( 0 )->GetRequestedRegion().GetSize() );
    typename RealImageType::RegionType::IndexType index;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      index[d] = 0;
      }
    region.SetIndex( index );

    typename RealImageType::Pointer priorImage = RealImageType::New();
    priorImage->SetRegions( region );
    priorImage->SetOrigin( this->GetIntensityImage( 0 )->GetOrigin() );
    priorImage->SetDirection( this->GetIntensityImage( 0 )->GetDirection() );
    priorImage->SetSpacing( this->GetIntensityImage( 0 )->GetSpacing() );
    priorImage->Allocate();
    priorImage->FillBuffer( 0 );

    typename SparseImageType::Pointer sparsePriorImage =
      this->m_PriorProbabilitySparseImages[whichClass - 1];

    typename SparseImageType::PointsContainer::ConstIterator It =
      sparsePriorImage->GetPoints()->Begin();
    typename SparseImageType::PointDataContainer::ConstIterator ItD =
      sparsePriorImage->GetPointData()->Begin();
    while( It != sparsePriorImage->GetPoints()->End() )
      {
      unsigned long number = static_cast<unsigned long>( It.Value()[0] );
      typename RealImageType::IndexType index =
        this->NumberToIndex( number, priorImage->GetRequestedRegion().GetSize() );
      priorImage->SetPixel( index, ItD.Value() );

      ++It;
      ++ItD;
      }

    return priorImage;
    }
  else
    {
    return this->m_PriorProbabilityImages[whichClass - 1];
    }
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::SetIntensityImage( unsigned int which, const ImageType * image )
{
  if( which > this->m_NumberOfIntensityImages - 1 )
    {
    this->m_NumberOfIntensityImages = which + 1;
    }
  this->SetNthInput( 2 + this->m_NumberOfClasses + which,
                     const_cast<ImageType *>( image ) );
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
const typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::ImageType
* AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GetIntensityImage( unsigned int which ) const
  {
  const ImageType *image;

  if( which == 0 )
    {
    image = dynamic_cast<const ImageType *>( this->ProcessObject::GetInput( 0 ) );
    }
  else if( which > 0 && which <= this->m_NumberOfIntensityImages )
    {
    image = dynamic_cast<const ImageType *>(
        this->ProcessObject::GetInput( 2 + this->m_NumberOfClasses + which ) );
    }
  else
    {
    itkExceptionMacro( "Image " << which << " is outside the range "
                                << "[1+m_NumberOfClasses...1+m_NumberOfClasses+m_NumberOfIntensityImages]." )
    }
  return image;
  }

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GenerateData()
{
  /**
   * Assign Gaussian likelihood functions if mixture model components are absent
   */
  typedef ants::Statistics::GaussianListSampleFunction
    <SampleType, float, float> LikelihoodType;
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    if( !this->GetLikelihoodFunction( n ) )
      {
      typename LikelihoodType::Pointer gaussianLikelihood =
        LikelihoodType::New();
      this->SetLikelihoodFunction( n, gaussianLikelihood );
      }
    }

  /**
   * Initialize the class labeling and the likelihood parameters
   */
  this->GenerateInitialClassLabeling();

  /**
   * Iterate until convergence or iterative exhaustion.
   */
  IterationReporter reporter( this, 0, 1 );

  bool     isConverged = false;
  RealType probabilityNew = 0.0;
  RealType probabilityOld = NumericTraits<RealType>::NonpositiveMin();

  unsigned int iteration = 0;
  while( !isConverged && iteration++ < this->m_MaximumNumberOfIterations )
    {
    /**
     * Clear the current posterior probability images to force
     * recalculation of the posterior probability images.
     */
    this->m_PosteriorProbabilityImages.clear();

    TimeProbe timer;
    timer.Start();
    probabilityNew = this->UpdateClassLabeling();
    timer.Stop();

    this->m_CurrentConvergenceMeasurement = probabilityNew - probabilityOld;

    if( this->m_CurrentConvergenceMeasurement < this->m_ConvergenceThreshold )
      {
      isConverged = true;
      }
    probabilityOld = probabilityNew;

    itkDebugMacro( "Iteration: " << probabilityNew );

    this->m_ElapsedIterations++;

    reporter.CompletedStep();
    }
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GenerateInitialClassLabeling()
{
  this->AllocateOutputs();
  this->GetOutput()->FillBuffer( NumericTraits<MaskLabelType>::Zero );

  switch( this->m_InitializationStrategy )
    {
    case Random:
      {
      typedef itk::Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
      typename GeneratorType::Pointer generator = GeneratorType::New();

      ImageRegionIterator<ClassifiedImageType> It( this->GetOutput(),
                                                   this->GetOutput()->GetRequestedRegion() );
      for( It.GoToBegin(); !It.IsAtEnd(); ++It )
        {
        LabelType label = generator->GetIntegerVariate(
            this->m_NumberOfClasses - 1 ) + 1;
        It.Set( label );
        }
      break;
      }
    case KMeans:  default:
      {
      this->GenerateInitialClassLabelingWithKMeansClustering();
      break;
      }
    case Otsu:
      {
      this->GenerateInitialClassLabelingWithOtsuThresholding();
      break;
      }
    case PriorProbabilityImages:
      {
      /**
       * Check for proper setting of prior probability images.
       */
      bool isOkay = true;
      for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
        {
        if( !this->GetPriorProbabilityImage( n + 1 ) )
          {
          isOkay = false;
          break;
          }
        }
      if( isOkay )
        {
        this->GenerateInitialClassLabelingWithPriorProbabilityImages();
        }
      else
        {
        itkWarningMacro( "The prior probability images were not set correctly."
                         << "Initializing with kmeans instead." );
        this->GenerateInitialClassLabelingWithKMeansClustering();
        }
      break;
      }
    case PriorLabelImage:
      {
      typedef ImageDuplicator<ClassifiedImageType> DuplicatorType;
      typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
      duplicator->SetInputImage( this->GetPriorLabelImage() );
      duplicator->Update();
      this->SetNthOutput( 0, duplicator->GetOutput() );
      break;
      }
    }

  /**
   * Calculate the initial parameters of the mixture model from the
   * initial labeling, i.e. the proportion, mean, and covariance for each label.
   */
  this->m_MixtureModelProportions.SetSize( this->m_NumberOfClasses );

  unsigned int totalSampleSize = 0;

  std::vector<typename SampleType::Pointer> samples;
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    typename SampleType::Pointer sample = SampleType::New();
    sample->SetMeasurementVectorSize( this->m_NumberOfIntensityImages );
    samples.push_back( sample );
    }

  /**
   * Accumulate the sample array for all labels.  Also accumulate the
   * prior probability weights, if applicable.
   */
  ImageRegionIteratorWithIndex<ClassifiedImageType> ItO( this->GetOutput(),
                                                         this->GetOutput()->GetRequestedRegion() );
  for( ItO.GoToBegin(); !ItO.IsAtEnd(); ++ItO )
    {
    LabelType label = ItO.Get();
    if( label == 0 )
      {
      continue;
      }
    typename SampleType::MeasurementVectorType measurement;
    ::itk::Statistics::MeasurementVectorTraits::SetLength( measurement,
                                                           this->m_NumberOfIntensityImages );
    for( unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++ )
      {
      measurement[i] =
        this->GetIntensityImage( i )->GetPixel( ItO.GetIndex() );
      }
    samples[label - 1]->PushBack( measurement );
    }

  /**
   * Create the weight array now that we know the sample sizes.
   */
  Array<unsigned int> count( this->m_NumberOfClasses );
  count.Fill( 0 );
  std::vector<typename LikelihoodFunctionType::WeightArrayType> weights;
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    totalSampleSize += samples[n]->Size();
    WeightArrayType weightArray( samples[n]->Size() );
    weightArray.Fill( 1.0 );
    weights.push_back( weightArray );
    }
  if( this->m_InitializationStrategy == PriorProbabilityImages )
    {
    for( ItO.GoToBegin(); !ItO.IsAtEnd(); ++ItO )
      {
      LabelType label = ItO.Get();
      if( label == 0 )
        {
        continue;
        }
      typename RealImageType::Pointer priorProbabilityImage =
        this->GetPriorProbabilityImage( label );
      weights[label - 1].SetElement( count[label - 1]++,
                                     priorProbabilityImage->GetPixel( ItO.GetIndex() ) );
      }
    }
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    this->m_MixtureModelComponents[n]->SetWeights( &weights[n] );
    this->m_MixtureModelComponents[n]->SetInputListSample( samples[n] );
    this->m_MixtureModelProportions[n] =
      static_cast<RealType>( samples[n]->Size() )
      / static_cast<RealType>( totalSampleSize );
    }
  for( unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++ )
    {
    ControlPointLatticeContainerType container;
    this->m_ControlPointLattices.push_back( container );
    for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
      {
      this->m_ControlPointLattices[i].push_back( NULL );
      }
    }
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GenerateInitialClassLabelingWithPriorProbabilityImages()
{
  // We first normalize prior probability images by keeping track of
  //   1. the sum of the prior probabilities at each pixel
  //   2. the max prior probability value
  //   3. which prior probability image corresponds to the max prior value

  typename RealImageType::Pointer sumPriorProbabilityImage = RealImageType::New();
  sumPriorProbabilityImage->SetOrigin( this->GetInput()->GetOrigin() );
  sumPriorProbabilityImage->SetSpacing( this->GetInput()->GetSpacing() );
  sumPriorProbabilityImage->SetRegions( this->GetInput()->GetRequestedRegion() );
  sumPriorProbabilityImage->SetDirection( this->GetInput()->GetDirection() );
  sumPriorProbabilityImage->Allocate();
  sumPriorProbabilityImage->FillBuffer( 0 );

  typename RealImageType::Pointer maxPriorProbabilityImage = RealImageType::New();
  maxPriorProbabilityImage->SetOrigin( this->GetInput()->GetOrigin() );
  maxPriorProbabilityImage->SetSpacing( this->GetInput()->GetSpacing() );
  maxPriorProbabilityImage->SetRegions( this->GetInput()->GetRequestedRegion() );
  maxPriorProbabilityImage->SetDirection( this->GetInput()->GetDirection() );
  maxPriorProbabilityImage->Allocate();
  maxPriorProbabilityImage->FillBuffer( 0 );
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    typename RealImageType::Pointer priorProbabilityImage =
      this->GetPriorProbabilityImage( n + 1 );

    ImageRegionIteratorWithIndex<ImageType> ItP( priorProbabilityImage,
                                                 priorProbabilityImage->GetRequestedRegion() );
    ImageRegionIterator<RealImageType> ItS( sumPriorProbabilityImage,
                                            sumPriorProbabilityImage->GetRequestedRegion() );
    ImageRegionIterator<RealImageType> ItM( maxPriorProbabilityImage,
                                            maxPriorProbabilityImage->GetRequestedRegion() );
    ImageRegionIterator<ClassifiedImageType> ItO( this->GetOutput(),
                                                  this->GetOutput()->GetRequestedRegion() );

    ItP.GoToBegin();
    ItS.GoToBegin();
    ItM.GoToBegin();
    ItO.GoToBegin();

    while( !ItP.IsAtEnd() )
      {
      if( !this->GetMaskImage() ||
          this->GetMaskImage()->GetPixel( ItP.GetIndex() ) == this->m_MaskLabel )
        {
        if( ItP.Get() > ItM.Get() )
          {
          ItM.Set( ItP.Get() );
          ItO.Set( n + 1 );
          }
        else if( ItP.Get() == ItM.Get() )
          {
          if( n == 0 )
            {
            ItO.Set( 1 );
            }
          else  // if maximum probabilities are the same, randomly select one
            {
            typedef itk::Statistics::MersenneTwisterRandomVariateGenerator
              GeneratorType;
            typename GeneratorType::Pointer generator = GeneratorType::New();
            if( generator->GetIntegerVariate( 1 ) )
              {
              ItO.Set( n + 1 );
              }
            }
          }
        }
      ++ItP;
      ++ItS;
      ++ItM;
      ++ItO;
      }
    }
  // Now we can normalize each prior probability image by dividing by the sum
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    typename RealImageType::Pointer priorProbabilityImage =
      this->GetPriorProbabilityImage( n + 1 );

    ImageRegionIteratorWithIndex<ImageType> ItP( priorProbabilityImage,
                                                 priorProbabilityImage->GetRequestedRegion() );
    ImageRegionIterator<RealImageType> ItS( sumPriorProbabilityImage,
                                            sumPriorProbabilityImage->GetRequestedRegion() );
    ImageRegionIterator<RealImageType> ItM( maxPriorProbabilityImage,
                                            maxPriorProbabilityImage->GetRequestedRegion() );
    ImageRegionIterator<ClassifiedImageType> ItO( this->GetOutput(),
                                                  this->GetOutput()->GetRequestedRegion() );

    ItP.GoToBegin();
    ItS.GoToBegin();
    ItM.GoToBegin();
    ItO.GoToBegin();

    while( !ItP.IsAtEnd() )
      {
      if( !this->GetMaskImage() ||
          this->GetMaskImage()->GetPixel( ItP.GetIndex() ) == this->m_MaskLabel )
        {
        if( ItM.Get() <= this->m_PriorProbabilityThreshold )
          {
          ItO.Set( NumericTraits<LabelType>::Zero );
          ItP.Set( NumericTraits<RealType>::Zero );
          }
        else
          {
          ItP.Set( ItP.Get() / ItS.Get() );
          }
        }
      ++ItP;
      ++ItS;
      ++ItM;
      ++ItO;
      }

    this->SetPriorProbabilityImage( n + 1, priorProbabilityImage );
    }

  /**
   * Set the initial output to be the prior label image.  This way we can
   * propogate the segmentation solution to regions of non-zero probability
   * but where the mask exists.
   */
  typedef ImageDuplicator<ClassifiedImageType> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage( this->GetOutput() );
  duplicator->Update();

  this->SetPriorLabelImage( duplicator->GetOutput() );
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GenerateInitialClassLabelingWithOtsuThresholding()
{
  RealType maxValue = NumericTraits<RealType>::min();
  RealType minValue = NumericTraits<RealType>::max();

  ImageRegionConstIteratorWithIndex<ImageType> ItI( this->GetInput(),
                                                    this->GetInput()->GetRequestedRegion() );
  for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
    {
    if( !this->GetMaskImage() ||
        this->GetMaskImage()->GetPixel( ItI.GetIndex() )
        == this->m_MaskLabel )
      {
      if( ItI.Get() < minValue )
        {
        minValue = ItI.Get();
        }
      else if( ItI.Get() > maxValue )
        {
        maxValue = ItI.Get();
        }
      }
    }

  typedef LabelStatisticsImageFilter<ImageType, MaskImageType> StatsType;
  typename StatsType::Pointer stats = StatsType::New();
  stats->SetInput( this->GetInput() );
  if( this->GetMaskImage() )
    {
    stats->SetLabelInput(
      const_cast<MaskImageType *>( this->GetMaskImage() ) );
    }
  else
    {
    typename MaskImageType::Pointer maskImage = MaskImageType::New();
    maskImage->SetOrigin( this->GetOutput()->GetOrigin() );
    maskImage->SetSpacing( this->GetOutput()->GetSpacing() );
    maskImage->SetRegions( this->GetOutput()->GetRequestedRegion() );
    maskImage->SetDirection( this->GetOutput()->GetDirection() );
    maskImage->Allocate();
    maskImage->FillBuffer( this->m_MaskLabel );
    stats->SetLabelInput( maskImage );
    }
  stats->UseHistogramsOn();
  stats->SetHistogramParameters( 200, minValue, maxValue );
  stats->Update();

  typedef OtsuMultipleThresholdsCalculator<typename StatsType::HistogramType>
    OtsuType;
  typename OtsuType::Pointer otsu = OtsuType::New();
  otsu->SetInputHistogram( stats->GetHistogram( this->m_MaskLabel ) );
  otsu->SetNumberOfThresholds( this->m_NumberOfClasses - 1 );
  otsu->Update();

  typename OtsuType::OutputType thresholds = otsu->GetOutput();

  ImageRegionIterator<ClassifiedImageType> ItO( this->GetOutput(),
                                                this->GetOutput()->GetRequestedRegion() );
  for( ItI.GoToBegin(), ItO.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItO )
    {
    LabelType label = NumericTraits<LabelType>::Zero;
    if( !this->GetMaskImage() ||
        this->GetMaskImage()->GetPixel( ItI.GetIndex() ) == this->m_MaskLabel )
      {
      if( ItI.Get() < thresholds[0] )
        {
        label = NumericTraits<LabelType>::One;
        }
      else
        {
        bool thresholdFound = false;
        for( unsigned int i = 1; i < thresholds.size(); i++ )
          {
          if( ItI.Get() >= thresholds[i - 1] && ItI.Get() <= thresholds[i] )
            {
            label = static_cast<LabelType>( i + 1 );
            thresholdFound = true;
            break;
            }
          }
        if( !thresholdFound )
          {
          label = static_cast<LabelType>( thresholds.size() + 1 );
          }
        }
      }
    ItO.Set( label );
    }
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GenerateInitialClassLabelingWithKMeansClustering()
{
  /**
   * We first perform kmeans on the first image and use the results to
   * seed the second run of kmeans using all the images.
   */
  typedef LabelStatisticsImageFilter<ImageType, MaskImageType> StatsType;
  typename StatsType::Pointer stats = StatsType::New();
  stats->SetInput( this->GetInput() );
  if( this->GetMaskImage() )
    {
    stats->SetLabelInput(
      const_cast<MaskImageType *>( this->GetMaskImage() ) );
    }
  else
    {
    typename MaskImageType::Pointer maskImage = MaskImageType::New();
    maskImage->SetOrigin( this->GetOutput()->GetOrigin() );
    maskImage->SetSpacing( this->GetOutput()->GetSpacing() );
    maskImage->SetRegions( this->GetOutput()->GetRequestedRegion() );
    maskImage->SetDirection( this->GetOutput()->GetDirection() );
    maskImage->Allocate();
    maskImage->FillBuffer( this->m_MaskLabel );
    stats->SetLabelInput( maskImage );
    }
  stats->UseHistogramsOff();
  stats->Update();

  RealType minValue = stats->GetMinimum( this->m_MaskLabel );
  RealType maxValue = stats->GetMaximum( this->m_MaskLabel );

  /**
   * The code below can be replaced by itkListSampleToImageFilter when we
   * migrate over to the Statistics classes current in the Review/ directory.
   */
  typename SampleType::Pointer sample = SampleType::New();
  sample->SetMeasurementVectorSize( 1 );

  ImageRegionConstIteratorWithIndex<ImageType> ItI( this->GetInput(),
                                                    this->GetInput()->GetRequestedRegion() );
  for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
    {
    if( !this->GetMaskImage() ||
        this->GetMaskImage()->GetPixel( ItI.GetIndex() ) == this->m_MaskLabel )
      {
      typename SampleType::MeasurementVectorType measurement;
      measurement.SetSize( 1 );
      measurement[0] = ItI.Get();
      sample->PushBack( measurement );
      }
    }

  typedef itk::Statistics::WeightedCentroidKdTreeGenerator<SampleType>
    TreeGeneratorType;
  typename TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
  treeGenerator->SetSample( sample );
  treeGenerator->SetBucketSize( 16 );
  treeGenerator->Update();

  /**
   * Guess initial class means by dividing the dynamic range of the first image
   *  into equal intervals.
   */

  typedef typename TreeGeneratorType::KdTreeType                TreeType;
  typedef itk::Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;
  typedef typename EstimatorType::ParametersType                ParametersType;

  ParametersType initialMeans( this->m_NumberOfClasses );
  initialMeans.Fill( 0.0 );
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    initialMeans[n] = minValue + ( maxValue - minValue )
      * ( static_cast<RealType>( n ) + 0.5 )
      / static_cast<RealType>( this->m_NumberOfClasses );
    }

  typename EstimatorType::Pointer estimator = EstimatorType::New();
  estimator->SetParameters( initialMeans );
  estimator->SetKdTree( treeGenerator->GetOutput() );
  estimator->SetMaximumIteration( 200 );
  estimator->SetCentroidPositionChangesThreshold( 0.0 );
  estimator->StartOptimization();

  /**
   * Classify the samples
   */
  typedef itk::Statistics::MinimumDecisionRule2 DecisionRuleType;
  typename DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();

  typedef itk::Statistics::SampleClassifierFilter<SampleType> ClassifierType;
  typename ClassifierType::Pointer classifier = ClassifierType::New();
  classifier->SetDecisionRule( decisionRule );
  classifier->SetInput( sample );
  classifier->SetNumberOfClasses( this->m_NumberOfClasses );

  typename ClassifierType::ClassLabelVectorObjectType::Pointer classLabels =
    ClassifierType::ClassLabelVectorObjectType::New();
  classifier->SetClassLabels( classLabels );
  typename ClassifierType::ClassLabelVectorType & classLabelVector =
    classLabels->Get();

  /**
   * Order the cluster means so that the lowest mean of the input image
   * corresponds to label '1', the second lowest to label '2', etc.
   */
  std::vector<RealType> estimatorParameters;
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    estimatorParameters.push_back( estimator->GetParameters()[n] );
    }
  std::sort( estimatorParameters.begin(), estimatorParameters.end() );

  typedef itk::Statistics::DistanceToCentroidMembershipFunction
    <MeasurementVectorType> MembershipFunctionType;
  typename ClassifierType::MembershipFunctionVectorObjectType::Pointer
  membershipFunctions = ClassifierType::MembershipFunctionVectorObjectType::New();
  typename ClassifierType::MembershipFunctionVectorType & membershipFunctionsVector =
    membershipFunctions->Get();

  classifier->SetMembershipFunctions( membershipFunctions );
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    typename MembershipFunctionType::Pointer
    membershipFunction = MembershipFunctionType::New();
    membershipFunction->SetMeasurementVectorSize(
      sample->GetMeasurementVectorSize() );
    typename MembershipFunctionType::CentroidType centroid;
    ::itk::Statistics::MeasurementVectorTraits::SetLength( centroid,
                                                           sample->GetMeasurementVectorSize() );
    centroid[0] = estimatorParameters[n];
    membershipFunction->SetCentroid( centroid );
    membershipFunctionsVector.push_back( membershipFunction.GetPointer() );

    classLabelVector.push_back(
      static_cast<typename ClassifierType::ClassLabelType>( n + 1 ) );
    }
  classifier->Update();

  /**
   * Classify the voxels
   */
  typedef typename ClassifierType::MembershipSampleType ClassifierOutputType;
  typedef typename ClassifierOutputType::ConstIterator  LabelIterator;

  ImageRegionIteratorWithIndex<ClassifiedImageType> ItO( this->GetOutput(),
                                                         this->GetOutput()->GetRequestedRegion() );
  ItO.GoToBegin();
  LabelIterator it = classifier->GetOutput()->Begin();
  while( it != classifier->GetOutput()->End() )
    {
    if( !this->GetMaskImage() ||
        this->GetMaskImage()->GetPixel( ItO.GetIndex() ) == this->m_MaskLabel )
      {
      ItO.Set( it.GetClassLabel() );
      ++it;
      }
    else
      {
      ItO.Set( NumericTraits<LabelType>::Zero );
      }
    ++ItO;
    }

  /**
   * If there are 1 or more auxiliary images, use the results from the first
   * kmeans grouping to seed a second invocation of the algorithm using both
   * the input image and the auxiliary images.
   */

  if( this->m_NumberOfIntensityImages == 0 )
    {
    return;
    }
  else
    {
    Array<RealType> minValues;
    minValues.SetSize( this->m_NumberOfIntensityImages );
    minValues.Fill( NumericTraits<RealType>::max() );
    Array<RealType> maxValues;
    maxValues.SetSize( this->m_NumberOfIntensityImages );
    maxValues.Fill( NumericTraits<RealType>::NonpositiveMin() );

    minValues[0] = minValue;
    maxValues[0] = maxValue;
    for( unsigned int i = 1; i < this->m_NumberOfIntensityImages; i++ )
      {
      typedef LabelStatisticsImageFilter<ImageType, ClassifiedImageType> StatsType;
      typename StatsType::Pointer stats = StatsType::New();
      stats->SetInput( this->GetIntensityImage( i ) );
      stats->SetLabelInput( this->GetOutput() );
      stats->UseHistogramsOff();
      stats->Update();

      minValues[i] = stats->GetMinimum( this->m_MaskLabel );
      maxValues[i] = stats->GetMaximum( this->m_MaskLabel );
      }

    typedef itk::Statistics::WeightedCentroidKdTreeGenerator<SampleType>
      TreeGeneratorType;
    typedef typename TreeGeneratorType::KdTreeType                TreeType;
    typedef itk::Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;
    typedef typename EstimatorType::ParametersType                ParametersType;

    /**
     * Guess initial class means by dividing the dynamic range of the first image
     *  into equal intervals.
     */
    typename EstimatorType::Pointer estimator = EstimatorType::New();
    ParametersType initialMeans( this->m_NumberOfClasses
                                 * ( this->m_NumberOfIntensityImages ) );
    initialMeans.Fill( 0.0 );
    for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
      {
      for( unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++ )
        {
        initialMeans[( this->m_NumberOfIntensityImages ) * n + i] = minValues[i]
          + ( maxValues[0] - minValues[0] ) * ( static_cast<RealType>( n ) + 0.5 )
          / static_cast<RealType>( this->m_NumberOfClasses );
        }
      }

    typename SampleType::Pointer sample = SampleType::New();
    sample->SetMeasurementVectorSize( this->m_NumberOfIntensityImages );
    for( ItO.GoToBegin(); !ItO.IsAtEnd(); ++ItO )
      {
      if( !this->GetMaskImage() ||
          this->GetMaskImage()->GetPixel( ItO.GetIndex() ) == this->m_MaskLabel )
        {
        typename SampleType::MeasurementVectorType measurement;
        measurement.SetSize( this->m_NumberOfIntensityImages );
        for( unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++ )
          {
          measurement[i] =
            this->GetIntensityImage( i )->GetPixel( ItO.GetIndex() );
          measurement[i] = minValues[0] + ( measurement[i] - minValues[i] )
            * ( maxValues[0] - minValues[0] ) / ( maxValues[i] - minValues[i] );
          }
        sample->PushBack( measurement );
        }
      }

    typename TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
    treeGenerator->SetSample( sample );
    treeGenerator->SetBucketSize( 16 );
    treeGenerator->Update();

    estimator->SetParameters( initialMeans );
    estimator->SetKdTree( treeGenerator->GetOutput() );
    estimator->SetMaximumIteration( 200 );
    estimator->SetCentroidPositionChangesThreshold( 0.0 );
    estimator->StartOptimization();

    /**
     * Classify the samples
     */
    typedef itk::Statistics::MinimumDecisionRule2 DecisionRuleType;
    typename DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();

    typedef itk::Statistics::SampleClassifierFilter<SampleType> ClassifierType;
    typename ClassifierType::Pointer classifier = ClassifierType::New();
    classifier->SetDecisionRule( decisionRule );
    classifier->SetInput( sample );
    classifier->SetNumberOfClasses( this->m_NumberOfClasses );

    typename ClassifierType::ClassLabelVectorObjectType::Pointer classLabels =
      ClassifierType::ClassLabelVectorObjectType::New();
    classifier->SetClassLabels( classLabels );
    typename ClassifierType::ClassLabelVectorType & classLabelVector =
      classLabels->Get();

    /**
     * Order the cluster means so that the lowest mean of the input image
     * corresponds to label '1', the second lowest to label '2', etc.
     */
    std::vector<RealType> estimatorParameters;
    for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
      {
      estimatorParameters.push_back( estimator->GetParameters()[n] );
      }
    std::sort( estimatorParameters.begin(), estimatorParameters.end() );

    typedef itk::Statistics::DistanceToCentroidMembershipFunction
      <MeasurementVectorType> MembershipFunctionType;
    typename ClassifierType::MembershipFunctionVectorObjectType::Pointer
    membershipFunctions = ClassifierType::MembershipFunctionVectorObjectType::New();
    typename ClassifierType::MembershipFunctionVectorType & membershipFunctionsVector =
      membershipFunctions->Get();

    classifier->SetMembershipFunctions( membershipFunctions );
    for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
      {
      typename MembershipFunctionType::Pointer
      membershipFunction = MembershipFunctionType::New();
      membershipFunction->SetMeasurementVectorSize(
        sample->GetMeasurementVectorSize() );
      typename MembershipFunctionType::CentroidType centroid;
      ::itk::Statistics::MeasurementVectorTraits::SetLength( centroid,
                                                             sample->GetMeasurementVectorSize() );
      for( unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++ )
        {
        centroid[i] = estimator->GetParameters()[
            ( this->m_NumberOfIntensityImages ) * n + i];
        }
      membershipFunction->SetCentroid( centroid );
      membershipFunctionsVector.push_back( membershipFunction.GetPointer() );

      classLabelVector.push_back(
        static_cast<typename ClassifierType::ClassLabelType>( n + 1 ) );
      }
    classifier->Update();

    /**
     * Classify the voxels
     */
    typedef typename ClassifierType::MembershipSampleType ClassifierOutputType;
    typedef typename ClassifierOutputType::ConstIterator  LabelIterator;

    ItO.GoToBegin();
    LabelIterator it = classifier->GetOutput()->Begin();
    while( it != classifier->GetOutput()->End() )
      {
      if( !this->GetMaskImage() ||
          this->GetMaskImage()->GetPixel( ItO.GetIndex() ) == this->m_MaskLabel )
        {
        ItO.Set( it.GetClassLabel() );
        ++it;
        }
      else
        {
        ItO.Set( NumericTraits<LabelType>::Zero );
        }
      ++ItO;
      }
    }
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::RealType
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::UpdateClassLabeling()
{
  typename RealImageType::Pointer maxPosteriorProbabilityImage =
    RealImageType::New();
  maxPosteriorProbabilityImage->SetOrigin( this->GetOutput()->GetOrigin() );
  maxPosteriorProbabilityImage->SetSpacing( this->GetOutput()->GetSpacing() );
  maxPosteriorProbabilityImage->SetDirection(
    this->GetOutput()->GetDirection() );
  maxPosteriorProbabilityImage->SetRegions(
    this->GetOutput()->GetRequestedRegion() );
  maxPosteriorProbabilityImage->Allocate();
  maxPosteriorProbabilityImage->FillBuffer( NumericTraits<RealType>::Zero );

  typename ClassifiedImageType::Pointer maxLabels =
    ClassifiedImageType::New();
  maxLabels->SetRegions( this->GetOutput()->GetRequestedRegion() );
  maxLabels->SetOrigin( this->GetOutput()->GetOrigin() );
  maxLabels->SetSpacing( this->GetOutput()->GetSpacing() );
  maxLabels->SetDirection( this->GetOutput()->GetDirection() );
  maxLabels->Allocate();
  maxLabels->FillBuffer( NumericTraits<LabelType>::Zero );

  typename RealImageType::Pointer weightedPriorProbabilityImage =
    RealImageType::New();
  weightedPriorProbabilityImage->SetRegions(
    this->GetOutput()->GetRequestedRegion() );
  weightedPriorProbabilityImage->SetOrigin( this->GetOutput()->GetOrigin() );
  weightedPriorProbabilityImage->SetSpacing( this->GetOutput()->GetSpacing() );
  weightedPriorProbabilityImage->SetDirection(
    this->GetOutput()->GetDirection() );
  weightedPriorProbabilityImage->Allocate();
  weightedPriorProbabilityImage->FillBuffer( NumericTraits<RealType>::Zero );

  /**
   * Accumulate the sample array for all labels.
   */
  std::vector<typename SampleType::Pointer> samples;
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    typename SampleType::Pointer sample = SampleType::New();
    samples.push_back( this->GetScalarSamples( n + 1 ) );
    }
  unsigned int totalSampleSize = samples[0]->Size();

  Array<unsigned int> count( this->m_NumberOfClasses );
  count.Fill( 0 );
  std::vector<WeightArrayType> weights;
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    WeightArrayType weightArray( samples[n]->Size() );
    weightArray.Fill( 1.0 );
    weights.push_back( weightArray );
    }

  Array<RealType> sumPosteriors( this->m_NumberOfClasses );
  sumPosteriors.Fill( 0.0 );
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    typename RealImageType::Pointer posteriorProbabilityImage
      = this->GetPosteriorProbabilityImage( n + 1 );

    ImageRegionIteratorWithIndex<ClassifiedImageType> ItO( maxLabels,
                                                           maxLabels->GetRequestedRegion() );
    ImageRegionConstIterator<RealImageType> ItP( posteriorProbabilityImage,
                                                 posteriorProbabilityImage->GetRequestedRegion() );
    ImageRegionIterator<RealImageType> ItM( maxPosteriorProbabilityImage,
                                            maxPosteriorProbabilityImage->GetRequestedRegion() );

    ItP.GoToBegin();
    ItM.GoToBegin();
    ItO.GoToBegin();
    while( !ItP.IsAtEnd() )
      {
      if( !this->GetMaskImage() ||
          this->GetMaskImage()->GetPixel( ItO.GetIndex() ) == this->m_MaskLabel )
        {
        RealType posteriorProbability = ItP.Get();
        weights[n].SetElement( count[n]++, posteriorProbability );

        if( posteriorProbability > 0.0 && posteriorProbability >= ItM.Get() )
          {
          ItM.Set( posteriorProbability );
          ItO.Set( static_cast<LabelType>( n + 1 ) );
          }
        sumPosteriors[n] += posteriorProbability;
        }
      ++ItP;
      ++ItM;
      ++ItO;
      }

    /**
     * Perform the following calculation as a preprocessing step to update the
     * class proportions.
     */
    typename RealImageType::Pointer distancePriorProbabilityImage =
      this->GetDistancePriorProbabilityImageFromPriorLabelImage( n + 1 );
    typename RealImageType::Pointer priorProbabilityImage =
      this->GetPriorProbabilityImage( n + 1 );

    ImageRegionIteratorWithIndex<RealImageType> ItW(
      weightedPriorProbabilityImage,
      weightedPriorProbabilityImage->GetRequestedRegion() );
    for( ItW.GoToBegin(); !ItW.IsAtEnd(); ++ItW )
      {
      if( !this->GetMaskImage() ||
          this->GetMaskImage()->GetPixel( ItW.GetIndex() ) == this->m_MaskLabel )
        {
        RealType priorProbability = 0.0;
        if( this->m_InitializationStrategy == PriorLabelImage ||
            this->m_InitializationStrategy == PriorProbabilityImages )
          {
          if( priorProbabilityImage )
            {
            priorProbability = priorProbabilityImage->GetPixel( ItW.GetIndex() );
            }
          if( priorProbability <= this->m_PriorProbabilityThreshold &&
              distancePriorProbabilityImage )
            {
            priorProbability =
              distancePriorProbabilityImage->GetPixel( ItW.GetIndex() );
            }
          }
        else
          {
          priorProbability = 1.0;
          }
        ItW.Set( ItW.Get() + this->m_MixtureModelProportions[n]
                 * priorProbability );
        }
      }
    }
  this->SetNthOutput( 0, maxLabels );
  /**
   * Update the class proportions
   */
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    RealType denominator = 0.0;

    typename RealImageType::Pointer distancePriorProbabilityImage =
      this->GetDistancePriorProbabilityImageFromPriorLabelImage( n + 1 );
    typename RealImageType::Pointer priorProbabilityImage =
      this->GetPriorProbabilityImage( n + 1 );

    ImageRegionIteratorWithIndex<RealImageType> ItW(
      weightedPriorProbabilityImage,
      weightedPriorProbabilityImage->GetRequestedRegion() );
    for( ItW.GoToBegin(); !ItW.IsAtEnd(); ++ItW )
      {
      if( !this->GetMaskImage() ||
          this->GetMaskImage()->GetPixel( ItW.GetIndex() ) == this->m_MaskLabel )
        {
        RealType priorProbability = 0.0;
        if( this->m_InitializationStrategy == PriorLabelImage ||
            this->m_InitializationStrategy == PriorProbabilityImages )
          {
          if( priorProbabilityImage )
            {
            priorProbability = priorProbabilityImage->GetPixel( ItW.GetIndex() );
            }
          if( priorProbability == 0.0 && distancePriorProbabilityImage )
            {
            priorProbability =
              distancePriorProbabilityImage->GetPixel( ItW.GetIndex() );
            }
          }
        else
          {
          priorProbability = 1.0;
          }
        if( ItW.Get() > 0.0 )
          {
          denominator += ( priorProbability / ItW.Get() );
          }
        }
      }
    if( denominator > 0.0 )
      {
      this->m_MixtureModelProportions[n] = sumPosteriors[n] / denominator;
      }
    else
      {
      this->m_MixtureModelProportions[n] = 0.0;
      }
    }
  /**
   * Calculate the parameters of the mixture model.
   */
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    this->m_MixtureModelComponents[n]->SetWeights( &weights[n] );
    this->m_MixtureModelComponents[n]->SetInputListSample( samples[n] );
    }

  /**
   * Calculate the maximum posterior probability sum over the region of
   * interest.  This quantity should increase at each iteration.
   */
  RealType maxPosteriorSum = 0.0;

  ImageRegionIteratorWithIndex<RealImageType> ItM( maxPosteriorProbabilityImage,
                                                   maxPosteriorProbabilityImage->GetRequestedRegion() );
  for( ItM.GoToBegin(); !ItM.IsAtEnd(); ++ItM )
    {
    if( !this->GetMaskImage() ||
        this->GetMaskImage()->GetPixel( ItM.GetIndex() ) == this->m_MaskLabel )
      {
      maxPosteriorSum += ItM.Get();
      }
    }

  return maxPosteriorSum / static_cast<RealType>( totalSampleSize );
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::SamplePointer
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GetScalarSamples( unsigned int whichClass )
{
  /**
   * This function returns a set of samples for each class such that each
   * measurement vector of the returned SampleType corresponds to a single
   * voxel across the set of auxiliary and input images.
   */

  std::vector<typename SampleType::Pointer> samples;

  /**
   * Accumulate the samples in individual SampleTypes.  This allows us to
   * "filter" the samples of each auxiliary/input image.  This filtering
   * could including outlier winsorization, log transformation, and/or
   * converting vector/tensor auxiliary images to scalar data for
   * modeling.
   */
  for( unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++ )
    {
    typename SampleType::Pointer sample = SampleType::New();
    sample->SetMeasurementVectorSize( 1 );
    samples.push_back( sample );
    }

  ImageRegionIteratorWithIndex<ClassifiedImageType> ItO( this->GetOutput(),
                                                         this->GetOutput()->GetRequestedRegion() );
  for( ItO.GoToBegin(); !ItO.IsAtEnd(); ++ItO )
    {
    if( !this->GetMaskImage() ||
        this->GetMaskImage()->GetPixel( ItO.GetIndex() ) == this->m_MaskLabel )
      {
      for( unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++ )
        {
        typename SampleType::MeasurementVectorType measurement;
        ::itk::Statistics::MeasurementVectorTraits::SetLength( measurement, 1 );
        measurement.SetSize( 1 );
        measurement[0] =
          this->GetIntensityImage( i )->GetPixel( ItO.GetIndex() );
        samples[i]->PushBack( measurement );
        }
      }
    }

  /**
   * Simultaneously filter the samples and accumulate for return.
   */
  typename SampleType::Pointer scalarSamples = SampleType::New();
  scalarSamples->SetMeasurementVectorSize( this->m_NumberOfIntensityImages );
  for( unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++ )
    {
    typename SampleType::Pointer univariateSamples = SampleType::New();
    univariateSamples->SetMeasurementVectorSize( 1 );
    if( this->m_OutlierHandlingFilter )
      {
      this->m_OutlierHandlingFilter->SetInput( samples[i] );
      this->m_OutlierHandlingFilter->Update();
      univariateSamples = this->m_OutlierHandlingFilter->GetOutput();
      }
    else
      {
      univariateSamples = samples[i];
      }

    if( i == 0 )
      {
      typename SampleType::ConstIterator It = univariateSamples->Begin();
      while( It != univariateSamples->End() )
        {
        typename SampleType::MeasurementVectorType measurement;
        measurement.SetSize( this->m_NumberOfIntensityImages );
        ::itk::Statistics::MeasurementVectorTraits::SetLength( measurement,
                                                               this->m_NumberOfIntensityImages );
        measurement[0] = It.GetMeasurementVector()[0];
        scalarSamples->PushBack( measurement );
        ++It;
        }
      }
    else
      {
      typename SampleType::Iterator ItS = scalarSamples->Begin();
      typename SampleType::ConstIterator It = univariateSamples->Begin();
      while( ItS != scalarSamples->End() )
        {
        scalarSamples->SetMeasurement( ItS.GetInstanceIdentifier(), i,
                                       It.GetMeasurementVector()[0] );
        ++It;
        ++ItS;
        }
      }
    }

  return scalarSamples;
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::RealImageType::Pointer
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GetPosteriorProbabilityImage( unsigned int whichClass )
{
  if( whichClass > this->m_NumberOfClasses )
    {
    itkExceptionMacro(
      "Requested class is greater than the number of classes." );
    }

  /**
   * If memory minimization is turned off and if the posterior probability
   * images have already been calculated, simply return the probability
   * image for the requested class.  Otherwise, calculate the probability
   * image.
   */
  if( whichClass <= this->m_PosteriorProbabilityImages.size() )
    {
    return this->m_PosteriorProbabilityImages[whichClass - 1];
    }
  else
    {
    /**
     * Here we assume that the calling function is invoked in order such
     * that GetPosteriorProbabilityImage( 1 ) is called before
     * GetPosteriorProbabilityImage( 2 ), etc.  As such, when this part of
     * the code is reached and the class requested is '1', we assume that
     * the sum of the posterior probability images needs to be calculated
     * for normalization purposes.  This sum is then saved for subsequent calls.
     */
    typename RealImageType::Pointer posteriorProbabilityImage =
      RealImageType::New();
    posteriorProbabilityImage->SetRegions(
      this->GetOutput()->GetRequestedRegion() );
    posteriorProbabilityImage->SetOrigin( this->GetOutput()->GetOrigin() );
    posteriorProbabilityImage->SetSpacing( this->GetOutput()->GetSpacing() );
    posteriorProbabilityImage->SetDirection( this->GetOutput()->GetDirection() );
    posteriorProbabilityImage->Allocate();
    posteriorProbabilityImage->FillBuffer( 0 );

    /**
     * Calculate the sum of the probability images.  Also, store the
     * posterior probability images if m_MinimizeMemoryUsage == false.
     */
    if( whichClass == 1 )
      {
      this->m_SumPosteriorProbabilityImage = RealImageType::New();
      this->m_SumPosteriorProbabilityImage->SetRegions(
        this->GetOutput()->GetRequestedRegion() );
      this->m_SumPosteriorProbabilityImage->SetOrigin(
        this->GetOutput()->GetOrigin() );
      this->m_SumPosteriorProbabilityImage->SetSpacing(
        this->GetOutput()->GetSpacing() );
      this->m_SumPosteriorProbabilityImage->SetDirection(
        this->GetOutput()->GetDirection() );
      this->m_SumPosteriorProbabilityImage->Allocate();
      this->m_SumPosteriorProbabilityImage->FillBuffer( 0 );
      for( unsigned int c = 0; c < this->m_NumberOfClasses; c++ )
        {
        std::vector<typename RealImageType::Pointer> smoothImages;

        if( this->m_InitializationStrategy == PriorProbabilityImages ||
            this->m_InitializationStrategy == PriorLabelImage )
          {
          for( unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++ )
            {
            if( this->m_AdaptiveSmoothingWeights.size() > i &&
                this->m_AdaptiveSmoothingWeights[i] > 0.0 )
              {
              smoothImages.push_back(
                this->CalculateSmoothIntensityImageFromPriorProbabilityImage( i,
                                                                              c + 1 ) );
              }
            else
              {
              smoothImages.push_back( NULL );
              }
            }
          }

        typename RealImageType::Pointer distancePriorProbabilityImage =
          this->GetDistancePriorProbabilityImageFromPriorLabelImage( c + 1 );
        typename RealImageType::Pointer priorProbabilityImage =
          this->GetPriorProbabilityImage( c + 1 );

        typename NeighborhoodIterator<ClassifiedImageType>::RadiusType radius;
        unsigned int neighborhoodSize = 1;
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          neighborhoodSize *= ( 2 * this->m_MRFRadius[d] + 1 );
          radius[d] = this->m_MRFRadius[d];
          }

        ConstNeighborhoodIterator<ClassifiedImageType> ItO( radius,
                                                            this->GetOutput(),
                                                            this->GetOutput()->GetRequestedRegion() );
        ImageRegionIterator<RealImageType> ItS(
          this->m_SumPosteriorProbabilityImage,
          this->m_SumPosteriorProbabilityImage->GetRequestedRegion() );
        for( ItO.GoToBegin(), ItS.GoToBegin(); !ItO.IsAtEnd(); ++ItO, ++ItS )
          {
          if( !this->GetMaskImage() ||
              this->GetMaskImage()->GetPixel( ItO.GetIndex() ) == this->m_MaskLabel )
            {
            RealType mrfPrior = 1.0;
            if( this->m_MRFSmoothingFactor > 0.0 && neighborhoodSize > 1 )
              {
              Array<RealType> weightedNumberOfClassNeighbors;
              weightedNumberOfClassNeighbors.SetSize( this->m_NumberOfClasses );
              weightedNumberOfClassNeighbors.Fill( 0.0 );
              for( unsigned int n = 0; n < neighborhoodSize; n++ )
                {
                if( n == static_cast<unsigned int>( 0.5 * neighborhoodSize ) )
                  {
                  continue;
                  }
                bool      isInBounds = false;
                LabelType label = ItO.GetPixel( n, isInBounds );
                if( !isInBounds || label == 0 )
                  {
                  continue;
                  }
                typename ClassifiedImageType::OffsetType offset
                  = ItO.GetOffset( n );

                double distance = 0.0;
                for( unsigned int d = 0; d < ImageDimension; d++ )
                  {
                  distance += vnl_math_sqr( offset[d]
                                            * this->GetOutput()->GetSpacing()[d] );
                  }
                distance = vcl_sqrt( distance );
                weightedNumberOfClassNeighbors[label - 1] += ( 1.0 / distance );
                }
              RealType numerator = vcl_exp( this->m_MRFSmoothingFactor
                                            * weightedNumberOfClassNeighbors[c] );
              RealType denominator = 0.0;
              for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
                {
                denominator += vcl_exp( this->m_MRFSmoothingFactor
                                        * weightedNumberOfClassNeighbors[n] );
                }
              mrfPrior = numerator / denominator;
              if( vnl_math_isinf( mrfPrior ) || vnl_math_isnan( mrfPrior ) )
                {
                mrfPrior = 1.0;
                }
              }

            RealType priorProbability = 0.0;
            if( this->m_InitializationStrategy == PriorLabelImage ||
                this->m_InitializationStrategy == PriorProbabilityImages )
              {
              if( priorProbabilityImage )
                {
                priorProbability = priorProbabilityImage->GetPixel( ItO.GetIndex() );
                }
              if( priorProbability <= this->m_PriorProbabilityThreshold &&
                  distancePriorProbabilityImage )
                {
                priorProbability = distancePriorProbabilityImage->GetPixel( ItO.GetIndex() );
                }
              }
            else
              {
              priorProbability = 1.0;
              }

            MeasurementVectorType measurement;
            measurement.SetSize( this->m_NumberOfIntensityImages );
            for( unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++ )
              {
              measurement[i] =
                this->GetIntensityImage( i )->GetPixel( ItO.GetIndex() );

              if( ( this->m_InitializationStrategy == PriorProbabilityImages ||
                    this->m_InitializationStrategy == PriorLabelImage ) &&
                  smoothImages[i] )
                {
                measurement[i] = ( 1.0 - this->m_AdaptiveSmoothingWeights[i] )
                  * measurement[i] + this->m_AdaptiveSmoothingWeights[i]
                  * smoothImages[i]->GetPixel( ItO.GetIndex() );
                }
              }
            RealType likelihood =
              this->m_MixtureModelComponents[c]->Evaluate( measurement );
            RealType posteriorProbability = this->m_PriorProbabilityWeight
              * this->m_MixtureModelProportions[c] * likelihood
              * mrfPrior * priorProbability + ( 1.0
                                                - this->m_PriorProbabilityWeight )
              * this->m_MixtureModelProportions[c] * likelihood * mrfPrior;

            if( vnl_math_isnan( posteriorProbability ) ||
                vnl_math_isinf( posteriorProbability ) )
              {
              posteriorProbability = 0.0;
              }

            ItS.Set( ItS.Get() + posteriorProbability  );
            if( ( c == 0 ) || !this->m_MinimizeMemoryUsage )
              {
              posteriorProbabilityImage->SetPixel( ItO.GetIndex(),
                                                   posteriorProbability );
              }
            }
          }
        if( !this->m_MinimizeMemoryUsage )
          {
          typedef ImageDuplicator<RealImageType> DuplicatorType;
          typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
          duplicator->SetInputImage( posteriorProbabilityImage );
          duplicator->Update();

          this->m_PosteriorProbabilityImages.push_back(
            duplicator->GetOutput() );
          }
        }

      /**
       * Normalize the posterior probability image(s).
       */
      ImageRegionIterator<RealImageType> ItS(
        this->m_SumPosteriorProbabilityImage,
        this->m_SumPosteriorProbabilityImage->GetRequestedRegion() );
      if( this->m_MinimizeMemoryUsage )
        {
        ImageRegionIterator<RealImageType> ItP( posteriorProbabilityImage,
                                                posteriorProbabilityImage->GetRequestedRegion() );
        for( ItP.GoToBegin(), ItS.GoToBegin(); !ItS.IsAtEnd(); ++ItP, ++ItS )
          {
          if( ItS.Get() > 0 )
            {
            ItP.Set( ItP.Get() / ItS.Get() );
            }
          }
        return posteriorProbabilityImage;
        }
      else
        {
        for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
          {
          ImageRegionIterator<RealImageType> ItP(
            this->m_PosteriorProbabilityImages[n],
            this->m_PosteriorProbabilityImages[n]->GetRequestedRegion() );
          for( ItP.GoToBegin(), ItS.GoToBegin(); !ItS.IsAtEnd(); ++ItP, ++ItS )
            {
            if( ItS.Get() > 0 )
              {
              ItP.Set( ItP.Get() / ItS.Get() );
              }
            }
          }
        return this->m_PosteriorProbabilityImages[0];
        }
      }
    else // whichClass > 1
      {
      std::vector<typename RealImageType::Pointer> smoothImages;

      if( this->m_InitializationStrategy == PriorProbabilityImages ||
          this->m_InitializationStrategy == PriorLabelImage )
        {
        for( unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++ )
          {
          if( this->m_AdaptiveSmoothingWeights.size() > i &&
              this->m_AdaptiveSmoothingWeights[i] > 0.0 )
            {
            smoothImages.push_back(
              this->CalculateSmoothIntensityImageFromPriorProbabilityImage( i,
                                                                            whichClass ) );
            }
          else
            {
            smoothImages.push_back( NULL );
            }
          }
        }

      typename RealImageType::Pointer distancePriorProbabilityImage =
        this->GetDistancePriorProbabilityImageFromPriorLabelImage( whichClass );
      typename RealImageType::Pointer priorProbabilityImage =
        this->GetPriorProbabilityImage( whichClass );

      typename NeighborhoodIterator<ClassifiedImageType>::RadiusType radius;
      unsigned int neighborhoodSize = 1;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        neighborhoodSize *= ( 2 * this->m_MRFRadius[d] + 1 );
        radius[d] = this->m_MRFRadius[d];
        }

      ConstNeighborhoodIterator<ClassifiedImageType> ItO( radius,
                                                          this->GetOutput(), this->GetOutput()->GetRequestedRegion() );
      for( ItO.GoToBegin(); !ItO.IsAtEnd(); ++ItO )
        {
        if( !this->GetMaskImage() ||
            this->GetMaskImage()->GetPixel( ItO.GetIndex() ) == this->m_MaskLabel )
          {
          RealType mrfPrior = 1.0;
          if( this->m_MRFSmoothingFactor > 0.0 && neighborhoodSize > 1 )
            {
            Array<RealType> weightedNumberOfClassNeighbors;
            weightedNumberOfClassNeighbors.SetSize( this->m_NumberOfClasses );
            weightedNumberOfClassNeighbors.Fill( 0.0 );
            for( unsigned int n = 0; n < neighborhoodSize; n++ )
              {
              if( n == static_cast<unsigned int>( 0.5 * neighborhoodSize ) )
                {
                continue;
                }
              bool      isInBounds = false;
              LabelType label = ItO.GetPixel( n, isInBounds );
              if( !isInBounds || label == 0 )
                {
                continue;
                }
              typename ClassifiedImageType::OffsetType offset
                = ItO.GetOffset( n );

              double distance = 0.0;
              for( unsigned int d = 0; d < ImageDimension; d++ )
                {
                distance += vnl_math_sqr( offset[d]
                                          * this->GetOutput()->GetSpacing()[d] );
                }
              distance = vcl_sqrt( distance );
              weightedNumberOfClassNeighbors[label - 1] += ( 1.0 / distance );
              }
            RealType numerator = vcl_exp( this->m_MRFSmoothingFactor
                                          * weightedNumberOfClassNeighbors[whichClass - 1] );
            RealType denominator = 0.0;
            for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
              {
              denominator += vcl_exp( this->m_MRFSmoothingFactor
                                      * weightedNumberOfClassNeighbors[n] );
              }
            mrfPrior = numerator / denominator;
            if( vnl_math_isinf( mrfPrior ) || vnl_math_isnan( mrfPrior ) )
              {
              mrfPrior = 1.0;
              }
            }

          RealType priorProbability = 0.0;
          if( this->m_InitializationStrategy == PriorLabelImage ||
              this->m_InitializationStrategy == PriorProbabilityImages )
            {
            if( priorProbabilityImage )
              {
              priorProbability = priorProbabilityImage->GetPixel( ItO.GetIndex() );
              }
            if( priorProbability <= this->m_PriorProbabilityThreshold &&
                distancePriorProbabilityImage )
              {
              priorProbability = distancePriorProbabilityImage->GetPixel( ItO.GetIndex() );
              }
            if( priorProbability <= this->m_PriorProbabilityThreshold )
              {
              priorProbability = 1.0;
              }
            }
          else
            {
            priorProbability = 1.0;
            }

          MeasurementVectorType measurement;
          measurement.SetSize( this->m_NumberOfIntensityImages );
          for( unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++ )
            {
            measurement[i] =
              this->GetIntensityImage( i )->GetPixel( ItO.GetIndex() );

            if( ( this->m_InitializationStrategy == PriorProbabilityImages ||
                  this->m_InitializationStrategy == PriorLabelImage ) &&
                smoothImages[i] )
              {
              measurement[i] = ( 1.0 - this->m_AdaptiveSmoothingWeights[i] )
                * measurement[i] + this->m_AdaptiveSmoothingWeights[i]
                * smoothImages[i]->GetPixel( ItO.GetIndex() );
              }
            }
          RealType likelihood =
            this->m_MixtureModelComponents[whichClass - 1]->Evaluate( measurement );
          RealType posteriorProbability = this->m_PriorProbabilityWeight
            * this->m_MixtureModelProportions[whichClass - 1] * likelihood
            * mrfPrior * priorProbability + ( 1.0
                                              - this->m_PriorProbabilityWeight )
            * this->m_MixtureModelProportions[whichClass - 1] * likelihood
            * mrfPrior;

          if( vnl_math_isnan( posteriorProbability ) ||
              vnl_math_isinf( posteriorProbability ) )
            {
            posteriorProbability = 0.0;
            }
          posteriorProbabilityImage->SetPixel( ItO.GetIndex(),
                                               posteriorProbability );
          }
        }

      /**
       * Normalize the posterior probability image(s).
       */
      ImageRegionIterator<RealImageType> ItS(
        this->m_SumPosteriorProbabilityImage,
        this->m_SumPosteriorProbabilityImage->GetRequestedRegion() );
      ImageRegionIterator<RealImageType> ItP( posteriorProbabilityImage,
                                              posteriorProbabilityImage->GetRequestedRegion() );
      for( ItP.GoToBegin(), ItS.GoToBegin(); !ItS.IsAtEnd(); ++ItP, ++ItS )
        {
        if( ItS.Get() > 0 )
          {
          ItP.Set( ItP.Get() / ItS.Get() );
          }
        }
      return posteriorProbabilityImage;
      }
    }
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::RealImageType::Pointer
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GetDistancePriorProbabilityImageFromPriorLabelImage( unsigned int whichClass )
{
  if( this->m_InitializationStrategy != PriorLabelImage &&
      this->m_InitializationStrategy != PriorProbabilityImages )
    {
    return NULL;
    }

  if( whichClass > this->m_NumberOfClasses )
    {
    itkExceptionMacro(
      "Requested class is greater than the number of classes." );
    }

  /**
   * If memory minimization is turned off and if the distance prior probability
   * images have already been calculated, simply return the probability
   * image for the requested class.  Otherwise, calculate the probability
   * image.
   */
  if( whichClass <= this->m_DistancePriorProbabilityImages.size() )
    {
    return this->m_DistancePriorProbabilityImages[whichClass - 1];
    }
  else
    {
    /**
     * Here we assume that the calling function is invoked in order such
     * that GetDistancePriorImage( 1 ) is called before
     * GetDistancePriorImage( 2 ), etc.  As such, when this part of
     * the code is reached and the class requested is '1', we assume that
     * the sum of the distance prior probability images needs to be calculated
     * for normalization purposes.  This sum is then saved for subsequent calls.
     */
    typename RealImageType::Pointer distancePriorProbabilityImage = NULL;

    /**
     * Calculate the sum of the distance probability images.  Also, store the
     * distance probability images if m_MinimizeMemoryUsage == false.
     */
    if( whichClass == 1 )
      {
      this->m_SumDistancePriorProbabilityImage = RealImageType::New();
      this->m_SumDistancePriorProbabilityImage->SetRegions(
        this->GetOutput()->GetRequestedRegion() );
      this->m_SumDistancePriorProbabilityImage->SetOrigin(
        this->GetOutput()->GetOrigin() );
      this->m_SumDistancePriorProbabilityImage->SetSpacing(
        this->GetOutput()->GetSpacing() );
      this->m_SumDistancePriorProbabilityImage->SetDirection(
        this->GetOutput()->GetDirection() );
      this->m_SumDistancePriorProbabilityImage->Allocate();
      this->m_SumDistancePriorProbabilityImage->FillBuffer( 0 );
      for( unsigned int c = 0; c < this->m_NumberOfClasses; c++ )
        {
        typedef BinaryThresholdImageFilter<ClassifiedImageType, RealImageType>
          ThresholderType;
        typename ThresholderType::Pointer thresholder = ThresholderType::New();
        if( this->m_InitializationStrategy == PriorLabelImage )
          {
          thresholder->SetInput( const_cast<ClassifiedImageType *>(
                                   this->GetPriorLabelImage() ) );
          }
        else
          {
          thresholder->SetInput( this->GetOutput() );
          }
        thresholder->SetInsideValue( 1 );
        thresholder->SetOutsideValue( 0 );
        thresholder->SetLowerThreshold( static_cast<LabelType>( c + 1 ) );
        thresholder->SetUpperThreshold( static_cast<LabelType>( c + 1 ) );
        thresholder->Update();

        typename RealImageType::Pointer distanceImage = RealImageType::New();

        if( this->m_UseEuclideanDistanceForPriorLabels )
          {
          typedef SignedMaurerDistanceMapImageFilter
            <RealImageType, RealImageType> DistancerType;
          typename DistancerType::Pointer distancer = DistancerType::New();
          distancer->SetInput( thresholder->GetOutput() );
          distancer->SetSquaredDistance( true );
          distancer->SetUseImageSpacing( true );
          distancer->SetInsideIsPositive( false );
          distancer->Update();

          distanceImage = distancer->GetOutput();
          }
        else
          {
          typedef BinaryContourImageFilter<RealImageType, RealImageType>
            ContourFilterType;
          typename ContourFilterType::Pointer contour = ContourFilterType::New();
          contour->SetInput( thresholder->GetOutput() );
          contour->FullyConnectedOff();
          contour->SetBackgroundValue( 0 );
          contour->SetForegroundValue( 1 );
          contour->Update();

          typedef FastMarchingImageFilter<RealImageType, RealImageType>
            FastMarchingFilterType;
          typename FastMarchingFilterType::Pointer fastMarching
            = FastMarchingFilterType::New();

          typedef CastImageFilter<MaskImageType, RealImageType> CasterType;
          typename CasterType::Pointer caster = CasterType::New();
          if( this->GetMaskImage() )
            {
            caster->SetInput( const_cast<MaskImageType *>( this->GetMaskImage() ) );
            caster->Update();
            fastMarching->SetInput( caster->GetOutput() );
            }
          else
            {
            fastMarching->SetSpeedConstant( 1.0 );
            fastMarching->SetOverrideOutputInformation( true );
            fastMarching->SetOutputOrigin( this->GetOutput()->GetOrigin() );
            fastMarching->SetOutputSpacing( this->GetOutput()->GetSpacing() );
            fastMarching->SetOutputRegion( this->GetOutput()->GetRequestedRegion() );
            fastMarching->SetOutputDirection( this->GetOutput()->GetDirection() );
            }

          typedef typename FastMarchingFilterType::NodeContainer NodeContainer;
          typedef typename FastMarchingFilterType::NodeType      NodeType;
          typename NodeContainer::Pointer trialPoints = NodeContainer::New();
          trialPoints->Initialize();

          unsigned long trialCount = 0;

          ImageRegionIteratorWithIndex<RealImageType> ItC(
            contour->GetOutput(), contour->GetOutput()->GetRequestedRegion() );
          for( ItC.GoToBegin(); !ItC.IsAtEnd(); ++ItC )
            {
            if( ItC.Get() == contour->GetForegroundValue() )
              {
              NodeType node;
              node.SetValue( 0.0 );
              node.SetIndex( ItC.GetIndex() );
              trialPoints->InsertElement( trialCount++, node );
              }
            }
          fastMarching->SetTrialPoints( trialPoints );
          fastMarching->SetStoppingValue( NumericTraits<RealType>::max() );
//           fastMarching->SetTopologyCheck( FastMarchingFilterType::None );
          fastMarching->Update();

          ImageRegionIterator<RealImageType> ItT( thresholder->GetOutput(),
                                                  thresholder->GetOutput()->GetRequestedRegion() );
          ImageRegionIterator<RealImageType> ItF( fastMarching->GetOutput(),
                                                  fastMarching->GetOutput()->GetRequestedRegion() );
          for( ItT.GoToBegin(), ItF.GoToBegin(); !ItT.IsAtEnd(); ++ItT, ++ItF )
            {
            RealType distance = ItF.Get();
            ItF.Set( distance * distance );
            if( ItT.Get() == 1 )
              {
              ItF.Set( -ItF.Get() );
              }
            }

          distanceImage = fastMarching->GetOutput();
          }

        RealType maximumInteriorDistance = 0.0;

        ImageRegionIterator<RealImageType> ItD( distanceImage,
                                                distanceImage->GetRequestedRegion() );
        for( ItD.GoToBegin(); !ItD.IsAtEnd(); ++ItD )
          {
          if( ItD.Get() < 0 &&
              maximumInteriorDistance < vnl_math_abs( ItD.Get() ) )
            {
            maximumInteriorDistance = vnl_math_abs( ItD.Get() );
            }
          }

        RealType labelSigma = 0.0;
        RealType labelBoundaryProbability = 1.0;

        typename LabelParameterMapType::iterator it =
          this->m_PriorLabelParameterMap.find( c + 1 );
        if( it != this->m_PriorLabelParameterMap.end() )
          {
          labelSigma = ( it->second ).first;
          labelBoundaryProbability = ( it->second ).second;
          }
        for( ItD.GoToBegin(); !ItD.IsAtEnd(); ++ItD )
          {
          if( labelSigma == 0 )
            {
            if( ItD.Get() <= 0 )
              {
              ItD.Set( labelBoundaryProbability );
              }
            else
              {
              ItD.Set( 0.0 );
              }
            }
          else if( ItD.Get() >= 0 )
            {
            ItD.Set( labelBoundaryProbability
                     * vcl_exp( -ItD.Get() / vnl_math_sqr( labelSigma ) ) );
            }
          else if( ItD.Get() < 0 )
            {
            ItD.Set( 1.0 - ( 1.0 - labelBoundaryProbability )
                     * ( maximumInteriorDistance - vnl_math_abs( ItD.Get() ) )
                     / ( maximumInteriorDistance ) );
            }
          }

        typedef AddImageFilter<RealImageType, RealImageType, RealImageType>
          AdderType;
        typename AdderType::Pointer adder = AdderType::New();
        adder->SetInput1( this->m_SumDistancePriorProbabilityImage );
        adder->SetInput2( distanceImage );
        adder->Update();

        this->m_SumDistancePriorProbabilityImage = adder->GetOutput();

        if( ( c == 0 ) && this->m_MinimizeMemoryUsage )
          {
          distancePriorProbabilityImage = distanceImage;
          }
        if( !this->m_MinimizeMemoryUsage )
          {
          this->m_DistancePriorProbabilityImages.push_back(
            distanceImage );
          }
        }

      /**
       * Normalize the distance prior probability image(s).
       */
      ImageRegionIterator<RealImageType> ItS(
        this->m_SumDistancePriorProbabilityImage,
        this->m_SumDistancePriorProbabilityImage->GetRequestedRegion() );
      if( this->m_MinimizeMemoryUsage )
        {
        ImageRegionIterator<RealImageType> ItD( distancePriorProbabilityImage,
                                                distancePriorProbabilityImage->GetRequestedRegion() );
        for( ItD.GoToBegin(), ItS.GoToBegin(); !ItS.IsAtEnd(); ++ItD, ++ItS )
          {
          ItD.Set( ItD.Get() - ( ItS.Get() - ItD.Get() ) );
          if( ItD.Get() < 0 )
            {
            ItD.Set( 0 );
            }
          }
        return distancePriorProbabilityImage;
        }
      else
        {
        for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
          {
          ImageRegionIterator<RealImageType> ItD(
            this->m_DistancePriorProbabilityImages[n],
            this->m_DistancePriorProbabilityImages[n]->GetRequestedRegion() );
          for( ItD.GoToBegin(), ItS.GoToBegin(); !ItS.IsAtEnd(); ++ItD, ++ItS )
            {
            ItD.Set( ItD.Get() - ( ItS.Get() - ItD.Get() ) );
            if( ItD.Get() < 0 )
              {
              ItD.Set( 0 );
              }
            }
          }
        return this->m_DistancePriorProbabilityImages[0];
        }
      }
    else // whichClass > 1
      {
      typedef BinaryThresholdImageFilter<ClassifiedImageType, RealImageType>
        ThresholderType;
      typename ThresholderType::Pointer thresholder = ThresholderType::New();
      if( this->m_InitializationStrategy == PriorLabelImage )
        {
        thresholder->SetInput( const_cast<ClassifiedImageType *>(
                                 this->GetPriorLabelImage() ) );
        }
      else
        {
        thresholder->SetInput( this->GetOutput() );
        }
      thresholder->SetInsideValue( 1 );
      thresholder->SetOutsideValue( 0 );
      thresholder->SetLowerThreshold( static_cast<LabelType>( whichClass ) );
      thresholder->SetUpperThreshold( static_cast<LabelType>( whichClass ) );
      thresholder->Update();

      typename RealImageType::Pointer distanceImage = RealImageType::New();

      if( this->m_UseEuclideanDistanceForPriorLabels )
        {
        typedef SignedMaurerDistanceMapImageFilter
          <RealImageType, RealImageType> DistancerType;
        typename DistancerType::Pointer distancer = DistancerType::New();
        distancer->SetInput( thresholder->GetOutput() );
        distancer->SetSquaredDistance( true );
        distancer->SetUseImageSpacing( true );
        distancer->SetInsideIsPositive( false );
        distancer->Update();

        distanceImage = distancer->GetOutput();
        }
      else
        {
        typedef BinaryContourImageFilter<RealImageType, RealImageType>
          ContourFilterType;
        typename ContourFilterType::Pointer contour = ContourFilterType::New();
        contour->SetInput( thresholder->GetOutput() );
        contour->FullyConnectedOff();
        contour->SetBackgroundValue( 0 );
        contour->SetForegroundValue( 1 );
        contour->Update();

        typedef FastMarchingImageFilter<RealImageType, RealImageType>
          FastMarchingFilterType;
        typename FastMarchingFilterType::Pointer fastMarching
          = FastMarchingFilterType::New();

        typedef CastImageFilter<MaskImageType, RealImageType> CasterType;
        typename CasterType::Pointer caster = CasterType::New();
        if( this->GetMaskImage() )
          {
          caster->SetInput( const_cast<MaskImageType *>( this->GetMaskImage() ) );
          caster->Update();
          fastMarching->SetInput( caster->GetOutput() );
          }
        else
          {
          fastMarching->SetSpeedConstant( 1.0 );
          fastMarching->SetOverrideOutputInformation( true );
          fastMarching->SetOutputOrigin( this->GetOutput()->GetOrigin() );
          fastMarching->SetOutputSpacing( this->GetOutput()->GetSpacing() );
          fastMarching->SetOutputRegion( this->GetOutput()->GetRequestedRegion() );
          fastMarching->SetOutputDirection( this->GetOutput()->GetDirection() );
          }

        typedef typename FastMarchingFilterType::NodeContainer NodeContainer;
        typedef typename FastMarchingFilterType::NodeType      NodeType;
        typename NodeContainer::Pointer trialPoints = NodeContainer::New();
        trialPoints->Initialize();

        unsigned long trialCount = 0;

        ImageRegionIteratorWithIndex<RealImageType> ItC(
          contour->GetOutput(), contour->GetOutput()->GetRequestedRegion() );
        for( ItC.GoToBegin(); !ItC.IsAtEnd(); ++ItC )
          {
          if( ItC.Get() == contour->GetForegroundValue() )
            {
            NodeType node;
            node.SetValue( 0.0 );
            node.SetIndex( ItC.GetIndex() );
            trialPoints->InsertElement( trialCount++, node );
            }
          }
        fastMarching->SetTrialPoints( trialPoints );
        fastMarching->SetStoppingValue( NumericTraits<RealType>::max() );
//           fastMarching->SetTopologyCheck( FastMarchingFilterType::None );
        fastMarching->Update();

        ImageRegionIterator<RealImageType> ItT( thresholder->GetOutput(),
                                                thresholder->GetOutput()->GetRequestedRegion() );
        ImageRegionIterator<RealImageType> ItF( fastMarching->GetOutput(),
                                                fastMarching->GetOutput()->GetRequestedRegion() );
        for( ItT.GoToBegin(), ItF.GoToBegin(); !ItT.IsAtEnd(); ++ItT, ++ItF )
          {
          RealType distance = ItF.Get();
          ItF.Set( distance * distance );
          if( ItT.Get() == 1 )
            {
            ItF.Set( -ItF.Get() );
            }
          }
        distanceImage = fastMarching->GetOutput();
        }

      distancePriorProbabilityImage = distanceImage;

      RealType maximumInteriorDistance = 0.0;

      ImageRegionIterator<RealImageType> ItD( distancePriorProbabilityImage,
                                              distancePriorProbabilityImage->GetRequestedRegion() );
      for( ItD.GoToBegin(); !ItD.IsAtEnd(); ++ItD )
        {
        if( ItD.Get() < 0 &&
            maximumInteriorDistance < vnl_math_abs( ItD.Get() ) )
          {
          maximumInteriorDistance = vnl_math_abs( ItD.Get() );
          }
        }

      RealType labelSigma = 0.0;
      RealType labelBoundaryProbability = 1.0;

      typename LabelParameterMapType::iterator it =
        this->m_PriorLabelParameterMap.find( whichClass );
      if( it != this->m_PriorLabelParameterMap.end() )
        {
        labelSigma = ( it->second ).first;
        labelBoundaryProbability = ( it->second ).second;
        }
      for( ItD.GoToBegin(); !ItD.IsAtEnd(); ++ItD )
        {
        if( labelSigma == 0 )
          {
          if( ItD.Get() <= 0 )
            {
            ItD.Set( labelBoundaryProbability );
            }
          else
            {
            ItD.Set( 0.0 );
            }
          }
        else if( ItD.Get() >= 0 )
          {
          ItD.Set( labelBoundaryProbability
                   * vcl_exp( -ItD.Get() / vnl_math_sqr( labelSigma ) ) );
          }
        else if( ItD.Get() < 0 )
          {
          ItD.Set( 1.0 - ( 1.0 - labelBoundaryProbability )
                   * ( maximumInteriorDistance - vnl_math_abs( ItD.Get() ) )
                   / ( maximumInteriorDistance ) );
          }
        }

      /**
       * Normalize the distance prior probability image(s).
       */
      ImageRegionIterator<RealImageType> ItS(
        this->m_SumDistancePriorProbabilityImage,
        this->m_SumDistancePriorProbabilityImage->GetRequestedRegion() );
      for( ItD.GoToBegin(), ItS.GoToBegin(); !ItS.IsAtEnd(); ++ItD, ++ItS )
        {
        ItD.Set( ItD.Get() - ( ItS.Get() - ItD.Get() ) );
        if( ItD.Get() < 0 )
          {
          ItD.Set( 0 );
          }
        }
      return distancePriorProbabilityImage;
      }
    }
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::RealImageType::Pointer
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::CalculateSmoothIntensityImageFromPriorProbabilityImage( unsigned int whichImage,
                                                          unsigned int whichClass )
{
  typename ScalarImageType::Pointer bsplineImage;

  if( this->m_ControlPointLattices[whichImage][whichClass - 1].GetPointer()
      != NULL )
    {
    typedef BSplineControlPointImageFilter<ControlPointLatticeType,
                                           ScalarImageType> BSplineReconstructorType;
    typename BSplineReconstructorType::Pointer bspliner
      = BSplineReconstructorType::New();

    bspliner->SetInput( this->m_ControlPointLattices[whichImage][whichClass - 1] );
    bspliner->SetSize( this->GetInput()->GetRequestedRegion().GetSize() );
    bspliner->SetSpacing( this->GetInput()->GetSpacing() );
    bspliner->SetOrigin( this->GetInput()->GetOrigin() );
    bspliner->SetDirection( this->GetInput()->GetDirection() );
    bspliner->SetSplineOrder( this->m_SplineOrder );
    bspliner->Update();

    bsplineImage = bspliner->GetOutput();
    }
  else
    {
    typename PointSetType::Pointer points = PointSetType::New();
    points->Initialize();

    typedef typename BSplineFilterType::WeightsContainerType WeightsType;
    typename WeightsType::Pointer weights = WeightsType::New();
    weights->Initialize();

    typename RealImageType::Pointer probabilityImage;
    if( this->m_InitializationStrategy == PriorProbabilityImages )
      {
      probabilityImage = this->GetPriorProbabilityImage( whichClass );
      }
    else
      {
      typedef BinaryThresholdImageFilter<ClassifiedImageType, RealImageType>
        ThresholderType;
      typename ThresholderType::Pointer thresholder = ThresholderType::New();
      thresholder->SetInput( const_cast<ClassifiedImageType *>(
                               this->GetPriorLabelImage() ) );
      thresholder->SetInsideValue( 1 );
      thresholder->SetOutsideValue( 0 );
      thresholder->SetLowerThreshold( static_cast<LabelType>( whichClass ) );
      thresholder->SetUpperThreshold( static_cast<LabelType>( whichClass ) );
      thresholder->Update();

      probabilityImage = thresholder->GetOutput();
      }

    typename RealImageType::DirectionType originalDirection
      = probabilityImage->GetDirection();
    typename RealImageType::DirectionType identity;
    identity.SetIdentity();
    probabilityImage->SetDirection( identity );

    unsigned long count = 0;

    ImageRegionIteratorWithIndex<RealImageType> ItP( probabilityImage,
                                                     probabilityImage->GetBufferedRegion() );
    for( ItP.GoToBegin(); !ItP.IsAtEnd(); ++ItP )
      {
      if( !this->GetMaskImage() ||
          this->GetMaskImage()->GetPixel( ItP.GetIndex() ) == this->m_MaskLabel )
        {
        if( ItP.Get() >= 0.5 )
          {
          typename RealImageType::PointType imagePoint;
          probabilityImage->TransformIndexToPhysicalPoint(
            ItP.GetIndex(), imagePoint );

          typename PointSetType::PointType bsplinePoint;
          bsplinePoint.CastFrom( imagePoint );

          ScalarType intensity;

          intensity[0] = this->GetIntensityImage(
              whichImage )->GetPixel( ItP.GetIndex() );

          points->SetPoint( count, bsplinePoint );
          points->SetPointData( count, intensity );
          weights->InsertElement( count, ItP.Get() );

          count++;
          }
        }
      }
    probabilityImage->SetDirection( originalDirection );

    typename BSplineFilterType::ArrayType numberOfControlPoints;
    typename BSplineFilterType::ArrayType numberOfLevels;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      numberOfControlPoints[d] = this->m_NumberOfControlPoints[d];
      numberOfLevels[d] = this->m_NumberOfLevels[d];
      }

    typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
    bspliner->SetInput( points );
    bspliner->SetPointWeights( weights );
    bspliner->SetNumberOfLevels( numberOfLevels );
    bspliner->SetSplineOrder( this->m_SplineOrder );
    bspliner->SetNumberOfControlPoints( numberOfControlPoints );
    bspliner->SetSize( this->GetOutput()->GetLargestPossibleRegion().GetSize() );
    bspliner->SetOrigin( this->GetOutput()->GetOrigin() );
    bspliner->SetDirection( this->GetOutput()->GetDirection() );
    bspliner->SetSpacing( this->GetOutput()->GetSpacing() );
    bspliner->SetGenerateOutputImage( true );
    bspliner->Update();

    bsplineImage = bspliner->GetOutput();

    this->m_ControlPointLattices[whichImage][whichClass - 1] =
      bspliner->GetPhiLattice();
    }

  typedef VectorIndexSelectionCastImageFilter
    <ScalarImageType, RealImageType> CasterType;
  typename CasterType::Pointer caster = CasterType::New();
  caster->SetInput( bsplineImage );
  caster->SetIndex( 0 );
  caster->Update();

  return caster->GetOutput();
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Maximum number of iterations: "
     << this->m_MaximumNumberOfIterations << std::endl;
  os << indent << "Convergence threshold: "
     << this->m_ConvergenceThreshold << std::endl;
  os << indent << "Mask label: "
     << static_cast<typename NumericTraits<LabelType>::PrintType>
    ( this->m_MaskLabel ) << std::endl;
  os << indent << "Number of classes: "
     << this->m_NumberOfClasses << std::endl;
  os << indent << "Minimize memory usage:";
  if( this->m_MinimizeMemoryUsage )
    {
    os << " true";
    if( this->m_MinimizeMemoryUsage &&
        this->m_InitializationStrategy == PriorProbabilityImages )
      {
      os << " (prior probability threshold = "
         << this->m_PriorProbabilityThreshold << ")" << std::endl;
      }
    os << std::endl;
    }
  else
    {
    os << " false" << std::endl;
    }

  os << indent << "Initialization strategy: ";

  switch( this->m_InitializationStrategy )
    {
    case Random:
      {
      os << "Random" << std::endl;
      break;
      }
    case KMeans:
      {
      os << "K means clustering" << std::endl;
      break;
      }
    case Otsu:
      {
      os << "Otsu thresholding" << std::endl;
      break;
      }
    case PriorProbabilityImages:
      {
      os << "Prior probability images" << std::endl;
      os << indent << "  Use Euclidean distance for prior labels:";
      if( this->m_UseEuclideanDistanceForPriorLabels )
        {
        os << " true" << std::endl;
        }
      else
        {
        os << " false" << std::endl;
        }
      if( this->m_PriorLabelParameterMap.size() > 0 )
        {
        os << indent << "  Specified prior label parameters:" << std::endl;
        typename LabelParameterMapType::const_iterator it;
        for( it = this->m_PriorLabelParameterMap.begin(); it !=
             this->m_PriorLabelParameterMap.end(); ++it )
          {
          RealType label = it->first;
          RealType sigma = ( it->second ).first;
          RealType boundaryProbability = ( it->second ).second;
          os << indent << "    Class " << label
             << ": sigma = " << sigma
             << ", boundary probability = " << boundaryProbability << std::endl;
          }
        }
      break;
      }
    case PriorLabelImage:
      {
      os << "Prior label image" << std::endl;
      os << indent << "  Use Euclidean distance for prior labels:";
      if( this->m_UseEuclideanDistanceForPriorLabels )
        {
        os << " true" << std::endl;
        }
      else
        {
        os << " false" << std::endl;
        }
      os << indent << "  Specified prior label parameters:" << std::endl;
      typename LabelParameterMapType::const_iterator it;
      for( it = this->m_PriorLabelParameterMap.begin(); it !=
           this->m_PriorLabelParameterMap.end(); ++it )
        {
        RealType label = it->first;
        RealType sigma = ( it->second ).first;
        RealType boundaryProbability = ( it->second ).second;
        os << indent << "    Class " << label
           << ": sigma = " << sigma
           << ", boundary probability = " << boundaryProbability << std::endl;
        }
      break;
      }
    }

  os << indent << "MRF parameters" << std::endl;
  os << indent << "  MRF smoothing factor: "
     << this->m_MRFSmoothingFactor << std::endl;
  os << indent << "  MRF radius: "
     << this->m_MRFRadius << std::endl;

  if( this->m_OutlierHandlingFilter )
    {
    os << indent << "Outlier handling " << std::endl;
    this->m_OutlierHandlingFilter->Print( os, indent.GetNextIndent() );
    }
  else
    {
    os << indent << "No outlier handling." << std::endl;
    }

  if( ( this->m_InitializationStrategy == PriorProbabilityImages ||
        this->m_InitializationStrategy == PriorLabelImage ) &&
      this->m_AdaptiveSmoothingWeights.size() > 0 )
    {
    os << indent << "Adaptive smoothing weights: [";
    for( unsigned int i = 0; i < this->m_AdaptiveSmoothingWeights.size() - 1; i++ )
      {
      os << this->m_AdaptiveSmoothingWeights[i] << ", " << std::endl;
      }
    os << this->m_AdaptiveSmoothingWeights[
      this->m_AdaptiveSmoothingWeights.size() - 1] << "]" << std::endl;
    os << indent << "BSpline smoothing" << std::endl;
    os << indent << "  Spline order: "
       << this->m_SplineOrder << std::endl;
    os << indent << "  Number of levels: "
       << this->m_NumberOfLevels << std::endl;
    os << indent << "  Number of initial control points: "
       << this->m_NumberOfControlPoints << std::endl;
    }
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    os << indent << "Class " << n + 1 << ": proportion = "
       << this->m_MixtureModelProportions[n] << std::endl;
    this->m_MixtureModelComponents[n]->Print( os, indent.GetNextIndent() );
    }
}
} // namespace ants
} // namespace itk

#endif
