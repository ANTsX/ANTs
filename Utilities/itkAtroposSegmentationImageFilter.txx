/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAtroposSegmentationImageFilter.txx,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkAtroposSegmentationImageFilter_txx
#define __itkAtroposSegmentationImageFilter_txx

#include "itkAtroposSegmentationImageFilter.h"

#include "itkAddImageFilter.h"
#include "itkAddConstantToImageFilter.h"
#include "itkBinaryContourImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkFastMarchingImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageToListGenerator.h"
#include "itkIterationReporter.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkEuclideanDistance.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkMinimumDecisionRule.h"
#include "itkOtsuMultipleThresholdsCalculator.h"
#include "itkSampleClassifier.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkWeightedCentroidKdTreeGenerator.h"

#include "itkTimeProbe.h"

#include "vnl/vnl_vector.h"

namespace itk
{
template <class TInputImage, class TMaskImage, class TClassifiedImage>
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::AtroposSegmentationImageFilter()
{
  this->ProcessObject::SetNumberOfRequiredInputs( 1 );

  this->m_NumberOfClasses = 3;
  this->m_NumberOfAuxiliaryImages = 0;
  this->m_MaximumNumberOfIterations = 5;
  this->m_ElapsedIterations = 0;
  this->m_ConvergenceThreshold = 0.001;

  this->m_MaskLabel = NumericTraits<LabelType>::One;

  this->m_InitializationStrategy = Otsu;

  this->m_AdaptiveSmoothingWeights.clear();
  this->m_PriorLabelParameterMap.clear();

  this->m_MRFSmoothingFactor = 0.3;
  this->m_MRFSigmoidAlpha = 0.0;
  this->m_MRFSigmoidBeta = 0.0;
  this->m_MRFRadius.Fill( 1 );

  this->m_SplineOrder = 3;
  this->m_NumberOfLevels.Fill( 6 );
  this->m_NumberOfControlPoints.Fill( this->m_SplineOrder + 1 );

  this->m_MinimizeMemoryUsage = false;
  this->m_PosteriorProbabilityImages.clear();
  this->m_DistancePriorProbabilityImages.clear();
  this->m_UseEuclideanDistanceForPriorLabels = false;
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
  unsigned int whichClass, const RealImageType * prior )
{
  if( whichClass < 1 || whichClass > this->m_NumberOfClasses )
    {
    itkExceptionMacro(
      "The prior probability images are inputs 3...3+m_NumberOfClasses-1.  "
      << "The requested image should be in the range [1, m_NumberOfClasses]" )
    }

  this->SetNthInput( 2 + whichClass, const_cast<RealImageType *>( prior ) );
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
const typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::RealImageType
* AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GetPriorProbabilityImage( unsigned int whichClass ) const
  {
  if( whichClass < 1 || whichClass > this->m_NumberOfClasses )
    {
    itkExceptionMacro(
      "The prior probability images are inputs 3...3+m_NumberOfClasses-1.  "
      << "The requested image should be in the range [1, m_NumberOfClasses]" )
    }

  const RealImageType *priorImage =
    dynamic_cast<const RealImageType *>(
      this->ProcessObject::GetInput( 2 + whichClass ) );

  return priorImage;
  }

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::SetAuxiliaryImage( unsigned int which, const ImageType * image )
{
  if( which > this->m_NumberOfAuxiliaryImages )
    {
    this->m_NumberOfAuxiliaryImages = which;
    }
  this->SetNthInput( 2 + this->m_NumberOfClasses + which,
                     const_cast<ImageType *>( image ) );
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
const typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::ImageType
* AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GetAuxiliaryImage( unsigned int which ) const
  {
  if( which < 1 || which > this->m_NumberOfAuxiliaryImages )
    {
    itkExceptionMacro( "The auxiliary images are inputs "
                       << "1+m_NumberOfClasses...1+m_NumberOfClasses+m_NumberOfAuxiliaryImages-1." )
    }

  const ImageType *image = dynamic_cast<const ImageType *>(
      this->ProcessObject::GetInput( 2 + this->m_NumberOfClasses + which ) );

  return image;
  }

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GenerateData()
{
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
    probabilityNew = this->UpdateClassParametersAndLabeling();
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
      typedef Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
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
    case KMeans:
      {
      this->GenerateInitialClassLabelingWithKMeansClustering();
      break;
      }
    case Otsu:  default:
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
        this->NormalizePriorProbabilityImages();
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
   * Calculate the initial parameters of the Gaussian mixture model from the
   * initial labeling, i.e. the proportion, mean, and covariance for each label.
   */
  this->m_GaussianMixtureModel.clear();
  this->m_GaussianMixtureModelProportions.SetSize( this->m_NumberOfClasses );

  unsigned int totalSampleSize = 0;

  std::vector<typename SampleType::Pointer> samples;
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    typename SampleType::Pointer sample = SampleType::New();
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
    measurement.SetSize( this->m_NumberOfAuxiliaryImages + 1 );
    for( unsigned int i = 0; i < this->m_NumberOfAuxiliaryImages + 1; i++ )
      {
      if( i == 0 )
        {
        measurement[i] = this->GetInput()->GetPixel( ItO.GetIndex() );
        }
      else
        {
        measurement[i] =
          this->GetAuxiliaryImage( i )->GetPixel( ItO.GetIndex() );
        }
      }
    samples[label - 1]->PushBack( measurement );
    }

  /**
   * If using prior probability images, create the weight array now that we
   * know the sample sizes.
   */
  Array<unsigned int> count( this->m_NumberOfClasses );
  count.Fill( 0 );
  std::vector<typename MeanCalculatorType::WeightArrayType> weights;
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    totalSampleSize += samples[n]->Size();
    typename MeanCalculatorType::WeightArrayType weightArray;
    weightArray.SetSize( samples[n]->Size() );
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
      weights[label - 1].SetElement( count[label - 1]++,
                                     this->GetPriorProbabilityImage( label )->GetPixel( ItO.GetIndex() ) );
      }
    }
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    typename GaussianType::Pointer gaussian = GaussianType::New();
    if( samples[n]->Size() > 0 )
      {
      typename MeanCalculatorType::Pointer meanCalculator =
        MeanCalculatorType::New();
      meanCalculator->SetWeights( &weights[n] );
      meanCalculator->SetInputSample( samples[n] );
      meanCalculator->Update();

      typename CovarianceCalculatorType::Pointer covarianceCalculator =
        CovarianceCalculatorType::New();
      covarianceCalculator->SetWeights( &weights[n] );
      covarianceCalculator->SetMean( meanCalculator->GetOutput() );
      covarianceCalculator->SetInputSample( samples[n] );
      covarianceCalculator->Update();

      gaussian->SetMean( *meanCalculator->GetOutput() );
      gaussian->SetCovariance( *covarianceCalculator->GetOutput() );
      }
    else
      {
      typename GaussianType::MeanType mean;
      mean.SetSize( this->m_NumberOfAuxiliaryImages + 1 );
      mean.Fill( 0.0 );
      typename GaussianType::CovarianceType covariance(
        this->m_NumberOfAuxiliaryImages + 1, this->m_NumberOfAuxiliaryImages + 1 );
      covariance.Fill( 0.0 );

      gaussian->SetMean( mean );
      gaussian->SetCovariance( covariance );
      }
    this->m_GaussianMixtureModel.push_back( gaussian );
    this->m_GaussianMixtureModelProportions[n] =
      static_cast<RealType>( samples[n]->Size() )
      / static_cast<RealType>( totalSampleSize );
    }
  for( unsigned int i = 0; i < this->m_NumberOfAuxiliaryImages + 1; i++ )
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
::NormalizePriorProbabilityImages()
{
  ImageRegionConstIteratorWithIndex<ImageType> ItI( this->GetInput(),
                                                    this->GetInput()->GetRequestedRegion() );
  for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
    {
    vnl_vector<RealType> priorProbabilities( this->m_NumberOfClasses );
    for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
      {
      priorProbabilities[n] =
        this->GetPriorProbabilityImage( n + 1 )->GetPixel( ItI.GetIndex() );
      }
    /**
     * We are assuming that each voxel in the mask is to be labeled even
     * if the sum of prior probabilities is zero at that point.  To change this
     * assumption, comment out the following conditional.
     */
    if( priorProbabilities.sum() == 0.0 )
      {
      priorProbabilities.fill( 1.0 );
      }

    RealType                  maxValue = priorProbabilities.max_value();
    std::vector<unsigned int> argMax;
    for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
      {
      if( maxValue == priorProbabilities[n] )
        {
        argMax.push_back( n );
        }
      }
    unsigned int whichArgIsMax = 0;
    if( argMax.size() > 1 )
      {
      typedef Statistics::MersenneTwisterRandomVariateGenerator GeneratorType;
      typename GeneratorType::Pointer generator = GeneratorType::New();
      whichArgIsMax = generator->GetIntegerVariate(
          argMax.size() - 1 ) + 1;
      }
    priorProbabilities[argMax[whichArgIsMax]] += 1e-6;
    priorProbabilities /= priorProbabilities.sum();
    for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
      {
      typename RealImageType::Pointer priorProbabilityImage
        = const_cast<RealImageType *>( this->GetPriorProbabilityImage( n + 1 ) );
      priorProbabilityImage->SetPixel( ItI.GetIndex(), priorProbabilities[n] );
      }
    }
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GenerateInitialClassLabelingWithPriorProbabilityImages()
{
  this->GetOutput()->FillBuffer( NumericTraits<LabelType>::Zero );

  ImageRegionIteratorWithIndex<ClassifiedImageType> ItO( this->GetOutput(),
                                                         this->GetOutput()->GetRequestedRegion() );
  for( ItO.GoToBegin(); !ItO.IsAtEnd(); ++ItO )
    {
    if( !this->GetMaskImage() ||
        this->GetMaskImage()->GetPixel( ItO.GetIndex() )
        == this->m_MaskLabel )
      {
      vnl_vector<RealType> priorProbabilities( this->m_NumberOfClasses );
      for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
        {
        priorProbabilities[n] =
          this->GetPriorProbabilityImage( n + 1 )->GetPixel( ItO.GetIndex() );
        }
      RealType     maxValue = priorProbabilities.max_value();
      unsigned int argMax = this->m_NumberOfClasses;
      for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
        {
        if( maxValue == priorProbabilities[n] )
          {
          argMax = n;
          break;
          }
        }
      LabelType maxLabel = static_cast<LabelType>( argMax + 1 );
      ItO.Set( maxLabel );
      }
    }
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GenerateInitialClassLabelingWithOtsuThresholding()
{
  RealType maxValue = itk::NumericTraits<RealType>::min();
  RealType minValue = itk::NumericTraits<RealType>::max();

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
    this->GetOutput()->FillBuffer( this->m_MaskLabel );
    stats->SetLabelInput( this->GetOutput() );
    }
  stats->UseHistogramsOn();
  stats->SetHistogramParameters( 200, minValue, maxValue );
  stats->Update();

  typedef itk::OtsuMultipleThresholdsCalculator<typename StatsType::HistogramType>
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
   * Acquire the samples but rescale so that they share the same dynamic range.
   * This (is kind of a hack that?) allows us to use the EuclideanDistance
   * membership function.
   */
  Array<RealType> minValues;
  minValues.SetSize( this->m_NumberOfAuxiliaryImages + 1 );
  minValues.Fill( NumericTraits<RealType>::max() );
  Array<RealType> maxValues;
  maxValues.SetSize( this->m_NumberOfAuxiliaryImages + 1 );
  maxValues.Fill( NumericTraits<RealType>::NonpositiveMin() );

  ImageRegionIteratorWithIndex<ClassifiedImageType> ItO( this->GetOutput(),
                                                         this->GetOutput()->GetRequestedRegion() );
  for( ItO.GoToBegin(); !ItO.IsAtEnd(); ++ItO )
    {
    if( !this->GetMaskImage() || this->GetMaskImage()->GetPixel( ItO.GetIndex() )
        == this->m_MaskLabel )
      {
      for( unsigned int i = 0; i < this->m_NumberOfAuxiliaryImages + 1; i++ )
        {
        RealType value = 0.0;
        if( i == 0 )
          {
          value = this->GetInput()->GetPixel( ItO.GetIndex() );
          }
        else
          {
          value = this->GetAuxiliaryImage( i )->GetPixel( ItO.GetIndex() );
          }
        if( value < minValues[i] )
          {
          minValues[i] = value;
          }
        else if( value > maxValues[i] )
          {
          maxValues[i] = value;
          }
        }
      }
    }

  typename SampleType::Pointer sample = SampleType::New();
  for( ItO.GoToBegin(); !ItO.IsAtEnd(); ++ItO )
    {
    if( !this->GetMaskImage() ||
        this->GetMaskImage()->GetPixel( ItO.GetIndex() ) == this->m_MaskLabel )
      {
      typename SampleType::MeasurementVectorType measurement;
      measurement.SetSize( this->m_NumberOfAuxiliaryImages + 1 );
      for( unsigned int i = 0; i < this->m_NumberOfAuxiliaryImages + 1; i++ )
        {
        if( i == 0 )
          {
          measurement[i] = this->GetInput()->GetPixel( ItO.GetIndex() );
          }
        else
          {
          measurement[i] =
            this->GetAuxiliaryImage( i )->GetPixel( ItO.GetIndex() );
          }
        measurement[i] = minValues[0] + ( measurement[i] - minValues[i] )
          * ( maxValues[0] - minValues[0] ) / ( maxValues[i] - minValues[i] );
        }
      sample->PushBack( measurement );
      }
    }

  typedef Statistics::WeightedCentroidKdTreeGenerator<SampleType>
    TreeGeneratorType;
  typedef typename TreeGeneratorType::KdTreeType           TreeType;
  typedef Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;
  typedef typename EstimatorType::ParametersType           ParametersType;

  typename TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
  treeGenerator->SetSample( sample );
  treeGenerator->SetBucketSize( 16 );
  treeGenerator->Update();

  /**
   * Guess initial class means by dividing the dynamic range
   *  into equal intervals.
   */
  typename EstimatorType::Pointer estimator = EstimatorType::New();
  ParametersType initialMeans( this->m_NumberOfClasses
                               * ( this->m_NumberOfAuxiliaryImages + 1 ) );
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    for( unsigned int i = 0; i < this->m_NumberOfAuxiliaryImages + 1; i++ )
      {
      initialMeans[( this->m_NumberOfAuxiliaryImages + 1 ) * n + i] = minValues[0]
        + ( maxValues[0] - minValues[0] ) * static_cast<RealType>( n + 1 )
        / static_cast<RealType>( this->m_NumberOfClasses + 1 );
      }
    }
  estimator->SetParameters( initialMeans );
  estimator->SetKdTree( treeGenerator->GetOutput() );
  estimator->SetMaximumIteration( 200 );
  estimator->SetCentroidPositionChangesThreshold( 0.0 );
  estimator->StartOptimization();

  /**
   * Classify the samples
   */
  typedef Statistics::SampleClassifier<SampleType> ClassifierType;
  typedef MinimumDecisionRule                      DecisionRuleType;
  typename DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();
  typename ClassifierType::Pointer classifier = ClassifierType::New();

  classifier->SetDecisionRule( decisionRule.GetPointer() );
  classifier->SetSample( sample );
  classifier->SetNumberOfClasses( this->m_NumberOfClasses );

  typedef std::vector<unsigned int> ClassLabelVectorType;
  ClassLabelVectorType classLabels;
  classLabels.resize( this->m_NumberOfClasses );

  /**
   * Order the cluster means so that the lowest mean of the input image
   * corresponds to label '1', the second lowest to label '2', etc.
   */

  std::vector<OrderingPair> estimatorParameters;
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    RealType measurementInputImage = estimator->GetParameters()[
        ( this->m_NumberOfAuxiliaryImages + 1 ) * n];
    MeasurementVectorType measurement;
    measurement.SetSize( this->m_NumberOfAuxiliaryImages + 1 );
    for( unsigned int i = 0; i < this->m_NumberOfAuxiliaryImages + 1; i++ )
      {
      measurement[i] = estimator->GetParameters()[
          ( this->m_NumberOfAuxiliaryImages + 1 ) * n + i];
      }
    OrderingPair measurementPair( measurementInputImage, measurement );
    estimatorParameters.push_back( measurementPair );
    }
  std::sort( estimatorParameters.begin(), estimatorParameters.end(),
             SortOrderedPair() );

  typedef itk::Statistics::EuclideanDistance<MeasurementVectorType>
    MembershipFunctionType;
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    classLabels[n] = n + 1;
    typename MembershipFunctionType::Pointer
    membershipFunction = MembershipFunctionType::New();
    typename MembershipFunctionType::OriginType origin(
      sample->GetMeasurementVectorSize() );
    for( unsigned int i = 0; i < this->m_NumberOfAuxiliaryImages + 1; i++ )
      {
      origin[i] = estimatorParameters[n].second[i];
      }
    membershipFunction->SetOrigin( origin );
    classifier->AddMembershipFunction( membershipFunction.GetPointer() );
    }
  classifier->SetMembershipFunctionClassLabels( classLabels );
  classifier->Update();

  /**
   * Classify the voxels
   */
  typedef typename ClassifierType::OutputType          ClassifierOutputType;
  typedef typename ClassifierOutputType::ConstIterator LabelIterator;

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

template <class TInputImage, class TMaskImage, class TClassifiedImage>
typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::RealType
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::UpdateClassParametersAndLabeling()
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
      if( !this->GetMaskImage() || this->GetMaskImage()->GetPixel(
            ItO.GetIndex() ) == this->m_MaskLabel )
        {
        RealType posteriorProbability = ItP.Get();
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
    typename RealImageType::Pointer distancePriorProbabilityImage = NULL;
    typename RealImageType::ConstPointer priorProbabilityImage = NULL;

    if( this->m_InitializationStrategy == PriorProbabilityImages )
      {
      priorProbabilityImage = const_cast<RealImageType *>(
          this->GetPriorProbabilityImage( n + 1 ) );
      }
    else if( this->m_InitializationStrategy == PriorLabelImage )
      {
      distancePriorProbabilityImage
        = this->GetDistancePriorProbabilityImageFromPriorLabelImage( n + 1 );
      }

    ImageRegionIteratorWithIndex<RealImageType> ItW(
      weightedPriorProbabilityImage,
      weightedPriorProbabilityImage->GetRequestedRegion() );
    for( ItW.GoToBegin(); !ItW.IsAtEnd(); ++ItW )
      {
      if( !this->GetMaskImage() || this->GetMaskImage()->GetPixel(
            ItW.GetIndex() ) == this->m_MaskLabel )
        {
        RealType priorProbability = 1.0;
        if( priorProbabilityImage )
          {
          priorProbability = priorProbabilityImage->GetPixel( ItW.GetIndex() );
          }
        else if( distancePriorProbabilityImage )
          {
          priorProbability =
            distancePriorProbabilityImage->GetPixel( ItW.GetIndex() );
          }
        ItW.Set( ItW.Get() + this->m_GaussianMixtureModelProportions[n]
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

    typename RealImageType::Pointer distancePriorProbabilityImage = NULL;
    typename RealImageType::ConstPointer priorProbabilityImage = NULL;

    if( this->m_InitializationStrategy == PriorProbabilityImages )
      {
      priorProbabilityImage = const_cast<RealImageType *>(
          this->GetPriorProbabilityImage( n + 1 ) );
      }
    else if( this->m_InitializationStrategy == PriorLabelImage )
      {
      distancePriorProbabilityImage
        = this->GetDistancePriorProbabilityImageFromPriorLabelImage( n + 1 );
      }

    ImageRegionIteratorWithIndex<RealImageType> ItW(
      weightedPriorProbabilityImage,
      weightedPriorProbabilityImage->GetRequestedRegion() );
    for( ItW.GoToBegin(); !ItW.IsAtEnd(); ++ItW )
      {
      if( !this->GetMaskImage() || this->GetMaskImage()->GetPixel(
            ItW.GetIndex() ) == this->m_MaskLabel )
        {
        RealType priorProbability = 1.0;
        if( priorProbabilityImage )
          {
          priorProbability = priorProbabilityImage->GetPixel( ItW.GetIndex() );
          }
        else if( distancePriorProbabilityImage )
          {
          priorProbability =
            distancePriorProbabilityImage->GetPixel( ItW.GetIndex() );
          }
        if( ItW.Get() > 0.0 )
          {
          denominator += ( priorProbability / ItW.Get() );
          }
        }
      }
    if( denominator > 0.0 )
      {
      this->m_GaussianMixtureModelProportions[n] =
        sumPosteriors[n] / denominator;
      }
    else
      {
      this->m_GaussianMixtureModelProportions[n] = 0.0;
      }
    }

  /**
   * Update the class means and variances
   */

  /**
   * Calculate the initial parameters of the Gaussian mixture model from the
   * initial labeling, i.e. the proportion, mean, and covariance for each label.
   */
  this->m_GaussianMixtureModel.clear();
  this->m_GaussianMixtureModelProportions.SetSize( this->m_NumberOfClasses );

  unsigned int totalSampleSize = 0;

  std::vector<typename SampleType::Pointer> samples;
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    typename SampleType::Pointer sample = SampleType::New();
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
    measurement.SetSize( this->m_NumberOfAuxiliaryImages + 1 );
    for( unsigned int i = 0; i < this->m_NumberOfAuxiliaryImages + 1; i++ )
      {
      if( i == 0 )
        {
        measurement[i] = this->GetInput()->GetPixel( ItO.GetIndex() );
        }
      else
        {
        measurement[i] =
          this->GetAuxiliaryImage( i )->GetPixel( ItO.GetIndex() );
        }
      }
    samples[label - 1]->PushBack( measurement );
    }

  Array<unsigned int> count( this->m_NumberOfClasses );
  count.Fill( 0 );
  std::vector<typename MeanCalculatorType::WeightArrayType> weights;
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    totalSampleSize += samples[n]->Size();
    typename MeanCalculatorType::WeightArrayType weightArray;
    weightArray.SetSize( samples[n]->Size() );
    weightArray.Fill( 1.0 );
    weights.push_back( weightArray );
    }
  for( ItO.GoToBegin(); !ItO.IsAtEnd(); ++ItO )
    {
    LabelType label = ItO.Get();
    if( label == 0 )
      {
      continue;
      }
    weights[label - 1].SetElement( count[label - 1]++,
                                   maxPosteriorProbabilityImage->GetPixel( ItO.GetIndex() ) );
    }
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    typename GaussianType::Pointer gaussian = GaussianType::New();
    if( samples[n]->Size() > 0 )
      {
      typename MeanCalculatorType::Pointer meanCalculator =
        MeanCalculatorType::New();
      meanCalculator->SetWeights( &weights[n] );
      meanCalculator->SetInputSample( samples[n] );
      meanCalculator->Update();

      typename CovarianceCalculatorType::Pointer covarianceCalculator =
        CovarianceCalculatorType::New();
      covarianceCalculator->SetWeights( &weights[n] );
      covarianceCalculator->SetMean( meanCalculator->GetOutput() );
      covarianceCalculator->SetInputSample( samples[n] );
      covarianceCalculator->Update();

      gaussian->SetMean( *meanCalculator->GetOutput() );
      gaussian->SetCovariance( *covarianceCalculator->GetOutput() );
      }
    else
      {
      typename GaussianType::MeanType mean;
      mean.SetSize( this->m_NumberOfAuxiliaryImages + 1 );
      mean.Fill( 0.0 );
      typename GaussianType::CovarianceType covariance(
        this->m_NumberOfAuxiliaryImages + 1, this->m_NumberOfAuxiliaryImages + 1 );
      covariance.Fill( 0.0 );

      gaussian->SetMean( mean );
      gaussian->SetCovariance( covariance );
      }
    this->m_GaussianMixtureModel.push_back( gaussian );
    this->m_GaussianMixtureModelProportions[n] =
      static_cast<RealType>( samples[n]->Size() )
      / static_cast<RealType>( totalSampleSize );

    sumPosteriors[n] = weights[n].sum();
    }

  return sumPosteriors.sum() / static_cast<RealType>( totalSampleSize );
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
     * Here we assume that the encompassing function is called in order such
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
        typename RealImageType::Pointer distancePriorProbabilityImage = NULL;
        typename RealImageType::ConstPointer priorProbabilityImage = NULL;

        if( this->m_InitializationStrategy == PriorProbabilityImages ||
            this->m_InitializationStrategy == PriorLabelImage )
          {
          for( unsigned int i = 0; i < this->m_NumberOfAuxiliaryImages + 1; i++ )
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

        if( this->m_InitializationStrategy == PriorProbabilityImages )
          {
          priorProbabilityImage = const_cast<RealImageType *>(
              this->GetPriorProbabilityImage( c + 1 ) );
          }
        else if( this->m_InitializationStrategy == PriorLabelImage )
          {
          distancePriorProbabilityImage
            = this->GetDistancePriorProbabilityImageFromPriorLabelImage( c + 1 );
          }

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
          if( !this->GetMaskImage() || this->GetMaskImage()->GetPixel(
                ItO.GetIndex() ) == this->m_MaskLabel )
            {
            RealType mrfPrior = 1.0;
            if( this->m_MRFSmoothingFactor > 0.0 && neighborhoodSize > 1 )
              {
              RealType weightedNumberOfClassNeighbors = 0.0;
              RealType weightedTotalNumberOfNeighbors = 0.0;
              for( unsigned int n = 0; n < neighborhoodSize; n++ )
                {
                if( n == static_cast<unsigned int>( 0.5 * neighborhoodSize ) )
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

                bool      isInBounds = false;
                LabelType label = ItO.GetPixel( n, isInBounds );
                if( isInBounds )
                  {
                  if( static_cast<unsigned int>( label ) == c + 1 )
                    {
                    weightedNumberOfClassNeighbors += ( 1.0 / distance );
                    }
                  weightedTotalNumberOfNeighbors += ( 1.0 / distance );
                  }
                }

              if( weightedTotalNumberOfNeighbors > 0.0 )
                {
                RealType ratio = weightedNumberOfClassNeighbors
                  / weightedTotalNumberOfNeighbors;
                mrfPrior = vcl_exp( -( 1.0 - ratio )
                                    / this->m_MRFSmoothingFactor );
                }
              }

            RealType prior = 1.0;
            if( priorProbabilityImage )
              {
              prior = priorProbabilityImage->GetPixel( ItO.GetIndex() );
              }
            else if( distancePriorProbabilityImage )
              {
              prior = distancePriorProbabilityImage->GetPixel( ItO.GetIndex() );
              }

            typename GaussianType::MeanType oldMean =
              this->m_GaussianMixtureModel[c]->GetMean();
            typename GaussianType::MeanType newMean = oldMean;
            for( unsigned int i = 0; i < this->m_NumberOfAuxiliaryImages + 1; i++ )
              {
              if( smoothImages.size() > i && smoothImages[i] )
                {
                newMean[i] = ( 1.0 - this->m_AdaptiveSmoothingWeights[i] )
                  * newMean[i] + this->m_AdaptiveSmoothingWeights[i]
                  * smoothImages[i]->GetPixel( ItO.GetIndex() );
                }
              }
            this->m_GaussianMixtureModel[c]->SetMean( newMean );

            typename GaussianType::MeasurementVectorType measurement;
            measurement.SetSize( this->m_NumberOfAuxiliaryImages + 1 );
            for( unsigned int i = 0; i < this->m_NumberOfAuxiliaryImages + 1; i++ )
              {
              if( i == 0 )
                {
                measurement[i] = this->GetInput()->GetPixel( ItO.GetIndex() );
                }
              else
                {
                measurement[i] =
                  this->GetAuxiliaryImage( i )->GetPixel( ItO.GetIndex() );
                }
              }

            RealType likelihood =
              this->m_GaussianMixtureModel[c]->Evaluate( measurement );
            RealType posteriorProbability = likelihood * mrfPrior * prior;

            this->m_GaussianMixtureModel[c]->SetMean( oldMean );

            if( this->m_MRFSigmoidAlpha > 0.0 )
              {
              posteriorProbability = 1.0 / ( 1.0 + vcl_exp(
                                               -( posteriorProbability - this->m_MRFSigmoidBeta )
                                               / this->m_MRFSigmoidAlpha ) );
              }

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
      typename RealImageType::Pointer distancePriorProbabilityImage = NULL;
      typename RealImageType::ConstPointer priorProbabilityImage = NULL;

      if( this->m_InitializationStrategy == PriorProbabilityImages ||
          this->m_InitializationStrategy == PriorLabelImage )
        {
        for( unsigned int i = 0; i < this->m_NumberOfAuxiliaryImages + 1; i++ )
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

      if( this->m_InitializationStrategy == PriorProbabilityImages )
        {
        priorProbabilityImage = const_cast<RealImageType *>(
            this->GetPriorProbabilityImage( whichClass ) );
        }
      else if( this->m_InitializationStrategy == PriorLabelImage )
        {
        distancePriorProbabilityImage =
          this->GetDistancePriorProbabilityImageFromPriorLabelImage( whichClass );
        }

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
        if( !this->GetMaskImage() || this->GetMaskImage()->GetPixel(
              ItO.GetIndex() ) == this->m_MaskLabel )
          {
          RealType mrfPrior = 1.0;
          if( this->m_MRFSmoothingFactor > 0.0 && neighborhoodSize > 1 )
            {
            RealType weightedNumberOfClassNeighbors = 0.0;
            RealType weightedTotalNumberOfNeighbors = 0.0;
            for( unsigned int n = 0; n < neighborhoodSize; n++ )
              {
              if( n == static_cast<unsigned int>( 0.5 * neighborhoodSize ) )
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

              bool      isInBounds = false;
              LabelType label = ItO.GetPixel( n, isInBounds );
              if( isInBounds )
                {
                if( static_cast<unsigned int>( label ) == whichClass )
                  {
                  weightedNumberOfClassNeighbors += ( 1.0 / distance );
                  }
                weightedTotalNumberOfNeighbors += ( 1.0 / distance );
                }
              }

            if( weightedTotalNumberOfNeighbors > 0.0 )
              {
              RealType ratio = weightedNumberOfClassNeighbors
                / weightedTotalNumberOfNeighbors;
              mrfPrior = vcl_exp( -( 1.0 - ratio )
                                  / this->m_MRFSmoothingFactor );
              }
            }

          RealType prior = 1.0;
          if( priorProbabilityImage )
            {
            prior = priorProbabilityImage->GetPixel( ItO.GetIndex() );
            }
          else if( distancePriorProbabilityImage )
            {
            prior = distancePriorProbabilityImage->GetPixel( ItO.GetIndex() );
            }

          typename GaussianType::MeanType oldMean =
            this->m_GaussianMixtureModel[whichClass - 1]->GetMean();
          typename GaussianType::MeanType newMean = oldMean;
          for( unsigned int i = 0; i < this->m_NumberOfAuxiliaryImages + 1; i++ )
            {
            if( smoothImages.size() > i && smoothImages[i] )
              {
              newMean[i] = ( 1.0 - this->m_AdaptiveSmoothingWeights[i] )
                * newMean[i] + this->m_AdaptiveSmoothingWeights[i]
                * smoothImages[i]->GetPixel( ItO.GetIndex() );
              }
            }
          this->m_GaussianMixtureModel[whichClass - 1]->SetMean( newMean );

          MeasurementVectorType measurement;
          measurement.SetSize( this->m_NumberOfAuxiliaryImages + 1 );
          for( unsigned int i = 0; i < this->m_NumberOfAuxiliaryImages + 1; i++ )
            {
            if( i == 0 )
              {
              measurement[i] = this->GetInput()->GetPixel( ItO.GetIndex() );
              }
            else
              {
              measurement[i] =
                this->GetAuxiliaryImage( i )->GetPixel( ItO.GetIndex() );
              }
            }
          RealType likelihood =
            this->m_GaussianMixtureModel[whichClass - 1]->Evaluate( measurement );
          RealType posteriorProbability = likelihood * mrfPrior * prior;

          this->m_GaussianMixtureModel[whichClass - 1]->SetMean( oldMean );

          if( this->m_MRFSigmoidAlpha > 0.0 )
            {
            posteriorProbability = 1.0 / ( 1.0 + vcl_exp(
                                             -( posteriorProbability - this->m_MRFSigmoidBeta )
                                             / this->m_MRFSigmoidAlpha ) );
            }

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
     * Here we assume that the encompassing function is called in order such
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
        thresholder->SetInput( const_cast<ClassifiedImageType *>(
                                 this->GetPriorLabelImage() ) );
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
          fastMarching->SetInput( contour->GetOutput() );
          fastMarching->SetStoppingValue( NumericTraits<RealType>::max() );
//           fastMarching->SetTopologyCheck( FastMarchingFilterType::None );
          fastMarching->Update();

          ImageRegionIterator<RealImageType> ItT( thresholder->GetOutput(),
                                                  thresholder->GetOutput()->GetRequestedRegion() );
          ImageRegionIterator<RealImageType> ItF( fastMarching->GetOutput(),
                                                  fastMarching->GetOutput()->GetRequestedRegion() );
          for( ItT.GoToBegin(), ItF.GoToBegin(); !ItT.IsAtEnd(); ++ItT, ++ItF )
            {
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

        RealType labelSigma = 0.1;
        RealType labelBoundaryProbability = 0.75;

        typename LabelParameterMapType::iterator it =
          this->m_PriorLabelParameterMap.find( c + 1 );
        if( it == this->m_PriorLabelParameterMap.end() )
          {
          itkWarningMacro( "The parameters for label \'" << c + 1
                                                         << "\' are not specified.  Using the default values of "
                                                         << "sigma = " << labelSigma << ", boundary probability = "
                                                         << labelBoundaryProbability );
          }
        else
          {
          labelSigma = ( it->second ).first;
          labelBoundaryProbability = ( it->second ).second;
          }
        for( ItD.GoToBegin(); !ItD.IsAtEnd(); ++ItD )
          {
          if( labelSigma == 0 )
            {
            ItD.Set( 0.0 );
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
      thresholder->SetInput( const_cast<ClassifiedImageType *>(
                               this->GetPriorLabelImage() ) );
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
        fastMarching->SetInput( contour->GetOutput() );
        fastMarching->SetStoppingValue( NumericTraits<RealType>::max() );
        // fastMarching->SetTopologyCheck( FastMarchingFilterType::None );
        fastMarching->Update();

        ImageRegionIterator<RealImageType> ItT( thresholder->GetOutput(),
                                                thresholder->GetOutput()->GetRequestedRegion() );
        ImageRegionIterator<RealImageType> ItF( fastMarching->GetOutput(),
                                                fastMarching->GetOutput()->GetRequestedRegion() );
        for( ItT.GoToBegin(), ItF.GoToBegin(); !ItT.IsAtEnd(); ++ItT, ++ItF )
          {
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

      RealType labelSigma = 0.1;
      RealType labelBoundaryProbability = 0.75;

      typename LabelParameterMapType::iterator it =
        this->m_PriorLabelParameterMap.find( whichClass );
      if( it == this->m_PriorLabelParameterMap.end() )
        {
        itkWarningMacro( "The parameters for label \'" << whichClass
                                                       << "\' are not specified.  Using the default values of "
                                                       << "sigma = " << labelSigma << ", boundary probability = "
                                                       << labelBoundaryProbability );
        }
      else
        {
        labelSigma = ( it->second ).first;
        labelBoundaryProbability = ( it->second ).second;
        }
      for( ItD.GoToBegin(); !ItD.IsAtEnd(); ++ItD )
        {
        if( labelSigma == 0 )
          {
          ItD.Set( 0.0 );
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
      probabilityImage = const_cast<RealImageType *>(
          this->GetPriorProbabilityImage( whichClass ) );
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

    ImageRegionConstIteratorWithIndex<RealImageType> ItP( probabilityImage,
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

          /**
           * Since the b-spline fitting routine is biased towards 0, we
           * subtract the mean from the intensity which is followed by fitting.
           * After fitting we add the mean to the sampled b-spline object
           * and the governing control point lattice.
           */
          typename ImageType::PixelType val;
          if( whichImage == 0 )
            {
            val = this->GetInput()->GetPixel( ItP.GetIndex() );
            }
          else
            {
            val = this->GetAuxiliaryImage( whichImage )->GetPixel(
                ItP.GetIndex() );
            }
          intensity[0] = val
            - this->m_GaussianMixtureModel[whichClass - 1]->GetMean()[whichImage];

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

    typedef AddConstantToImageFilter<ScalarImageType,
                                     ScalarType, ScalarImageType> AdderType;

    /**
     * Now we need to add the mean to both the sampled b-spline object
     * and to the control points.  This latter step is valid since
     * b-spline objects are affine invariant, i.e. an affine transformation
     * of the control points followed by reconstruction is equivalent
     * to a direct affine transformation of the reconstructed b-spline object.
     *
     */
    ScalarType classMean;
    classMean[0] =
      this->m_GaussianMixtureModel[whichClass - 1]->GetMean()[whichImage];

    typename AdderType::Pointer adder1 = AdderType::New();
    adder1->SetInput( bsplineImage );
    adder1->SetConstant( classMean );
    adder1->Update();
    bsplineImage = adder1->GetOutput();

    typename AdderType::Pointer adder2 = AdderType::New();
    adder2->SetInput( bspliner->GetPhiLattice() );
    adder2->SetConstant( classMean );
    adder2->Update();

    this->m_ControlPointLattices[whichImage][whichClass - 1] =
      adder2->GetOutput();
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
      break;
      }
    case PriorLabelImage:
      {
      os << "Prior label image" << std::endl;
      os << indent << "  Use Euclidean distance for Prior Labels:" << std::endl;
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
  os << indent << "  MRF sigmoid alpha: "
     << this->m_MRFSigmoidAlpha << std::endl;
  os << indent << "  MRF sigmoid beta: "
     << this->m_MRFSigmoidBeta << std::endl;

  if( this->m_InitializationStrategy == PriorProbabilityImages ||
      this->m_InitializationStrategy == PriorLabelImage &&
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

  os << indent << "Class parameters" << std::endl;
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    os << indent << "  Class " << n + 1 << ": ";
    os << "proportion = " << this->m_GaussianMixtureModelProportions[n]
       << ", mean = " << this->m_GaussianMixtureModel[n]->GetMean() << ", ";
    typename GaussianType::CovarianceType covariance =
      this->m_GaussianMixtureModel[n]->GetCovariance();
    os << "covariance = [";
    for( unsigned int r = 0; r < covariance.Rows(); r++ )
      {
      for( unsigned int c = 0; c < covariance.Cols() - 1; c++ )
        {
        os << covariance( r, c ) << ", ";
        }
      if( r == covariance.Rows() - 1 )
        {
        os << covariance( r, covariance.Cols() - 1 ) << "]" << std::endl;
        }
      else
        {
        os << covariance( r, covariance.Cols() - 1 ) << "; ";
        }
      }
    }
}
} // namespace itk

#endif
