/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkApocritaSegmentationImageFilter.txx,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkApocritaSegmentationImageFilter_txx
#define __itkApocritaSegmentationImageFilter_txx

#include "itkApocritaSegmentationImageFilter.h"

#include "itkAddImageFilter.h"
// #include "itkAddConstantToImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkEuclideanDistance.h"
#include "itkGaussianMixtureModelComponent.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageToListGenerator.h"
#include "itkIterationReporter.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkMinimumDecisionRule.h"
#include "itkOtsuMultipleThresholdsCalculator.h"
#include "itkSampleClassifier.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkWeightedCentroidKdTreeGenerator.h"

#include "itkTimeProbe.h"

#include "vnl/vnl_vector.h"

#include <algorithm>

namespace itk
{
template <class TInputImage, class TMaskImage, class TClassifiedImage>
ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::ApocritaSegmentationImageFilter()
{
  this->ProcessObject::SetNumberOfRequiredInputs( 1 );

  this->m_NumberOfClasses = 3;
  this->m_MaximumNumberOfIterations = 5;
  this->m_ElapsedIterations = 0;
  this->m_ConvergenceThreshold = 0.001;

  this->m_MaskLabel = NumericTraits<LabelType>::One;

  this->m_InitializationStrategy = Otsu;

  this->m_PriorProbabilityWeighting = 0.0;
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
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::~ApocritaSegmentationImageFilter()
{
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::SetMaskImage( const MaskImageType * mask )
{
  this->SetNthInput( 1, const_cast<MaskImageType *>( mask ) );
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
const typename ApocritaSegmentationImageFilter
<TInputImage, TMaskImage, TClassifiedImage>::MaskImageType
* ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GetMaskImage() const
  {
  const MaskImageType * maskImage =
    dynamic_cast<const MaskImageType *>( this->ProcessObject::GetInput( 1 ) );

  return maskImage;
  }

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::SetPriorLabelImage( const ClassifiedImageType * prior )
{
  this->m_InitializationStrategy = PriorLabelImage;
  this->SetNthInput( 2, const_cast<ClassifiedImageType *>( prior ) );
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
const typename ApocritaSegmentationImageFilter
<TInputImage, TMaskImage, TClassifiedImage>::ClassifiedImageType
* ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GetPriorLabelImage() const
  {
  const ClassifiedImageType * prior =
    dynamic_cast<const ClassifiedImageType *>(
      this->ProcessObject::GetInput( 2 ) );

  return prior;
  }

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::SetPriorProbabilityImage(
  unsigned int whichClass, const RealImageType * prior )
{
  if( whichClass < 1 || whichClass > this->m_NumberOfClasses )
    {
    itkExceptionMacro(
      "The prior probability images are inputs 2...2+m_NumberOfClasses-1.  "
      << "The requested image should be in the range [1, m_NumberOfClasses]" )
    }

  this->SetNthInput( 1 + whichClass, const_cast<RealImageType *>( prior ) );
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
const typename ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::RealImageType
* ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GetPriorProbabilityImage( unsigned int whichClass ) const
  {
  if( whichClass < 1 || whichClass > this->m_NumberOfClasses )
    {
    itkExceptionMacro(
      "The prior probability images are inputs 2...2+m_NumberOfClasses-1.  "
      << "The requested image should be in the range [1, m_NumberOfClasses]" )
    }

  const RealImageType *priorImage =
    dynamic_cast<const RealImageType *>(
      this->ProcessObject::GetInput( 1 + whichClass ) );

  return priorImage;
  }

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
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

    std::cout << "Elapsed time: " << timer.GetMeanTime() << std::endl;

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
ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
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
      this->m_PriorProbabilityWeighting = 0.0;
      break;
      }
    case Otsu:  default:
      {
      this->GenerateInitialClassLabelingWithOtsuThresholding();
      this->m_PriorProbabilityWeighting = 0.0;
      break;
      }
    case PriorProbabilityImages:
      {
      // Check for proper setting of prior probability images.
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
        this->m_PriorProbabilityWeighting = 0.0;
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

  typedef LabelStatisticsImageFilter<ImageType, MaskImageType> StatsType;
  typename StatsType::Pointer initialStats = StatsType::New();
  initialStats->SetInput( this->GetInput() );
  initialStats->SetLabelInput( this->GetOutput() );
  initialStats->UseHistogramsOff();
  initialStats->Update();

  this->m_CurrentClassParameters.clear();

  RealType sumCount = 0.0;
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    ParametersType params( 3 );

    params[0] = initialStats->GetMean( static_cast<LabelType>( n + 1 ) );
    params[1] = vnl_math_sqr(
        initialStats->GetSigma( static_cast<LabelType>( n + 1 ) ) );
    params[2] = static_cast<RealType>(
        initialStats->GetCount( static_cast<LabelType>( n + 1 ) ) );
    sumCount += params[2];

    this->m_CurrentClassParameters.push_back( params );
    }
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    this->m_CurrentClassParameters[n][2] /= sumCount;
    }

  if( this->m_InitializationStrategy == PriorProbabilityImages ||
      this->m_InitializationStrategy == PriorLabelImage )
    {
    for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
      {
      this->m_ControlPointLattices.push_back( NULL );
      }
    }
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
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
    if( priorProbabilities.sum() > 0.0 )
      {
      priorProbabilities /= priorProbabilities.sum();
      RealType maxValue = priorProbabilities.max_value();
      if( maxValue < 0.5 )
        {
        unsigned int argMax = 0;
        for( unsigned int i = 0; i < priorProbabilities.size(); i++ )
          {
          if( priorProbabilities[i] == maxValue )
            {
            argMax = i;
            break;
            }
          }
        RealType probabilityDifference = 0.5 - priorProbabilities[argMax];
        for( unsigned int i = 0; i < priorProbabilities.size(); i++ )
          {
          if( i == argMax )
            {
            continue;
            }
          priorProbabilities[i] += ( probabilityDifference
                                     / static_cast<RealType>( priorProbabilities.size() - 1 ) );
          }
        priorProbabilities[argMax] = 0.5;
        }
      }
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
ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GenerateInitialClassLabelingWithPriorProbabilityImages()
{
  this->GetOutput()->FillBuffer( NumericTraits<LabelType>::Zero );
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    typedef BinaryThresholdImageFilter<RealImageType, ClassifiedImageType>
      ThresholderType;
    typename ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetInput( this->GetPriorProbabilityImage( n + 1 ) );
    thresholder->SetInsideValue( n + 1 );
    thresholder->SetOutsideValue( 0 );
    thresholder->SetLowerThreshold( 0.5 );
    thresholder->SetUpperThreshold( 1.0 );

    typedef AddImageFilter<ClassifiedImageType, ClassifiedImageType,
                           ClassifiedImageType> AdderType;
    typename AdderType::Pointer adder = AdderType::New();
    adder->SetInput1( this->GetOutput() );
    adder->SetInput2( thresholder->GetOutput() );
    adder->Update();

    this->SetNthOutput( 0, adder->GetOutput() );
    }

  if( this->GetMaskImage() )
    {
    typedef MaskImageFilter<ClassifiedImageType, MaskImageType> MaskerType;
    typename MaskerType::Pointer masker = MaskerType::New();
    masker->SetInput1( this->GetOutput() );
    masker->SetInput2( this->GetMaskImage() );
    masker->Update();

    this->SetNthOutput( 0, masker->GetOutput() );
    }
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
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
ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GenerateInitialClassLabelingWithKMeansClustering()
{
  typedef typename Statistics::ImageToListGenerator<ImageType, MaskImageType>
    ListSampleGeneratorType;
  typename ListSampleGeneratorType::Pointer sampler
    = ListSampleGeneratorType::New();
  sampler->SetInput( this->GetInput() );
  if( this->GetMaskImage() )
    {
    sampler->SetMaskImage( this->GetMaskImage() );
    sampler->SetMaskValue( this->m_MaskLabel );
    }
  sampler->Update();

  typedef typename ListSampleGeneratorType::ListSampleType ListSampleType;
  typedef typename ListSampleGeneratorType::MeasurementVectorType
    MeasurementVectorType;
  typedef Statistics::WeightedCentroidKdTreeGenerator
    <ListSampleType> TreeGeneratorType;
  typedef typename TreeGeneratorType::KdTreeType           TreeType;
  typedef Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;
  typedef typename EstimatorType::ParametersType           ParametersType;

  typename TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
  treeGenerator->SetSample( const_cast<ListSampleType *>(
                              sampler->GetListSample() ) );
  treeGenerator->SetBucketSize( 16 );
  treeGenerator->Update();

  /**
   * Guess initial class means by dividing the dynamic range
   *  into equal intervals.
   */

  RealType maxValue = itk::NumericTraits<RealType>::min();
  RealType minValue = itk::NumericTraits<RealType>::max();

  ImageRegionConstIteratorWithIndex<ImageType> ItI( this->GetInput(),
                                                    this->GetInput()->GetRequestedRegion() );
  for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
    {
    if( !this->GetMaskImage() || this->GetMaskImage()->GetPixel( ItI.GetIndex() )
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

  typename EstimatorType::Pointer estimator = EstimatorType::New();
  ParametersType initialMeans( this->m_NumberOfClasses );
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    initialMeans[n] = minValue + ( maxValue - minValue )
      * static_cast<RealType>( n + 1 )
      / static_cast<RealType>( this->m_NumberOfClasses + 1 );
    }
  estimator->SetParameters( initialMeans );

  estimator->SetKdTree( treeGenerator->GetOutput() );
  estimator->SetMaximumIteration( 200 );
  estimator->SetCentroidPositionChangesThreshold( 0.0 );
  estimator->StartOptimization();

  // Now classify the samples

  typedef Statistics::SampleClassifier<ListSampleType> ClassifierType;
  typedef MinimumDecisionRule                          DecisionRuleType;
  typename DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();
  typename ClassifierType::Pointer classifier = ClassifierType::New();

  classifier->SetDecisionRule( decisionRule.GetPointer() );
  classifier->SetSample( const_cast<ListSampleType *>(
                           sampler->GetListSample() ) );
  classifier->SetNumberOfClasses( this->m_NumberOfClasses  );

  typedef itk::Statistics::EuclideanDistance<MeasurementVectorType>
    MembershipFunctionType;

  typedef std::vector<unsigned int> ClassLabelVectorType;
  ClassLabelVectorType classLabels;
  classLabels.resize( this->m_NumberOfClasses );

  // Order the cluster means so that the lowest mean corresponds to label '1',
  //  the second lowest to label '2', etc.
  std::vector<RealType> estimatorParameters;
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    estimatorParameters.push_back( estimator->GetParameters()[n] );
    }
  std::sort( estimatorParameters.begin(), estimatorParameters.end() );
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    classLabels[n] = n + 1;
    typename MembershipFunctionType::Pointer
    membershipFunction = MembershipFunctionType::New();
    typename MembershipFunctionType::OriginType origin(
      sampler->GetMeasurementVectorSize() );
    origin[0] = estimatorParameters[n];
    membershipFunction->SetOrigin( origin );
    classifier->AddMembershipFunction( membershipFunction.GetPointer() );
    }
  classifier->SetMembershipFunctionClassLabels( classLabels );
  classifier->Update();

  // Now classify the pixels

  ImageRegionIteratorWithIndex<ClassifiedImageType> ItO( this->GetOutput(),
                                                         this->GetOutput()->GetRequestedRegion() );
  typedef typename ClassifierType::OutputType          ClassifierOutputType;
  typedef typename ClassifierOutputType::ConstIterator LabelIterator;

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
typename ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::RealType
ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::UpdateClassParametersAndLabeling()
{
  typename RealImageType::Pointer maxProbabilityImage =
    RealImageType::New();
  maxProbabilityImage->SetRegions( this->GetOutput()->GetRequestedRegion() );
  maxProbabilityImage->SetOrigin( this->GetOutput()->GetOrigin() );
  maxProbabilityImage->SetSpacing( this->GetOutput()->GetSpacing() );
  maxProbabilityImage->SetDirection( this->GetOutput()->GetDirection() );
  maxProbabilityImage->Allocate();
  maxProbabilityImage->FillBuffer( NumericTraits<RealType>::Zero );

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

  ParametersType sumPosteriors( this->m_NumberOfClasses );
  sumPosteriors.Fill( 0.0 );
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    typename RealImageType::Pointer posteriorProbabilityImage
      = this->GetPosteriorProbabilityImage( n + 1 );

    ImageRegionIteratorWithIndex<ClassifiedImageType> ItO( maxLabels,
                                                           maxLabels->GetRequestedRegion() );
    ImageRegionConstIterator<RealImageType> ItP( posteriorProbabilityImage,
                                                 posteriorProbabilityImage->GetRequestedRegion() );
    ImageRegionIterator<RealImageType> ItM( maxProbabilityImage,
                                            maxProbabilityImage->GetRequestedRegion() );

    ItP.GoToBegin();
    ItM.GoToBegin();
    ItO.GoToBegin();
    while( !ItP.IsAtEnd() )
      {
      if( !this->GetMaskImage() || this->GetMaskImage()->GetPixel(
            ItO.GetIndex() ) == this->m_MaskLabel )
        {
        if( ItP.Get() >= ItM.Get() )
          {
          ItM.Set( ItP.Get() );
          ItO.Set( static_cast<LabelType>( n + 1 ) );
          }
        sumPosteriors[n] += ItP.Get();
        }
      ++ItP;
      ++ItM;
      ++ItO;
      }

    // Perform the following calculation to update the class proportions
    typename RealImageType::Pointer distancePriorProbabilityImage = NULL;
    typename RealImageType::ConstPointer priorProbabilityImage = NULL;

    if( this->m_PriorProbabilityWeighting > 0.0 &&
        this->m_InitializationStrategy == PriorProbabilityImages )
      {
      priorProbabilityImage = const_cast<RealImageType *>(
          this->GetPriorProbabilityImage( n + 1 ) );
      }
    else if(  this->m_PriorProbabilityWeighting > 0.0 &&
              this->m_InitializationStrategy == PriorLabelImage )
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
        ItW.Set( ItW.Get() + this->m_CurrentClassParameters[n][2]
                 * priorProbability );
        }
      }
    }
  // Update the class proportions
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    RealType denominator = 0.0;

    typename RealImageType::Pointer distancePriorProbabilityImage = NULL;
    typename RealImageType::ConstPointer priorProbabilityImage = NULL;

    if( this->m_PriorProbabilityWeighting > 0.0 &&
        this->m_InitializationStrategy == PriorProbabilityImages )
      {
      priorProbabilityImage = const_cast<RealImageType *>(
          this->GetPriorProbabilityImage( n + 1 ) );
      }
    else if(  this->m_PriorProbabilityWeighting > 0.0 &&
              this->m_InitializationStrategy == PriorLabelImage )
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
        denominator += ( priorProbability / ItW.Get() );
        }
      }
    this->m_CurrentClassParameters[n][2] = sumPosteriors[n] / denominator;
    }

  // now update the class means and variances

  vnl_vector<RealType> N( this->m_NumberOfClasses );
  N.fill( 0.0 );

  ImageRegionIteratorWithIndex<ClassifiedImageType> ItO( maxLabels,
                                                         maxLabels->GetRequestedRegion() );
  ImageRegionIterator<RealImageType> ItM( maxProbabilityImage,
                                          maxProbabilityImage->GetRequestedRegion() );
  ImageRegionConstIterator<ImageType> ItI( this->GetInput(),
                                           this->GetInput()->GetRequestedRegion() );

  ItI.GoToBegin();
  ItM.GoToBegin();
  ItO.GoToBegin();

  unsigned long voxelCount = 0;
/**/
// below -- use the maximum prob to get the mean / variance
//
  while( !ItO.IsAtEnd() )
    {
    if( !this->GetMaskImage() || this->GetMaskImage()->GetPixel(
          ItO.GetIndex() ) == this->m_MaskLabel )
      {
      RealType intensity = static_cast<RealType>( ItI.Get() );
      RealType weight = ItM.Get();

      voxelCount++;

      // running weighted mean and variance formulation

      unsigned int n = ItO.Get() - 1;

      N[n] += weight;
      this->m_CurrentClassParameters[n][0] = ( ( N[n] - weight )
                                               * this->m_CurrentClassParameters[n][0] + weight * intensity ) / N[n];

      if( N[n] > weight )
        {
        this->m_CurrentClassParameters[n][1]
          = this->m_CurrentClassParameters[n][1] * ( N[n] - weight )
            / N[n] + vnl_math_sqr( intensity
                                   - this->m_CurrentClassParameters[n][0] ) * weight / ( N[n] - weight );
        }
      }
    ++ItI;
    ++ItM;
    ++ItO;
    }

/*

// here get mean / var over all classes
   for (unsigned int n=0;  n < this->m_NumberOfClasses; n++) {

  ItI.GoToBegin();
  ItM.GoToBegin();
  ItO.GoToBegin();

    typename RealImageType::Pointer posteriorProbabilityImage
      = this->GetPosteriorProbabilityImage( n + 1 );

   while( !ItI.IsAtEnd() )
    {
    if( !this->GetMaskImage() || this->GetMaskImage()->GetPixel(
      ItI.GetIndex() ) == this->m_MaskLabel )
      {
      RealType intensity = static_cast<RealType>( ItI.Get() );
      RealType weight =   posteriorProbabilityImage->GetPixel(ItI.GetIndex());
      voxelCount++;

      // running weighted mean and variance formulation


      N[n] += weight;
      this->m_CurrentClassParameters[n][0] = ( ( N[n] - weight ) *
        this->m_CurrentClassParameters[n][0] + weight * intensity ) / N[n];

      if( N[n] > weight )
        {
        this->m_CurrentClassParameters[n][1]
          = this->m_CurrentClassParameters[n][1] * ( N[n] - weight )
          / N[n] + vnl_math_sqr( intensity -
          this->m_CurrentClassParameters[n][0] ) * weight / ( N[n] - weight );
        }
      }
    ++ItI;
    }
   }

*/

  this->SetNthOutput( 0, maxLabels );

  return N.sum() / static_cast<RealType>( voxelCount );
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
typename ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::RealImageType::Pointer
ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
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
        typename RealImageType::Pointer smoothImage = NULL;
        typename RealImageType::Pointer distancePriorProbabilityImage = NULL;
        typename RealImageType::ConstPointer priorProbabilityImage = NULL;

        if( this->m_PriorProbabilityWeighting > 0.0 )
          {
          smoothImage =
            this->CalculateSmoothIntensityImageFromPriorProbabilityImage(
              c + 1 );
          }

        if( this->m_PriorProbabilityWeighting > 0.0 &&
            this->m_InitializationStrategy == PriorProbabilityImages )
          {
          priorProbabilityImage = const_cast<RealImageType *>(
              this->GetPriorProbabilityImage( c + 1 ) );
          }
        else if(  this->m_PriorProbabilityWeighting > 0.0 &&
                  this->m_InitializationStrategy == PriorLabelImage )
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

        ImageRegionConstIterator<ImageType> ItI( this->GetInput(),
                                                 this->GetInput()->GetRequestedRegion() );
        ConstNeighborhoodIterator<ClassifiedImageType> ItO( radius,
                                                            this->GetOutput(),
                                                            this->GetOutput()->GetRequestedRegion() );
        ImageRegionIterator<RealImageType> ItS(
          this->m_SumPosteriorProbabilityImage,
          this->m_SumPosteriorProbabilityImage->GetRequestedRegion() );
        for( ItI.GoToBegin(), ItO.GoToBegin(), ItS.GoToBegin(); !ItI.IsAtEnd();
             ++ItI, ++ItO, ++ItS )
          {
          if( !this->GetMaskImage() || this->GetMaskImage()->GetPixel(
                ItO.GetIndex() ) == this->m_MaskLabel )
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
            RealType ratio = weightedNumberOfClassNeighbors
              / weightedTotalNumberOfNeighbors;

            RealType mrfPrior = 1.0;
            if( this->m_MRFSmoothingFactor > 0.0 )
              {
              mrfPrior = vcl_exp( -( 1.0 - ratio )
                                  / this->m_MRFSmoothingFactor );
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

            RealType mu = this->m_CurrentClassParameters[c][0];
            if( smoothImage )
              {
              mu = ( 1.0 - this->m_PriorProbabilityWeighting ) * mu
                + this->m_PriorProbabilityWeighting
                * smoothImage->GetPixel( ItO.GetIndex() );
              }
            RealType likelihood = 1.0 / vcl_sqrt( 2.0 * vnl_math::pi
                                                  * this->m_CurrentClassParameters[c][1] )
              * vcl_exp( -0.5 * vnl_math_sqr( ItI.Get() - mu )
                         / this->m_CurrentClassParameters[c][1] );

            RealType posteriorProbability = likelihood * mrfPrior * prior
              * this->m_CurrentClassParameters[c][2];

            if( this->m_MRFSigmoidAlpha > 0.0 )
              {
              posteriorProbability = 1.0 / ( 1.0 + vcl_exp(
                                               -( posteriorProbability - this->m_MRFSigmoidBeta )
                                               / this->m_MRFSigmoidAlpha ) );
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
      typename RealImageType::Pointer smoothImage = NULL;
      typename RealImageType::Pointer distancePriorProbabilityImage = NULL;
      typename RealImageType::ConstPointer priorProbabilityImage = NULL;

      if( this->m_PriorProbabilityWeighting > 0.0 )
        {
        smoothImage =
          this->CalculateSmoothIntensityImageFromPriorProbabilityImage(
            whichClass );
        }

      if( this->m_PriorProbabilityWeighting > 0.0 &&
          this->m_InitializationStrategy == PriorProbabilityImages )
        {
        priorProbabilityImage = const_cast<RealImageType *>(
            this->GetPriorProbabilityImage( whichClass ) );
        }
      else if(  this->m_PriorProbabilityWeighting > 0.0 &&
                this->m_InitializationStrategy == PriorLabelImage )
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

      ImageRegionConstIterator<ImageType> ItI( this->GetInput(),
                                               this->GetInput()->GetRequestedRegion() );
      ConstNeighborhoodIterator<ClassifiedImageType> ItO( radius,
                                                          this->GetOutput(), this->GetOutput()->GetRequestedRegion() );
      for( ItI.GoToBegin(), ItO.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItO )
        {
        if( !this->GetMaskImage() || this->GetMaskImage()->GetPixel(
              ItO.GetIndex() ) == this->m_MaskLabel )
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
          RealType ratio = weightedNumberOfClassNeighbors
            / weightedTotalNumberOfNeighbors;

          RealType mrfPrior = 1.0;
          if( this->m_MRFSmoothingFactor > 0.0 )
            {
            mrfPrior = vcl_exp( -( 1.0 - ratio )
                                / this->m_MRFSmoothingFactor );
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

          RealType mu = this->m_CurrentClassParameters[whichClass - 1][0];
          if( smoothImage )
            {
            mu = ( 1.0 - this->m_PriorProbabilityWeighting ) * mu
              + this->m_PriorProbabilityWeighting
              * smoothImage->GetPixel( ItO.GetIndex() );
            }
          RealType likelihood = 1.0 / vcl_sqrt( 2.0 * vnl_math::pi
                                                * this->m_CurrentClassParameters[whichClass - 1][1] )
            * vcl_exp( -0.5 * vnl_math_sqr( ItI.Get() - mu )
                       / this->m_CurrentClassParameters[whichClass - 1][1] );

          RealType posteriorProbability = likelihood * mrfPrior * prior
            * this->m_CurrentClassParameters[whichClass - 1][2];

          if( this->m_MRFSigmoidAlpha > 0.0 )
            {
            posteriorProbability = 1.0 / ( 1.0 + vcl_exp(
                                             -( posteriorProbability - this->m_MRFSigmoidBeta )
                                             / this->m_MRFSigmoidAlpha ) );
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
typename ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::RealImageType::Pointer
ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
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

        typedef SignedMaurerDistanceMapImageFilter
          <RealImageType, RealImageType> DistancerType;
        typename DistancerType::Pointer distancer = DistancerType::New();
        distancer->SetInput( thresholder->GetOutput() );
        distancer->SetSquaredDistance( true );
        distancer->SetUseImageSpacing( true );
        distancer->SetInsideIsPositive( false );
        distancer->Update();

        RealType maximumInteriorDistance = 0.0;

        ImageRegionIterator<RealImageType> ItD( distancer->GetOutput(),
                                                distancer->GetOutput()->GetRequestedRegion() );
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
        adder->SetInput2( distancer->GetOutput() );
        adder->Update();

        this->m_SumDistancePriorProbabilityImage = adder->GetOutput();

        if( ( c == 0 ) && this->m_MinimizeMemoryUsage )
          {
          distancePriorProbabilityImage = distancer->GetOutput();
          }
        if( !this->m_MinimizeMemoryUsage )
          {
          this->m_DistancePriorProbabilityImages.push_back(
            distancer->GetOutput() );
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

      typedef SignedMaurerDistanceMapImageFilter
        <RealImageType, RealImageType> DistancerType;
      typename DistancerType::Pointer distancer = DistancerType::New();
      distancer->SetInput( thresholder->GetOutput() );
      distancer->SetSquaredDistance( true );
      distancer->SetUseImageSpacing( true );
      distancer->SetInsideIsPositive( false );
      distancer->Update();

      distancePriorProbabilityImage = distancer->GetOutput();

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
typename ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::RealImageType::Pointer
ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::CalculateSmoothIntensityImageFromPriorProbabilityImage( unsigned int whichClass )
{
  typename ScalarImageType::Pointer bsplineImage;

  if( this->m_ControlPointLattices[whichClass - 1].GetPointer() != NULL )
    {
    typedef BSplineControlPointImageFilter<ControlPointLatticeType,
                                           ScalarImageType> BSplineReconstructorType;
    typename BSplineReconstructorType::Pointer bspliner
      = BSplineReconstructorType::New();

    bspliner->SetInput( this->m_ControlPointLattices[whichClass - 1] );
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

    ImageRegionConstIterator<ImageType> ItI( this->GetInput(),
                                             this->GetInput()->GetRequestedRegion() );
    ImageRegionConstIteratorWithIndex<RealImageType> ItP( probabilityImage,
                                                          probabilityImage->GetBufferedRegion() );
    for( ItI.GoToBegin(), ItP.GoToBegin(); !ItP.IsAtEnd(); ++ItI, ++ItP )
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
          intensity[0] = ItI.Get()
            - this->m_CurrentClassParameters[whichClass - 1][0];

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

    std::cout << " DISABLED DUE TO LACK OF AddConstantToImageFilter<ScalarImageType> FILTER in ITK " << std::endl;
//    typedef AddConstantToImageFilter<ScalarImageType,
//    ScalarType, ScalarImageType> AdderType;

    /**
     * Now we need to add the mean to both the sampled b-spline object
     * and to the control points.  This latter step is valid since
     * b-spline objects are affine invariant, i.e. an affine transformation
     * of the control points followed by reconstruction is equivalent
     * to a direct affine transformation of the reconstructed b-spline object.
     *
     */
/*
    ScalarType classMean;
    classMean[0] = this->m_CurrentClassParameters[whichClass-1][0];

    typename AdderType::Pointer adder1 = AdderType::New();
    adder1->SetInput( bsplineImage );
    adder1->SetConstant( classMean );
    adder1->Update();
    bsplineImage = adder1->GetOutput();

    typename AdderType::Pointer adder2 = AdderType::New();
    adder2->SetInput( bspliner->GetPhiLattice() );
    adder2->SetConstant( classMean );
    adder2->Update();

    this->m_ControlPointLattices[whichClass-1] = adder2->GetOutput();
*/
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
ApocritaSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
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
      os << indent << "  Prior probability weighting: "
         << this->m_PriorProbabilityWeighting << std::endl;
      break;
      }
    case PriorLabelImage:
      {
      os << "Prior label image" << std::endl;
      os << indent << "  Prior probability weighting: "
         << this->m_PriorProbabilityWeighting << std::endl;
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

  if( this->m_PriorProbabilityWeighting > 0.0 )
    {
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
    os << "mean = " << this->m_CurrentClassParameters[n][0] << ", ";
    os << "variance = " << this->m_CurrentClassParameters[n][1] << ", ";
    os << "proportion = " << this->m_CurrentClassParameters[n][2] << "."
       << std::endl;
    }
}
} // namespace itk

#endif
