/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkWASPSegmentationImageFilter.txx,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkWASPSegmentationImageFilter_txx
#define __itkWASPSegmentationImageFilter_txx
#include "ReadWriteImage.h"
#include "itkWASPSegmentationImageFilter.h"
#include "itkSurfaceImageCurvature.h"
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
#include "itkMinimumDecisionRule.h"
#include "itkOtsuMultipleThresholdsCalculator.h"
#include "itkSampleClassifier.h"
#include "itkSigmoidImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkWeightedCentroidKdTreeGenerator.h"

#include "itkTimeProbe.h"
#include "itkImageFileWriter.h"

#include "vnl/vnl_vector.h"

#include <algorithm>

namespace itk
{
template <class TInputImage, class TMaskImage, class TClassifiedImage>
WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::WASPSegmentationImageFilter()
{
  this->ProcessObject::SetNumberOfRequiredInputs( 1 );

  this->m_NumberOfClasses = 3;
  this->m_MaximumNumberOfIterations = 5;
  this->m_ElapsedIterations = 0;
  this->m_ConvergenceThreshold = 0.001;

  this->m_MaskLabel = NumericTraits<LabelType>::One;

  this->m_InitializationStrategy = Otsu;

  this->m_PriorProbabilityWeighting = 0.0;

  this->m_MRFSmoothingFactor = 0.3;
  this->m_MRFSigmoidAlpha = 0.1;
  this->m_MRFSigmoidBeta = 0.25;
  this->m_MRFRadius.Fill( 1 );
  this->m_SumProbabilityImage = NULL;
  this->m_SplineOrder = 3;
  this->m_NumberOfLevels.Fill( 8 );
  this->m_NumberOfControlPoints.Fill( this->m_SplineOrder + 1 );
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::~WASPSegmentationImageFilter()
{
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::SetMaskImage( const MaskImageType * mask )
{
  this->SetNthInput( 1, const_cast<MaskImageType *>( mask ) );
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
const typename WASPSegmentationImageFilter
<TInputImage, TMaskImage, TClassifiedImage>::MaskImageType
* WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GetMaskImage() const
  {
  const MaskImageType * maskImage =
    dynamic_cast<const MaskImageType *>( this->ProcessObject::GetInput( 1 ) );

  return maskImage;
  }

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::SetPriorLabelImage( const ClassifiedImageType * prior )
{
  this->m_InitializationStrategy = PriorLabelImage;
  this->SetNthInput( 2, const_cast<ClassifiedImageType *>( prior ) );
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
const typename WASPSegmentationImageFilter
<TInputImage, TMaskImage, TClassifiedImage>::ClassifiedImageType
* WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GetPriorLabelImage() const
  {
  const ClassifiedImageType * prior =
    dynamic_cast<const ClassifiedImageType *>(
      this->ProcessObject::GetInput( 2 ) );

  return prior;
  }

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::SetPriorProbabilityImage(
  unsigned int whichClass, const RealImageType * prior )
{
  if( whichClass < 1 || whichClass > this->m_NumberOfClasses )
    {
    itkExceptionMacro(
      "The prior probability images are inputs 2...2+m_NumberOfClasses-1.  "
      << "The requested image should be in the range [1, m_NumberOfClasses]" )
    }

  this->SetNthInput( 2 + whichClass, const_cast<RealImageType *>( prior ) );
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
const typename WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::RealImageType
* WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
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
      this->ProcessObject::GetInput( 2 + whichClass ) );

  return priorImage;
  }

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
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
    TimeProbe timer;
    timer.Start();
//    probabilityNew = this->StraightUpdateClassParametersAndLabeling();
    probabilityNew = this->RecursiveUpdateClassParametersAndLabeling();
    timer.Stop();

    std::cout << "Elapsed time: " << timer.GetMeanTime() << std::endl;

    this->m_CurrentConvergenceMeasurement = probabilityNew - probabilityOld;

    if( this->m_CurrentConvergenceMeasurement < this->m_ConvergenceThreshold && this->m_ConvergenceThreshold < 1 )
      {
      isConverged = true;
      }
    probabilityOld = probabilityNew;

    itkDebugMacro( "Iteration: " << probabilityNew );
    std::cout << " pNew " << probabilityNew << std::endl;
    this->m_ElapsedIterations++;

    reporter.CompletedStep();
    }
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GenerateInitialClassLabeling()
{
  this->AllocateOutputs();
  this->GetOutput()->FillBuffer( NumericTraits<MaskLabelType>::Zero );

  switch( this->m_InitializationStrategy )
    {
    case KMeans: default:
      {
      this->GenerateInitialClassLabelingWithKMeansClustering();
      this->m_PriorProbabilityWeighting = 0.0;
      break;
      }
    case Otsu:
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
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    ParametersType params( 2 );

    params[0] = initialStats->GetMean( static_cast<LabelType>( n + 1 ) );
    params[1] = vnl_math_sqr(
        initialStats->GetSigma( static_cast<LabelType>( n + 1 ) ) );

    this->m_CurrentClassParameters.push_back( params );
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
WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
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
WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GenerateInitialClassLabelingWithPriorProbabilityImages()
{
  this->GetOutput()->FillBuffer( NumericTraits<LabelType>::Zero );
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
//    std::cout << " bin1 " << std::endl;
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
WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
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

  this->SetNthInput( 2, const_cast<ClassifiedImageType *>( this->GetOutput() ) );
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::GenerateInitialClassLabelingWithKMeansClustering()
{
  std::cout << " KMeans " << std::endl;
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

  this->SetNthInput( 2, const_cast<ClassifiedImageType *>( this->GetOutput() ) );
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
typename WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::RealType
WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::RecursiveUpdateClassParametersAndLabeling()
{
  typename RealImageType::Pointer maxProbabilityImage =
    RealImageType::New();
  maxProbabilityImage->SetRegions( this->GetOutput()->GetRequestedRegion() );
  maxProbabilityImage->SetOrigin( this->GetOutput()->GetOrigin() );
  maxProbabilityImage->SetSpacing( this->GetOutput()->GetSpacing() );
  maxProbabilityImage->SetDirection( this->GetOutput()->GetDirection() );
  maxProbabilityImage->Allocate();
  maxProbabilityImage->FillBuffer( NumericTraits<RealType>::Zero );
  vnl_vector<double> oldmean( this->m_NumberOfClasses, 0 );
  vnl_vector<double> oldvar( this->m_NumberOfClasses, 0);

  this->m_SumProbabilityImage =  RealImageType::New();
  this->m_SumProbabilityImage->SetRegions( this->GetOutput()->GetRequestedRegion() );
  this->m_SumProbabilityImage->SetOrigin( this->GetOutput()->GetOrigin() );
  this->m_SumProbabilityImage->SetSpacing( this->GetOutput()->GetSpacing() );
  this->m_SumProbabilityImage->SetDirection( this->GetOutput()->GetDirection() );
  this->m_SumProbabilityImage->Allocate();
  this->m_SumProbabilityImage->FillBuffer( NumericTraits<RealType>::Zero );

  typename ClassifiedImageType::Pointer maxLabels =
    ClassifiedImageType::New();
  maxLabels->SetRegions( this->GetOutput()->GetRequestedRegion() );
  maxLabels->SetOrigin( this->GetOutput()->GetOrigin() );
  maxLabels->SetSpacing( this->GetOutput()->GetSpacing() );
  maxLabels->SetDirection( this->GetOutput()->GetDirection() );
  maxLabels->Allocate();
  maxLabels->FillBuffer( NumericTraits<LabelType>::Zero );
  ImageRegionIteratorWithIndex<ClassifiedImageType> ItO( maxLabels,
                                                         maxLabels->GetRequestedRegion() );

  this->AmassDistancePriors();

// this is the E-step  in the EM algorithm
  this->m_PosteriorImages.clear();
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    typename RealImageType::Pointer probabilityImage
      = this->CalculatePosteriorProbabilityImage( n + 1, true );
    ImageRegionIterator<RealImageType> ItP( probabilityImage,
                                            probabilityImage->GetRequestedRegion() );
    ImageRegionIterator<RealImageType> ItM( maxProbabilityImage,
                                            maxProbabilityImage->GetRequestedRegion() );
    ImageRegionIterator<RealImageType> ItS( this->m_SumProbabilityImage,
                                            this->m_SumProbabilityImage->GetRequestedRegion() );

    ItP.GoToBegin();
    ItM.GoToBegin();
    ItO.GoToBegin();
    ItS.GoToBegin();

    while( !ItP.IsAtEnd() )
      {
      if( !this->GetMaskImage() || this->GetMaskImage()->GetPixel(
            ItO.GetIndex() ) == this->m_MaskLabel )
        {
        if( this->m_DistanceImages.size() == this->m_NumberOfClasses  )
          {
          float spprob = 1.0;
          if( !this->m_DistanceImages[n]  )
            {
//  get local sum over all  classes with spatial prior
            for( unsigned int sn = 0; sn < this->m_NumberOfClasses; sn++ )
              {
              if( this->m_DistanceImages[sn] )
                {
                spprob = spprob - this->m_DistanceImages[sn]->GetPixel( ItO.GetIndex() );
                }
              }
            if( spprob < 0 )
              {
              spprob = 0;
              }
            ItP.Set(ItP.Get() * spprob);
            }
          }

        if( ItP.Get() >= ItM.Get() )
          {
          ItM.Set( ItP.Get() );
          ItO.Set( static_cast<LabelType>( n + 1 ) );
          }
        ItS.Set( ItS.Get() + ItP.Get() );
        }
      ++ItP;
      ++ItM;
      ++ItO;
      ++ItS;
      }

    oldmean[n] = this->m_CurrentClassParameters[n][0];
    oldvar[n] = this->m_CurrentClassParameters[n][1];
    this->m_PosteriorImages.push_back( probabilityImage );
    }

  // now update the class means and variances

  vnl_vector<double> N( this->m_NumberOfClasses );
  N.fill( 0.0 );
  double        Psum = 0;
  unsigned long ct = 0;
  unsigned int  n = 0;
// this is the M-step  in the EM algorithm
//  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
//   unsigned int n = this->m_ElapsedIterations % this->m_NumberOfClasses;
    {
    ImageRegionIterator<RealImageType> ItS( this->m_SumProbabilityImage,
                                            this->m_SumProbabilityImage->GetRequestedRegion() );
    ImageRegionIterator<RealImageType>  ItP( maxProbabilityImage, maxProbabilityImage->GetRequestedRegion() );
    ImageRegionConstIterator<ImageType> ItI( this->GetInput(), this->GetInput()->GetRequestedRegion() );

    ItI.GoToBegin();
    ItP.GoToBegin();
    ItS.GoToBegin();
    ItO.GoToBegin();

    while( !ItS.IsAtEnd() )
      {
      if( !this->GetMaskImage() || this->GetMaskImage()->GetPixel(ItS.GetIndex() ) == this->m_MaskLabel )
        {
        RealType intensity = static_cast<RealType>( ItI.Get() );
        float    sum = ItS.Get();
        if( sum == 0 )
          {
          sum = 1;
          }
/** Normalize probability images by total probability */
        for( unsigned int nn = 0; nn < this->m_NumberOfClasses; nn++ )
          {
          this->m_PosteriorImages[nn]->SetPixel(ItS.GetIndex(), this->m_PosteriorImages[nn]->GetPixel(
                                                  ItS.GetIndex() ) / sum);
          }
        n = ItO.Get() - 1;
        RealType weight = this->m_PosteriorImages[n]->GetPixel(ItS.GetIndex() ); // ItP.Get()/sum;
        // running weighted mean and variance formulation
        if(  weight  > 0. )
          {
          Psum += (weight);
          ct++;
          }
        if(  weight > 0. )
          {
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
          } // robust part
        }   // mask label
      ++ItI;
      ++ItP;
      ++ItS;
      ++ItO;
      }
    }
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    this->m_CurrentClassParameters[n][0] = this->m_CurrentClassParameters[n][0];
    this->m_CurrentClassParameters[n][1] = this->m_CurrentClassParameters[n][1];
    std::cout << "  Class " << n + 1 << ": ";
    std::cout << "mean = " << this->m_CurrentClassParameters[n][0] << ", ";
    std::cout << "variance = " << this->m_CurrentClassParameters[n][1] << "."  << std::endl;
    }
  this->SetNthOutput( 0, maxLabels );
  return Psum / (ct);
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
typename WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::RealType
WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::StraightUpdateClassParametersAndLabeling()
{
//	std::cout << " begin class param and label " << std::endl;
  vnl_vector<double> N( this->m_NumberOfClasses );
  vnl_vector<double> WT( this->m_NumberOfClasses );
  vnl_vector<double> INTENSbyWT( this->m_NumberOfClasses );
  vnl_vector<double> VARbyWT( this->m_NumberOfClasses );
  N.fill(0);
  WT.fill(0);
  INTENSbyWT.fill(0);
  VARbyWT.fill(0);

  typename RealImageType::Pointer maxProbabilityImage =
    RealImageType::New();
  maxProbabilityImage->SetRegions( this->GetOutput()->GetLargestPossibleRegion() );
  maxProbabilityImage->SetOrigin( this->GetOutput()->GetOrigin() );
  maxProbabilityImage->SetSpacing( this->GetOutput()->GetSpacing() );
  maxProbabilityImage->SetDirection( this->GetOutput()->GetDirection() );
  maxProbabilityImage->Allocate();
  maxProbabilityImage->FillBuffer( NumericTraits<RealType>::Zero );
  vnl_vector<double> oldmean( this->m_NumberOfClasses, 0 );
  vnl_vector<double> oldvar( this->m_NumberOfClasses, 0);

  this->m_SumProbabilityImage =  RealImageType::New();
  this->m_SumProbabilityImage->SetRegions( this->GetOutput()->GetLargestPossibleRegion() );
  this->m_SumProbabilityImage->SetOrigin( this->GetOutput()->GetOrigin() );
  this->m_SumProbabilityImage->SetSpacing( this->GetOutput()->GetSpacing() );
  this->m_SumProbabilityImage->SetDirection( this->GetOutput()->GetDirection() );
  this->m_SumProbabilityImage->Allocate();
  this->m_SumProbabilityImage->FillBuffer( NumericTraits<RealType>::Zero );

  typename ClassifiedImageType::Pointer maxLabels =
    ClassifiedImageType::New();
  maxLabels->SetRegions( this->GetOutput()->GetLargestPossibleRegion() );
  maxLabels->SetOrigin( this->GetOutput()->GetOrigin() );
  maxLabels->SetSpacing( this->GetOutput()->GetSpacing() );
  maxLabels->SetDirection( this->GetOutput()->GetDirection() );
  maxLabels->Allocate();
  maxLabels->FillBuffer( NumericTraits<LabelType>::Zero );

  ImageRegionIterator<RealImageType> ItS( this->m_SumProbabilityImage,
                                          this->m_SumProbabilityImage->GetLargestPossibleRegion() );
  ImageRegionIterator<RealImageType> ItP( maxProbabilityImage,
                                          maxProbabilityImage->GetLargestPossibleRegion() );
  ImageRegionConstIterator<ImageType>               ItI( this->GetInput(), this->GetInput()->GetLargestPossibleRegion() );
  ImageRegionIteratorWithIndex<ClassifiedImageType> ItO( maxLabels, maxLabels->GetLargestPossibleRegion() );

// first step -- calculate the distance maps, if they are required
//  std::vector<typename RealImageType::Pointer>  this->m_DistanceImages;
  if(      this->m_DistanceImages.size() == 0 )
    {
    for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
      {
      typename RealImageType::Pointer distanceImage = NULL;
      if( n <= this->m_PriorLabelSigmas.size() &&
          this->m_PriorLabelSigmas[n] > 0.0 )
        {
        typedef BinaryThresholdImageFilter<ClassifiedImageType, RealImageType>
        ThresholderType;
        typename ThresholderType::Pointer thresholder = ThresholderType::New();
        thresholder->SetInput( const_cast<ClassifiedImageType *>(
                                 this->GetPriorLabelImage() ) );
        thresholder->SetInsideValue( 1 );
        thresholder->SetOutsideValue( 0 );
        thresholder->SetLowerThreshold( static_cast<LabelType>( n + 1 ) );
        thresholder->SetUpperThreshold( static_cast<LabelType>( n + 1 ) );
        thresholder->Update();

        typedef SignedMaurerDistanceMapImageFilter
        <RealImageType, RealImageType> DistancerType;
        typename DistancerType::Pointer distancer = DistancerType::New();
        distancer->SetInput( thresholder->GetOutput() );
        distancer->SetSquaredDistance( true );
        distancer->SetUseImageSpacing( true );
        distancer->SetInsideIsPositive( false );
        distancer->Update();
        distanceImage = distancer->GetOutput();

        ImageRegionIterator<RealImageType> ItD( distanceImage,
                                                distanceImage->GetRequestedRegion() );
// get max dist
        float maxdist = 0;
        for( ItD.GoToBegin(); !ItD.IsAtEnd(); ++ItD )
          {
          if( ItD.Get() < 0.0 )
            {
            if( fabs(ItD.Get() ) > maxdist )
              {
              maxdist = fabs(ItD.Get() );
              }
            }
          }
        for( ItD.GoToBegin(); !ItD.IsAtEnd(); ++ItD )
          {
          float distancePrior = 1;
          float dist = ItD.Get();
          float critval = 0.;                              // the probability at the boundary will be 1.0-critval
          float delta = (maxdist - fabs(dist) ) / maxdist; // in range of zero to one
          // below, the value at the boundary (D=0) is 1-critval and reduces away from the boundary
          if( dist >= 0 )
            {
            distancePrior = vcl_exp( -1.0 * ItD.Get() / vnl_math_sqr( this->m_PriorLabelSigmas[n] ) ) * (1.0 - critval);
            }
          // below, the value inside the object (D>0) increases from 1-crtival to 1
          else
            {
            distancePrior = 1.0 - critval * delta;
            }
          ItD.Set(distancePrior);
          }
        } // end if for spatial prior
      this->m_DistanceImages.push_back( distanceImage );
      }
    }

  float sumpriorsigma = 0;
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    if(   this->m_PriorLabelSigmas.size() == this->m_NumberOfClasses )
      {
      sumpriorsigma += this->m_PriorLabelSigmas[n];
//      std::cout <<" sig " << this->m_PriorLabelSigmas[n] << std::endl;
      }
    }
  // std::cout <<" prior sigma size " << this->m_PriorLabelSigmas.size() << " : " << sumpriorsigma << std::endl;

// this is the E-step  in the EM algorithm
  this->m_PosteriorImages.clear();
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    typename RealImageType::Pointer probabilityImage
      = this->CalculatePosteriorProbabilityImage( n + 1, false );
    ImageRegionIteratorWithIndex<RealImageType> ItT( probabilityImage, probabilityImage->GetLargestPossibleRegion() );

    ItP.GoToBegin();
    ItO.GoToBegin();
    ItS.GoToBegin();
    ItT.GoToBegin();
    while( !ItP.IsAtEnd() )
      {
      if( !this->GetMaskImage() || this->GetMaskImage()->GetPixel(
            ItO.GetIndex() ) == this->m_MaskLabel )
        {
        float spprob = 1;
        //  adjust the probability by the spatial priors .. if any ...
        // first deal with the case when we do not have a spatial prior for this class, n.
        // in this case, the tissue priors are all value 1
        if( sumpriorsigma > 0  )
          {
          if( !this->m_DistanceImages[n]  )
            {
//  get local sum over all  classes with spatial prior
            for( unsigned int sn = 0; sn < this->m_NumberOfClasses; sn++ )
              {
              if( this->m_DistanceImages[sn] )
                {
                spprob = spprob - this->m_DistanceImages[sn]->GetPixel( ItO.GetIndex() );
                }
              }
            if( spprob < 0 )
              {
              spprob = 0;
              }
//	   if ( this->m_DistanceImages[3]->GetPixel( ItO.GetIndex() ) > 0.5  && n == 2 ) std::cout << " sprob " <<
//  this->m_DistanceImages[3]->GetPixel( ItO.GetIndex() )  << " gprob " << ItT.Get()  << " spprob " << spprob <<
// std::endl;
            }
          else
            {
            spprob = this->m_DistanceImages[n]->GetPixel( ItO.GetIndex() );
            }
          ItT.Set(ItT.Get() * spprob);
          }
        if( ItT.Get() >= ItP.Get() )
          {
          ItP.Set( ItT.Get() );
          ItO.Set( static_cast<LabelType>( n + 1 ) );
          }
        ItS.Set( ItS.Get() + ItT.Get() );
        }
      ++ItP;
      ++ItO;
      ++ItS;
      ++ItT;
      }

    // oldmean[n]= this->m_CurrentClassParameters[n][0] ;
    // oldvar[n]= this->m_CurrentClassParameters[n][1] ;
    this->m_PosteriorImages.push_back( probabilityImage );
    }
  std::cout << " done with posteriors " << std::endl;

/** Normalize probability images by total probability */
  ItO.GoToBegin();
  while( !ItO.IsAtEnd() )
    {
    if( !this->GetMaskImage() || this->GetMaskImage()->GetPixel(ItO.GetIndex() ) == this->m_MaskLabel )
      {
      for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
        {
        double sum = this->m_SumProbabilityImage->GetPixel( ItO.GetIndex() );
        if( sum > 0 )
          {
          this->m_PosteriorImages[n]->SetPixel(ItO.GetIndex(), this->m_PosteriorImages[n]->GetPixel(
                                                 ItO.GetIndex() ) / sum);
          }
        }
      float wt = this->m_PosteriorImages[ItO.Get() - 1]->GetPixel(ItO.GetIndex() );
      WT[ItO.Get() - 1] += wt;
      INTENSbyWT[ItO.Get() - 1] += wt * this->GetInput()->GetPixel(ItO.GetIndex() );
      }
    ++ItO;
    }

  std::cout << " done normalizing " << std::endl;
  // now update the class means and variances
  double        Psum = 0;
  unsigned long ct = 0;
  unsigned int  n = 0;
  ItI.GoToBegin();    ItP.GoToBegin();    ItS.GoToBegin();    ItO.GoToBegin();
  while( !ItS.IsAtEnd() )
    {
    if( !this->GetMaskImage() || this->GetMaskImage()->GetPixel(ItS.GetIndex() ) == this->m_MaskLabel )
      {
      RealType intensity = static_cast<RealType>( ItI.Get() );
      n = ItO.Get() - 1;
      RealType weight = this->m_PosteriorImages[n]->GetPixel(ItO.GetIndex() );
      Psum += (weight);
      ct++;
      VARbyWT[n] += weight * vnl_math_sqr(  intensity - INTENSbyWT[n] / WT[n]  );
      }
    ++ItI;
    ++ItP;
    ++ItS;
    ++ItO;
    }

  std::cout << " done with var " << std::endl;
  for( unsigned int n = 0; n < this->m_NumberOfClasses; n++ )
    {
    this->m_CurrentClassParameters[n][0] = INTENSbyWT[n] / WT[n]; // this->m_CurrentClassParameters[n][0]*wt2+oldmean[n]*wt1;
    this->m_CurrentClassParameters[n][1] = VARbyWT[n] / WT[n];    // this->m_CurrentClassParameters[n][1]*wt2+oldvar[n]*wt1;
    std::cout << "  Class " << n + 1 << ": ";
    std::cout << "mean = " << this->m_CurrentClassParameters[n][0] << ", ";
    std::cout << "variance = " << this->m_CurrentClassParameters[n][1] << "."  << std::endl;
    }
  std::cout << " output " << std::endl;
  this->SetNthOutput( 0, maxLabels );
  return Psum / (ct);
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
typename WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::RealImageType::Pointer
WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::CalculatePosteriorProbabilityImage( unsigned int whichClass, bool calcdist )
{
  if( whichClass > this->m_NumberOfClasses )
    {
    itkExceptionMacro( "Requested class is greater than the number of classes." );
    }

  typename RealImageType::Pointer smoothImage = NULL;
  typename RealImageType::Pointer distanceImage = NULL;
  typename RealImageType::ConstPointer priorProbabilityImage = NULL;
  if( this->m_PriorProbabilityWeighting > 0.0 )
    {
    smoothImage = this->CalculateSmoothIntensityImageFromPriorProbabilityImage( whichClass );
    }
// FIXME -- this looks wrong to transition from the prior to the thresholding below
//  only thresholding uses the distance map ...
  if( this->m_InitializationStrategy == PriorProbabilityImages )
    {
    priorProbabilityImage = const_cast<RealImageType *>(
        this->GetPriorProbabilityImage( whichClass ) );
    }
  else
    {
//  BA FIXME -- moving this to the calculation of the normalized across classes posteriors
    if(  whichClass <= this->m_PriorLabelSigmas.size() && this->m_PriorLabelSigmas[whichClass - 1] > 0.0
         && calcdist && this->m_DistanceImages.size() == this->m_PriorLabelSigmas.size() )
      {
      distanceImage = this->m_DistanceImages[whichClass - 1];
      }
    } // end if for spatial prior

//  std::cout <<"  FIXME -- do we need distance map probabilities to be normalized ? " << std::endl;
//  std::cout << " e.g. P_dist(x) = exp( - Dist(x) ) for the class of interest and 1-P_dist elsewhere " << std::endl;
//  std::cout <<" the issue is that if we are 'inside' the tissue , should we not exclude other tissue probabilities? "
// << std::endl;
  typename RealImageType::Pointer posteriorProbabilityImage =
    RealImageType::New();
  posteriorProbabilityImage->SetRegions(
    this->GetOutput()->GetRequestedRegion() );
  posteriorProbabilityImage->SetOrigin( this->GetOutput()->GetOrigin() );
  posteriorProbabilityImage->SetSpacing( this->GetOutput()->GetSpacing() );
  posteriorProbabilityImage->SetDirection( this->GetOutput()->GetDirection() );
  posteriorProbabilityImage->Allocate();
  posteriorProbabilityImage->FillBuffer( 0 );

  if( this->m_CurrentClassParameters[whichClass - 1][1] == 0 )
    {
    return posteriorProbabilityImage;
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
  ImageRegionIterator<RealImageType> ItP( posteriorProbabilityImage,
                                          posteriorProbabilityImage->GetRequestedRegion() );
  ConstNeighborhoodIterator<ClassifiedImageType> ItO( radius, this->GetOutput(),
                                                      this->GetOutput()->GetRequestedRegion() );

  ItI.GoToBegin();
  ItP.GoToBegin();
  ItO.GoToBegin();

  while( !ItI.IsAtEnd() )
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
        typename ClassifiedImageType::OffsetType offset = ItO.GetOffset( n );

        double distance = 0.0;
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          distance += vnl_math_sqr( this->GetOutput()->GetSpacing()[d] );
          }
        distance = vcl_sqrt( distance );
//        distance = 1.0;

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

      RealType mrfPrior = vcl_exp( -( 1.0 - ratio ) / this->m_MRFSmoothingFactor );

      RealType prior = 1.0;
      if( priorProbabilityImage )
        {
        prior = priorProbabilityImage->GetPixel( ItO.GetIndex() );
        }
      else if( distanceImage )
        {
        prior = distanceImage->GetPixel( ItO.GetIndex() );
        }

      RealType mu = this->m_CurrentClassParameters[whichClass - 1][0];
      if( smoothImage )
        {
        mu = ( 1.0 - this->m_PriorProbabilityWeighting ) * mu
          + this->m_PriorProbabilityWeighting
          * smoothImage-> GetPixel( ItO.GetIndex() );
        }
      RealType likelihood =
        vcl_exp( -0.5 * vnl_math_sqr( ItI.Get() - mu )
                 / this->m_CurrentClassParameters[whichClass - 1][1] );

      double finalprob = likelihood * mrfPrior * prior;
      if( this->m_MRFSigmoidAlpha > 0.0 )
        {
        finalprob = 1. / (1. + exp(-1.0 * ( finalprob - this->m_MRFSigmoidBeta) / this->m_MRFSigmoidAlpha ) ); // a
                                                                                                               // decision
                                                                                                               // function
        }
      ItP.Set( finalprob );
      }
    else
      {
      ItP.Set(0);
      }
    ++ItI;
    ++ItP;
    ++ItO;
    }

/*
  if( this->GetMaskImage() )
    {
    typedef MaskImageFilter
      <RealImageType, MaskImageType, RealImageType> MaskerType;
    typename MaskerType::Pointer masker = MaskerType::New();
    masker->SetInput1( posteriorProbabilityImage );
    masker->SetInput2( this->GetMaskImage() );
    masker->SetOutsideValue( 0 );
    masker->Update();
    posteriorProbabilityImage = masker->GetOutput();
    }
*/

  return posteriorProbabilityImage;
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::AmassDistancePriors()
{
  typename RealImageType::Pointer smoothImage = NULL;
  typename RealImageType::Pointer distanceImage = NULL;
  typename RealImageType::ConstPointer priorProbabilityImage = NULL;
  std::cout << " D-Img size " <<  this->m_DistanceImages.size()  << std::endl;
  if( this->m_DistanceImages.size()  == this->m_NumberOfClasses )
    {
    return;
    }
  for( unsigned int whichClass = 1;  whichClass <= this->m_PriorLabelSigmas.size();  whichClass++ )
    {
    if( whichClass <= this->m_PriorLabelSigmas.size() && this->m_PriorLabelSigmas[whichClass - 1] > 0.0
        && this->m_DistanceImages.size() < this->m_PriorLabelSigmas.size() )
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
      distanceImage = distancer->GetOutput();

      ImageRegionIterator<RealImageType> ItD( distanceImage,
                                              distanceImage->GetRequestedRegion() );
// get max dist
      float maxdist = 0;
      for( ItD.GoToBegin(); !ItD.IsAtEnd(); ++ItD )
        {
        if( ItD.Get() < 0.0 )
          {
          if( fabs(ItD.Get() ) > maxdist )
            {
            maxdist = fabs(ItD.Get() );
            }
          }
        }
      for( ItD.GoToBegin(); !ItD.IsAtEnd(); ++ItD )
        {
        float distancePrior = 1;
        float dist = ItD.Get();
        float boundaryprob = 0.75;
        float critval = 1.0 - boundaryprob;              // the probability at the boundary will be 1.0-critval
        float delta = (maxdist - fabs(dist) ) / maxdist; // in range of zero to one
        // below, the value at the boundary (D=0) is 1-critval and reduces away from the boundary
        if( dist >= 0 )
          {
          distancePrior =
            vcl_exp( -1.0 * ItD.Get() * ItD.Get()
                     / vnl_math_sqr( this->m_PriorLabelSigmas[whichClass - 1] ) ) * boundaryprob;
          }
        // below, the value inside the object (D>0) increases from 1-crtival to 1
        else
          {
          distancePrior = 1.0 - critval * delta;
          }
        ItD.Set(distancePrior);
        }
      std::cout << "  adding dist image " << whichClass << std::endl;
      }
    this->m_DistanceImages.push_back( distanceImage );
    }
  // end if for spatial prior
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
typename WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::RealImageType::Pointer
WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
::CalculateSmoothIntensityImageFromPriorProbabilityImage( unsigned int whichClass )
{
  typename RealImageType::Pointer smimg = RealImageType::New();
  smimg->SetRegions( this->GetOutput()->GetRequestedRegion() );
  smimg->SetOrigin( this->GetOutput()->GetOrigin() );
  smimg->SetSpacing( this->GetOutput()->GetSpacing() );
  smimg->SetDirection( this->GetOutput()->GetDirection() );
  smimg->Allocate();
  smimg->FillBuffer( NumericTraits<RealType>::Zero );
  typename RealImageType::Pointer smimg2 = RealImageType::New();
  smimg2->SetRegions( this->GetOutput()->GetRequestedRegion() );
  smimg2->SetOrigin( this->GetOutput()->GetOrigin() );
  smimg2->SetSpacing( this->GetOutput()->GetSpacing() );
  smimg2->SetDirection( this->GetOutput()->GetDirection() );
  smimg2->Allocate();
  smimg2->FillBuffer( NumericTraits<RealType>::Zero );

  unsigned int numrepeats = this->m_SplineOrder;

// set smimg equal to the original input image * threshold output
  float                                   meana = 0;
  unsigned long                           cta = 0;
  ImageRegionIteratorWithIndex<ImageType> ItI(smimg, smimg->GetRequestedRegion() );
  for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
    {
    ItI.Set(  this->m_CurrentClassParameters[whichClass - 1][0] );
    if( this->GetOutput()->GetPixel(ItI.GetIndex() ) == whichClass )
      {
      ItI.Set(this->GetInput()->GetPixel(ItI.GetIndex() ) );
      }
    }

  typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> dgf;
  typename dgf::Pointer filter = dgf::New();
  filter->SetVariance( this->m_SplineOrder );
  filter->SetUseImageSpacingOn();
  filter->SetMaximumError(.01f);
  filter->SetInput(smimg);
  filter->Update();
  return filter->GetOutput();

  std::cout << " starting mean for label " << meana / cta << " label " << whichClass;
  return smimg;
  if( this->m_ElapsedIterations == 0 )
    {
    return smimg;
    }
  typename NeighborhoodIterator<RealImageType>::RadiusType radius;
  unsigned int neighborhoodSize = 1;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    neighborhoodSize *= ( 2 * this->m_MRFRadius[d] + 1 );
    radius[d] = this->m_MRFRadius[d];
    }
  for( unsigned int nr = 0; nr < numrepeats; nr++ )
    {
    float                               smmean = 0;
    unsigned long                       smct = 0;
    NeighborhoodIterator<RealImageType> ItN( radius, smimg, smimg->GetLargestPossibleRegion() );
    ItN.GoToBegin();
    while( !ItN.IsAtEnd() )
      {
      if( !this->GetMaskImage() || this->GetMaskImage()->GetPixel(ItN.GetIndex() ) == this->m_MaskLabel )
        {
        unsigned long ct = 0;
        float         mean = 0;
        for( unsigned int n = 0; n < neighborhoodSize; n++ )
          {
          if( n == static_cast<unsigned int>( 0.5 * neighborhoodSize ) )
            {
            continue;
            }
          bool isInBounds = false;
          if( this->GetOutput()->GetPixel(  ItN.GetIndex(n) ) == whichClass )
            {
            isInBounds = true;
            }
          if( isInBounds )
            {
            float intensity = smimg->GetPixel( ItN.GetIndex(n) );
            mean += intensity;
            ct++;
            }
          }
        if( ct > 0 )
          {
          smct++;  smmean += mean / ct; smimg2->SetPixel(ItN.GetIndex(), mean / ct );
          }
        }
      ++ItN;
      }
    for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
      {
      ItI.Set( smimg2->GetPixel(ItI.GetIndex() ) );
      }

    std::cout << " smmean " << smmean / smct << " mean " << this->m_CurrentClassParameters[whichClass - 1][0];
    }
  std::cout << std::endl;
  WriteImage<RealImageType>( smimg, "temp.nii.gz");
  return smimg;

/*
     typedef itk::SurfaceImageCurvature<RealImageType>  ParamType;
     typename ParamType::Pointer Parameterizer=ParamType::New();
     Parameterizer->SetInput( thresholder->GetOutput() );
     Parameterizer->SetFunctionImage( const_cast<RealImageType *>(this->GetInput()));
     Parameterizer->SetNeighborhoodRadius( sig );
     Parameterizer->SetSigma(sig);
     Parameterizer->SetUseGeodesicNeighborhood(false);
     Parameterizer->SetUseLabel(true);
     Parameterizer->SetThreshold(0.5);
     Parameterizer->IntegrateFunctionOverSurface(true);
     for (unsigned int i=0; i<numrepeats; i++)
        Parameterizer->IntegrateFunctionOverSurface(true);
//    std::cout <<" end integration  " << std::endl;
    return Parameterizer->GetFunctionImage();

  typename ScalarImageType::Pointer bsplineImage;
  std::cout <<" Nulling the BSpline and fitting to current label set " << std::endl;
  this->m_ControlPointLattices[whichClass-1]=NULL;// BA test FIXME
  if( this->m_ControlPointLattices[whichClass-1].GetPointer() != NULL )
    {
    typedef BSplineControlPointImageFilter<ControlPointLatticeType,
      ScalarImageType> BSplineReconstructorType;
    typename BSplineReconstructorType::Pointer bspliner = BSplineReconstructorType::New();
    bspliner->SetInput( this->m_ControlPointLattices[whichClass-1] );
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

    typedef typename BSplineFilterType::WeightsContainerType  WeightsType;
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
//      std::cout << " bin3 " << std::endl;
      typedef BinaryThresholdImageFilter<ClassifiedImageType, RealImageType>
        ThresholderType;
      typename ThresholderType::Pointer thresholder = ThresholderType::New();
      thresholder->SetInput( const_cast<ClassifiedImageType *>(
      this->GetOutput() ) ); // BA test FIXME
//        this->GetPriorLabelImage() ) ); // BA test FIXME
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
          intensity[0] = ItI.Get() - this->m_CurrentClassParameters[whichClass-1][0];

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

    this->m_ControlPointLattices[whichClass-1] = bspliner->GetPhiLattice();
    }

  typedef VectorIndexSelectionCastImageFilter
    <ScalarImageType, RealImageType> CasterType;
  typename CasterType::Pointer caster = CasterType::New();
  caster->SetInput( bsplineImage );
  caster->SetIndex( 0 );
  caster->Update();
  typename RealImageType::Pointer realimg=caster->GetOutput();

// make average bspline intensity match the class mean
    ImageRegionIteratorWithIndex<RealImageType> ItB( realimg,
      realimg->GetBufferedRegion() );
    float bmean=0;
    double ct=1.e-9;
    for(  ItB.GoToBegin(); !ItB.IsAtEnd(); ++ItB )
     {
      if( !this->GetMaskImage() ||
        this->GetMaskImage()->GetPixel( ItB.GetIndex() ) == this->m_MaskLabel )
        {
    ItB.Set(ItB.Get()+this->m_CurrentClassParameters[whichClass-1][0]);
    if ( this->GetOutput()->GetPixel(ItB.GetIndex()) == whichClass  )  {
            bmean+=ItB.Get();
            ct++;
          }
        }
     }
    bmean/=ct;
    float bscale= this->m_CurrentClassParameters[whichClass-1][0]/bmean;
//    std::cout << " bscale " << bscale << " bmean " << bmean << " mean " <<  this->m_CurrentClassParameters[whichClass-1][0] << std::endl;
    for(  ItB.GoToBegin(); !ItB.IsAtEnd(); ++ItB )
     {
      if( !this->GetMaskImage() ||
        this->GetMaskImage()->GetPixel( ItB.GetIndex() ) == this->m_MaskLabel )
        {
          ItB.Set( ItB.Get()*bscale );
        }
     }

     return realimg;
*/
}

template <class TInputImage, class TMaskImage, class TClassifiedImage>
void
WASPSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>
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
      os << indent << "  Prior label sigmas" << std::endl;
      for( unsigned int n = 0; n < this->m_PriorLabelSigmas.size(); n++ )
        {
        os << indent << "    Class " << n + 1 << ": sigma = "
           << this->m_PriorLabelSigmas[n] << std::endl;
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
    os << "variance = " << this->m_CurrentClassParameters[n][1] << "."
       << std::endl;
    }
}
} // namespace itk

#endif
