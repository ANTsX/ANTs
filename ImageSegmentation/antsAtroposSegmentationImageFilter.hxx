/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
  https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt
  for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __antsAtroposSegmentationImageFilter_hxx
#define __antsAtroposSegmentationImageFilter_hxx

#include "antsGaussianListSampleFunction.h"
#include "antsAllocImage.h"
#include "itkAddImageFilter.h"
#include "itkAddConstantToImageFilter.h"
#include "itkBinaryContourImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkContinuousIndex.h"
#include "itkDistanceToCentroidMembershipFunction.h"
#include "itkFastMarchingImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkIterationReporter.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkMinimumDecisionRule.h"
#include "itkMultiplyImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkOtsuMultipleThresholdsCalculator.h"
#include "itkSampleClassifierFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkVariableSizeMatrix.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkWeightedCentroidKdTreeGenerator.h"

#include "vnl/vnl_vector.h"

namespace itk
{
namespace ants
{
template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::AtroposSegmentationImageFilter()
{
  this->ProcessObject::SetNumberOfRequiredInputs(1);
  this->m_NumberOfIntensityImages = 1;

  this->m_NumberOfTissueClasses = 3;
  this->m_NumberOfPartialVolumeClasses = 0;
  this->m_PartialVolumeClasses.clear();
  this->m_UsePartialVolumeLikelihoods = false;

  this->m_MaximumNumberOfIterations = 5;
  this->m_ElapsedIterations = 0;
  this->m_CurrentPosteriorProbability = 0.0;
  this->m_ConvergenceThreshold = 0.001;

  this->m_InitializationStrategy = KMeans;
  this->m_InitialKMeansParameters.SetSize(0);

  this->m_PosteriorProbabilityFormulation = Socrates;
  this->m_UseMixtureModelProportions = true;

  this->m_PriorProbabilityWeight = 0.0;
  this->m_AdaptiveSmoothingWeights.clear();
  this->m_PriorLabelParameterMap.clear();
  this->m_ProbabilityThreshold = 0.0;
  this->m_PriorProbabilityImages.clear();
  this->m_PriorProbabilitySparseImages.clear();

  this->m_IntensityImages.clear();
  this->m_PriorLabelImage = nullptr;
  this->m_MaskImage = nullptr;

  this->m_MRFSmoothingFactor = 0.3;
  this->m_MRFRadius.Fill(1);

  this->m_SplineOrder = 3;
  this->m_NumberOfLevels.Fill(6);
  this->m_NumberOfControlPoints.Fill(this->m_SplineOrder + 1);

  this->m_MinimizeMemoryUsage = false;

  this->m_UseEuclideanDistanceForPriorLabels = false;
  this->m_PosteriorProbabilityImages.clear();
  this->m_DistancePriorProbabilityImages.clear();

  this->m_OutlierHandlingFilter = nullptr;

  this->m_RandomizerInitializationSeed = std::numeric_limits<RandomizerSeedType>::quiet_NaN();
  this->m_Randomizer = RandomizerType::New();
  this->m_Randomizer->Initialize();

  this->m_MaximumICMCode = 0;
  this->m_InitialAnnealingTemperature = 1.0;
  this->m_AnnealingRate = 1.0;
  this->m_MinimumAnnealingTemperature = 0.1;
  this->m_ICMCodeImage = nullptr;
  this->m_UseAsynchronousUpdating = true;
  this->m_MaximumNumberOfICMIterations = 1;
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::~AtroposSegmentationImageFilter() = default;

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::SetRandomizerInitializationSeed(
  const RandomizerSeedType seed)
{
  if (seed != this->m_RandomizerInitializationSeed)
  {
    this->m_RandomizerInitializationSeed = seed;
    this->m_Randomizer->Initialize(this->m_RandomizerInitializationSeed);
    this->Modified();
  }
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::SetMaskImage(const MaskImageType * mask)
{
  this->SetNthInput(1, const_cast<MaskImageType *>(mask));
  this->m_MaskImage = mask;
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
const typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::MaskImageType *
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::GetMaskImage() const
{
  //  const MaskImageType * maskImage =
  //    dynamic_cast<const MaskImageType *>( this->ProcessObject::GetInput( 1 ) );
  //
  //  return maskImage;

  return this->m_MaskImage;
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::SetPriorLabelImage(
  const ClassifiedImageType * prior)
{
  this->SetNthInput(2, const_cast<ClassifiedImageType *>(prior));
  this->m_PriorLabelImage = prior;
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
const typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::ClassifiedImageType *
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::GetPriorLabelImage() const
{
  //  const ClassifiedImageType * prior =
  //    dynamic_cast<const ClassifiedImageType *>(
  //    this->ProcessObject::GetInput( 2 ) );
  //
  //  return prior;

  return this->m_PriorLabelImage;
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::SetPriorProbabilityImage(
  unsigned int    whichClass,
  RealImageType * priorImage)
{
  if (whichClass < 1 || whichClass > this->m_NumberOfTissueClasses)
  {
    itkExceptionMacro("The requested prior probability image = " << whichClass << " should be in the range [1, "
                                                                 << this->m_NumberOfTissueClasses << "]");
  }
  if (this->m_MinimizeMemoryUsage)
  {
    // To make matters simpler, we force the index to be zero
    //   for each priorImage image.

    typename RealImageType::IndexType startIndex = priorImage->GetRequestedRegion().GetIndex();

    typename SparseImageType::Pointer sparsePriorImage = SparseImageType::New();
    sparsePriorImage->Initialize();

    unsigned long count = 0;

    ImageRegionConstIteratorWithIndex<RealImageType> It(priorImage, priorImage->GetRequestedRegion());
    for (It.GoToBegin(); !It.IsAtEnd(); ++It)
    {
      if (It.Get() > this->m_ProbabilityThreshold)
      {
        typename RealImageType::IndexType index = It.GetIndex();
        for (unsigned int d = 0; d < ImageDimension; d++)
        {
          index[d] -= startIndex[d];
        }

        unsigned long number = this->IndexToNumber(index, priorImage->GetRequestedRegion().GetSize());

        typename SparseImageType::PointType imageNumberIndex;
        imageNumberIndex[0] = number;

        sparsePriorImage->SetPoint(count, imageNumberIndex);
        sparsePriorImage->SetPointData(count, It.Get());
        count++;
      }
    }
    if (this->m_PriorProbabilitySparseImages.size() < whichClass)
    {
      this->m_PriorProbabilitySparseImages.resize(whichClass);
    }
    this->m_PriorProbabilitySparseImages[whichClass - 1] = sparsePriorImage;
  }
  else
  {
    if (this->m_PriorProbabilityImages.size() < whichClass)
    {
      this->m_PriorProbabilityImages.resize(whichClass);
    }
    this->m_PriorProbabilityImages[whichClass - 1] = priorImage;
  }
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::RealImagePointer
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::GetPriorProbabilityImage(
  unsigned int whichClass) const
{
  if (this->m_InitializationStrategy != PriorProbabilityImages)
  {
    return nullptr;
  }
  if (this->m_NumberOfPartialVolumeClasses == 0 && whichClass > this->m_NumberOfTissueClasses)
  {
    itkExceptionMacro("The requested prior probability image = " << whichClass << " should be in the range [1, "
                                                                 << this->m_NumberOfTissueClasses << "]");
  }
  else if (whichClass > this->m_NumberOfTissueClasses)
  {
    return nullptr;
  }

  if ((this->m_MinimizeMemoryUsage && this->m_PriorProbabilitySparseImages.size() != this->m_NumberOfTissueClasses) ||
      (!this->m_MinimizeMemoryUsage && this->m_PriorProbabilityImages.size() != this->m_NumberOfTissueClasses))
  {
    itkExceptionMacro("The number of prior probability images does not "
                      << "equal the number of classes.");
  }

  if (this->m_MinimizeMemoryUsage)
  {
    // To make matters simpler, we forced the prior probability images
    //   to all have indices of zeros.
    typename RealImageType::RegionType region;
    region.SetSize(this->GetIntensityImage(0)->GetRequestedRegion().GetSize());
    typename RealImageType::RegionType::IndexType index;
    for (unsigned int d = 0; d < ImageDimension; d++)
    {
      index[d] = 0;
    }
    region.SetIndex(index);

    RealImagePointer priorImage = AllocImage<RealImageType>(this->GetIntensityImage(0), 0);

    typename SparseImageType::Pointer sparsePriorImage = this->m_PriorProbabilitySparseImages[whichClass - 1];

    typename SparseImageType::PointsContainer::ConstIterator    It = sparsePriorImage->GetPoints()->Begin();
    typename SparseImageType::PointDataContainer::ConstIterator ItD = sparsePriorImage->GetPointData()->Begin();
    while (It != sparsePriorImage->GetPoints()->End())
    {
      unsigned long                     number = static_cast<unsigned long>(It.Value()[0]);
      typename RealImageType::IndexType index2 =
        this->NumberToIndex(number, priorImage->GetRequestedRegion().GetSize());
      priorImage->SetPixel(index2, ItD.Value());

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

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::SetIntensityImage(unsigned int      which,
                                                                                             const ImageType * image)
{
  if (which == 0)
  {
    this->SetInput(image);
  }
  else if (which + 1 > this->m_NumberOfIntensityImages)
  {
    this->m_NumberOfIntensityImages = which + 1;
    this->SetNthInput(2 + which, const_cast<ImageType *>(image));
  }

  // Since we need fast access to these inputs, we maintain a separate
  // pointer array.

  if (which + 1 > this->m_IntensityImages.size())
  {
    this->m_IntensityImages.resize(which + 1);
  }
  this->m_IntensityImages[which] = image;
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
const typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::ImageType *
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::GetIntensityImage(unsigned int which) const
{
  const ImageType * image;

  //  if( which == 0 )
  //    {
  //    image = dynamic_cast<const ImageType *>( this->ProcessObject::GetInput( 0 ) );
  //    }
  //  else if( which > 0 && which < this->m_NumberOfIntensityImages )
  //    {
  //    image = dynamic_cast<const ImageType *>(
  //      this->ProcessObject::GetInput( 2 + which ) );
  //    }

  if (which < this->m_NumberOfIntensityImages)
  {
    image = dynamic_cast<const ImageType *>(this->m_IntensityImages[which]);
  }
  else
  {
    itkExceptionMacro("Image " << which << " is outside the range "
                               << "[1,...,1 + m_NumberOfIntensityImages].");
  }
  return image;
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::AddPartialVolumeLabelSet(
  PartialVolumeLabelSetType labelSet)
{
  //
  // Three checks:
  //   1.  need to see if the labels are in {1, NumberOfTissueClasses}
  //   2.  need to determine if labelSet is a duplicate
  //   3.  check if each label is unique
  //
  typename PartialVolumeLabelSetType::const_iterator it;
  for (it = labelSet.begin(); it != labelSet.end(); ++it)
  {
    if (*it < 1 || *it > this->m_NumberOfTissueClasses)
    {
      itkWarningMacro("The label " << *it << " is outside the specified "
                                   << "range of the specified tissue class labels.");
      return;
    }
  }

  typename PartialVolumeClassesType::const_iterator itp;
  for (itp = this->m_PartialVolumeClasses.begin(); itp != this->m_PartialVolumeClasses.end(); ++itp)
  {
    bool isDuplicate = true;
    if (labelSet.size() == itp->size())
    {
      typename PartialVolumeLabelSetType::const_iterator itc;
      typename PartialVolumeLabelSetType::const_iterator itl;
      for (itc = itp->begin(), itl = labelSet.begin(); itc != itp->end(); ++itc, ++itl)
      {
        if (*itl != *itc)
        {
          isDuplicate = false;
          break;
        }
      }
    }
    if (isDuplicate)
    {
      itkWarningMacro("Duplicate label set.");
      return;
    }
  }
  for (LabelType l = 1; l <= this->m_NumberOfTissueClasses; l++)
  {
    unsigned int cardinality = std::count(labelSet.begin(), labelSet.end(), l);
    if (cardinality > 1)
    {
      itkWarningMacro("Duplicate label " << l);
      return;
    }
  }

  this->m_PartialVolumeClasses.push_back(labelSet);
  this->Modified();
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::GenerateData()
{
  //
  // Assign Gaussian likelihood functions if mixture model components are absent
  //
  typedef ants::Statistics::GaussianListSampleFunction<SampleType, RealType, RealType> LikelihoodType;
  for (unsigned int n = 0; n < this->m_NumberOfTissueClasses; n++)
  {
    if (!this->GetLikelihoodFunction(n))
    {
      typename LikelihoodType::Pointer gaussianLikelihood = LikelihoodType::New();
      this->SetLikelihoodFunction(n, gaussianLikelihood);
    }
  }

  if (this->m_UsePartialVolumeLikelihoods)
  {
    this->m_NumberOfPartialVolumeClasses = this->m_PartialVolumeClasses.size();
  }
  else
  {
    this->m_NumberOfPartialVolumeClasses = 0;
  }

  //
  // Initialize the class labeling and the likelihood models
  //
  this->GenerateInitialClassLabeling();

  //
  // Iterate until convergence or iterative exhaustion.
  //
  IterationReporter reporter(this, 0, 1);

  //
  // Get spacing for mrf neighborhood evaluation loop
  //
  this->m_ImageSpacing = this->GetOutput()->GetSpacing();

  bool isConverged = false;
  this->m_CurrentPosteriorProbability = 0.0;
  RealType probabilityOld = NumericTraits<RealType>::NonpositiveMin();

  unsigned int iteration = 0;
  while (!isConverged && iteration++ < this->m_MaximumNumberOfIterations)
  {
    reporter.CompletedStep();

    this->m_CurrentPosteriorProbability = this->UpdateClassLabeling();

    if (this->m_CurrentPosteriorProbability - probabilityOld < this->m_ConvergenceThreshold)
    {
      isConverged = true;
    }
    probabilityOld = this->m_CurrentPosteriorProbability;

    itkDebugMacro("Iteration: " << this->m_ElapsedIterations
                                << ", "
                                   "current posterior probability = "
                                << this->m_CurrentPosteriorProbability);

    this->m_ElapsedIterations++;

    //
    // Clear the current posterior probability images to force
    // recalculation of the posterior probability images.
    //
    this->m_PosteriorProbabilityImages.clear();
  }
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::GenerateInitialClassLabeling()
{
  this->AllocateOutputs();
  this->GetOutput()->FillBuffer(NumericTraits<MaskLabelType>::ZeroValue());

  switch (this->m_InitializationStrategy)
  {
    case Random:
    {
      ImageRegionIterator<ClassifiedImageType> It(this->GetOutput(), this->GetOutput()->GetRequestedRegion());
      for (It.GoToBegin(); !It.IsAtEnd(); ++It)
      {
        LabelType label = this->m_Randomizer->GetIntegerVariate(this->m_NumberOfTissueClasses - 1) + 1;
        It.Set(label);
      }
      break;
    }
    case KMeans:
    default:
    {
      this->GenerateInitialClassLabelingWithKMeansClustering();
    }
    break;
    case Otsu:
    {
      this->GenerateInitialClassLabelingWithOtsuThresholding();
    }
    break;
    case PriorProbabilityImages:
    {
      this->GenerateInitialClassLabelingWithPriorProbabilityImages();
    }
    break;
    case PriorLabelImage:
    {
      if (this->GetMaskImage())
      {
        typedef MaskImageFilter<ClassifiedImageType, MaskImageType, ClassifiedImageType> MaskerType;
        typename MaskerType::Pointer                                                     masker = MaskerType::New();
        masker->SetInput1(this->GetPriorLabelImage());
        masker->SetInput2(this->GetMaskImage());
        masker->Update();

        this->SetNthOutput(0, masker->GetOutput());
      }
      else
      {
        typedef ImageDuplicator<ClassifiedImageType> DuplicatorType;
        typename DuplicatorType::Pointer             duplicator = DuplicatorType::New();
        duplicator->SetInputImage(this->GetPriorLabelImage());
        duplicator->Update();

        this->SetNthOutput(0, duplicator->GetOutput());
      }
    }
    break;
  }

  //
  // Calculate the initial parameters of the mixture model from the
  // initial labeling, i.e. the proportion, mean, and covariance for each label.
  //
  unsigned int totalNumberOfClasses = this->m_NumberOfTissueClasses + this->m_NumberOfPartialVolumeClasses;

  this->m_MixtureModelProportions.SetSize(totalNumberOfClasses);

  this->m_LabelVolumes.SetSize(totalNumberOfClasses);

  unsigned int totalSampleSize = 0;

  std::vector<typename SampleType::Pointer> samples;
  for (unsigned int n = 0; n < totalNumberOfClasses; n++)
  {
    typename SampleType::Pointer sample = SampleType::New();
    sample->SetMeasurementVectorSize(this->m_NumberOfIntensityImages);
    samples.push_back(sample);
  }

  //
  // Accumulate the sample array for all labels.  Also accumulate the
  // prior probability weights, if applicable.
  //
  ImageRegionIteratorWithIndex<ClassifiedImageType> ItO(this->GetOutput(), this->GetOutput()->GetRequestedRegion());
  for (ItO.GoToBegin(); !ItO.IsAtEnd(); ++ItO)
  {
    LabelType label = ItO.Get();
    if (label == 0)
    {
      continue;
    }
    typename SampleType::MeasurementVectorType measurement;
    NumericTraits<MeasurementVectorType>::SetLength(measurement, this->m_NumberOfIntensityImages);
    for (unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++)
    {
      measurement[i] = this->GetIntensityImage(i)->GetPixel(ItO.GetIndex());
    }
    samples[label - 1]->PushBack(measurement);
  }

  //
  // Create the weight array now that we know the sample sizes.
  //
  Array<unsigned int> count(totalNumberOfClasses);
  count.Fill(0);
  std::vector<WeightArrayType> weights;
  for (unsigned int n = 0; n < totalNumberOfClasses; n++)
  {
    totalSampleSize += samples[n]->Size();
    WeightArrayType weightArray(samples[n]->Size());
    weightArray.Fill(1.0);
    weights.push_back(weightArray);

    this->m_LabelVolumes[n] = samples[n]->Size();
  }
  if (this->m_InitializationStrategy == PriorProbabilityImages)
  {
    for (ItO.GoToBegin(); !ItO.IsAtEnd(); ++ItO)
    {
      LabelType label = ItO.Get();
      if (label == 0 || label > this->m_NumberOfTissueClasses)
      {
        continue;
      }
      RealImagePointer priorProbabilityImage = this->GetPriorProbabilityImage(label);
      weights[label - 1].SetElement(count[label - 1]++, priorProbabilityImage->GetPixel(ItO.GetIndex()));
    }
  }

  for (unsigned int n = 0; n < totalNumberOfClasses; n++)
  {
    if (n < this->m_NumberOfTissueClasses)
    {
      this->m_MixtureModelComponents[n]->SetListSampleWeights(&weights[n]);
      this->m_MixtureModelComponents[n]->SetInputListSample(samples[n]);
      // this->m_MixtureModelComponents[n]->ClearInputListSample();

      if (this->m_UseMixtureModelProportions)
      {
        this->m_MixtureModelProportions[n] =
          static_cast<RealType>(samples[n]->Size()) / static_cast<RealType>(totalSampleSize);
      }
      else
      {
        this->m_MixtureModelProportions[n] =
          NumericTraits<RealType>::OneValue() / static_cast<RealType>(totalNumberOfClasses);
      }
    }
    else
    {
      PartialVolumeLabelSetType labelSet = this->m_PartialVolumeClasses[n - this->m_NumberOfTissueClasses];
      for (unsigned d = 0; d < labelSet.size(); d++)
      {
        this->m_MixtureModelComponents[n]->SetListSampleWeights(d, &weights[labelSet[d] - 1]);
        this->m_MixtureModelComponents[n]->SetIndexedInputListSample(d, samples[labelSet[d] - 1]);
        // this->m_MixtureModelComponents[n]->ClearInputListSample(d);
      }

      this->m_MixtureModelProportions[n] = 0.0;
    }
  }
  for (unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++)
  {
    ControlPointLatticeContainerType container;
    this->m_ControlPointLattices.push_back(container);
    for (unsigned int n = 0; n < this->m_NumberOfTissueClasses; n++)
    {
      this->m_ControlPointLattices[i].push_back(nullptr);
    }
  }
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::
  GenerateInitialClassLabelingWithPriorProbabilityImages()
{
  // We first normalize prior probability images by keeping track of
  //   1. the sum of the prior probabilities at each pixel
  //   2. the max prior probability value
  //   3. which prior probability image corresponds to the max prior value

  RealImagePointer sumPriorProbabilityImage = AllocImage<RealImageType>(this->GetInput(), 0);

  RealImagePointer maxPriorProbabilityImage = AllocImage<RealImageType>(this->GetInput(), 0);

  for (unsigned int n = 0; n < this->m_NumberOfTissueClasses; n++)
  {
    RealImagePointer priorProbabilityImage = this->GetPriorProbabilityImage(n + 1);

    typedef AddImageFilter<RealImageType, RealImageType, RealImageType> AdderType;
    typename AdderType::Pointer                                         adder = AdderType::New();
    adder->SetInput1(sumPriorProbabilityImage);
    adder->SetInput2(priorProbabilityImage);
    adder->Update();

    sumPriorProbabilityImage = adder->GetOutput();

    ImageRegionIteratorWithIndex<RealImageType> ItP(priorProbabilityImage, priorProbabilityImage->GetRequestedRegion());
    ImageRegionIterator<RealImageType> ItM(maxPriorProbabilityImage, maxPriorProbabilityImage->GetRequestedRegion());
    ImageRegionIterator<ClassifiedImageType> ItO(this->GetOutput(), this->GetOutput()->GetRequestedRegion());

    ItP.GoToBegin();
    ItM.GoToBegin();
    ItO.GoToBegin();
    while (!ItP.IsAtEnd())
    {
      if (!this->GetMaskImage() ||
          this->GetMaskImage()->GetPixel(ItP.GetIndex()) != NumericTraits<MaskLabelType>::ZeroValue())
      {
        if (ItP.Get() > ItM.Get())
        {
          ItM.Set(ItP.Get());
          ItO.Set(n + 1);
        }
        else if (Math::FloatAlmostEqual(ItP.Get(), ItM.Get()))
        {
          if (n == 0)
          {
            ItO.Set(1);
          }
          else // if maximum probabilities are the same, randomly select one
          {
            if (this->m_Randomizer->GetIntegerVariate(1))
            {
              ItO.Set(n + 1);
            }
          }
        }
      }
      ++ItP;
      ++ItM;
      ++ItO;
    }
  }
  // Now we can normalize each prior probability image by dividing by the sum
  for (unsigned int n = 0; n < this->m_NumberOfTissueClasses; n++)
  {
    RealImagePointer priorProbabilityImage = this->GetPriorProbabilityImage(n + 1);

    ImageRegionIteratorWithIndex<RealImageType> ItP(priorProbabilityImage, priorProbabilityImage->GetRequestedRegion());
    ImageRegionIterator<RealImageType> ItS(sumPriorProbabilityImage, sumPriorProbabilityImage->GetRequestedRegion());
    ImageRegionIterator<RealImageType> ItM(maxPriorProbabilityImage, maxPriorProbabilityImage->GetRequestedRegion());
    ImageRegionIterator<ClassifiedImageType> ItO(this->GetOutput(), this->GetOutput()->GetRequestedRegion());

    ItP.GoToBegin();
    ItS.GoToBegin();
    ItM.GoToBegin();
    ItO.GoToBegin();

    while (!ItP.IsAtEnd())
    {
      if (!this->GetMaskImage() ||
          this->GetMaskImage()->GetPixel(ItP.GetIndex()) != NumericTraits<MaskLabelType>::ZeroValue())
      {
        if (ItM.Get() <= this->m_ProbabilityThreshold ||
            Math::FloatAlmostEqual(ItS.Get(), NumericTraits<typename ImageType::PixelType>::ZeroValue()))
        {
          ItO.Set(NumericTraits<LabelType>::ZeroValue());
          ItP.Set(NumericTraits<RealType>::ZeroValue());
        }
        else
        {
          ItP.Set(ItP.Get() / ItS.Get());
        }
      }
      ++ItP;
      ++ItS;
      ++ItM;
      ++ItO;
    }

    this->SetPriorProbabilityImage(n + 1, priorProbabilityImage);
  }

  //
  // Set the initial output to be the prior label image.  This way we can
  // propogate the segmentation solution to regions of non-zero probability
  // but where the mask exists.
  //
  typedef ImageDuplicator<ClassifiedImageType> DuplicatorType;
  typename DuplicatorType::Pointer             duplicator = DuplicatorType::New();
  duplicator->SetInputImage(this->GetOutput());
  duplicator->Update();

  this->SetPriorLabelImage(duplicator->GetOutput());
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::
  GenerateInitialClassLabelingWithOtsuThresholding()
{
  RealType maxValue = NumericTraits<RealType>::min();
  RealType minValue = NumericTraits<RealType>::max();

  ImageRegionConstIteratorWithIndex<ImageType> ItI(this->GetInput(), this->GetInput()->GetRequestedRegion());
  for (ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI)
  {
    if (!this->GetMaskImage() ||
        this->GetMaskImage()->GetPixel(ItI.GetIndex()) != NumericTraits<MaskLabelType>::ZeroValue())
    {
      if (ItI.Get() < minValue)
      {
        minValue = ItI.Get();
      }
      else if (ItI.Get() > maxValue)
      {
        maxValue = ItI.Get();
      }
    }
  }

  typedef LabelStatisticsImageFilter<ImageType, MaskImageType> StatsType;
  typename StatsType::Pointer                                  stats = StatsType::New();
  stats->SetInput(this->GetInput());
  if (this->GetMaskImage())
  {
    stats->SetLabelInput(const_cast<MaskImageType *>(this->GetMaskImage()));
  }
  else
  {
    typename MaskImageType::Pointer maskImage =
      AllocImage<MaskImageType>(this->GetOutput(), NumericTraits<LabelType>::OneValue());
    stats->SetLabelInput(maskImage);
  }
  stats->UseHistogramsOn();
  stats->SetHistogramParameters(200, minValue, maxValue);
  stats->Update();

  typedef OtsuMultipleThresholdsCalculator<typename StatsType::HistogramType> OtsuType;
  typename OtsuType::Pointer                                                  otsu = OtsuType::New();
  otsu->SetInputHistogram(stats->GetHistogram(NumericTraits<LabelType>::OneValue()));
  otsu->SetNumberOfThresholds(this->m_NumberOfTissueClasses - 1);
  otsu->Compute();

  typename OtsuType::OutputType thresholds = otsu->GetOutput();

  ImageRegionIterator<ClassifiedImageType> ItO(this->GetOutput(), this->GetOutput()->GetRequestedRegion());
  for (ItI.GoToBegin(), ItO.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItO)
  {
    LabelType label = NumericTraits<LabelType>::ZeroValue();
    if (!this->GetMaskImage() ||
        this->GetMaskImage()->GetPixel(ItI.GetIndex()) != NumericTraits<MaskLabelType>::ZeroValue())
    {
      if (ItI.Get() < static_cast<typename ImageType::PixelType>(thresholds[0]))
      {
        label = NumericTraits<LabelType>::OneValue();
      }
      else
      {
        bool thresholdFound = false;
        for (unsigned int i = 1; i < thresholds.size(); i++)
        {
          if (ItI.Get() >= static_cast<typename ImageType::PixelType>(thresholds[i - 1]) &&
              ItI.Get() <= static_cast<typename ImageType::PixelType>(thresholds[i]))
          {
            label = static_cast<LabelType>(i + 1);
            thresholdFound = true;
            break;
          }
        }
        if (!thresholdFound)
        {
          label = static_cast<LabelType>(thresholds.size() + 1);
        }
      }
    }
    ItO.Set(label);
  }
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::
  GenerateInitialClassLabelingWithKMeansClustering()
{
  //
  // We first perform kmeans on the first image and use the results to
  // seed the second run of kmeans using all the images.
  //
  typedef LabelStatisticsImageFilter<ImageType, MaskImageType> StatsType;
  typename StatsType::Pointer                                  stats = StatsType::New();
  stats->SetInput(this->GetInput());
  if (this->GetMaskImage())
  {
    stats->SetLabelInput(const_cast<MaskImageType *>(this->GetMaskImage()));
  }
  else
  {
    typename MaskImageType::Pointer maskImage =
      AllocImage<MaskImageType>(this->GetOutput(), NumericTraits<LabelType>::OneValue());
    stats->SetLabelInput(maskImage);
  }
  stats->UseHistogramsOff();
  stats->Update();

  const RealType minValue = stats->GetMinimum(NumericTraits<LabelType>::OneValue());
  const RealType maxValue = stats->GetMaximum(NumericTraits<LabelType>::OneValue());

  //
  // The code below can be replaced by itkListSampleToImageFilter when we
  // migrate over to the Statistics classes current in the Review/ directory.
  //
  typename SampleType::Pointer sample = SampleType::New();
  sample->SetMeasurementVectorSize(1);

  ImageRegionConstIteratorWithIndex<ImageType> ItI(this->GetInput(), this->GetInput()->GetRequestedRegion());
  for (ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI)
  {
    if (!this->GetMaskImage() ||
        this->GetMaskImage()->GetPixel(ItI.GetIndex()) != NumericTraits<MaskLabelType>::ZeroValue())
    {
      typename SampleType::MeasurementVectorType measurement;
      measurement.SetSize(1);
      measurement[0] = ItI.Get();
      sample->PushBack(measurement);
    }
  }

  typedef itk::Statistics::WeightedCentroidKdTreeGenerator<SampleType> TreeGeneratorType;
  typename TreeGeneratorType::Pointer                                  treeGenerator = TreeGeneratorType::New();
  treeGenerator->SetSample(sample);
  treeGenerator->SetBucketSize(16);
  treeGenerator->Update();

  typedef typename TreeGeneratorType::KdTreeType                TreeType;
  typedef itk::Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;
  typename EstimatorType::Pointer                               estimator = EstimatorType::New();
  estimator->SetKdTree(treeGenerator->GetOutput());
  estimator->SetMaximumIteration(200);
  estimator->SetCentroidPositionChangesThreshold(0.0);

  typename EstimatorType::ParametersType initialMeans(this->m_NumberOfTissueClasses);

  //
  // If the initial KMeans parameters are not set, guess initial class means by
  // dividing the dynamic range of the first image into equal intervals.
  //
  if (this->m_InitialKMeansParameters.Size() == this->m_NumberOfTissueClasses)
  {
    for (unsigned int n = 0; n < this->m_NumberOfTissueClasses; n++)
    {
      initialMeans[n] = this->m_InitialKMeansParameters[n];
    }
  }
  else
  {
    for (unsigned int n = 0; n < this->m_NumberOfTissueClasses; n++)
    {
      initialMeans[n] = minValue + (maxValue - minValue) * (static_cast<RealType>(n) + static_cast<RealType>(0.5)) /
                                     static_cast<RealType>(this->m_NumberOfTissueClasses);
    }
  }
  estimator->SetParameters(initialMeans);
  estimator->StartOptimization();

  //
  // Classify the samples
  //
  typedef itk::Statistics::MinimumDecisionRule DecisionRuleType;
  typename DecisionRuleType::Pointer           decisionRule = DecisionRuleType::New();

  typedef itk::Statistics::SampleClassifierFilter<SampleType> ClassifierType;
  typename ClassifierType::Pointer                            classifier = ClassifierType::New();
  classifier->SetDecisionRule(decisionRule);
  classifier->SetInput(sample);
  classifier->SetNumberOfClasses(this->m_NumberOfTissueClasses);

  typename ClassifierType::ClassLabelVectorObjectType::Pointer classLabels =
    ClassifierType::ClassLabelVectorObjectType::New();
  classifier->SetClassLabels(classLabels);
  typename ClassifierType::ClassLabelVectorType & classLabelVector = classLabels->Get();

  //
  // Order the cluster means so that the lowest mean of the input image
  // corresponds to label '1', the second lowest to label '2', etc.
  //
  std::vector<RealType> estimatorParameters;
  for (unsigned int n = 0; n < this->m_NumberOfTissueClasses; n++)
  {
    estimatorParameters.push_back(estimator->GetParameters()[n]);
  }
  std::sort(estimatorParameters.begin(), estimatorParameters.end());

  typedef itk::Statistics::DistanceToCentroidMembershipFunction<MeasurementVectorType> MembershipFunctionType;
  typename ClassifierType::MembershipFunctionVectorObjectType::Pointer                 membershipFunctions =
    ClassifierType::MembershipFunctionVectorObjectType::New();
  typename ClassifierType::MembershipFunctionVectorType & membershipFunctionsVector = membershipFunctions->Get();

  classifier->SetMembershipFunctions(membershipFunctions);
  for (unsigned int n = 0; n < this->m_NumberOfTissueClasses; n++)
  {
    typename MembershipFunctionType::Pointer membershipFunction = MembershipFunctionType::New();
    membershipFunction->SetMeasurementVectorSize(sample->GetMeasurementVectorSize());
    typename MembershipFunctionType::CentroidType centroid;
    NumericTraits<typename MembershipFunctionType::CentroidType>::SetLength(centroid,
                                                                            sample->GetMeasurementVectorSize());
    centroid[0] = estimatorParameters[n];
    membershipFunction->SetCentroid(centroid);
    membershipFunctionsVector.push_back(membershipFunction.GetPointer());

    classLabelVector.push_back(static_cast<typename ClassifierType::ClassLabelType>(n + 1));
  }
  classifier->Update();

  //
  // Classify the voxels
  //
  typedef typename ClassifierType::MembershipSampleType ClassifierOutputType;
  typedef typename ClassifierOutputType::ConstIterator  LabelIterator;

  ImageRegionIteratorWithIndex<ClassifiedImageType> ItO(this->GetOutput(), this->GetOutput()->GetRequestedRegion());
  ItO.GoToBegin();
  LabelIterator it = classifier->GetOutput()->Begin();
  while (it != classifier->GetOutput()->End())
  {
    if (!this->GetMaskImage() ||
        this->GetMaskImage()->GetPixel(ItO.GetIndex()) != NumericTraits<MaskLabelType>::ZeroValue())
    {
      ItO.Set(it.GetClassLabel());
      ++it;
    }
    else
    {
      ItO.Set(NumericTraits<LabelType>::ZeroValue());
    }
    ++ItO;
  }

  //
  // If there are more than one intensity images, use the results from the first
  // kmeans grouping to seed a second invoking of the algorithm using both
  // the input image and the auxiliary images.
  //
  if (this->m_NumberOfIntensityImages > 1)
  {
    VariableSizeMatrix<RealType> classMeanValues;
    classMeanValues.SetSize(this->m_NumberOfIntensityImages, this->m_NumberOfTissueClasses);
    for (unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++)
    {
      typedef LabelStatisticsImageFilter<ImageType, ClassifiedImageType> ClassStatsType;
      typename ClassStatsType::Pointer                                   stats2 = ClassStatsType::New();
      stats2->SetInput(this->GetIntensityImage(i));
      stats2->SetLabelInput(this->GetOutput());
      stats2->UseHistogramsOff();
      stats2->Update();
      for (unsigned int j = 0; j < this->m_NumberOfTissueClasses; j++)
      {
        classMeanValues(i, j) = stats->GetMean(j + 1);
      }
    }

    typename EstimatorType::Pointer        estimator2 = EstimatorType::New();
    typename EstimatorType::ParametersType initialMeans2(this->m_NumberOfTissueClasses *
                                                         (this->m_NumberOfIntensityImages));
    initialMeans2.Fill(0.0);
    for (unsigned int n = 0; n < this->m_NumberOfTissueClasses; n++)
    {
      for (unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++)
      {
        initialMeans2[this->m_NumberOfIntensityImages * n + i] = classMeanValues(i, n);
      }
    }

    typename SampleType::Pointer sample2 = SampleType::New();
    sample2->SetMeasurementVectorSize(this->m_NumberOfIntensityImages);
    for (ItO.GoToBegin(); !ItO.IsAtEnd(); ++ItO)
    {
      if (!this->GetMaskImage() ||
          this->GetMaskImage()->GetPixel(ItO.GetIndex()) != NumericTraits<LabelType>::ZeroValue())
      {
        typename SampleType::MeasurementVectorType measurement;
        measurement.SetSize(this->m_NumberOfIntensityImages);
        for (unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++)
        {
          measurement[i] = this->GetIntensityImage(i)->GetPixel(ItO.GetIndex());
        }
        sample2->PushBack(measurement);
      }
    }

    typename TreeGeneratorType::Pointer treeGenerator2 = TreeGeneratorType::New();
    treeGenerator2->SetSample(sample2);
    treeGenerator2->SetBucketSize(16);
    treeGenerator2->Update();

    estimator2->SetParameters(initialMeans2);
    estimator2->SetKdTree(treeGenerator2->GetOutput());
    estimator2->SetMaximumIteration(200);
    estimator2->SetCentroidPositionChangesThreshold(0.0);
    estimator2->StartOptimization();

    //
    // Classify the samples
    //
    typename DecisionRuleType::Pointer decisionRule2 = DecisionRuleType::New();

    typename ClassifierType::Pointer classifier2 = ClassifierType::New();
    classifier2->SetDecisionRule(decisionRule2);
    classifier2->SetInput(sample2);
    classifier2->SetNumberOfClasses(this->m_NumberOfTissueClasses);

    typename ClassifierType::ClassLabelVectorObjectType::Pointer classLabels2 =
      ClassifierType::ClassLabelVectorObjectType::New();
    classifier2->SetClassLabels(classLabels2);
    typename ClassifierType::ClassLabelVectorType & classLabelVector2 = classLabels2->Get();

    //
    // Order the cluster means so that the lowest mean of the input image
    // corresponds to label '1', the second lowest to label '2', etc.
    //
    std::vector<RealType> estimatorParameters2;
    for (unsigned int n = 0; n < this->m_NumberOfTissueClasses; n++)
    {
      estimatorParameters2.push_back(estimator2->GetParameters()[n]);
    }
    std::sort(estimatorParameters2.begin(), estimatorParameters2.end());

    typename ClassifierType::MembershipFunctionVectorObjectType::Pointer membershipFunctions2 =
      ClassifierType::MembershipFunctionVectorObjectType::New();
    typename ClassifierType::MembershipFunctionVectorType & membershipFunctionsVector2 = membershipFunctions2->Get();

    classifier2->SetMembershipFunctions(membershipFunctions2);
    for (unsigned int n = 0; n < this->m_NumberOfTissueClasses; n++)
    {
      typename MembershipFunctionType::Pointer membershipFunction2 = MembershipFunctionType::New();
      membershipFunction2->SetMeasurementVectorSize(sample2->GetMeasurementVectorSize());
      typename MembershipFunctionType::CentroidType centroid;
      NumericTraits<typename MembershipFunctionType::CentroidType>::SetLength(centroid,
                                                                              sample2->GetMeasurementVectorSize());
      for (unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++)
      {
        centroid[i] = estimator2->GetParameters()[(this->m_NumberOfIntensityImages) * n + i];
      }
      membershipFunction2->SetCentroid(centroid);
      membershipFunctionsVector2.push_back(membershipFunction2.GetPointer());

      classLabelVector2.push_back(static_cast<typename ClassifierType::ClassLabelType>(n + 1));
    }
    classifier2->Update();

    //
    // Classify the voxels
    //
    ItO.GoToBegin();
    LabelIterator it2 = classifier->GetOutput()->Begin();
    while (it2 != classifier->GetOutput()->End())
    {
      if (!this->GetMaskImage() ||
          this->GetMaskImage()->GetPixel(ItO.GetIndex()) != NumericTraits<LabelType>::ZeroValue())
      {
        ItO.Set(it2.GetClassLabel());
        ++it2;
      }
      else
      {
        ItO.Set(NumericTraits<LabelType>::ZeroValue());
      }
      ++ItO;
    }
  }
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::RealType
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::UpdateClassLabeling()
{
  double maxPosteriorSum = 0.0;

  if (this->m_UseAsynchronousUpdating)
  {
    if (this->m_MaximumICMCode == 0)
    {
      this->ComputeICMCodeImage();
    }

    typename NeighborhoodIterator<ClassifiedImageType>::RadiusType radius;
    for (unsigned int d = 0; d < ImageDimension; d++)
    {
      radius[d] = this->m_MRFRadius[d];
    }
    NeighborhoodIterator<ClassifiedImageType> ItO(radius, this->GetOutput(), this->GetOutput()->GetRequestedRegion());

    maxPosteriorSum = 0.0;
    double     oldMaxPosteriorSum = -1.0;
    unsigned int numberOfIterations = 0;
    while (maxPosteriorSum > oldMaxPosteriorSum && numberOfIterations++ < this->m_MaximumNumberOfICMIterations)
    {
      itkDebugMacro("ICM iteration: " << numberOfIterations);

      oldMaxPosteriorSum = maxPosteriorSum;

      maxPosteriorSum = 0.0;

      // Iterate randomly through the ICM codes so we don't visit the same codes
      // in the same order every iteration.  We use the Fisher-Yates shuffle to
      // create a random shuffling between 1 and m_MaximumICMCode.

      Array<LabelType> icmCodeSet(this->m_MaximumICMCode);
      for (unsigned int n = 1; n <= this->m_MaximumICMCode; n++)
      {
        icmCodeSet[n - 1] = static_cast<LabelType>(n);
      }
      for (int i = icmCodeSet.Size() - 1; i > 0; i--)
      {
        unsigned int j = this->m_Randomizer->GetIntegerVariate(i);
        LabelType    tmp = icmCodeSet[i];
        icmCodeSet[i] = icmCodeSet[j];
        icmCodeSet[j] = tmp;
      }
      for (unsigned int n = 0; n < icmCodeSet.Size(); n++)
      {
        for (ItO.GoToBegin(); !ItO.IsAtEnd(); ++ItO)
        {
          if (this->m_ICMCodeImage->GetPixel(ItO.GetIndex()) == icmCodeSet[n])
          {
            maxPosteriorSum += this->PerformLocalLabelingUpdate(ItO);
          }
        }
      }
      itkDebugMacro("ICM posterior probability sum: " << maxPosteriorSum);
    }
  }

  RealImagePointer maxPosteriorProbabilityImage =
    AllocImage<RealImageType>(this->GetOutput(), NumericTraits<RealType>::ZeroValue());

  typename ClassifiedImageType::Pointer maxLabels =
    AllocImage<ClassifiedImageType>(this->GetOutput(), NumericTraits<LabelType>::ZeroValue());

  unsigned int totalNumberOfClasses = this->m_NumberOfTissueClasses + this->m_NumberOfPartialVolumeClasses;

  Array<double> sumPosteriors(totalNumberOfClasses);
  sumPosteriors.Fill(0.0);

  typename SampleType::Pointer sample = SampleType::New();
  sample = this->GetScalarSamples();
  unsigned long totalSampleSize = sample->Size();
  for (unsigned int n = 0; n < totalNumberOfClasses; n++)
  {
    RealImagePointer posteriorProbabilityImage = this->GetPosteriorProbabilityImage(n + 1);

    ImageRegionIteratorWithIndex<ClassifiedImageType> ItO(maxLabels, maxLabels->GetRequestedRegion());
    ImageRegionConstIterator<RealImageType>           ItP(posteriorProbabilityImage,
                                                posteriorProbabilityImage->GetRequestedRegion());
    ImageRegionIterator<RealImageType>                ItM(maxPosteriorProbabilityImage,
                                           maxPosteriorProbabilityImage->GetRequestedRegion());

    WeightArrayType weights(totalSampleSize);

    unsigned long count = 0;

    ItP.GoToBegin();
    ItM.GoToBegin();
    ItO.GoToBegin();
    while (!ItP.IsAtEnd())
    {
      if (!this->GetMaskImage() ||
          this->GetMaskImage()->GetPixel(ItO.GetIndex()) != NumericTraits<LabelType>::ZeroValue())
      {
        RealType posteriorProbability = ItP.Get();
        weights.SetElement(count++, posteriorProbability);

        // The following commented lines enforce "hard EM" as opposed to "soft EM"
        // which uses probabilities.
        // if( this->GetOutput()->GetPixel( ItP.GetIndex() ) == n + 1 )
        //   {
        //   posteriorProbability = 1.0;
        //   }
        // else
        //   {
        //   posteriorProbability = 0.0;
        //   }

        if (posteriorProbability > ItM.Get())
        {
          ItM.Set(posteriorProbability);
          ItO.Set(static_cast<LabelType>(n + 1));
        }
        else if (Math::FloatAlmostEqual(posteriorProbability, ItM.Get()))
        {
          LabelType currentLabel = ItO.Get();
          if (currentLabel >= 1 && this->m_LabelVolumes[n] < this->m_LabelVolumes[currentLabel - 1])
          {
            ItO.Set(static_cast<LabelType>(n + 1));
          }
        }
        sumPosteriors[n] += static_cast<double>(posteriorProbability);
      }
      ++ItP;
      ++ItM;
      ++ItO;
    }

    if (n < this->m_NumberOfTissueClasses)
    {
      this->m_MixtureModelComponents[n]->SetListSampleWeights(&weights);
      this->m_MixtureModelComponents[n]->SetInputListSample(sample);
      // this->m_MixtureModelComponents[n]->ClearInputListSample();
    }
    else
    {
      PartialVolumeLabelSetType labelSet = this->m_PartialVolumeClasses[n - this->m_NumberOfTissueClasses];
      for (unsigned d = 0; d < labelSet.size(); d++)
      {
        if (n == labelSet[d] - 1)
        {
          this->m_MixtureModelComponents[n]->SetListSampleWeights(d, &weights);
          this->m_MixtureModelComponents[n]->SetIndexedInputListSample(d, sample);
          // this->m_MixtureModelComponents[n]->ClearInputListSample(d);
        }
      }
      this->m_MixtureModelProportions[n] = 0.0;
    }

    if (this->m_UseMixtureModelProportions)
    {
      this->m_MixtureModelProportions[n] = static_cast<RealType>(sumPosteriors[n] / static_cast<double>(totalSampleSize));
    }
    else
    {
      this->m_MixtureModelProportions[n] =
        NumericTraits<RealType>::OneValue() / static_cast<RealType>(totalNumberOfClasses);
    }
  }
  // Check to see if the resulting number of classes is equal to the total number of classes.
  using StatisticsImageFilterType = itk::StatisticsImageFilter<ClassifiedImageType>;
  auto statisticsImageFilter = StatisticsImageFilterType::New();
  statisticsImageFilter->SetInput(maxLabels);
  statisticsImageFilter->Update();

  if (statisticsImageFilter->GetMaximum() != totalNumberOfClasses)
  {
    itkExceptionMacro("Updating the class labeling resulted in " << statisticsImageFilter->GetMaximum() <<
                      " non-zero label sets but " << totalNumberOfClasses << " is requested.");
  }

  auto labelMapFilter = itk::LabelImageToShapeLabelMapFilter<ClassifiedImageType>::New();
  labelMapFilter->SetInput(maxLabels);
  labelMapFilter->SetComputeOrientedBoundingBox(false);
  labelMapFilter->SetComputePerimeter(false);
  labelMapFilter->SetComputeFeretDiameter(false);
  labelMapFilter->SetInput(maxLabels);
  labelMapFilter->Update();

  for (unsigned int n = 0; n < totalNumberOfClasses; n++)
  {
    auto labelObject = labelMapFilter->GetOutput()->GetNthLabelObject(n);
    this->m_LabelVolumes[n] = labelObject->GetNumberOfPixels();
  }

  this->SetNthOutput(0, maxLabels);

  //  The commented code below is used to calculate the mixture model proportions
  //  according to the formulae given in Ashburner et al., "Unified Segmentation",
  //  Neuroimage, 2005 Jul 1;26(3):839-51.
  //
  //    if( this->m_UseMixtureModelProportions )
  //      {
  //      //
  //      // Perform the following calculation as a preprocessing step to update the
  //      // class proportions.
  //      //
  //      RealImagePointer distancePriorProbabilityImage =
  //        this->GetDistancePriorProbabilityImage( n + 1 );
  //      RealImagePointer priorProbabilityImage =
  //        this->GetPriorProbabilityImage( n + 1 );
  //
  //      ImageRegionIteratorWithIndex<RealImageType> ItW(
  //        weightedPriorProbabilityImage,
  //        weightedPriorProbabilityImage->GetRequestedRegion() );
  //      for( ItW.GoToBegin(); !ItW.IsAtEnd(); ++ItW )
  //        {
  //        if( !this->GetMaskImage() ||
  //            this->GetMaskImage()->GetPixel( ItW.GetIndex() ) != NumericTraits<LabelType>::ZeroValue() )
  //          {
  //          RealType priorProbability = 0.0;
  //          if( this->m_InitializationStrategy == PriorLabelImage ||
  //            this->m_InitializationStrategy == PriorProbabilityImages )
  //            {
  //            if( priorProbabilityImage )
  //              {
  //              priorProbability = priorProbabilityImage->GetPixel( ItW.GetIndex() );
  //              }
  //            if( priorProbability <= this->m_ProbabilityThreshold &&
  //              distancePriorProbabilityImage )
  //              {
  //              priorProbability =
  //                distancePriorProbabilityImage->GetPixel( ItW.GetIndex() );
  //              }
  //               if( this->GetPriorLabelImage() )
  //              {
  //              if( priorProbability == 0.0 )
  //                {
  //                priorProbability = 1.0 / static_cast<RealType>(
  //                  this->m_NumberOfTissueClasses );
  //                }
  //                 else
  //                   {
  //                   priorProbability = 1.0;
  //                   }
  //              }
  //            }
  //          else
  //            {
  //            priorProbability = 1.0;
  //            }
  //          ItW.Set( ItW.Get() + this->m_MixtureModelProportions[n] *
  //            priorProbability );
  //          }
  //        }
  //      }
  //    }
  //  this->SetNthOutput( 0, maxLabels );
  //
  //  if( this->m_UseMixtureModelProportions )
  //    {
  //    //
  //    // Update the class proportions
  //    //
  //    for( unsigned int n = 0; n < this->m_NumberOfTissueClasses; n++ )
  //      {
  //      RealType denominator = 0.0;
  //
  //      RealImagePointer distancePriorProbabilityImage =
  //        this->GetDistancePriorProbabilityImage( n + 1 );
  //      RealImagePointer priorProbabilityImage =
  //        this->GetPriorProbabilityImage( n + 1 );
  //
  //      ImageRegionIteratorWithIndex<RealImageType> ItW(
  //        weightedPriorProbabilityImage,
  //        weightedPriorProbabilityImage->GetRequestedRegion() );
  //      for( ItW.GoToBegin(); !ItW.IsAtEnd(); ++ItW )
  //        {
  //        if( !this->GetMaskImage() || this->GetMaskImage()->GetPixel( ItW.GetIndex() ) !=
  //        NumericTraits<MaskLabelType>::ZeroValue() )
  //          {
  //          RealType priorProbability = 0.0;
  //          if( this->m_InitializationStrategy == PriorLabelImage ||
  //            this->m_InitializationStrategy == PriorProbabilityImages )
  //            {
  //            if( priorProbabilityImage )
  //              {
  //              priorProbability = priorProbabilityImage->GetPixel( ItW.GetIndex() );
  //              }
  //            if( priorProbability <= this->m_ProbabilityThreshold &&
  //              distancePriorProbabilityImage )
  //              {
  //              priorProbability =
  //                distancePriorProbabilityImage->GetPixel( ItW.GetIndex() );
  //              }
  //            if( this->GetPriorLabelImage() )
  //              {
  //              if( priorProbability == 0 )
  //                {
  //                priorProbability= 1.0 /
  //                  static_cast<RealType>( this->m_NumberOfTissueClasses );
  //                }
  //              else
  //                {
  //                   priorProbability = 1.0;
  //                   }
  //              }
  //            }
  //          else
  //            {
  //            priorProbability = 1.0;
  //            }
  //          if( ItW.Get() > 0.0 )
  //            {
  //            denominator += ( priorProbability / ItW.Get() );
  //            }
  //          }
  //        }
  //      if( denominator > 0.0 )
  //        {
  //        this->m_MixtureModelProportions[n] = sumPosteriors[n] / denominator;
  //        }
  //      else
  //        {
  //        this->m_MixtureModelProportions[n] = 0.0;
  //        }
  //      }
  //    }

  //
  // Calculate the maximum posterior probability sum over the region of
  // interest.  This quantity should increase at each iteration.
  //

  if (!this->m_UseAsynchronousUpdating)
  {
    maxPosteriorSum = 0.0;
    ImageRegionConstIteratorWithIndex<RealImageType> ItM(maxPosteriorProbabilityImage,
                                                         maxPosteriorProbabilityImage->GetRequestedRegion());
    for (ItM.GoToBegin(); !ItM.IsAtEnd(); ++ItM)
    {
      if (!this->GetMaskImage() ||
          this->GetMaskImage()->GetPixel(ItM.GetIndex()) != NumericTraits<MaskLabelType>::ZeroValue())
      {
        maxPosteriorSum += static_cast<double>(ItM.Get());
      }
    }
  }

  return static_cast<RealType>(maxPosteriorSum / static_cast<double>(totalSampleSize));
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::RealType
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::PerformLocalLabelingUpdate(
  NeighborhoodIterator<ClassifiedImageType> It)
{
  MeasurementVectorType measurement;

  measurement.SetSize(this->m_NumberOfIntensityImages);
  for (unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++)
  {
    measurement[i] = this->GetIntensityImage(i)->GetPixel(It.GetIndex());
  }

  RealType mrfSmoothingFactor = this->m_MRFSmoothingFactor;
  if (this->m_MRFCoefficientImage)
  {
    mrfSmoothingFactor = this->m_MRFCoefficientImage->GetPixel(It.GetIndex());
  }

  Array<RealType> mrfNeighborhoodWeights;
  this->EvaluateMRFNeighborhoodWeights(It, mrfNeighborhoodWeights);

  LabelType maxLabel = this->m_Randomizer->GetIntegerVariate(this->m_NumberOfTissueClasses - 1) + 1;
  RealType  maxPosteriorProbability = 0.0;
  RealType  sumPosteriorProbability = 0.0;

  unsigned int totalNumberOfClasses = this->m_NumberOfTissueClasses + this->m_NumberOfPartialVolumeClasses;
  for (unsigned int k = 0; k < totalNumberOfClasses; k++)
  {
    // Calculate likelihood probability

    RealType likelihood = this->m_MixtureModelComponents[k]->Evaluate(measurement);

    // Calculate the mrf prior probability

    RealType mrfPriorProbability = 1.0;
    if (mrfSmoothingFactor > NumericTraits<RealType>::ZeroValue() && (It.GetNeighborhood()).Size() > 1)
    {
      RealType numerator = std::exp(-mrfSmoothingFactor * mrfNeighborhoodWeights[k]);
      RealType denominator = 0.0;
      for (unsigned int n = 0; n < this->m_NumberOfTissueClasses; n++)
      {
        denominator += std::exp(-mrfSmoothingFactor * mrfNeighborhoodWeights[n]);
      }
      if (denominator > NumericTraits<RealType>::ZeroValue())
      {
        mrfPriorProbability = numerator / denominator;
      }
    }

    // Get the spatial prior probability

    RealType priorProbability = 1.0;
    if (this->GetPriorProbabilityImage(k + 1))
    {
      priorProbability = this->GetPriorProbabilityImage(k + 1)->GetPixel(It.GetIndex());
    }

    //
    // Calculate the local posterior probability.  Given that the
    // algorithm is meant to maximize the posterior probability of the
    // labeling configuration, this is the critical energy minimization
    // equation.
    //

    RealType posteriorProbability = this->CalculateLocalPosteriorProbability(
      this->m_MixtureModelProportions[k], priorProbability, 1, mrfPriorProbability, likelihood, It.GetIndex(), k + 1);

    if (std::isnan(posteriorProbability) || std::isinf(posteriorProbability))
    {
      posteriorProbability = 0.0;
    }
    sumPosteriorProbability += posteriorProbability;

    if (posteriorProbability > maxPosteriorProbability)
    {
      maxPosteriorProbability = posteriorProbability;
      maxLabel = static_cast<LabelType>(k + 1);
    }
  }
  It.SetCenterPixel(maxLabel);

  maxPosteriorProbability /= sumPosteriorProbability;

  if (std::isnan(maxPosteriorProbability) || std::isinf(maxPosteriorProbability))
  {
    maxPosteriorProbability = 0.0;
  }

  return maxPosteriorProbability;
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::SamplePointer
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::GetScalarSamples()
{
  //
  // This function returns a set of samples for each class such that each
  // measurement vector of the returned SampleType corresponds to a single
  // voxel across the set of auxiliary and input images.
  //

  std::vector<typename SampleType::Pointer> samples;

  //
  // Accumulate the samples in individual SampleTypes.  This allows us to
  // "filter" the samples of each auxiliary/input image.  This filtering
  // could including outlier winsorization, log transformation, and/or
  // converting vector/tensor auxiliary images to scalar data for
  // modeling.
  //
  for (unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++)
  {
    typename SampleType::Pointer sample = SampleType::New();
    sample->SetMeasurementVectorSize(1);
    samples.push_back(sample);
  }

  ImageRegionIteratorWithIndex<ClassifiedImageType> ItO(this->GetOutput(), this->GetOutput()->GetRequestedRegion());
  for (ItO.GoToBegin(); !ItO.IsAtEnd(); ++ItO)
  {
    if (!this->GetMaskImage() ||
        this->GetMaskImage()->GetPixel(ItO.GetIndex()) != NumericTraits<MaskLabelType>::ZeroValue())
    {
      for (unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++)
      {
        typename SampleType::MeasurementVectorType measurement;
        NumericTraits<MeasurementVectorType>::SetLength(measurement, 1);
        measurement.SetSize(1);
        measurement[0] = this->GetIntensityImage(i)->GetPixel(ItO.GetIndex());
        samples[i]->PushBack(measurement);
      }
    }
  }

  //
  // Simultaneously filter the samples and accumulate for return.
  //
  typename SampleType::Pointer scalarSamples = SampleType::New();
  scalarSamples->SetMeasurementVectorSize(this->m_NumberOfIntensityImages);
  for (unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++)
  {
    typename SampleType::Pointer univariateSamples = SampleType::New();
    univariateSamples->SetMeasurementVectorSize(1);
    if (this->m_OutlierHandlingFilter)
    {
      this->m_OutlierHandlingFilter->SetInputListSample(samples[i]);
      this->m_OutlierHandlingFilter->Update();
      univariateSamples = this->m_OutlierHandlingFilter->GetOutput();
    }
    else
    {
      univariateSamples = samples[i];
    }

    if (i == 0)
    {
      typename SampleType::ConstIterator It = univariateSamples->Begin();
      while (It != univariateSamples->End())
      {
        typename SampleType::MeasurementVectorType measurement;
        measurement.SetSize(this->m_NumberOfIntensityImages);
        NumericTraits<MeasurementVectorType>::SetLength(measurement, this->m_NumberOfIntensityImages);
        measurement[0] = It.GetMeasurementVector()[0];
        scalarSamples->PushBack(measurement);
        ++It;
      }
    }
    else
    {
      typename SampleType::Iterator      ItS = scalarSamples->Begin();
      typename SampleType::ConstIterator It = univariateSamples->Begin();
      while (ItS != scalarSamples->End())
      {
        scalarSamples->SetMeasurement(ItS.GetInstanceIdentifier(), i, It.GetMeasurementVector()[0]);
        ++It;
        ++ItS;
      }
    }
  }

  return scalarSamples;
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::RealType
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::CalculateLocalPosteriorProbability(
  RealType     itkNotUsed(mixtureModelProportion),
  RealType     spatialPriorProbability,
  RealType     distancePriorProbability,
  RealType     mrfPriorProbability,
  RealType     likelihood,
  IndexType    index,
  unsigned int whichClass)
{
  // Ensure that the probabilities are at least and epsilon > 0.0 due to numerical issues

  RealType probabilityEpsilon = 1.0e-10;

  spatialPriorProbability = std::max(static_cast<RealType>(probabilityEpsilon), spatialPriorProbability);
  distancePriorProbability = std::max(static_cast<RealType>(probabilityEpsilon), distancePriorProbability);
  mrfPriorProbability = std::max(static_cast<RealType>(probabilityEpsilon), mrfPriorProbability);
  likelihood = std::max(static_cast<RealType>(probabilityEpsilon), likelihood);

  RealType posteriorProbability = 0.0;

  switch (this->m_PosteriorProbabilityFormulation)
  {
    case Socrates:
    default:
    {
      /** R example
        library(ANTsR)
        post<-antsImageRead('testBrainSegmentationPosteriors2.nii.gz',2)
        prior<-antsImageRead('testBrainSegmentationPriorWarped2.nii.gz',2)
        mask<-antsImageRead('testBrainExtractionMask.nii.gz',2)
        pwt<-rev( c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)  )
        gz0<-mask > 0
        socrates<-function( likelihood, prior , priorwt  ) {likelihood^(1-priorwt) * prior^priorwt }
        png('atropossocrates.png',width=1280,height=1280)
        par( mfrow=c(3,4) )
        plotANTsImage( prior )
        for ( i in c(1:length(pwt)) )
        {
        tit<-paste( "Prior-Weight = ", as.character( pwt[i] ) )
        temp<-antsImageClone( prior )
        temp[ gz0 ]<-socrates( post[gz0], prior[gz0], pwt[i] )
        plotANTsImage( temp )
        }
        dev.off()
      */
      posteriorProbability =
        std::pow(static_cast<RealType>(spatialPriorProbability),
                 static_cast<RealType>(this->m_PriorProbabilityWeight)) *
        std::pow(static_cast<RealType>(likelihood * mrfPriorProbability),
                 static_cast<RealType>(NumericTraits<RealType>::OneValue() - this->m_PriorProbabilityWeight));
    }
    break;
    case Plato:
    {
      if (this->m_InitializationStrategy == PriorLabelImage && this->GetPriorLabelImage() &&
          this->GetPriorLabelImage()->GetPixel(index) == whichClass)
      {
        spatialPriorProbability = 1.0;
        likelihood = 1.0;
        mrfPriorProbability = 1.0;
      }
      else if (this->GetPriorProbabilityImage(whichClass) &&
               Math::FloatAlmostEqual(this->GetPriorProbabilityImage(whichClass)->GetPixel(index),
                                      NumericTraits<RealType>::OneValue()))
      {
        spatialPriorProbability = 1.0;
        likelihood = 1.0;
        mrfPriorProbability = 1.0;
      }
      posteriorProbability =
        std::pow(static_cast<RealType>(spatialPriorProbability),
                 static_cast<RealType>(this->m_PriorProbabilityWeight)) *
        std::pow(static_cast<RealType>(likelihood * mrfPriorProbability),
                 static_cast<RealType>(NumericTraits<RealType>::OneValue() - this->m_PriorProbabilityWeight));
    }
    break;
    case Aristotle:
    {
      posteriorProbability =
        std::pow(static_cast<RealType>(spatialPriorProbability * distancePriorProbability),
                 static_cast<RealType>(this->m_PriorProbabilityWeight)) *
        std::pow(static_cast<RealType>(likelihood * mrfPriorProbability),
                 static_cast<RealType>(NumericTraits<RealType>::OneValue() - this->m_PriorProbabilityWeight));
    }
    break;
    case Sigmoid:
    {
      /** some R code to show the effect
        pwt<-c(0,2^c(0:7))
        sig<-function( x , priorwt, offset=0.5 ) { 1.0 / ( 1 + exp(-(x-offset)*priorwt )  ) }
        probs<-c(0:100)/100
        par( mfrow=c(3,3) )
        for ( i in c(1:length(pwt)) )
        {
        tit<-paste("Prior-Weight = ",as.character(pwt[i]))
        plot(main=tit,sig(probs,pwt[i]),type='l',ylim=c(0,1)); points(probs,type='l',col='red')
        }
      */
      RealType totalNumberOfClasses =
        static_cast<RealType>(this->m_NumberOfTissueClasses + this->m_NumberOfPartialVolumeClasses);
      RealType offset = NumericTraits<RealType>::OneValue() / totalNumberOfClasses;
      spatialPriorProbability = NumericTraits<RealType>::OneValue() /
                                (NumericTraits<RealType>::OneValue() +
                                 std::exp(-NumericTraits<RealType>::OneValue() * (spatialPriorProbability - offset) *
                                          this->m_PriorProbabilityWeight));
      posteriorProbability =
        static_cast<RealType>(spatialPriorProbability) * static_cast<RealType>(likelihood * mrfPriorProbability);
    }
    break;
  }

  if (!Math::FloatAlmostEqual(this->m_InitialAnnealingTemperature, NumericTraits<RealType>::OneValue()))
  {
    RealType annealingTemperature = this->m_InitialAnnealingTemperature *
                                    std::pow(this->m_AnnealingRate, static_cast<RealType>(this->m_ElapsedIterations));

    annealingTemperature = std::max(annealingTemperature, this->m_MinimumAnnealingTemperature);

    posteriorProbability = std::pow(static_cast<RealType>(posteriorProbability),
                                    static_cast<RealType>(NumericTraits<RealType>::OneValue() / annealingTemperature));
  }
  if (std::isnan(posteriorProbability) || std::isinf(posteriorProbability))
  {
    posteriorProbability = 0.0;
  }

  return posteriorProbability;
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::EvaluateMRFNeighborhoodWeights(
  ConstNeighborhoodIterator<TClassifiedImage> It,
  Array<RealType> &                           mrfNeighborhoodWeights)
{
  unsigned int totalNumberOfClasses = this->m_NumberOfTissueClasses + this->m_NumberOfPartialVolumeClasses;

  RealType mrfSmoothingFactor = this->m_MRFSmoothingFactor;

  if (this->m_MRFCoefficientImage)
  {
    mrfSmoothingFactor = this->m_MRFCoefficientImage->GetPixel(It.GetIndex());
  }

  mrfNeighborhoodWeights.SetSize(totalNumberOfClasses);
  mrfNeighborhoodWeights.Fill(0.0);

  unsigned int neighborhoodSize = (It.GetNeighborhood()).Size();

  if (mrfSmoothingFactor > NumericTraits<RealType>::ZeroValue() && neighborhoodSize > 1)
  {
    for (unsigned int label = 1; label <= totalNumberOfClasses; label++)
    {
      for (unsigned int n = 0; n < neighborhoodSize; n++)
      {
        if (n == static_cast<unsigned int>(0.5 * neighborhoodSize))
        {
          continue;
        }
        bool      isInBounds = false;
        LabelType neighborLabel = It.GetPixel(n, isInBounds);
        if (!isInBounds || neighborLabel == 0)
        {
          continue;
        }
        typename ClassifiedImageType::OffsetType offset = It.GetOffset(n);

        RealType distance = 0.0;
        for (unsigned int d = 0; d < ImageDimension; d++)
        {
          distance += static_cast<RealType>(Math::sqr(offset[d] * this->m_ImageSpacing[d]));
        }
        distance = std::sqrt(distance);

        RealType delta = 0.0;
        if (label == neighborLabel)
        {
          if (this->m_NumberOfPartialVolumeClasses > 0)
          {
            delta = -2.0;
          }
          else
          {
            delta = -1.0;
          }
        }
        else
        {
          bool                                              isCommonTissue = false;
          typename PartialVolumeClassesType::const_iterator it;
          for (it = this->m_PartialVolumeClasses.begin(); it != this->m_PartialVolumeClasses.end(); ++it)
          {
            if (std::find(it->begin(), it->end(), static_cast<LabelType>(label)) != it->end() &&
                std::find(it->begin(), it->end(), static_cast<LabelType>(neighborLabel)) != it->end())
            {
              isCommonTissue = true;
              break;
            }
          }
          if (isCommonTissue)
          {
            delta = -1.0;
          }
          else if (this->m_NumberOfPartialVolumeClasses > 0)
          {
            delta = 1.0;
          }
          else
          {
            delta = 0.0;
          }
        }
        mrfNeighborhoodWeights[label - 1] += (delta / distance);
      }
    }
  }
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::RealImagePointer
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::GetPosteriorProbabilityImage(
  unsigned int whichClass)
{
  unsigned int totalNumberOfClasses = this->m_NumberOfTissueClasses + this->m_NumberOfPartialVolumeClasses;

  if (whichClass > totalNumberOfClasses)
  {
    itkExceptionMacro("Requested class is greater than the number of classes.");
  }

  //
  // If memory minimization is turned off and if the posterior probability
  // images have already been calculated, simply return the probability
  // image for the requested class.  Otherwise, calculate the probability
  // image.
  //
  if (whichClass <= this->m_PosteriorProbabilityImages.size())
  {
    return this->m_PosteriorProbabilityImages[whichClass - 1];
  }
  else
  {
    //
    // Here we assume that the calling function is invoked in order such
    // that GetPosteriorProbabilityImage( 1 ) is called before
    // GetPosteriorProbabilityImage( 2 ), etc.  As such, when this part of
    // the code is reached and the class requested is '1', we assume that
    // the sum of the posterior probability images needs to be calculated
    // for normalization purposes.  This sum is then saved for subsequent calls.
    //
    RealImagePointer posteriorProbabilityImage = AllocImage<RealImageType>(this->GetOutput(), 0);

    //
    // Calculate the sum of the probability images.  Also, store the
    // posterior probability images if m_MinimizeMemoryUsage == false.
    //
    if (whichClass == 1)
    {
      this->m_SumPosteriorProbabilityImage = AllocImage<RealImageType>(this->GetOutput(), 0);

      RealImagePointer sumPriorProbabilityImage = nullptr;

      if (this->m_InitializationStrategy == PriorLabelImage || this->m_InitializationStrategy == PriorProbabilityImages)
      {
        sumPriorProbabilityImage = AllocImage<RealImageType>(this->GetOutput(), 0);
        for (unsigned int c = 0; c < totalNumberOfClasses; c++)
        {
          RealImagePointer priorProbabilityImage = this->GetPriorProbabilityImage(c + 1);
          if (!priorProbabilityImage)
          {
            continue;
          }

          ImageRegionIteratorWithIndex<RealImageType> ItS(sumPriorProbabilityImage,
                                                          sumPriorProbabilityImage->GetLargestPossibleRegion());
          for (ItS.GoToBegin(); !ItS.IsAtEnd(); ++ItS)
          {
            RealType priorProbability = 0.0;
            if (priorProbabilityImage)
            {
              priorProbability = priorProbabilityImage->GetPixel(ItS.GetIndex());
            }
            else if (this->GetPriorLabelImage())
            {
              if (Math::FloatAlmostEqual(priorProbability, NumericTraits<RealType>::ZeroValue()))
              {
                priorProbability = NumericTraits<RealType>::OneValue() / static_cast<RealType>(totalNumberOfClasses);
              }
              else
              {
                priorProbability = NumericTraits<RealType>::OneValue();
              }
            }
            ItS.Set(ItS.Get() + this->m_MixtureModelProportions[c] * priorProbability);
          }
        }
      }
      for (unsigned int c = 0; c < totalNumberOfClasses; c++)
      {
        std::vector<RealImagePointer> smoothImages;

        if (this->m_InitializationStrategy == PriorProbabilityImages ||
            (this->m_InitializationStrategy == PriorLabelImage && c < this->m_NumberOfTissueClasses))
        {
          for (unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++)
          {
            if (this->m_AdaptiveSmoothingWeights.size() > i &&
                this->m_AdaptiveSmoothingWeights[i] > NumericTraits<RealType>::ZeroValue())
            {
              smoothImages.push_back(this->GetSmoothIntensityImageFromPriorImage(i, c + 1));
            }
            else
            {
              smoothImages.push_back(nullptr);
            }
          }
        }

        RealImagePointer distancePriorProbabilityImage = this->GetDistancePriorProbabilityImage(c + 1);
        RealImagePointer priorProbabilityImage = this->GetPriorProbabilityImage(c + 1);

        typename NeighborhoodIterator<ClassifiedImageType>::RadiusType radius;
        unsigned int                                                   neighborhoodSize = 1;
        for (unsigned int d = 0; d < ImageDimension; d++)
        {
          neighborhoodSize *= (2 * this->m_MRFRadius[d] + 1);
          radius[d] = this->m_MRFRadius[d];
        }

    auto multiThreader = this->GetMultiThreader();
    using RegionType = typename RealImageType::RegionType;
    multiThreader->template ParallelizeImageRegion<ImageDimension>(
      this->GetOutput()->GetRequestedRegion(),
      [this, &radius, &c, &totalNumberOfClasses, &distancePriorProbabilityImage, &priorProbabilityImage, &sumPriorProbabilityImage, &smoothImages, &posteriorProbabilityImage](const RegionType & outputRegionForThread) {

        ConstNeighborhoodIterator<ClassifiedImageType> ItO(
          radius, this->GetOutput(), outputRegionForThread);
        ImageRegionIterator<RealImageType> ItS(this->m_SumPosteriorProbabilityImage,
                                               outputRegionForThread);
        for (ItO.GoToBegin(), ItS.GoToBegin(); !ItO.IsAtEnd(); ++ItO, ++ItS)
        {
          if (!this->GetMaskImage() ||
              this->GetMaskImage()->GetPixel(ItO.GetIndex()) != NumericTraits<MaskLabelType>::ZeroValue())
          {
            RealType mrfSmoothingFactor = this->m_MRFSmoothingFactor;
            if (this->m_MRFCoefficientImage)
            {
              mrfSmoothingFactor = this->m_MRFCoefficientImage->GetPixel(ItO.GetIndex());
            }
            //
            // Perform mrf prior calculation
            //
            RealType mrfPriorProbability = NumericTraits<RealType>::OneValue();
            if (mrfSmoothingFactor > NumericTraits<RealType>::ZeroValue() && (ItO.GetNeighborhood()).Size() > 1)
            {
              Array<RealType> mrfNeighborhoodWeights;
              this->EvaluateMRFNeighborhoodWeights(ItO, mrfNeighborhoodWeights);

              RealType numerator = std::exp(-mrfSmoothingFactor * mrfNeighborhoodWeights[c]);
              RealType denominator = NumericTraits<RealType>::ZeroValue();
              for (unsigned int n = 0; n < totalNumberOfClasses; n++)
              {
                denominator += std::exp(-mrfSmoothingFactor * mrfNeighborhoodWeights[n]);
              }
              if (denominator > NumericTraits<RealType>::ZeroValue())
              {
                mrfPriorProbability = numerator / denominator;
              }
            }

            //
            // Perform prior calculation using both the mixing proportions
            // and template-based prior images (if available)
            //
            RealType priorProbability = 0.0;
            RealType distancePriorProbability = 0.0;
            if (this->m_InitializationStrategy == PriorLabelImage ||
                this->m_InitializationStrategy == PriorProbabilityImages)
            {
              if (distancePriorProbabilityImage)
              {
                distancePriorProbability = distancePriorProbabilityImage->GetPixel(ItO.GetIndex());
              }
              if (priorProbabilityImage)
              {
                priorProbability = priorProbabilityImage->GetPixel(ItO.GetIndex());
              }
              else if (this->GetPriorLabelImage())
              {
                if (Math::FloatAlmostEqual(priorProbability, NumericTraits<RealType>::ZeroValue()))
                {
                  priorProbability = NumericTraits<RealType>::OneValue() / static_cast<RealType>(totalNumberOfClasses);
                }
                else
                {
                  priorProbability = NumericTraits<RealType>::OneValue();
                }
              }
              RealType sumPriorProbability = sumPriorProbabilityImage->GetPixel(ItO.GetIndex());
              if (sumPriorProbability > this->m_ProbabilityThreshold)
              {
                priorProbability *= (this->m_MixtureModelProportions[c] / sumPriorProbability);
              }
              else if (distancePriorProbabilityImage)
              {
                priorProbability = distancePriorProbability;
              }
            }

            MeasurementVectorType measurement;
            measurement.SetSize(this->m_NumberOfIntensityImages);
            for (unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++)
            {
              measurement[i] = this->GetIntensityImage(i)->GetPixel(ItO.GetIndex());

              if ((this->m_InitializationStrategy == PriorProbabilityImages ||
                   this->m_InitializationStrategy == PriorLabelImage) &&
                  smoothImages[i])
              {
                measurement[i] =
                  (NumericTraits<RealType>::OneValue() - this->m_AdaptiveSmoothingWeights[i]) * measurement[i] +
                  this->m_AdaptiveSmoothingWeights[i] * smoothImages[i]->GetPixel(ItO.GetIndex());
              }
            }

            //
            // Calculate likelihood probability from the model
            //
            RealType likelihood = this->m_MixtureModelComponents[c]->Evaluate(measurement);

            //
            // Calculate the local posterior probability.  Given that the
            // algorithm is meant to maximize the posterior probability of the
            // labeling configuration, this is the critical energy minimization
            // equation.
            //
            RealType posteriorProbability = this->CalculateLocalPosteriorProbability(this->m_MixtureModelProportions[c],
                                                                                     priorProbability,
                                                                                     distancePriorProbability,
                                                                                     mrfPriorProbability,
                                                                                     likelihood,
                                                                                     ItO.GetIndex(),
                                                                                     c + 1);

            if (std::isnan(posteriorProbability) || std::isinf(posteriorProbability))
            {
              posteriorProbability = 0.0;
            }

            if ((c == 0) || !this->m_MinimizeMemoryUsage)
            {
              posteriorProbabilityImage->SetPixel(ItO.GetIndex(), posteriorProbability);
            }

            //
            // Calculate a running total of the posterior probabilities.
            //
            ItS.Set(ItS.Get() + posteriorProbability);
          }
        }
      }, nullptr); // end ParallelizeImageRegion
        if (!this->m_MinimizeMemoryUsage)
        {
          typedef ImageDuplicator<RealImageType> DuplicatorType;
          typename DuplicatorType::Pointer       duplicator = DuplicatorType::New();
          duplicator->SetInputImage(posteriorProbabilityImage);
          duplicator->Update();

          this->m_PosteriorProbabilityImages.push_back(duplicator->GetOutput());
        }
      }

      //
      // Normalize the posterior probability image(s).
      //
      ImageRegionIterator<RealImageType> ItS(this->m_SumPosteriorProbabilityImage,
                                             this->m_SumPosteriorProbabilityImage->GetRequestedRegion());
      if (this->m_MinimizeMemoryUsage)
      {
        ImageRegionIterator<RealImageType> ItP(posteriorProbabilityImage,
                                               posteriorProbabilityImage->GetRequestedRegion());
        for (ItP.GoToBegin(), ItS.GoToBegin(); !ItS.IsAtEnd(); ++ItP, ++ItS)
        {
          if (ItS.Get() > 0)
          {
            ItP.Set(ItP.Get() / ItS.Get());
          }
        }
        return posteriorProbabilityImage;
      }
      else
      {
        for (unsigned int n = 0; n < totalNumberOfClasses; n++)
        {
          ImageRegionIterator<RealImageType> ItP(this->m_PosteriorProbabilityImages[n],
                                                 this->m_PosteriorProbabilityImages[n]->GetRequestedRegion());
          for (ItP.GoToBegin(), ItS.GoToBegin(); !ItS.IsAtEnd(); ++ItP, ++ItS)
          {
            if (ItS.Get() > 0)
            {
              ItP.Set(ItP.Get() / ItS.Get());
            }
          }
        }
        return this->m_PosteriorProbabilityImages[0];
      }
    }
    else // whichClass > 1
    {
      RealImagePointer sumPriorProbabilityImage = nullptr;

      if (this->m_InitializationStrategy == PriorLabelImage || this->m_InitializationStrategy == PriorProbabilityImages)
      {
        sumPriorProbabilityImage = RealImageType::New();
        sumPriorProbabilityImage->CopyInformation(this->GetOutput());
        sumPriorProbabilityImage->SetRegions(this->GetOutput()->GetRequestedRegion());
        sumPriorProbabilityImage->AllocateInitialized();
        for (unsigned int c = 0; c < totalNumberOfClasses; c++)
        {
          RealImagePointer priorProbabilityImage = this->GetPriorProbabilityImage(c + 1);
          if (!priorProbabilityImage)
          {
            continue;
          }

          ImageRegionIteratorWithIndex<RealImageType> ItS(sumPriorProbabilityImage,
                                                          sumPriorProbabilityImage->GetLargestPossibleRegion());
          for (ItS.GoToBegin(); !ItS.IsAtEnd(); ++ItS)
          {
            RealType priorProbability = 0.0;
            if (priorProbabilityImage)
            {
              priorProbability = priorProbabilityImage->GetPixel(ItS.GetIndex());
            }
            else if (this->GetPriorLabelImage())
            {
              if (Math::FloatAlmostEqual(priorProbability, NumericTraits<RealType>::ZeroValue()))
              {
                priorProbability = NumericTraits<RealType>::OneValue() / static_cast<RealType>(totalNumberOfClasses);
              }
              else
              {
                priorProbability = NumericTraits<RealType>::OneValue();
              }
            }
            ItS.Set(ItS.Get() + this->m_MixtureModelProportions[c] * priorProbability);
          }
        }
      }

      std::vector<RealImagePointer> smoothImages;

      if (this->m_InitializationStrategy == PriorProbabilityImages || this->m_InitializationStrategy == PriorLabelImage)
      {
        for (unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++)
        {
          if (this->m_AdaptiveSmoothingWeights.size() > i &&
              this->m_AdaptiveSmoothingWeights[i] > NumericTraits<RealType>::ZeroValue())
          {
            smoothImages.push_back(this->GetSmoothIntensityImageFromPriorImage(i, whichClass));
          }
          else
          {
            smoothImages.push_back(nullptr);
          }
        }
      }

      RealImagePointer distancePriorProbabilityImage = this->GetDistancePriorProbabilityImage(whichClass);
      RealImagePointer priorProbabilityImage = this->GetPriorProbabilityImage(whichClass);

      typename NeighborhoodIterator<ClassifiedImageType>::RadiusType radius;
      for (unsigned int d = 0; d < ImageDimension; d++)
      {
        radius[d] = this->m_MRFRadius[d];
      }

      ConstNeighborhoodIterator<ClassifiedImageType> ItO(
        radius, this->GetOutput(), this->GetOutput()->GetRequestedRegion());
      for (ItO.GoToBegin(); !ItO.IsAtEnd(); ++ItO)
      {
        if (!this->GetMaskImage() ||
            this->GetMaskImage()->GetPixel(ItO.GetIndex()) != NumericTraits<MaskLabelType>::ZeroValue())
        {
          RealType mrfSmoothingFactor = this->m_MRFSmoothingFactor;
          if (this->m_MRFCoefficientImage)
          {
            mrfSmoothingFactor = this->m_MRFCoefficientImage->GetPixel(ItO.GetIndex());
          }
          //
          // Perform mrf prior calculation
          //
          RealType mrfPriorProbability = NumericTraits<RealType>::OneValue();
          if (mrfSmoothingFactor > NumericTraits<RealType>::ZeroValue() && (ItO.GetNeighborhood()).Size() > 1)
          {
            Array<RealType> mrfNeighborhoodWeights;
            this->EvaluateMRFNeighborhoodWeights(ItO, mrfNeighborhoodWeights);

            RealType numerator = std::exp(-mrfSmoothingFactor * mrfNeighborhoodWeights[whichClass - 1]);
            RealType denominator = NumericTraits<RealType>::ZeroValue();
            for (unsigned int n = 0; n < totalNumberOfClasses; n++)
            {
              denominator += std::exp(-mrfSmoothingFactor * mrfNeighborhoodWeights[n]);
            }
            if (denominator > NumericTraits<RealType>::ZeroValue())
            {
              mrfPriorProbability = numerator / denominator;
            }
          }

          //
          // Perform prior calculation using both the mixing proportions
          // and template-based prior images (if available)
          //
          RealType priorProbability = NumericTraits<RealType>::ZeroValue();
          RealType distancePriorProbability = NumericTraits<RealType>::ZeroValue();
          if (this->m_InitializationStrategy == PriorLabelImage ||
              this->m_InitializationStrategy == PriorProbabilityImages)
          {
            if (distancePriorProbabilityImage)
            {
              distancePriorProbability = distancePriorProbabilityImage->GetPixel(ItO.GetIndex());
            }
            if (priorProbabilityImage)
            {
              priorProbability = priorProbabilityImage->GetPixel(ItO.GetIndex());
            }
            RealType sumPriorProbability = sumPriorProbabilityImage->GetPixel(ItO.GetIndex());
            if (sumPriorProbability > this->m_ProbabilityThreshold)
            {
              priorProbability *= (this->m_MixtureModelProportions[whichClass - 1] / sumPriorProbability);
            }
            else if (distancePriorProbabilityImage)
            {
              priorProbability = distancePriorProbability;
            }
          }

          MeasurementVectorType measurement;
          measurement.SetSize(this->m_NumberOfIntensityImages);
          for (unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++)
          {
            measurement[i] = this->GetIntensityImage(i)->GetPixel(ItO.GetIndex());

            if ((this->m_InitializationStrategy == PriorProbabilityImages ||
                 this->m_InitializationStrategy == PriorLabelImage) &&
                smoothImages[i])
            {
              measurement[i] =
                (NumericTraits<RealType>::OneValue() - this->m_AdaptiveSmoothingWeights[i]) * measurement[i] +
                this->m_AdaptiveSmoothingWeights[i] * smoothImages[i]->GetPixel(ItO.GetIndex());
            }
          }

          //
          // Calculate likelihood probability from the model
          //
          RealType likelihood = this->m_MixtureModelComponents[whichClass - 1]->Evaluate(measurement);

          //
          // Calculate the local posterior probability.  Given that the
          // algorithm is meant to maximize the posterior probability of the
          // labeling configuration, this is the critical energy minimization
          // equation.
          //
          RealType posteriorProbability =
            this->CalculateLocalPosteriorProbability(this->m_MixtureModelProportions[whichClass - 1],
                                                     priorProbability,
                                                     distancePriorProbability,
                                                     mrfPriorProbability,
                                                     likelihood,
                                                     ItO.GetIndex(),
                                                     whichClass);

          if (std::isnan(posteriorProbability) || std::isinf(posteriorProbability))
          {
            posteriorProbability = NumericTraits<RealType>::ZeroValue();
          }

          posteriorProbabilityImage->SetPixel(ItO.GetIndex(), posteriorProbability);
        }
      }

      //
      // Normalize the posterior probability image(s).
      //
      ImageRegionIterator<RealImageType> ItS(this->m_SumPosteriorProbabilityImage,
                                             this->m_SumPosteriorProbabilityImage->GetRequestedRegion());
      ImageRegionIterator<RealImageType> ItP(posteriorProbabilityImage,
                                             posteriorProbabilityImage->GetRequestedRegion());
      for (ItP.GoToBegin(), ItS.GoToBegin(); !ItS.IsAtEnd(); ++ItP, ++ItS)
      {
        if (ItS.Get() > 0)
        {
          ItP.Set(ItP.Get() / ItS.Get());
        }
      }
      return posteriorProbabilityImage;
    }
  }
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::RealImagePointer
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::GetDistancePriorProbabilityImage(
  unsigned int whichClass)
{
  unsigned int totalNumberOfClasses = this->m_NumberOfTissueClasses + this->m_NumberOfPartialVolumeClasses;

  if ((this->m_InitializationStrategy != PriorLabelImage && this->m_InitializationStrategy != PriorProbabilityImages) ||
      (whichClass > this->m_NumberOfTissueClasses && whichClass <= totalNumberOfClasses))
  {
    return nullptr;
  }
  if (this->m_NumberOfPartialVolumeClasses == 0 && whichClass > this->m_NumberOfTissueClasses)
  {
    itkExceptionMacro("The requested distance prior probability image = "
                      << whichClass << " should be in the range [1, " << this->m_NumberOfTissueClasses << "]");
  }
  else if (whichClass > this->m_NumberOfTissueClasses)
  {
    return nullptr;
  }

  //
  // If memory minimization is turned off and if the distance prior probability
  // images have already been calculated, simply return the probability
  // image for the requested class.  Otherwise, calculate the probability
  // image.
  //
  if (whichClass <= this->m_DistancePriorProbabilityImages.size())
  {
    return this->m_DistancePriorProbabilityImages[whichClass - 1];
  }
  else
  {
    //
    // Here we assume that the calling function is invoked in order such
    // that GetDistancePriorImage( 1 ) is called before
    // GetDistancePriorImage( 2 ), etc.  As such, when this part of
    // the code is reached and the class requested is '1', we assume that
    // the sum of the distance prior probability images needs to be calculated
    // for normalization purposes.  This sum is then saved for subsequent calls.
    //
    RealImagePointer distancePriorProbabilityImage = nullptr;

    //
    // Calculate the sum of the distance probability images.  Also, store the
    // distance probability images if m_MinimizeMemoryUsage == false.
    //
    if (whichClass == 1)
    {
      this->m_SumDistancePriorProbabilityImage = RealImageType::New();
      this->m_SumDistancePriorProbabilityImage->CopyInformation(this->GetOutput());
      this->m_SumDistancePriorProbabilityImage->SetRegions(this->GetOutput()->GetRequestedRegion());
      this->m_SumDistancePriorProbabilityImage->AllocateInitialized();
      for (unsigned int c = 0; c < this->m_NumberOfTissueClasses; c++)
      {
        typedef BinaryThresholdImageFilter<ClassifiedImageType, RealImageType> ThresholderType;
        typename ThresholderType::Pointer                                      thresholder = ThresholderType::New();
        if (this->m_InitializationStrategy == PriorLabelImage)
        {
          thresholder->SetInput(const_cast<ClassifiedImageType *>(this->GetPriorLabelImage()));
        }
        else
        {
          thresholder->SetInput(this->GetOutput());
        }
        thresholder->SetInsideValue(1);
        thresholder->SetOutsideValue(0);
        thresholder->SetLowerThreshold(static_cast<LabelType>(c + 1));
        thresholder->SetUpperThreshold(static_cast<LabelType>(c + 1));
        thresholder->Update();

        RealImagePointer distanceImage = RealImageType::New();

        if (this->m_UseEuclideanDistanceForPriorLabels)
        {
          typedef SignedMaurerDistanceMapImageFilter<RealImageType, RealImageType> DistancerType;
          typename DistancerType::Pointer                                          distancer = DistancerType::New();
          distancer->SetInput(thresholder->GetOutput());
          distancer->SetSquaredDistance(false);
          distancer->SetUseImageSpacing(true);
          distancer->SetInsideIsPositive(false);
          distancer->Update();

          distanceImage = distancer->GetOutput();
        }
        else
        {
          typedef BinaryContourImageFilter<RealImageType, RealImageType> ContourFilterType;
          typename ContourFilterType::Pointer                            contour = ContourFilterType::New();
          contour->SetInput(thresholder->GetOutput());
          contour->FullyConnectedOff();
          contour->SetBackgroundValue(0);
          contour->SetForegroundValue(1);
          contour->Update();

          typedef FastMarchingImageFilter<RealImageType, RealImageType> FastMarchingFilterType;
          typename FastMarchingFilterType::Pointer                      fastMarching = FastMarchingFilterType::New();

          typedef CastImageFilter<MaskImageType, RealImageType> CasterType;
          typename CasterType::Pointer                          caster = CasterType::New();
          if (this->GetMaskImage())
          {
            caster->SetInput(const_cast<MaskImageType *>(this->GetMaskImage()));
            caster->Update();
            fastMarching->SetInput(caster->GetOutput());
          }
          else
          {
            fastMarching->SetSpeedConstant(1.0);
            fastMarching->SetOverrideOutputInformation(true);
            fastMarching->SetOutputOrigin(this->GetOutput()->GetOrigin());
            fastMarching->SetOutputSpacing(this->GetOutput()->GetSpacing());
            fastMarching->SetOutputRegion(this->GetOutput()->GetRequestedRegion());
            fastMarching->SetOutputDirection(this->GetOutput()->GetDirection());
          }

          typedef typename FastMarchingFilterType::NodeContainer NodeContainer;
          typedef typename FastMarchingFilterType::NodeType      NodeType;
          typename NodeContainer::Pointer                        trialPoints = NodeContainer::New();
          trialPoints->Initialize();

          unsigned long trialCount = 0;

          ImageRegionIteratorWithIndex<RealImageType> ItC(contour->GetOutput(),
                                                          contour->GetOutput()->GetRequestedRegion());
          for (ItC.GoToBegin(); !ItC.IsAtEnd(); ++ItC)
          {
            if (Math::FloatAlmostEqual(ItC.Get(), contour->GetForegroundValue()))
            {
              NodeType node;
              node.SetValue(0.0);
              node.SetIndex(ItC.GetIndex());
              trialPoints->InsertElement(trialCount++, node);
            }
          }
          fastMarching->SetTrialPoints(trialPoints);
          fastMarching->SetStoppingValue(NumericTraits<RealType>::max());
          //           fastMarching->SetTopologyCheck( FastMarchingFilterType::None );
          fastMarching->Update();

          ImageRegionIterator<RealImageType> ItT(thresholder->GetOutput(),
                                                 thresholder->GetOutput()->GetRequestedRegion());
          ImageRegionIterator<RealImageType> ItF(fastMarching->GetOutput(),
                                                 fastMarching->GetOutput()->GetRequestedRegion());
          for (ItT.GoToBegin(), ItF.GoToBegin(); !ItT.IsAtEnd(); ++ItT, ++ItF)
          {
            if (Math::FloatAlmostEqual(ItT.Get(), NumericTraits<float>::OneValue()))
            {
              ItF.Set(-ItF.Get());
            }
          }

          distanceImage = fastMarching->GetOutput();
        }

        RealType maximumInteriorDistance = 0.0;

        ImageRegionIterator<RealImageType> ItD(distanceImage, distanceImage->GetRequestedRegion());
        for (ItD.GoToBegin(); !ItD.IsAtEnd(); ++ItD)
        {
          if (ItD.Get() < 0 && maximumInteriorDistance < Math::abs(ItD.Get()))
          {
            maximumInteriorDistance = Math::abs(ItD.Get());
          }
        }

        RealType labelLambda = 0.0;
        RealType labelBoundaryProbability = 1.0;

        typename LabelParameterMapType::iterator it = this->m_PriorLabelParameterMap.find(c + 1);
        if (it != this->m_PriorLabelParameterMap.end())
        {
          labelLambda = (it->second).first;
          labelBoundaryProbability = (it->second).second;
        }
        for (ItD.GoToBegin(); !ItD.IsAtEnd(); ++ItD)
        {
          if (Math::FloatAlmostEqual(labelLambda, NumericTraits<RealType>::ZeroValue()))
          {
            if (ItD.Get() <= NumericTraits<float>::ZeroValue())
            {
              ItD.Set(labelBoundaryProbability);
            }
            else
            {
              ItD.Set(NumericTraits<RealType>::ZeroValue());
            }
          }
          else if (ItD.Get() >= NumericTraits<float>::ZeroValue())
          {
            ItD.Set(labelBoundaryProbability * std::exp(-labelLambda * ItD.Get()));
          }
          else if (ItD.Get() < NumericTraits<float>::ZeroValue())
          {
            ItD.Set(NumericTraits<RealType>::OneValue() -
                    (NumericTraits<RealType>::OneValue() - labelBoundaryProbability) *
                      (maximumInteriorDistance - Math::abs(ItD.Get())) / (maximumInteriorDistance));
          }
        }

        typedef AddImageFilter<RealImageType, RealImageType, RealImageType> AdderType;
        typename AdderType::Pointer                                         adder = AdderType::New();
        adder->SetInput1(this->m_SumDistancePriorProbabilityImage);
        adder->SetInput2(distanceImage);
        adder->Update();

        this->m_SumDistancePriorProbabilityImage = adder->GetOutput();

        if ((c == 0) && this->m_MinimizeMemoryUsage)
        {
          distancePriorProbabilityImage = distanceImage;
        }
        if (!this->m_MinimizeMemoryUsage)
        {
          this->m_DistancePriorProbabilityImages.push_back(distanceImage);
        }
      }

      //
      // Normalize the distance prior probability image(s).
      //
      ImageRegionIterator<RealImageType> ItS(this->m_SumDistancePriorProbabilityImage,
                                             this->m_SumDistancePriorProbabilityImage->GetRequestedRegion());
      if (this->m_MinimizeMemoryUsage)
      {
        ImageRegionIteratorWithIndex<RealImageType> ItD(distancePriorProbabilityImage,
                                                        distancePriorProbabilityImage->GetRequestedRegion());
        for (ItD.GoToBegin(), ItS.GoToBegin(); !ItS.IsAtEnd(); ++ItD, ++ItS)
        {
          if (!this->GetMaskImage() ||
              this->GetMaskImage()->GetPixel(ItD.GetIndex()) != NumericTraits<MaskLabelType>::ZeroValue())
          {
            if (ItS.Get() <= this->m_ProbabilityThreshold)
            {
              ItD.Set(NumericTraits<RealType>::ZeroValue());
            }
            else
            {
              ItD.Set(ItD.Get() / ItS.Get());
            }
          }
        }
        return distancePriorProbabilityImage;
      }
      else
      {
        for (unsigned int c = 0; c < this->m_NumberOfTissueClasses; c++)
        {
          ImageRegionIteratorWithIndex<RealImageType> ItD(
            this->m_DistancePriorProbabilityImages[c], this->m_DistancePriorProbabilityImages[c]->GetRequestedRegion());
          for (ItD.GoToBegin(), ItS.GoToBegin(); !ItS.IsAtEnd(); ++ItD, ++ItS)
          {
            if (!this->GetMaskImage() ||
                this->GetMaskImage()->GetPixel(ItD.GetIndex()) != NumericTraits<MaskLabelType>::ZeroValue())
            {
              if (ItS.Get() <= this->m_ProbabilityThreshold)
              {
                ItD.Set(NumericTraits<RealType>::ZeroValue());
              }
              else
              {
                ItD.Set(ItD.Get() / ItS.Get());
              }
            }
          }
        }
        return this->m_DistancePriorProbabilityImages[0];
      }
    }
    else // whichClass > 1
    {
      typedef BinaryThresholdImageFilter<ClassifiedImageType, RealImageType> ThresholderType;
      typename ThresholderType::Pointer                                      thresholder = ThresholderType::New();
      if (this->m_InitializationStrategy == PriorLabelImage)
      {
        thresholder->SetInput(const_cast<ClassifiedImageType *>(this->GetPriorLabelImage()));
      }
      else
      {
        thresholder->SetInput(this->GetOutput());
      }
      thresholder->SetInsideValue(1);
      thresholder->SetOutsideValue(0);
      thresholder->SetLowerThreshold(static_cast<LabelType>(whichClass));
      thresholder->SetUpperThreshold(static_cast<LabelType>(whichClass));
      thresholder->Update();

      RealImagePointer distanceImage = RealImageType::New();

      if (this->m_UseEuclideanDistanceForPriorLabels)
      {
        typedef SignedMaurerDistanceMapImageFilter<RealImageType, RealImageType> DistancerType;
        typename DistancerType::Pointer                                          distancer = DistancerType::New();
        distancer->SetInput(thresholder->GetOutput());
        distancer->SetSquaredDistance(false);
        distancer->SetUseImageSpacing(true);
        distancer->SetInsideIsPositive(false);
        distancer->Update();

        distanceImage = distancer->GetOutput();
      }
      else
      {
        typedef BinaryContourImageFilter<RealImageType, RealImageType> ContourFilterType;
        typename ContourFilterType::Pointer                            contour = ContourFilterType::New();
        contour->SetInput(thresholder->GetOutput());
        contour->FullyConnectedOff();
        contour->SetBackgroundValue(0);
        contour->SetForegroundValue(1);
        contour->Update();

        typedef FastMarchingImageFilter<RealImageType, RealImageType> FastMarchingFilterType;
        typename FastMarchingFilterType::Pointer                      fastMarching = FastMarchingFilterType::New();

        typedef CastImageFilter<MaskImageType, RealImageType> CasterType;
        typename CasterType::Pointer                          caster = CasterType::New();
        if (this->GetMaskImage())
        {
          caster->SetInput(const_cast<MaskImageType *>(this->GetMaskImage()));
          caster->Update();
          fastMarching->SetInput(caster->GetOutput());
        }
        else
        {
          fastMarching->SetSpeedConstant(1.0);
          fastMarching->SetOverrideOutputInformation(true);
          fastMarching->SetOutputOrigin(this->GetOutput()->GetOrigin());
          fastMarching->SetOutputSpacing(this->GetOutput()->GetSpacing());
          fastMarching->SetOutputRegion(this->GetOutput()->GetRequestedRegion());
          fastMarching->SetOutputDirection(this->GetOutput()->GetDirection());
        }

        typedef typename FastMarchingFilterType::NodeContainer NodeContainer;
        typedef typename FastMarchingFilterType::NodeType      NodeType;
        typename NodeContainer::Pointer                        trialPoints = NodeContainer::New();
        trialPoints->Initialize();

        unsigned long trialCount = 0;

        ImageRegionIteratorWithIndex<RealImageType> ItC(contour->GetOutput(),
                                                        contour->GetOutput()->GetRequestedRegion());
        for (ItC.GoToBegin(); !ItC.IsAtEnd(); ++ItC)
        {
          if (Math::FloatAlmostEqual(ItC.Get(), contour->GetForegroundValue()))
          {
            NodeType node;
            node.SetValue(0.0);
            node.SetIndex(ItC.GetIndex());
            trialPoints->InsertElement(trialCount++, node);
          }
        }
        fastMarching->SetTrialPoints(trialPoints);
        fastMarching->SetStoppingValue(NumericTraits<RealType>::max());
        //         fastMarching->SetTopologyCheck( FastMarchingFilterType::None );
        fastMarching->Update();

        ImageRegionIterator<RealImageType> ItT(thresholder->GetOutput(),
                                               thresholder->GetOutput()->GetRequestedRegion());
        ImageRegionIterator<RealImageType> ItF(fastMarching->GetOutput(),
                                               fastMarching->GetOutput()->GetRequestedRegion());
        for (ItT.GoToBegin(), ItF.GoToBegin(); !ItT.IsAtEnd(); ++ItT, ++ItF)
        {
          if (Math::FloatAlmostEqual(ItT.Get(), itk::NumericTraits<float>::OneValue()))
          {
            ItF.Set(-ItF.Get());
          }
        }
        distanceImage = fastMarching->GetOutput();
      }

      distancePriorProbabilityImage = distanceImage;

      RealType maximumInteriorDistance = 0.0;

      ImageRegionIterator<RealImageType> ItD(distancePriorProbabilityImage,
                                             distancePriorProbabilityImage->GetRequestedRegion());
      for (ItD.GoToBegin(); !ItD.IsAtEnd(); ++ItD)
      {
        if (ItD.Get() < 0 && maximumInteriorDistance < Math::abs(ItD.Get()))
        {
          maximumInteriorDistance = Math::abs(ItD.Get());
        }
      }

      RealType labelLambda = 0.0;
      RealType labelBoundaryProbability = 1.0;

      typename LabelParameterMapType::iterator it = this->m_PriorLabelParameterMap.find(whichClass);
      if (it != this->m_PriorLabelParameterMap.end())
      {
        labelLambda = (it->second).first;
        labelBoundaryProbability = (it->second).second;
      }
      for (ItD.GoToBegin(); !ItD.IsAtEnd(); ++ItD)
      {
        if (Math::FloatAlmostEqual(labelLambda, NumericTraits<RealType>::ZeroValue()))
        {
          if (ItD.Get() <= NumericTraits<RealType>::ZeroValue())
          {
            ItD.Set(labelBoundaryProbability);
          }
          else
          {
            ItD.Set(NumericTraits<RealType>::ZeroValue());
          }
        }
        else if (ItD.Get() >= NumericTraits<RealType>::ZeroValue())
        {
          ItD.Set(labelBoundaryProbability * std::exp(-labelLambda * ItD.Get()));
        }
        else if (ItD.Get() < 0)
        {
          ItD.Set(NumericTraits<RealType>::OneValue() -
                  (NumericTraits<RealType>::OneValue() - labelBoundaryProbability) *
                    (maximumInteriorDistance - Math::abs(ItD.Get())) / (maximumInteriorDistance));
        }
      }

      //
      // Normalize the distance prior probability image(s).
      //
      ImageRegionIterator<RealImageType> ItS(this->m_SumDistancePriorProbabilityImage,
                                             this->m_SumDistancePriorProbabilityImage->GetRequestedRegion());
      for (ItD.GoToBegin(), ItS.GoToBegin(); !ItS.IsAtEnd(); ++ItD, ++ItS)
      {
        if (!this->GetMaskImage() ||
            this->GetMaskImage()->GetPixel(ItD.GetIndex()) != NumericTraits<MaskLabelType>::ZeroValue())
        {
          if (ItS.Get() <= this->m_ProbabilityThreshold)
          {
            ItD.Set(NumericTraits<RealType>::ZeroValue());
          }
          else
          {
            ItD.Set(ItD.Get() / ItS.Get());
          }
        }
      }
      return distancePriorProbabilityImage;
    }
  }
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::RealImagePointer
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::GetSmoothIntensityImageFromPriorImage(
  unsigned int whichImage,
  unsigned int whichClass)
{
  typename ScalarImageType::Pointer bsplineImage;

  if (this->m_ControlPointLattices[whichImage][whichClass - 1].GetPointer() != nullptr)
  {
    typedef BSplineControlPointImageFilter<ControlPointLatticeType, ScalarImageType> BSplineReconstructorType;
    typename BSplineReconstructorType::Pointer bspliner = BSplineReconstructorType::New();

    bspliner->SetInput(this->m_ControlPointLattices[whichImage][whichClass - 1]);
    bspliner->SetSize(this->GetInput()->GetRequestedRegion().GetSize());
    bspliner->SetSpacing(this->GetInput()->GetSpacing());
    bspliner->SetOrigin(this->GetInput()->GetOrigin());
    bspliner->SetDirection(this->GetInput()->GetDirection());
    bspliner->SetSplineOrder(this->m_SplineOrder);
    bspliner->Update();

    bsplineImage = bspliner->GetOutput();
  }
  else
  {
    typename PointSetType::Pointer points = PointSetType::New();
    points->Initialize();

    typedef typename BSplineFilterType::WeightsContainerType WeightsType;
    typename WeightsType::Pointer                            weights = WeightsType::New();
    weights->Initialize();

    RealImagePointer probabilityImage;
    if (this->m_InitializationStrategy == PriorProbabilityImages)
    {
      probabilityImage = this->GetPriorProbabilityImage(whichClass);
    }
    else
    {
      typedef BinaryThresholdImageFilter<ClassifiedImageType, RealImageType> ThresholderType;
      typename ThresholderType::Pointer                                      thresholder = ThresholderType::New();
      thresholder->SetInput(const_cast<ClassifiedImageType *>(this->GetPriorLabelImage()));
      thresholder->SetInsideValue(1);
      thresholder->SetOutsideValue(0);
      thresholder->SetLowerThreshold(static_cast<LabelType>(whichClass));
      thresholder->SetUpperThreshold(static_cast<LabelType>(whichClass));
      thresholder->Update();

      probabilityImage = thresholder->GetOutput();
    }

    typename RealImageType::DirectionType originalDirection = probabilityImage->GetDirection();
    typename RealImageType::DirectionType identity;
    identity.SetIdentity();
    probabilityImage->SetDirection(identity);

    unsigned long count = 0;

    ImageRegionIteratorWithIndex<RealImageType> ItP(probabilityImage, probabilityImage->GetBufferedRegion());
    for (ItP.GoToBegin(); !ItP.IsAtEnd(); ++ItP)
    {
      if (!this->GetMaskImage() ||
          this->GetMaskImage()->GetPixel(ItP.GetIndex()) != NumericTraits<MaskLabelType>::ZeroValue())
      {
        if (ItP.Get() >= static_cast<RealType>(0.5))
        {
          typename RealImageType::PointType imagePoint;
          probabilityImage->TransformIndexToPhysicalPoint(ItP.GetIndex(), imagePoint);

          typename PointSetType::PointType bsplinePoint;
          bsplinePoint.CastFrom(imagePoint);

          ScalarType intensity;

          intensity[0] = this->GetIntensityImage(whichImage)->GetPixel(ItP.GetIndex());

          points->SetPoint(count, bsplinePoint);
          points->SetPointData(count, intensity);
          weights->InsertElement(count, ItP.Get());

          count++;
        }
      }
    }
    probabilityImage->SetDirection(originalDirection);

    typename BSplineFilterType::ArrayType numberOfControlPoints;
    typename BSplineFilterType::ArrayType numberOfLevels;
    for (unsigned int d = 0; d < ImageDimension; d++)
    {
      numberOfControlPoints[d] = this->m_NumberOfControlPoints[d];
      numberOfLevels[d] = this->m_NumberOfLevels[d];
    }

    typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
    bspliner->SetInput(points);
    bspliner->SetPointWeights(weights);
    bspliner->SetNumberOfLevels(numberOfLevels);
    bspliner->SetSplineOrder(this->m_SplineOrder);
    bspliner->SetNumberOfControlPoints(numberOfControlPoints);
    bspliner->SetSize(this->GetOutput()->GetLargestPossibleRegion().GetSize());
    bspliner->SetOrigin(this->GetOutput()->GetOrigin());
    bspliner->SetDirection(this->GetOutput()->GetDirection());
    bspliner->SetSpacing(this->GetOutput()->GetSpacing());
    bspliner->SetGenerateOutputImage(true);
    bspliner->Update();

    bsplineImage = bspliner->GetOutput();

    this->m_ControlPointLattices[whichImage][whichClass - 1] = bspliner->GetPhiLattice();
  }

  typedef VectorIndexSelectionCastImageFilter<ScalarImageType, RealImageType> CasterType;
  typename CasterType::Pointer                                                caster = CasterType::New();
  caster->SetInput(bsplineImage);
  caster->SetIndex(0);
  caster->Update();

  return caster->GetOutput();
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
typename AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::RealImagePointer
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::GetLikelihoodImage(unsigned int whichClass)
{
  RealImagePointer likelihoodImage = RealImageType::New();

  likelihoodImage->CopyInformation(this->GetInput());
  likelihoodImage->SetRegions(this->GetInput()->GetRequestedRegion());
  likelihoodImage->AllocateInitialized();

  std::vector<RealImagePointer> smoothImages;
  if (this->m_InitializationStrategy == PriorProbabilityImages || this->m_InitializationStrategy == PriorLabelImage)
  {
    for (unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++)
    {
      if (this->m_AdaptiveSmoothingWeights.size() > i &&
          this->m_AdaptiveSmoothingWeights[i] > NumericTraits<RealType>::ZeroValue())
      {
        smoothImages.push_back(this->GetSmoothIntensityImageFromPriorImage(i, whichClass));
      }
      else
      {
        smoothImages.push_back(nullptr);
      }
    }
  }

  ImageRegionIteratorWithIndex<RealImageType> It(likelihoodImage, likelihoodImage->GetRequestedRegion());
  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    if (!this->GetMaskImage() ||
        this->GetMaskImage()->GetPixel(It.GetIndex()) != NumericTraits<MaskLabelType>::ZeroValue())
    {
      MeasurementVectorType measurement;
      measurement.SetSize(this->m_NumberOfIntensityImages);
      for (unsigned int i = 0; i < this->m_NumberOfIntensityImages; i++)
      {
        measurement[i] = this->GetIntensityImage(i)->GetPixel(It.GetIndex());

        if ((this->m_InitializationStrategy == PriorProbabilityImages ||
             this->m_InitializationStrategy == PriorLabelImage) &&
            smoothImages[i])
        {
          measurement[i] =
            (NumericTraits<RealType>::OneValue() - this->m_AdaptiveSmoothingWeights[i]) * measurement[i] +
            this->m_AdaptiveSmoothingWeights[i] * smoothImages[i]->GetPixel(It.GetIndex());
        }
      }
      RealType likelihood = this->m_MixtureModelComponents[whichClass - 1]->Evaluate(measurement);
      It.Set(likelihood);
    }
  }

  return likelihoodImage;
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::ComputeICMCodeImage()
{
  if (!this->m_ICMCodeImage)
  {
    this->m_ICMCodeImage = ClassifiedImageType::New();
    this->m_ICMCodeImage->CopyInformation(this->GetInput());
    this->m_ICMCodeImage->SetRegions(this->GetInput()->GetRequestedRegion());
    this->m_ICMCodeImage->AllocateInitialized();
  }

  typename NeighborhoodIterator<ClassifiedImageType>::RadiusType radius;
  unsigned int                                                   neighborhoodSize = 1;
  for (unsigned int d = 0; d < ImageDimension; d++)
  {
    neighborhoodSize *= (2 * this->m_MRFRadius[d] + 1);
    radius[d] = this->m_MRFRadius[d];
  }

  NeighborhoodIterator<ClassifiedImageType> It(
    radius, this->m_ICMCodeImage, this->m_ICMCodeImage->GetRequestedRegion());
  this->m_MaximumICMCode = 1;

  bool codingComplete = false;
  while (!codingComplete)
  {
    codingComplete = true;

    It.GoToBegin();
    while (!It.IsAtEnd())
    {
      if (this->m_ICMCodeImage->GetPixel(It.GetIndex()) == 0 &&
          (!this->GetMaskImage() ||
           this->GetMaskImage()->GetPixel(It.GetIndex()) != NumericTraits<MaskLabelType>::ZeroValue()))
      {
        bool hasCode = false;
        for (unsigned int n = 0; n < neighborhoodSize; n++)
        {
          bool      isInBounds = false;
          LabelType label = It.GetPixel(n, isInBounds);
          if (isInBounds && label == this->m_MaximumICMCode)
          {
            hasCode = true;
            break;
          }
        }
        if (!hasCode)
        {
          It.SetCenterPixel(this->m_MaximumICMCode);
          codingComplete = false;
        }
      }
      ++It;
    }

    if (codingComplete == false)
    {
      this->m_MaximumICMCode++;
    }
  }

  this->m_MaximumICMCode--;
}

template <typename TInputImage, typename TMaskImage, typename TClassifiedImage>
void
AtroposSegmentationImageFilter<TInputImage, TMaskImage, TClassifiedImage>::PrintSelf(std::ostream & os,
                                                                                     Indent         indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Maximum number of iterations: " << this->m_MaximumNumberOfIterations << std::endl;
  os << indent << "Convergence threshold: " << this->m_ConvergenceThreshold << std::endl;
  os << indent << "Number of tissue classes: " << this->m_NumberOfTissueClasses << std::endl;
  os << indent << "Number of partial volume classes: " << this->m_NumberOfPartialVolumeClasses << std::endl;
  os << indent << "Minimize memory usage:";
  if (this->m_MinimizeMemoryUsage)
  {
    os << " true";
    if (this->m_MinimizeMemoryUsage && this->m_InitializationStrategy == PriorProbabilityImages)
    {
      os << " (prior probability threshold = " << this->m_ProbabilityThreshold << ")" << std::endl;
    }
    os << std::endl;
  }
  else
  {
    os << " false" << std::endl;
  }

  os << indent << "Initialization strategy: ";

  switch (this->m_InitializationStrategy)
  {
    case Random:
    {
      os << "Random" << std::endl;
    }
    break;
    case KMeans:
    {
      os << "K means clustering" << std::endl;
    }
    break;
    case Otsu:
    {
      os << "Otsu thresholding" << std::endl;
    }
    break;
    case PriorProbabilityImages:
    {
      os << "Prior probability images" << std::endl;
      os << indent << "  Use Euclidean distance for prior labels:";
      if (this->m_UseEuclideanDistanceForPriorLabels)
      {
        os << " true" << std::endl;
      }
      else
      {
        os << " false" << std::endl;
      }
      if (this->m_PriorLabelParameterMap.size() > 0)
      {
        os << indent << "  Specified prior label parameters:" << std::endl;
        typename LabelParameterMapType::const_iterator it;
        for (it = this->m_PriorLabelParameterMap.begin(); it != this->m_PriorLabelParameterMap.end(); ++it)
        {
          RealType label = it->first;
          RealType lambda = (it->second).first;
          RealType boundaryProbability = (it->second).second;
          os << indent << "    Class " << label << ": lambda = " << lambda
             << ", boundary probability = " << boundaryProbability << std::endl;
        }
      }
    }
    break;
    case PriorLabelImage:
    {
      os << "Prior label image" << std::endl;
      os << indent << "  Use Euclidean distance for prior labels:";
      if (this->m_UseEuclideanDistanceForPriorLabels)
      {
        os << " true" << std::endl;
      }
      else
      {
        os << " false" << std::endl;
      }
      os << indent << "  Specified prior label parameters:" << std::endl;
      typename LabelParameterMapType::const_iterator it;
      for (it = this->m_PriorLabelParameterMap.begin(); it != this->m_PriorLabelParameterMap.end(); ++it)
      {
        RealType label = it->first;
        RealType lambda = (it->second).first;
        RealType boundaryProbability = (it->second).second;
        os << indent << "    Class " << label << ": lambda = " << lambda
           << ", boundary probability = " << boundaryProbability << std::endl;
      }
    }
    break;
  }
  os << indent << "Posterior probability formulation: ";

  switch (this->m_PosteriorProbabilityFormulation)
  {
    case Socrates:
    {
      os << "Socrates" << std::endl;
    }
    break;
    case Plato:
    {
      os << "Plato" << std::endl;
    }
    break;
    case Aristotle:
    {
      os << "Aristotle" << std::endl;
    }
    break;
    case Sigmoid:
    {
      os << "Sigmoid" << std::endl;
    }
    break;
  }
  os << indent << "  initial annealing temperature = " << this->m_InitialAnnealingTemperature << std::endl;
  os << indent << "  annealing rate = " << this->m_AnnealingRate << std::endl;
  os << indent << "  minimum annealing temperature = " << this->m_MinimumAnnealingTemperature << std::endl;

  os << indent << "MRF parameters" << std::endl;
  if (this->m_MRFCoefficientImage)
  {
    os << indent << "  Using MRF coefficient image" << std::endl;
  }
  else
  {
    os << indent << "  MRF smoothing factor = " << this->m_MRFSmoothingFactor << std::endl;
  }
  os << indent << "  MRF radius = " << this->m_MRFRadius << std::endl;

  if (this->m_UseAsynchronousUpdating)
  {
    os << indent << "Use asynchronous updating of the labels." << std::endl;
    os << indent << "  ICM parameters" << std::endl;
    os << indent << "    maximum ICM code = " << this->m_MaximumICMCode << std::endl;
    os << indent << "    maximum number of ICM iterations = " << this->m_MaximumNumberOfICMIterations << std::endl;
  }
  else
  {
    os << indent << "Use synchronous updating of the labels." << std::endl;
  }

  if (this->m_OutlierHandlingFilter)
  {
    os << indent << "Outlier handling " << std::endl;
    this->m_OutlierHandlingFilter->Print(os, indent.GetNextIndent());
  }
  else
  {
    os << indent << "No outlier handling." << std::endl;
  }

  if ((this->m_InitializationStrategy == PriorProbabilityImages || this->m_InitializationStrategy == PriorLabelImage) &&
      this->m_AdaptiveSmoothingWeights.size() > 0)
  {
    os << indent << "Adaptive smoothing weights: [";
    for (unsigned int i = 0; i < this->m_AdaptiveSmoothingWeights.size() - 1; i++)
    {
      os << this->m_AdaptiveSmoothingWeights[i] << ", ";
    }
    os << this->m_AdaptiveSmoothingWeights[this->m_AdaptiveSmoothingWeights.size() - 1] << "]" << std::endl;
    os << indent << "B-spline smoothing" << std::endl;
    os << indent << "  spline order = " << this->m_SplineOrder << std::endl;
    os << indent << "  number of levels = " << this->m_NumberOfLevels << std::endl;
    os << indent << "  number of initial control points = " << this->m_NumberOfControlPoints << std::endl;
  }
  for (unsigned int n = 0; n < this->m_NumberOfTissueClasses; n++)
  {
    if (this->m_MixtureModelProportions.size() > n)
    {
      os << indent << "Tissue class " << n + 1 << ": proportion = " << this->m_MixtureModelProportions[n] << std::endl;
    }
    if (this->m_MixtureModelComponents.size() > n)
    {
      this->m_MixtureModelComponents[n]->Print(os, indent.GetNextIndent());
    }
  }

  if (this->m_UsePartialVolumeLikelihoods)
  {
    unsigned int                                      n = this->m_NumberOfTissueClasses;
    typename PartialVolumeClassesType::const_iterator it;
    for (it = this->m_PartialVolumeClasses.begin(); it != this->m_PartialVolumeClasses.end(); ++it)
    {
      os << indent << "Partial volume class " << n + 1 << " [labels: ";
      for (unsigned int l = 0; l < it->size() - 1; l++)
      {
        os << (*it)[l] << 'x';
      }
      os << (*it)[it->size() - 1] << "]: proportion = " << this->m_MixtureModelProportions[n] << std::endl;

      this->m_MixtureModelComponents[n++]->Print(os, indent.GetNextIndent());
    }
  }
}
} // namespace ants
} // namespace itk

#endif
