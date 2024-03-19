/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDiReCTImageFilter_hxx
#define __itkDiReCTImageFilter_hxx


#include "itkAddImageFilter.h"
#include "itkBinaryContourImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkComposeDisplacementFieldsImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkDisplacementFieldToBSplineImageFilter.h"
#include "itkGaussianOperator.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImportImageFilter.h"
#include "itkInvertDisplacementFieldImageFilter.h"
#include "itkIterationReporter.h"
#include "itkMaskedSmoothingImageFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkPasteImageFilter.h"
#include "itkMultiplyByConstantImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkVectorNeighborhoodOperatorImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkWindowConvergenceMonitoringFunction.h"

namespace itk
{
template <typename TInputImage, typename TOutputImage>
DiReCTImageFilter<TInputImage, TOutputImage>::DiReCTImageFilter()
  : m_ThicknessPriorEstimate(10.0)
  , m_SmoothingVariance(1.0)
  , m_SmoothingVelocityFieldVariance(1.5)
  , m_BSplineSmoothingIsotropicMeshSpacing(5.75)
  , m_InitialGradientStep(0.025)
  , m_CurrentGradientStep(0.025)
  , m_NumberOfIntegrationPoints(10)
  , m_SparseImageNeighborhoodRadius(2)
  , m_GrayMatterLabel(2)
  , m_WhiteMatterLabel(3)
  , m_MaximumNumberOfIterations(45)
  , m_MaximumNumberOfInvertDisplacementFieldIterations(20)
  , m_CurrentEnergy(NumericTraits<RealType>::max())
  , m_ConvergenceThreshold(0.0)
  , m_ConvergenceWindowSize(10)
  , m_UseBSplineSmoothing(false)
  , m_UseMaskedSmoothing(false)
  , m_RestrictDeformation(false)
  , m_IncludeCumulativeVelocityFields(false)
  , m_TimeSmoothingVariance(1.0)
{
  this->m_ThicknessPriorImage = nullptr;
  this->m_SparseMatrixIndexImage = nullptr;
  this->m_TimePoints.clear();

  this->SetNumberOfRequiredOutputs(3);

  // thickness image
  this->SetNthOutput(0, this->MakeOutput(0));

  // forward/inverse cumulative velocity fields
  this->SetNthOutput(1, this->MakeOutput(1));
  this->SetNthOutput(2, this->MakeOutput(2));
}

template <typename TInputImage, typename TOutputImage>
DiReCTImageFilter<TInputImage, TOutputImage>::~DiReCTImageFilter() = default;

template <typename TInputImage, typename TOutputImage>
auto
DiReCTImageFilter<TInputImage, TOutputImage>::MakeOutput(
  DataObjectPointerArraySizeType idx) -> DataObjectPointer
{
  if (idx > 0)
  {
    return CumulativeVelocityFieldType::New().GetPointer();
  }
  return Superclass::MakeOutput(idx);
} 

template <typename TInputImage, typename TOutputImage>
auto
DiReCTImageFilter<TInputImage, TOutputImage>::GetThicknessImage() -> OutputImageType *
{
  return dynamic_cast<OutputImageType *>(this->ProcessObject::GetOutput(0));
}

template <typename TInputImage, typename TOutputImage>
auto
DiReCTImageFilter<TInputImage, TOutputImage>::GetForwardCumulativeVelocityField() -> CumulativeVelocityFieldType *
{
  return dynamic_cast<CumulativeVelocityFieldType *>(this->ProcessObject::GetOutput(1));
}

template <typename TInputImage, typename TOutputImage>
auto
DiReCTImageFilter<TInputImage, TOutputImage>::GetInverseCumulativeVelocityField() -> CumulativeVelocityFieldType *
{
  return dynamic_cast<CumulativeVelocityFieldType *>(this->ProcessObject::GetOutput(2));
}

template <typename TInputImage, typename TOutputImage>
void
DiReCTImageFilter<TInputImage, TOutputImage>::GenerateOutputInformation()
{
  RealImagePointer corticalThicknessImage = this->GetThicknessImage();
  corticalThicknessImage->CopyInformation(this->GetInput());
  corticalThicknessImage->SetRegions(this->GetInput()->GetLargestPossibleRegion());

  if (this->m_IncludeCumulativeVelocityFields)
  {
    CumulativeVelocityFieldPointer forwardCumulativeVelocityField = this->GetForwardCumulativeVelocityField();
    CumulativeVelocityFieldPointer inverseCumulativeVelocityField = this->GetInverseCumulativeVelocityField();

    typename CumulativeVelocityFieldType::PointType origin;
    typename CumulativeVelocityFieldType::DirectionType direction;
    typename CumulativeVelocityFieldType::SpacingType spacing;
    typename CumulativeVelocityFieldType::RegionType region;
    typename CumulativeVelocityFieldType::RegionType::SizeType size;
    
    for (unsigned int d = 0; d < ImageDimension; d++ )
    {
      origin[d] = this->GetInput()->GetOrigin()[d];
      spacing[d] = this->GetInput()->GetSpacing()[d];
      size[d] = this->GetInput()->GetLargestPossibleRegion().GetSize()[d];
    }
    origin[ImageDimension] = 0.0;
    spacing[ImageDimension] = 1.0;
    direction.SetIdentity();
    size[ImageDimension] = this->m_NumberOfIntegrationPoints;

    inverseCumulativeVelocityField->SetSpacing(spacing);
    inverseCumulativeVelocityField->SetOrigin(origin);
    inverseCumulativeVelocityField->SetDirection(direction);
    inverseCumulativeVelocityField->SetRegions(size);

    forwardCumulativeVelocityField->SetSpacing(spacing);
    forwardCumulativeVelocityField->SetOrigin(origin);
    forwardCumulativeVelocityField->SetDirection(direction);
    forwardCumulativeVelocityField->SetRegions(size);
  } 
}

template <typename TInputImage, typename TOutputImage>
void
DiReCTImageFilter<TInputImage, TOutputImage>::GenerateData()
{
  this->m_CurrentGradientStep = this->m_InitialGradientStep;

  // Convert all input direction matrices to identities saving the original
  // directions to put back at the end of filter processing. We do this
  // because the white and gray matters reside in the same space and the
  // assumption simplifies the underlying registration code.

  const bool filterHandlesMemory = false;

  typename InputImageType::DirectionType identity;
  identity.SetIdentity();

  using SegmentationImageImporterType = ImportImageFilter<InputPixelType, ImageDimension>;
  typename SegmentationImageImporterType::Pointer segmentationImageImporter = SegmentationImageImporterType::New();
  segmentationImageImporter->SetImportPointer(
    const_cast<InputPixelType *>(this->GetSegmentationImage()->GetBufferPointer()),
    (this->GetSegmentationImage()->GetBufferedRegion()).GetNumberOfPixels(),
    filterHandlesMemory);
  segmentationImageImporter->SetRegion(this->GetSegmentationImage()->GetBufferedRegion());
  segmentationImageImporter->SetOrigin(this->GetSegmentationImage()->GetOrigin());
  segmentationImageImporter->SetSpacing(this->GetSegmentationImage()->GetSpacing());
  segmentationImageImporter->SetDirection(identity);

  InputImagePointer segmentationImage = segmentationImageImporter->GetOutput();
  segmentationImage->Update();
  segmentationImage->DisconnectPipeline();

  using ProbablilityImageImporterType = ImportImageFilter<RealType, ImageDimension>;

  typename ProbablilityImageImporterType::Pointer grayMatterProbabilityImageImporter =
    ProbablilityImageImporterType::New();
  grayMatterProbabilityImageImporter->SetImportPointer(
    const_cast<RealType *>(this->GetGrayMatterProbabilityImage()->GetBufferPointer()),
    (this->GetGrayMatterProbabilityImage()->GetBufferedRegion()).GetNumberOfPixels(),
    filterHandlesMemory);
  grayMatterProbabilityImageImporter->SetRegion(this->GetGrayMatterProbabilityImage()->GetBufferedRegion());
  grayMatterProbabilityImageImporter->SetOrigin(this->GetGrayMatterProbabilityImage()->GetOrigin());
  grayMatterProbabilityImageImporter->SetSpacing(this->GetGrayMatterProbabilityImage()->GetSpacing());
  grayMatterProbabilityImageImporter->SetDirection(identity);

  RealImagePointer grayMatterProbabilityImage = grayMatterProbabilityImageImporter->GetOutput();
  grayMatterProbabilityImage->Update();
  grayMatterProbabilityImage->DisconnectPipeline();

  typename ProbablilityImageImporterType::Pointer whiteMatterProbabilityImageImporter =
    ProbablilityImageImporterType::New();
  whiteMatterProbabilityImageImporter->SetImportPointer(
    const_cast<RealType *>(this->GetWhiteMatterProbabilityImage()->GetBufferPointer()),
    (this->GetWhiteMatterProbabilityImage()->GetBufferedRegion()).GetNumberOfPixels(),
    filterHandlesMemory);
  whiteMatterProbabilityImageImporter->SetRegion(this->GetWhiteMatterProbabilityImage()->GetBufferedRegion());
  whiteMatterProbabilityImageImporter->SetOrigin(this->GetWhiteMatterProbabilityImage()->GetOrigin());
  whiteMatterProbabilityImageImporter->SetSpacing(this->GetWhiteMatterProbabilityImage()->GetSpacing());
  whiteMatterProbabilityImageImporter->SetDirection(identity);

  RealImagePointer whiteMatterProbabilityImage = whiteMatterProbabilityImageImporter->GetOutput();
  whiteMatterProbabilityImage->Update();
  whiteMatterProbabilityImage->DisconnectPipeline();

  // Extract the gray and white matter segmentations and combine to form the
  // gm/wm region.  Dilate the latter region by 1 voxel.

  InputImagePointer grayMatter = this->ExtractRegion(segmentationImage, this->m_GrayMatterLabel);
  InputImagePointer whiteMatter = this->ExtractRegion(segmentationImage, this->m_WhiteMatterLabel);

  using AdderType = AddImageFilter<InputImageType, InputImageType, InputImageType>;
  typename AdderType::Pointer adder = AdderType::New();
  adder->SetInput1(grayMatter);
  adder->SetInput2(whiteMatter);
  adder->Update();

  InputImagePointer thresholdedRegion = this->ExtractRegion(const_cast<const InputImageType *>(adder->GetOutput()), 1);

  // Extract the white and gm/wm matter contours

  InputImagePointer matterContours = this->ExtractRegionalContours(thresholdedRegion, 1);
  InputImagePointer whiteMatterContoursTmp = this->ExtractRegionalContours(segmentationImage, this->m_WhiteMatterLabel);

  using CasterType = CastImageFilter<InputImageType, RealImageType>;
  typename CasterType::Pointer caster = CasterType::New();
  caster->SetInput(whiteMatterContoursTmp);
  caster->Update();
  RealImagePointer whiteMatterContours = caster->GetOutput();

  // Initialize fields and images.
  VectorType zeroVector(0.0);

  RealImagePointer corticalThicknessImage = this->GetThicknessImage();
  corticalThicknessImage->SetDirection(segmentationImage->GetDirection());
  corticalThicknessImage->AllocateInitialized();

  if (this->m_IncludeCumulativeVelocityFields)
  {
    CumulativeVelocityFieldPointer forwardCumulativeVelocityField = this->GetForwardCumulativeVelocityField();
    forwardCumulativeVelocityField->AllocateInitialized();
    CumulativeVelocityFieldPointer inverseCumulativeVelocityField = this->GetInverseCumulativeVelocityField();
    inverseCumulativeVelocityField->AllocateInitialized();
  }

  DisplacementFieldPointer forwardIncrementalField = DisplacementFieldType::New();
  forwardIncrementalField->CopyInformation(segmentationImage);
  forwardIncrementalField->SetRegions(segmentationImage->GetRequestedRegion());
  forwardIncrementalField->AllocateInitialized();

  RealImagePointer hitImage = RealImageType::New();
  hitImage->CopyInformation(segmentationImage);
  hitImage->SetRegions(segmentationImage->GetRequestedRegion());
  hitImage->AllocateInitialized();

  DisplacementFieldPointer integratedField = DisplacementFieldType::New();
  integratedField->CopyInformation(segmentationImage);
  integratedField->SetRegions(segmentationImage->GetRequestedRegion());
  integratedField->AllocateInitialized();

  DisplacementFieldPointer inverseField = DisplacementFieldType::New();
  inverseField->CopyInformation(segmentationImage);
  inverseField->SetRegions(segmentationImage->GetRequestedRegion());
  inverseField->AllocateInitialized();

  DisplacementFieldPointer inverseIncrementalField = DisplacementFieldType::New();
  inverseIncrementalField->CopyInformation(segmentationImage);
  inverseIncrementalField->SetRegions(segmentationImage->GetRequestedRegion());
  inverseIncrementalField->AllocateInitialized();

  RealImagePointer speedImage = RealImageType::New();
  speedImage->CopyInformation(segmentationImage);
  speedImage->SetRegions(segmentationImage->GetRequestedRegion());
  speedImage->AllocateInitialized();

  RealImagePointer thicknessImage = RealImageType::New();
  thicknessImage->CopyInformation(segmentationImage);
  thicknessImage->SetRegions(segmentationImage->GetRequestedRegion());
  thicknessImage->AllocateInitialized();

  RealImagePointer totalImage = RealImageType::New();
  totalImage->CopyInformation(segmentationImage);
  totalImage->SetRegions(segmentationImage->GetRequestedRegion());
  totalImage->AllocateInitialized();

  DisplacementFieldPointer velocityField = DisplacementFieldType::New();
  velocityField->CopyInformation(segmentationImage);
  velocityField->SetRegions(segmentationImage->GetRequestedRegion());
  velocityField->AllocateInitialized();

  // Instantiate most of the iterators all in one place

  ImageRegionIterator<RealImageType>         ItCorticalThicknessImage(corticalThicknessImage,
                                                                      corticalThicknessImage->GetRequestedRegion());
  ImageRegionConstIterator<RealImageType>    ItGrayMatterProbabilityMap(grayMatterProbabilityImage,
                                                                        grayMatterProbabilityImage->GetRequestedRegion());
  ImageRegionIterator<RealImageType>         ItHitImage(hitImage, hitImage->GetRequestedRegion());
  ImageRegionIterator<DisplacementFieldType> ItForwardIncrementalField(forwardIncrementalField,
                                                                       forwardIncrementalField->GetRequestedRegion());
  ImageRegionConstIterator<InputImageType>   ItMatterContours(matterContours, matterContours->GetRequestedRegion());
  ImageRegionIterator<DisplacementFieldType> ItInverseIncrementalField(inverseIncrementalField,
                                                                       inverseIncrementalField->GetRequestedRegion());
  ImageRegionConstIteratorWithIndex<InputImageType> ItSegmentationImage(segmentationImage,
                                                                        segmentationImage->GetRequestedRegion());
  ImageRegionIterator<RealImageType>                ItSpeedImage(speedImage, speedImage->GetRequestedRegion());
  ImageRegionIterator<RealImageType>      ItThicknessImage(thicknessImage, thicknessImage->GetRequestedRegion());
  ImageRegionIterator<RealImageType>      ItTotalImage(totalImage, totalImage->GetRequestedRegion());
  ImageRegionConstIterator<RealImageType> ItWhiteMatterContours(whiteMatterContours,
                                                                whiteMatterContours->GetRequestedRegion());

  // Monitor the convergence
  using ConvergenceMonitoringType = typename Function::WindowConvergenceMonitoringFunction<double>;
  ConvergenceMonitoringType::Pointer convergenceMonitoring = ConvergenceMonitoringType::New();
  convergenceMonitoring->SetWindowSize(this->m_ConvergenceWindowSize);

  // Instantiate the progress reporter

  IterationReporter reporter(this, 0, 1);

  bool isConverged = false;
  this->m_CurrentConvergenceMeasurement = NumericTraits<RealType>::max();
  this->m_ElapsedIterations = 0;
  while (this->m_ElapsedIterations++ < this->m_MaximumNumberOfIterations && isConverged == false)
  {
    this->MakeThicknessImage(hitImage, totalImage, segmentationImage, thicknessImage);

    RealType      priorEnergy = 0;
    unsigned long priorEnergyCount = 0;

    RealType currentEnergy = 0.0;
    RealType numberOfGrayMatterVoxels = 0.0;

    forwardIncrementalField->FillBuffer(zeroVector);
    inverseField->FillBuffer(zeroVector);
    inverseIncrementalField->FillBuffer(zeroVector);

    hitImage->FillBuffer(0.0);
    totalImage->FillBuffer(0.0);

    ImageRegionIterator<DisplacementFieldType> ItVelocityField(velocityField, velocityField->GetRequestedRegion());

    unsigned int integrationPoint = 0;
    while (integrationPoint++ < this->m_NumberOfIntegrationPoints)
    {
      using ComposerType = ComposeDisplacementFieldsImageFilter<DisplacementFieldType>;
      typename ComposerType::Pointer composer = ComposerType::New();
      composer->SetDisplacementField(inverseIncrementalField);
      composer->SetWarpingField(inverseField);

      inverseField = composer->GetOutput();
      inverseField->Update();
      inverseField->DisconnectPipeline();

      RealImagePointer warpedWhiteMatterProbabilityImage = this->WarpImage(whiteMatterProbabilityImage, inverseField);
      RealImagePointer warpedWhiteMatterContours = this->WarpImage(whiteMatterContours, inverseField);
      RealImagePointer warpedThicknessImage = this->WarpImage(thicknessImage, inverseField);

      using GradientImageFilterType = GradientRecursiveGaussianImageFilter<RealImageType, DisplacementFieldType>;
      typename GradientImageFilterType::Pointer gradientFilter = GradientImageFilterType::New();
      gradientFilter->SetInput(warpedWhiteMatterProbabilityImage);
      gradientFilter->SetSigma(this->m_SmoothingVariance);
      gradientFilter->Update();

      DisplacementFieldPointer gradientImage = gradientFilter->GetOutput();

      // Instantiate the iterators all in one place

      ImageRegionIterator<DisplacementFieldType> ItGradientImage(gradientImage, gradientImage->GetRequestedRegion());
      ImageRegionIterator<RealImageType>         ItWarpedWhiteMatterProbabilityMap(
        warpedWhiteMatterProbabilityImage, warpedWhiteMatterProbabilityImage->GetRequestedRegion());
      ImageRegionIterator<RealImageType>         ItWarpedWhiteMatterContours(warpedWhiteMatterContours,
                                                                     warpedWhiteMatterContours->GetRequestedRegion());
      ImageRegionIterator<DisplacementFieldType> ItInverseField(inverseField, inverseField->GetRequestedRegion());
      ImageRegionIterator<DisplacementFieldType> ItIntegratedField(integratedField,
                                                                   integratedField->GetRequestedRegion());

      // Generate speed image

      speedImage->FillBuffer(NumericTraits<RealType>::Zero);

      ItGradientImage.GoToBegin();
      ItGrayMatterProbabilityMap.GoToBegin();
      ItSegmentationImage.GoToBegin();
      ItSpeedImage.GoToBegin();
      ItWarpedWhiteMatterProbabilityMap.GoToBegin();

      const typename InputImageType::PixelType grayMatterPixel = static_cast<typename InputImageType::PixelType>(
                                                 this->m_GrayMatterLabel),
                                               whiteMatterPixel = static_cast<typename InputImageType::PixelType>(
                                                 this->m_WhiteMatterLabel);
      while (!ItSegmentationImage.IsAtEnd())
      {
        if (ItSegmentationImage.Get() == grayMatterPixel)
        {
          RealType norm = (ItGradientImage.Get()).GetNorm();
          if (norm > static_cast<RealType>(1e-3) && !std::isnan(norm) && !std::isinf(norm))
          {
            ItGradientImage.Set(ItGradientImage.Get() / norm);
          }
          else
          {
            ItGradientImage.Set(zeroVector);
          }
          RealType delta = (ItWarpedWhiteMatterProbabilityMap.Get() - ItGrayMatterProbabilityMap.Get());
          currentEnergy += Math::abs(delta);
          numberOfGrayMatterVoxels++;
          RealType speedValue = -delta * ItGrayMatterProbabilityMap.Get() * this->m_CurrentGradientStep;
          if (std::isnan(speedValue) || std::isinf(speedValue))
          {
            speedValue = NumericTraits<RealType>::Zero;
          }
          ItSpeedImage.Set(speedValue);
        }
        else
        {
          ItSpeedImage.Set(NumericTraits<RealType>::Zero);
        }

        ++ItGradientImage;
        ++ItGrayMatterProbabilityMap;
        ++ItSegmentationImage;
        ++ItSpeedImage;
        ++ItWarpedWhiteMatterProbabilityMap;
      }

      // Calculate objective function value

      ItForwardIncrementalField.GoToBegin();
      ItGradientImage.GoToBegin();
      ItSegmentationImage.GoToBegin();
      ItSpeedImage.GoToBegin();
      ItWhiteMatterContours.GoToBegin();

      while (!ItSegmentationImage.IsAtEnd())
      {
        typename InputImageType::IndexType index = ItSegmentationImage.GetIndex();
        typename InputImageType::PixelType segmentationValue = ItSegmentationImage.Get();

        ItForwardIncrementalField.Set(ItForwardIncrementalField.Get() + ItGradientImage.Get() * ItSpeedImage.Get());
        if (segmentationValue == grayMatterPixel || segmentationValue == whiteMatterPixel)
        {
          if (integrationPoint == 1)
          {
            typename InputImageType::PixelType whiteMatterContoursValue =
              static_cast<typename InputImageType::PixelType>(ItWhiteMatterContours.Get());
            hitImage->SetPixel(index, whiteMatterContoursValue);

            VectorType vector = integratedField->GetPixel(index);
            RealType   weightedNorm = vector.GetNorm() * whiteMatterContoursValue;

            thicknessImage->SetPixel(index, weightedNorm);
            totalImage->SetPixel(index, weightedNorm);
          }
          else if (segmentationValue == grayMatterPixel)
          {
            hitImage->SetPixel(index, hitImage->GetPixel(index) + warpedWhiteMatterContours->GetPixel(index));
            totalImage->SetPixel(index, totalImage->GetPixel(index) + warpedThicknessImage->GetPixel(index));
          }
        }

        ++ItForwardIncrementalField;
        ++ItGradientImage;
        ++ItSegmentationImage;
        ++ItSpeedImage;
        ++ItWhiteMatterContours;
      }

      ItSegmentationImage.GoToBegin();
      ItMatterContours.GoToBegin();
      ItWhiteMatterContours.GoToBegin();
      ItVelocityField.GoToBegin();
      ItInverseIncrementalField.GoToBegin();
      ItIntegratedField.GoToBegin();
      ItInverseField.GoToBegin();
      while (!ItSegmentationImage.IsAtEnd())
      {
        typename InputImageType::PixelType segmentationValue = ItSegmentationImage.Get();
        typename InputImageType::PixelType whiteMatterContoursValue =
          static_cast<typename InputImageType::PixelType>(ItWhiteMatterContours.Get());
        typename InputImageType::PixelType matterContoursValue = ItMatterContours.Get();

        if (segmentationValue == 0 ||
            (whiteMatterContoursValue == 0 && matterContoursValue == 0 && segmentationValue != this->m_GrayMatterLabel))
        {
          ItInverseField.Set(zeroVector);
          ItVelocityField.Set(zeroVector);
          ItIntegratedField.Set(zeroVector);
        }

        ItInverseIncrementalField.Set(ItVelocityField.Get());

        ++ItSegmentationImage;
        ++ItMatterContours;
        ++ItWhiteMatterContours;
        ++ItVelocityField;
        ++ItInverseIncrementalField;
        ++ItInverseField;
        ++ItIntegratedField;
      }

      if (integrationPoint == 1)
      {
        integratedField->FillBuffer(zeroVector);
      }

      using InverterType = InvertDisplacementFieldImageFilter<DisplacementFieldType>;

      typename InverterType::Pointer inverter1 = InverterType::New();
      inverter1->SetInput(inverseField);
      inverter1->SetInverseFieldInitialEstimate(integratedField);
      inverter1->SetMaximumNumberOfIterations(this->m_MaximumNumberOfInvertDisplacementFieldIterations);
      inverter1->SetMeanErrorToleranceThreshold(0.001);
      inverter1->SetMaxErrorToleranceThreshold(0.1);
      if (this->m_UseMaskedSmoothing)
      {
        inverter1->SetEnforceBoundaryCondition(false);
      }
      inverter1->Update();

      integratedField = inverter1->GetOutput();
      integratedField->DisconnectPipeline();

      typename InverterType::Pointer inverter2 = InverterType::New();
      inverter2->SetInput(integratedField);
      inverter2->SetInverseFieldInitialEstimate(inverseField);
      inverter2->SetMaximumNumberOfIterations(this->m_MaximumNumberOfInvertDisplacementFieldIterations);
      inverter2->SetMeanErrorToleranceThreshold(0.001);
      inverter2->SetMaxErrorToleranceThreshold(0.1);
      if (this->m_UseMaskedSmoothing)
      {
        inverter2->SetEnforceBoundaryCondition(false);
      }
      inverter2->Update();

      inverseField = inverter2->GetOutput();
      inverseField->DisconnectPipeline();

      if (this->m_IncludeCumulativeVelocityFields)
      {
        CumulativeVelocityFieldPointer forwardCumulativeVelocityField = this->GetForwardCumulativeVelocityField();
        CumulativeVelocityFieldPointer inverseCumulativeVelocityField = this->GetInverseCumulativeVelocityField();

        using PasterType = PasteImageFilter<CumulativeVelocityFieldType, DisplacementFieldType>;

        typename CumulativeVelocityFieldType::IndexType destinationIndex;
        destinationIndex.Fill(NumericTraits<typename IndexType::IndexValueType>::Zero);
        destinationIndex[ImageDimension] = static_cast<typename IndexType::IndexValueType>(integrationPoint);

        typename PasterType::Pointer paster1 = PasterType::New();
        paster1->SetSourceImage(integratedField);
        paster1->SetSourceRegion(integratedField->GetLargestPossibleRegion());
        paster1->SetDestinationImage(forwardCumulativeVelocityField);
        paster1->SetDestinationIndex(destinationIndex);
        paster1->Update();

        this->ProcessObject::SetNthOutput(1, paster1->GetOutput());

        typename PasterType::Pointer paster2 = PasterType::New();
        paster2->SetSourceImage(inverseField);
        paster2->SetSourceRegion(inverseField->GetLargestPossibleRegion());
        paster2->SetDestinationImage(inverseCumulativeVelocityField);
        paster2->SetDestinationIndex(destinationIndex);
        paster2->Update();

        this->ProcessObject::SetNthOutput(2, paster2->GetOutput());
      }
    }

    // calculate the size of the solution to allow us to adjust the
    // gradient step length.

    //     RealType maxNorm = 0.0;
    //
    //     typename InputImageType::SpacingType spacing = grayMatter->GetSpacing();
    //
    //     ImageRegionIterator<DisplacementFieldType> ItIntegratedField2(
    //       integratedField,
    //       integratedField->GetRequestedRegion() );
    //
    //     ItIntegratedField2.GoToBegin();
    //     for( ItIntegratedField2.GoToBegin(); !ItIntegratedField2.IsAtEnd();
    //          ++ItIntegratedField2 )
    //       {
    //       VectorType vector = ItIntegratedField2.Get();
    //       for( unsigned int d = 0; d < ImageDimension; d++ )
    //         {
    //         vector[d] = vector[d] / spacing[d];
    //         }
    //       RealType norm = vector.GetNorm();
    //       if( norm > maxNorm )
    //         {
    //         maxNorm = norm;
    //         }
    //       }
    //     itkDebugMacro( "   MaxNorm = " << maxNorm );
    //
    //     if( this->m_ElapsedIterations == 2 )
    //       {
    //       this->m_CurrentGradientStep = this->m_CurrentGradientStep * 1.0 / maxNorm;
    //       velocityField->FillBuffer( zeroVector );
    //       }

    RealImagePointer smoothHitImage;
    RealImagePointer smoothTotalImage;
    if (this->m_SmoothingVariance > NumericTraits<RealType>::ZeroValue())
    {
      smoothHitImage = this->SmoothImage(hitImage, this->m_SmoothingVariance);
      smoothTotalImage = this->SmoothImage(totalImage, this->m_SmoothingVariance);
    }
    else
    {
      smoothHitImage = hitImage;
      smoothTotalImage = totalImage;
    }

    ImageRegionConstIterator<RealImageType> ItSmoothHitImage(smoothHitImage, smoothHitImage->GetRequestedRegion());
    ImageRegionConstIterator<RealImageType> ItSmoothTotalImage(smoothTotalImage,
                                                               smoothTotalImage->GetRequestedRegion());

    ItCorticalThicknessImage.GoToBegin();
    ItForwardIncrementalField.GoToBegin();
    ItSegmentationImage.GoToBegin();
    ItSmoothHitImage.GoToBegin();
    ItSmoothTotalImage.GoToBegin();
    ItVelocityField.GoToBegin();

    while (!ItSegmentationImage.IsAtEnd())
    {
      ItVelocityField.Set(ItVelocityField.Get() + ItForwardIncrementalField.Get());
      const typename InputImageType::PixelType grayMatterPixel =
        static_cast<typename InputImageType::PixelType>(this->m_GrayMatterLabel);
      if (ItSegmentationImage.Get() == grayMatterPixel)
      {
        RealType thicknessValue = NumericTraits<RealType>::ZeroValue();
        if (ItSmoothHitImage.Get() > static_cast<RealType>(0.001))
        {
          thicknessValue = ItSmoothTotalImage.Get() / ItSmoothHitImage.Get();
          if (thicknessValue < NumericTraits<RealType>::ZeroValue())
          {
            thicknessValue = NumericTraits<RealType>::ZeroValue();
          }
          if (!this->m_ThicknessPriorImage && (thicknessValue > this->m_ThicknessPriorEstimate))
          {
            RealType fraction = this->m_ThicknessPriorEstimate / thicknessValue;
            ItVelocityField.Set(ItVelocityField.Get() * itk::Math::sqr(fraction));
          }
          else if (this->m_ThicknessPriorImage)
          {
            typename RealImageType::IndexType index = ItSegmentationImage.GetIndex();
            RealType                          thicknessPrior = this->m_ThicknessPriorImage->GetPixel(index);
            if ((thicknessPrior > NumericTraits<RealType>::ZeroValue()) && (thicknessValue > thicknessPrior))
            {
              priorEnergy += itk::Math::abs(thicknessPrior - thicknessValue);
              priorEnergyCount++;

              RealType fraction = thicknessPrior / thicknessValue;
              ItVelocityField.Set(ItVelocityField.Get() * itk::Math::sqr(fraction));
            }
          }
        }
        ItCorticalThicknessImage.Set(thicknessValue);
      }

      ++ItCorticalThicknessImage;
      ++ItForwardIncrementalField;
      ++ItSmoothHitImage;
      ++ItSegmentationImage;
      ++ItSmoothTotalImage;
      ++ItVelocityField;
    }

    if (this->m_UseBSplineSmoothing)
    {
      velocityField = this->BSplineSmoothDisplacementField(velocityField, this->m_BSplineSmoothingIsotropicMeshSpacing);
    }
    else if (this->m_UseMaskedSmoothing)
    {
      using MaskedSmootherType = MaskedSmoothingImageFilter<DisplacementFieldType, InputImageType, DisplacementFieldType>;
      typename MaskedSmootherType::Pointer maskedSmoother = MaskedSmootherType::New();
      maskedSmoother->SetInput(velocityField);
      maskedSmoother->SetMaskImage(thresholdedRegion);
      maskedSmoother->SetSmoothingVariance(this->m_SmoothingVariance);
      maskedSmoother->SetSparseImageNeighborhoodRadius(this->m_SparseImageNeighborhoodRadius);
      if (this->m_TimePoints.size() > 0)
      {
        maskedSmoother->SetTimeSmoothingVariance(this->m_TimeSmoothingVariance);
        maskedSmoother->SetTimePoints(this->m_TimePoints);
      }
      maskedSmoother->Update();

      velocityField = maskedSmoother->GetOutput();
    }
    else
    {
      velocityField = this->GaussianSmoothDisplacementField(velocityField, this->m_SmoothingVelocityFieldVariance);
    }

    // Calculate current energy and current convergence measurement

    currentEnergy /= numberOfGrayMatterVoxels;
    priorEnergy /= priorEnergyCount;

    if (this->m_ThicknessPriorImage)
    {
      itkDebugMacro("   PriorEnergy = " << priorEnergy);
    }
    this->m_CurrentEnergy = currentEnergy;

    convergenceMonitoring->AddEnergyValue(this->m_CurrentEnergy);
    this->m_CurrentConvergenceMeasurement = convergenceMonitoring->GetConvergenceValue();

    if (this->m_CurrentConvergenceMeasurement < this->m_ConvergenceThreshold)
    {
      isConverged = true;
    }
    reporter.CompletedStep();
  }

  // Replace the identity direction with the original direction in the outputs

  corticalThicknessImage->SetDirection(this->GetSegmentationImage()->GetDirection());

  if (this->m_IncludeCumulativeVelocityFields)
  {
    typename CumulativeVelocityFieldType::DirectionType direction;
    direction.SetIdentity();

    for (unsigned int d = 0; d < ImageDimension; d++)
    {
      for (unsigned int e = 0; e < ImageDimension; e++) 
      { 
        direction(d, e) = this->GetSegmentationImage()->GetDirection()(d, e);
      }
    }
   
    this->GetForwardCumulativeVelocityField()->SetDirection(direction);
    this->GetInverseCumulativeVelocityField()->SetDirection(direction);
  }  
}

template <typename TInputImage, typename TOutputImage>
typename DiReCTImageFilter<TInputImage, TOutputImage>::InputImagePointer
DiReCTImageFilter<TInputImage, TOutputImage>::ExtractRegion(
  const InputImageType *                                           segmentationImage,
  typename DiReCTImageFilter<TInputImage, TOutputImage>::LabelType whichRegion)
{
  using ThresholderType = BinaryThresholdImageFilter<InputImageType, InputImageType>;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput(segmentationImage);
  thresholder->SetLowerThreshold(whichRegion);
  thresholder->SetUpperThreshold(whichRegion);
  thresholder->SetInsideValue(1);
  thresholder->SetOutsideValue(0);
  thresholder->Update();

  InputImagePointer thresholdRegion = thresholder->GetOutput();
  thresholdRegion->Update();
  thresholdRegion->DisconnectPipeline();

  return thresholdRegion;
}

template <typename TInputImage, typename TOutputImage>
typename DiReCTImageFilter<TInputImage, TOutputImage>::InputImagePointer
DiReCTImageFilter<TInputImage, TOutputImage>::ExtractRegionalContours(
  const InputImageType *                                           segmentationImage,
  typename DiReCTImageFilter<TInputImage, TOutputImage>::LabelType whichRegion)
{
  InputImagePointer thresholdedRegion = this->ExtractRegion(segmentationImage, whichRegion);

  using ContourFilterType = BinaryContourImageFilter<InputImageType, InputImageType>;
  typename ContourFilterType::Pointer contourFilter = ContourFilterType::New();
  contourFilter->SetInput(thresholdedRegion);
  contourFilter->SetFullyConnected(true);
  contourFilter->SetBackgroundValue(0);
  contourFilter->SetForegroundValue(1);

  InputImagePointer contours = contourFilter->GetOutput();
  contours->Update();
  contours->DisconnectPipeline();
  contours->SetRegions(segmentationImage->GetRequestedRegion());

  return contours;
}


template <typename TInputImage, typename TOutputImage>
void
DiReCTImageFilter<TInputImage, TOutputImage>::MakeThicknessImage(RealImagePointer  hitImage,
                                                                 RealImagePointer  totalImage,
                                                                 InputImagePointer segmentationImage,
                                                                 RealImagePointer  corticalThicknessImage)
{
  RealImagePointer smoothHitImage;
  RealImagePointer smoothTotalImage;
  if (this->m_SmoothingVariance > NumericTraits<RealType>::ZeroValue())
  {
    smoothHitImage = this->SmoothImage(hitImage, this->m_SmoothingVariance);
    smoothTotalImage = this->SmoothImage(totalImage, this->m_SmoothingVariance);
  }
  else
  {
    smoothHitImage = hitImage;
    smoothTotalImage = totalImage;
  }

  ImageRegionIterator<RealImageType>                ItCorticalThicknessImage(corticalThicknessImage,
                                                              corticalThicknessImage->GetRequestedRegion());
  ImageRegionConstIteratorWithIndex<InputImageType> ItSegmentationImage(segmentationImage,
                                                                        segmentationImage->GetRequestedRegion());

  ImageRegionConstIterator<RealImageType> ItSmoothHitImage(smoothHitImage, smoothHitImage->GetRequestedRegion());
  ImageRegionConstIterator<RealImageType> ItSmoothTotalImage(smoothTotalImage, smoothTotalImage->GetRequestedRegion());

  ItCorticalThicknessImage.GoToBegin();
  ItSegmentationImage.GoToBegin();
  ItSmoothHitImage.GoToBegin();
  ItSmoothTotalImage.GoToBegin();

  RealType      meanThickness = 0;
  while (!ItSegmentationImage.IsAtEnd())
  {
    const typename InputImageType::PixelType grayMatterPixel =
      static_cast<typename InputImageType::PixelType>(this->m_GrayMatterLabel);
    if (ItSegmentationImage.Get() == grayMatterPixel)
    {
      RealType thicknessValue = NumericTraits<RealType>::ZeroValue();
      if (ItSmoothHitImage.Get() > static_cast<RealType>(0.001))
      {
        thicknessValue = ItSmoothTotalImage.Get() / ItSmoothHitImage.Get();
        meanThickness += thicknessValue;
        if (thicknessValue < NumericTraits<RealType>::ZeroValue())
        {
          thicknessValue = NumericTraits<RealType>::ZeroValue();
        }
      }
      ItCorticalThicknessImage.Set(thicknessValue);
    }
    ++ItCorticalThicknessImage;
    ++ItSmoothHitImage;
    ++ItSegmentationImage;
    ++ItSmoothTotalImage;
  }
}

template <typename TInputImage, typename TOutputImage>
typename DiReCTImageFilter<TInputImage, TOutputImage>::RealImagePointer
DiReCTImageFilter<TInputImage, TOutputImage>::WarpImage(const RealImageType *         inputImage,
                                                        const DisplacementFieldType * displacementField)
{
  using InterpolatorType = itk::LinearInterpolateImageFunction<RealImageType, double>;
  auto interpolator = InterpolatorType::New();

  using WarperType = WarpImageFilter<RealImageType, RealImageType, DisplacementFieldType>;
  typename WarperType::Pointer warper = WarperType::New();
  warper->SetInput(inputImage);
  warper->SetInterpolator(interpolator);
  warper->SetEdgePaddingValue(0);
  warper->SetOutputSpacing(displacementField->GetSpacing());
  warper->SetOutputOrigin(displacementField->GetOrigin());
  warper->SetOutputDirection(displacementField->GetDirection());
  warper->SetDisplacementField(displacementField);

  RealImagePointer warpedImage = warper->GetOutput();
  warpedImage->Update();
  warpedImage->DisconnectPipeline();

  return warpedImage;
}

template <typename TInputImage, typename TOutputImage>
typename DiReCTImageFilter<TInputImage, TOutputImage>::DisplacementFieldPointer
DiReCTImageFilter<TInputImage, TOutputImage>::GaussianSmoothDisplacementField(const DisplacementFieldType * inputField,
                                                                              const RealType                variance)
{
  using DuplicatorType = ImageDuplicator<DisplacementFieldType>;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(inputField);
  duplicator->Update();
  DisplacementFieldPointer outputField = duplicator->GetOutput();

  using SmootherType = VectorNeighborhoodOperatorImageFilter<DisplacementFieldType, DisplacementFieldType>;
  typename SmootherType::Pointer smoother = SmootherType::New();

  using GaussianType = GaussianOperator<VectorValueType, ImageDimension>;
  GaussianType gaussian;
  gaussian.SetVariance(variance);
  gaussian.SetMaximumError(0.001);
  for (unsigned int d = 0; d < ImageDimension; d++)
  {
    gaussian.SetDirection(d);
    gaussian.SetMaximumKernelWidth(outputField->GetRequestedRegion().GetSize()[d]);
    gaussian.CreateDirectional();

    smoother->SetOperator(gaussian);
    smoother->SetInput(outputField);

    outputField = smoother->GetOutput();
    outputField->Update();
    outputField->DisconnectPipeline();
  }

  // Ensure zero motion on the boundary

  RealType weight1 = NumericTraits<RealType>::OneValue();
  if (variance < static_cast<RealType>(0.5))
  {
    weight1 = NumericTraits<RealType>::OneValue() -
              NumericTraits<RealType>::OneValue() * (variance / static_cast<RealType>(0.5));
  }
  RealType weight2 = NumericTraits<RealType>::OneValue() - weight1;

  using MultiplierType = MultiplyByConstantImageFilter<DisplacementFieldType, RealType, DisplacementFieldType>;

  typename MultiplierType::Pointer multiplier1 = MultiplierType::New();
  multiplier1->SetConstant2(weight1);
  multiplier1->SetInput1(outputField);

  typename MultiplierType::Pointer multiplier2 = MultiplierType::New();
  multiplier2->SetConstant2(weight2);
  multiplier2->SetInput1(inputField);

  using AdderType = AddImageFilter<DisplacementFieldType, DisplacementFieldType, DisplacementFieldType>;
  typename AdderType::Pointer adder = AdderType::New();
  adder->SetInput1(multiplier1->GetOutput());
  adder->SetInput2(multiplier2->GetOutput());

  outputField = adder->GetOutput();
  outputField->Update();
  outputField->DisconnectPipeline();

  VectorType zeroVector(0.0);

  ImageLinearIteratorWithIndex<DisplacementFieldType> It(outputField, outputField->GetRequestedRegion());
  for (unsigned int d = 0; d < ImageDimension; d++)
  {
    It.SetDirection(d);
    It.GoToBegin();
    while (!It.IsAtEnd())
    {
      It.GoToBeginOfLine();
      It.Set(zeroVector);
      It.GoToEndOfLine();
      --It;
      It.Set(zeroVector);

      It.NextLine();
    }
  }

  return outputField;
}

template <typename TInputImage, typename TOutputImage>
typename DiReCTImageFilter<TInputImage, TOutputImage>::DisplacementFieldPointer
DiReCTImageFilter<TInputImage, TOutputImage>::BSplineSmoothDisplacementField(const DisplacementFieldType * inputField,
                                                                             const RealType isotropicMeshSpacing)
{
  using BSplineFilterType = DisplacementFieldToBSplineImageFilter<DisplacementFieldType>;

  // calculate the number of control points based on the isotropic mesh spacing

  typename BSplineFilterType::ArrayType ncps;

  for (unsigned int d = 0; d < ImageDimension; d++)
  {
    RealType domain = static_cast<RealType>(inputField->GetLargestPossibleRegion().GetSize()[d] - 1) *
                      static_cast<RealType>(inputField->GetSpacing()[d]);
    ncps[d] = static_cast<unsigned int>(std::ceil(domain / isotropicMeshSpacing));
  }

  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
  bspliner->SetDisplacementField(inputField);
  bspliner->SetNumberOfControlPoints(ncps);
  bspliner->SetSplineOrder(3);
  bspliner->SetNumberOfFittingLevels(1);
  bspliner->SetEnforceStationaryBoundary(true);
  bspliner->SetEstimateInverse(false);
  bspliner->Update();

  DisplacementFieldPointer outputField = bspliner->GetOutput();
  outputField->DisconnectPipeline();

  return outputField;
}

template <typename TInputImage, typename TOutputImage>
typename DiReCTImageFilter<TInputImage, TOutputImage>::RealImagePointer
DiReCTImageFilter<TInputImage, TOutputImage>::SmoothImage(const RealImageType * inputImage, const RealType variance)
{
  using SmootherType = DiscreteGaussianImageFilter<RealImageType, RealImageType>;
  typename SmootherType::Pointer smoother = SmootherType::New();
  smoother->SetVariance(variance);
  smoother->SetUseImageSpacing(false);
  smoother->SetMaximumError(0.01);
  smoother->SetInput(inputImage);

  typename RealImageType::Pointer smoothImage = smoother->GetOutput();
  smoothImage->Update();
  smoothImage->DisconnectPipeline();

  return smoothImage;
}

/**
 * Standard "PrintSelf" method
 */
template <typename TInputImage, typename TOutputImage>
void
DiReCTImageFilter<TInputImage, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Gray matter label = " << this->m_GrayMatterLabel << std::endl;
  os << indent << "White matter label = " << this->m_WhiteMatterLabel << std::endl;
  os << indent << "Maximum number of iterations = " << this->m_MaximumNumberOfIterations << std::endl;
  os << indent << "Thickness prior estimate = " << this->m_ThicknessPriorEstimate << std::endl;
  os << indent << "Smoothing variance = " << this->m_SmoothingVariance << std::endl;
  if (this->m_UseBSplineSmoothing)
  {
    os << indent << "B-spline smoothing isotropic mesh spacing = " << this->m_BSplineSmoothingIsotropicMeshSpacing
       << std::endl;
  }
  else
  {
    os << indent << "Smoothing velocity field variance = " << this->m_SmoothingVelocityFieldVariance << std::endl;
  }
  os << indent << "Use masked smoothing = " << static_cast<int>(this->m_UseMaskedSmoothing) << std::endl;
  os << indent << "Number of integration points = " << this->m_NumberOfIntegrationPoints << std::endl;
  os << indent << "Maximum number of invert displacement field iterations = "
     << this->m_MaximumNumberOfInvertDisplacementFieldIterations << std::endl;
  os << indent << "Initial gradient step = " << this->m_InitialGradientStep << std::endl;
  os << indent << "Current gradient step = " << this->m_CurrentGradientStep << std::endl;
  os << indent << "Convergence threshold = " << this->m_ConvergenceThreshold << std::endl;
  os << indent << "Convergence window size = " << this->m_ConvergenceWindowSize << std::endl;

  os << indent << "Time smoothing variance = " << this->m_TimeSmoothingVariance << std::endl;
}
} // end namespace itk

#endif
