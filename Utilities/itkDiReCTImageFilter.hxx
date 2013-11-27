/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkLabelOverlapMeasuresImageFilter.hxx,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDiReCTImageFilter_hxx
#define __itkDiReCTImageFilter_hxx

#include "itkDiReCTImageFilter.h"

#include "itkAddImageFilter.h"
// #include "itkAndImageFilter.h"
#include "itkBinaryContourImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkComposeDisplacementFieldsImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkGaussianOperator.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImportImageFilter.h"
#include "itkInvertDisplacementFieldImageFilter.h"
#include "itkIterationReporter.h"
#include "itkMaximumImageFilter.h"
#include "itkMultiplyByConstantImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkVectorNeighborhoodOperatorImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkWindowConvergenceMonitoringFunction.h"


#include "itkImageFileWriter.h"
#include "ReadWriteImage.h"

namespace itk
{
template <class TInputImage, class TOutputImage>
DiReCTImageFilter<TInputImage, TOutputImage>
::DiReCTImageFilter() :
  m_ThicknessPriorEstimate( 10.0 ),
  m_SmoothingVariance( 1.0 ),
  m_SmoothingVelocityFieldVariance( 1.5 ),
  m_InitialGradientStep( 0.025 ),
  m_CurrentGradientStep( 0.025 ),
  m_NumberOfIntegrationPoints( 10 ),
  m_GrayMatterLabel( 2 ),
  m_WhiteMatterLabel( 3 ),
  m_MaximumNumberOfIterations( 50 ),
  m_MaximumNumberOfInvertDisplacementFieldIterations( 20 ),
  m_CurrentEnergy( NumericTraits<RealType>::max() ),
  m_ConvergenceThreshold( 0.001 ),
  m_ConvergenceWindowSize( 10 )
{
  this->m_ThicknessPriorImage = NULL;
  this->SetNumberOfRequiredInputs( 3 );
}

template <class TInputImage, class TOutputImage>
DiReCTImageFilter<TInputImage, TOutputImage>
::~DiReCTImageFilter()
{
}

template <class TInputImage, class TOutputImage>
void
DiReCTImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  if ( this->m_ThicknessPriorImage )
    {
    std::cout << "Using prior thickness image." << std::endl;
    }
  this->m_CurrentGradientStep = this->m_InitialGradientStep;

  // Convert all input direction matrices to identities saving the original
  // directions to put back at the end of filter processing. We do this
  // because the white and gray matters reside in the same space and the
  // assumption simplifies the underlying registration code.

  const bool filterHandlesMemory = false;

  typename InputImageType::DirectionType identity;
  identity.SetIdentity();

  typedef ImportImageFilter<InputPixelType, ImageDimension> SegmentationImageImporterType;
  typename SegmentationImageImporterType::Pointer segmentationImageImporter = SegmentationImageImporterType::New();
  segmentationImageImporter->SetImportPointer( const_cast<InputPixelType *>(
                                                 this->GetSegmentationImage()->GetBufferPointer() ),
                                               ( this->GetSegmentationImage()->GetBufferedRegion() ).GetNumberOfPixels(),
                                               filterHandlesMemory );
  segmentationImageImporter->SetRegion( this->GetSegmentationImage()->GetBufferedRegion() );
  segmentationImageImporter->SetOrigin( this->GetSegmentationImage()->GetOrigin() );
  segmentationImageImporter->SetSpacing( this->GetSegmentationImage()->GetSpacing() );
  segmentationImageImporter->SetDirection( identity );

  InputImagePointer segmentationImage = segmentationImageImporter->GetOutput();
  segmentationImage->Update();
  segmentationImage->DisconnectPipeline();

  typedef ImportImageFilter<RealType, ImageDimension> ProbablilityImageImporterType;

  typename ProbablilityImageImporterType::Pointer grayMatterProbabilityImageImporter = ProbablilityImageImporterType::New();
  grayMatterProbabilityImageImporter->SetImportPointer( const_cast<RealType *>(
                                                          this->GetGrayMatterProbabilityImage()->GetBufferPointer() ),
                                                        ( this->GetGrayMatterProbabilityImage()->GetBufferedRegion() ).
                                                        GetNumberOfPixels(), filterHandlesMemory );
  grayMatterProbabilityImageImporter->SetRegion( this->GetGrayMatterProbabilityImage()->GetBufferedRegion() );
  grayMatterProbabilityImageImporter->SetOrigin( this->GetGrayMatterProbabilityImage()->GetOrigin() );
  grayMatterProbabilityImageImporter->SetSpacing( this->GetGrayMatterProbabilityImage()->GetSpacing() );
  grayMatterProbabilityImageImporter->SetDirection( identity );

  RealImagePointer grayMatterProbabilityImage = grayMatterProbabilityImageImporter->GetOutput();
  grayMatterProbabilityImage->Update();
  grayMatterProbabilityImage->DisconnectPipeline();

  typename ProbablilityImageImporterType::Pointer whiteMatterProbabilityImageImporter = ProbablilityImageImporterType::New();
  whiteMatterProbabilityImageImporter->SetImportPointer( const_cast<RealType *>(
                                                           this->GetWhiteMatterProbabilityImage()->GetBufferPointer() ),
                                                         ( this->GetWhiteMatterProbabilityImage()->GetBufferedRegion() )
                                                         .GetNumberOfPixels(), filterHandlesMemory );
  whiteMatterProbabilityImageImporter->SetRegion( this->GetWhiteMatterProbabilityImage()->GetBufferedRegion() );
  whiteMatterProbabilityImageImporter->SetOrigin( this->GetWhiteMatterProbabilityImage()->GetOrigin() );
  whiteMatterProbabilityImageImporter->SetSpacing( this->GetWhiteMatterProbabilityImage()->GetSpacing() );
  whiteMatterProbabilityImageImporter->SetDirection( identity );

  RealImagePointer whiteMatterProbabilityImage = whiteMatterProbabilityImageImporter->GetOutput();
  whiteMatterProbabilityImage->Update();
  whiteMatterProbabilityImage->DisconnectPipeline();

  // Extract the gray and white matter segmentations and combine to form the
  // gm/wm region.  Dilate the latter region by 1 voxel.

  InputImagePointer grayMatter = this->ExtractRegion( segmentationImage, this->m_GrayMatterLabel );
  InputImagePointer whiteMatter = this->ExtractRegion( segmentationImage, this->m_WhiteMatterLabel );

  typedef AddImageFilter<InputImageType, InputImageType, InputImageType> AdderType;
  typename AdderType::Pointer adder = AdderType::New();
  adder->SetInput1( grayMatter );
  adder->SetInput2( whiteMatter );
  adder->Update();

  InputImagePointer thresholdedRegion = this->ExtractRegion( const_cast<const InputImageType *>( adder->GetOutput() ), 1 );

  // Extract the white and gm/wm matter contours

  InputImagePointer matterContours = this->ExtractRegionalContours( thresholdedRegion, 1 );
  InputImagePointer whiteMatterContoursTmp = this->ExtractRegionalContours( segmentationImage, this->m_WhiteMatterLabel );

  typedef CastImageFilter<InputImageType, RealImageType> CasterType;
  typename CasterType::Pointer caster = CasterType::New();
  caster->SetInput( whiteMatterContoursTmp );
  caster->Update();
  RealImagePointer whiteMatterContours = caster->GetOutput();

  // Initialize fields and images.

  VectorType zeroVector( 0.0 );

  RealImagePointer corticalThicknessImage = RealImageType::New();
  corticalThicknessImage->CopyInformation( segmentationImage );
  corticalThicknessImage->SetRegions( segmentationImage->GetRequestedRegion() );
  corticalThicknessImage->Allocate();
  corticalThicknessImage->FillBuffer( 0.0 );

  DisplacementFieldPointer forwardIncrementalField = DisplacementFieldType::New();
  forwardIncrementalField->CopyInformation( segmentationImage );
  forwardIncrementalField->SetRegions( segmentationImage->GetRequestedRegion() );
  forwardIncrementalField->Allocate();

  RealImagePointer hitImage = RealImageType::New();
  hitImage->CopyInformation( segmentationImage );
  hitImage->SetRegions( segmentationImage->GetRequestedRegion() );
  hitImage->Allocate();

  DisplacementFieldPointer integratedField = DisplacementFieldType::New();
  integratedField->CopyInformation( segmentationImage );
  integratedField->SetRegions( segmentationImage->GetRequestedRegion() );
  integratedField->Allocate();
  integratedField->FillBuffer( zeroVector );

  DisplacementFieldPointer inverseField = DisplacementFieldType::New();
  inverseField->CopyInformation( segmentationImage );
  inverseField->SetRegions( segmentationImage->GetRequestedRegion() );
  inverseField->Allocate();

  DisplacementFieldPointer inverseIncrementalField = DisplacementFieldType::New();
  inverseIncrementalField->CopyInformation( segmentationImage );
  inverseIncrementalField->SetRegions( segmentationImage->GetRequestedRegion() );
  inverseIncrementalField->Allocate();

  RealImagePointer speedImage = RealImageType::New();
  speedImage->CopyInformation( segmentationImage );
  speedImage->SetRegions( segmentationImage->GetRequestedRegion() );
  speedImage->Allocate();

  RealImagePointer thicknessImage = RealImageType::New();
  thicknessImage->CopyInformation( segmentationImage );
  thicknessImage->SetRegions( segmentationImage->GetRequestedRegion() );
  thicknessImage->Allocate();

  RealImagePointer totalImage = RealImageType::New();
  totalImage->CopyInformation( segmentationImage );
  totalImage->SetRegions( segmentationImage->GetRequestedRegion() );
  totalImage->Allocate();

  DisplacementFieldPointer velocityField = DisplacementFieldType::New();
  velocityField->CopyInformation( segmentationImage );
  velocityField->SetRegions( segmentationImage->GetRequestedRegion() );
  velocityField->Allocate();
  velocityField->FillBuffer( zeroVector );

  // Instantiate most of the iterators all in one place

  ImageRegionIterator<RealImageType> ItCorticalThicknessImage(
    corticalThicknessImage,
    corticalThicknessImage->GetRequestedRegion() );
  ImageRegionConstIterator<RealImageType> ItGrayMatterProbabilityMap(
    grayMatterProbabilityImage,
    grayMatterProbabilityImage->GetRequestedRegion() );
  ImageRegionIterator<RealImageType> ItHitImage(
    hitImage,
    hitImage->GetRequestedRegion() );
  ImageRegionIterator<DisplacementFieldType> ItForwardIncrementalField(
    forwardIncrementalField,
    forwardIncrementalField->GetRequestedRegion() );
  ImageRegionConstIterator<InputImageType> ItMatterContours(
    matterContours,
    matterContours->GetRequestedRegion() );
  ImageRegionIterator<DisplacementFieldType> ItInverseIncrementalField(
    inverseIncrementalField,
    inverseIncrementalField->GetRequestedRegion() );
  ImageRegionConstIteratorWithIndex<InputImageType> ItSegmentationImage(
    segmentationImage,
    segmentationImage->GetRequestedRegion() );
  ImageRegionIterator<RealImageType> ItSpeedImage(
    speedImage,
    speedImage->GetRequestedRegion() );
  ImageRegionIterator<RealImageType> ItThicknessImage(
    thicknessImage,
    thicknessImage->GetRequestedRegion() );
  ImageRegionIterator<RealImageType> ItTotalImage(
    totalImage,
    totalImage->GetRequestedRegion() );
  ImageRegionConstIterator<RealImageType> ItWhiteMatterContours(
    whiteMatterContours,
    whiteMatterContours->GetRequestedRegion() );

  // Monitor the convergence
  typedef Function::WindowConvergenceMonitoringFunction<double> ConvergenceMonitoringType;
  ConvergenceMonitoringType::Pointer convergenceMonitoring = ConvergenceMonitoringType::New();
  convergenceMonitoring->SetWindowSize( this->m_ConvergenceWindowSize );

  // Instantiate the progress reporter

  IterationReporter reporter( this, 0, 1 );

  bool isConverged = false;
  this->m_CurrentConvergenceMeasurement = NumericTraits<RealType>::max();
  this->m_ElapsedIterations = 0;
  while( this->m_ElapsedIterations++ < this->m_MaximumNumberOfIterations && isConverged == false )
    {
    this->MakeThicknessImage( hitImage, totalImage, segmentationImage, thicknessImage );

    RealType priorEnergy = 0;
    unsigned long priorEnergyCount = 0;

    RealType currentEnergy = 0.0;
    RealType numberOfGrayMatterVoxels = 0.0;

    forwardIncrementalField->FillBuffer( zeroVector );
    inverseField->FillBuffer( zeroVector );
    inverseIncrementalField->FillBuffer( zeroVector );

    hitImage->FillBuffer( 0.0 );
    totalImage->FillBuffer( 0.0 );

    ImageRegionIterator<DisplacementFieldType> ItVelocityField(
      velocityField,
      velocityField->GetRequestedRegion() );

    unsigned int integrationPoint = 0;
    while( integrationPoint++ < this->m_NumberOfIntegrationPoints )
      {
      typedef ComposeDisplacementFieldsImageFilter<DisplacementFieldType> ComposerType;
      typename ComposerType::Pointer composer = ComposerType::New();
      composer->SetDisplacementField( inverseIncrementalField );
      composer->SetWarpingField( inverseField );

      inverseField = composer->GetOutput();
      inverseField->Update();
      inverseField->DisconnectPipeline();

      RealImagePointer warpedWhiteMatterProbabilityImage = this->WarpImage( whiteMatterProbabilityImage, inverseField );
      RealImagePointer warpedWhiteMatterContours = this->WarpImage( whiteMatterContours, inverseField );
      RealImagePointer warpedThicknessImage = this->WarpImage( thicknessImage, inverseField );

      typedef GradientRecursiveGaussianImageFilter<RealImageType, DisplacementFieldType>
        GradientImageFilterType;
      typename GradientImageFilterType::Pointer gradientFilter =
        GradientImageFilterType::New();
      gradientFilter->SetInput( warpedWhiteMatterProbabilityImage );
      gradientFilter->SetSigma( this->m_SmoothingVariance );
      gradientFilter->Update();

      DisplacementFieldPointer gradientImage = gradientFilter->GetOutput();

      // Instantiate the iterators all in one place

      ImageRegionIterator<DisplacementFieldType> ItGradientImage(
        gradientImage,
        gradientImage->GetRequestedRegion() );
      ImageRegionIterator<RealImageType> ItWarpedWhiteMatterProbabilityMap(
        warpedWhiteMatterProbabilityImage,
        warpedWhiteMatterProbabilityImage->GetRequestedRegion() );
      ImageRegionIterator<RealImageType> ItWarpedWhiteMatterContours(
        warpedWhiteMatterContours,
        warpedWhiteMatterContours->GetRequestedRegion() );
      ImageRegionIterator<DisplacementFieldType> ItInverseField(
        inverseField,
        inverseField->GetRequestedRegion() );
      ImageRegionIterator<DisplacementFieldType> ItIntegratedField(
        integratedField,
        integratedField->GetRequestedRegion() );

      // Generate speed image

      speedImage->FillBuffer( 0.0 );

      ItGradientImage.GoToBegin();
      ItGrayMatterProbabilityMap.GoToBegin();
      ItSegmentationImage.GoToBegin();
      ItSpeedImage.GoToBegin();
      ItWarpedWhiteMatterProbabilityMap.GoToBegin();

      const typename InputImageType::PixelType grayMatterPixel =
        static_cast<typename InputImageType::PixelType>( this->m_GrayMatterLabel ),
      whiteMatterPixel =
        static_cast<typename InputImageType::PixelType>( this->m_WhiteMatterLabel );
      while( !ItSegmentationImage.IsAtEnd() )
        {
        if( ItSegmentationImage.Get() == grayMatterPixel )
          {
          RealType norm = ( ItGradientImage.Get() ).GetNorm();
          if( norm > 1e-3 && !vnl_math_isnan( norm ) && !vnl_math_isinf( norm ) )
            {
	           if( ! this->m_ThicknessPriorImage  )
	             {
	             ItGradientImage.Set( ItGradientImage.Get() / norm );
	             }
            }
          else
            {
            ItGradientImage.Set( zeroVector );
            }
          RealType delta = ( ItWarpedWhiteMatterProbabilityMap.Get() - ItGrayMatterProbabilityMap.Get() );
          currentEnergy += vnl_math_abs( delta );

										if( this->m_ThicknessPriorImage )
												{
												typename RealImageType::IndexType index = ItGrayMatterProbabilityMap.GetIndex();
												RealType thicknessPrior = this->m_ThicknessPriorImage->GetPixel( index );
												RealType thicknessValue = thicknessImage->GetPixel( index );
												if( thicknessValue > thicknessPrior )
												  {
												  priorEnergy += vnl_math_abs( thicknessPrior - thicknessValue );
												  priorEnergyCount++;
												  }
												if( ( thicknessPrior > NumericTraits<RealType>::min() ) &&
										      ( thicknessValue > thicknessPrior ) )
														{
														RealType downScale = vnl_math_abs( thicknessPrior / thicknessValue ) * this->m_ThicknessPriorEstimate;
														if( integrationPoint > 0 )
														  {
																velocityField->SetPixel( index , velocityField->GetPixel( index ) * downScale );
																}
														}
												}
          numberOfGrayMatterVoxels++;
          RealType speedValue = -1.0 * delta * ItGrayMatterProbabilityMap.Get() * this->m_CurrentGradientStep;
          if( vnl_math_isnan( speedValue ) || vnl_math_isinf( speedValue ) )
            {
            speedValue = 0.0;
            }
          ItSpeedImage.Set( speedValue );
          }
        else
          {
          ItSpeedImage.Set( 0.0 );
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

      while( !ItSegmentationImage.IsAtEnd() )
        {
        typename InputImageType::IndexType index =
          ItSegmentationImage.GetIndex();
        typename InputImageType::PixelType segmentationValue =
          ItSegmentationImage.Get();

        ItForwardIncrementalField.Set( ItForwardIncrementalField.Get()
                                       + ItGradientImage.Get() * ItSpeedImage.Get() );
        if( segmentationValue == grayMatterPixel || segmentationValue == whiteMatterPixel )
          {
          if( integrationPoint == 1 )
            {
            typename InputImageType::PixelType whiteMatterContoursValue =
              static_cast<typename InputImageType::PixelType>( ItWhiteMatterContours.Get() );
            hitImage->SetPixel( index, whiteMatterContoursValue );

            VectorType vector = integratedField->GetPixel( index );
            RealType   weightedNorm = vector.GetNorm() * whiteMatterContoursValue;

            thicknessImage->SetPixel( index, weightedNorm );
            totalImage->SetPixel( index, weightedNorm );
            }
          else if( segmentationValue == grayMatterPixel )
            {
            hitImage->SetPixel( index, hitImage->GetPixel( index )
                                + warpedWhiteMatterContours->GetPixel( index ) );
            totalImage->SetPixel( index, totalImage->GetPixel( index )
                                  + warpedThicknessImage->GetPixel( index ) );
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
      while( !ItSegmentationImage.IsAtEnd() )
        {
        typename InputImageType::PixelType segmentationValue =
          ItSegmentationImage.Get();
        typename InputImageType::PixelType whiteMatterContoursValue =
          static_cast<typename InputImageType::PixelType>( ItWhiteMatterContours.Get() );
        typename InputImageType::PixelType matterContoursValue = ItMatterContours.Get();

        if( segmentationValue == 0 ||
          ( whiteMatterContoursValue == 0 && matterContoursValue == 0 && segmentationValue != this->m_GrayMatterLabel ) )
          {
          ItInverseField.Set( zeroVector );
          ItVelocityField.Set( zeroVector );
          ItIntegratedField.Set( zeroVector );
          }

        ItInverseIncrementalField.Set( ItVelocityField.Get() );

        ++ItSegmentationImage;
        ++ItMatterContours;
        ++ItWhiteMatterContours;
        ++ItVelocityField;
        ++ItInverseIncrementalField;
        ++ItInverseField;
        ++ItIntegratedField;
        }

      if( integrationPoint == 1 )
        {
        integratedField->FillBuffer( zeroVector );
        }

      typedef InvertDisplacementFieldImageFilter<DisplacementFieldType> InverterType;

      typename InverterType::Pointer inverter1 = InverterType::New();
      inverter1->SetInput( inverseField );
      inverter1->SetInverseFieldInitialEstimate( integratedField );
      inverter1->SetMaximumNumberOfIterations( this->m_MaximumNumberOfInvertDisplacementFieldIterations );
      inverter1->SetMeanErrorToleranceThreshold( 0.001 );
      inverter1->SetMaxErrorToleranceThreshold( 0.1 );
      inverter1->Update();

      integratedField = inverter1->GetOutput();
      integratedField->DisconnectPipeline();

      typename InverterType::Pointer inverter2 = InverterType::New();
      inverter2->SetInput( integratedField );
      inverter2->SetInverseFieldInitialEstimate( inverseField );
      inverter2->SetMaximumNumberOfIterations( this->m_MaximumNumberOfInvertDisplacementFieldIterations );
      inverter2->SetMeanErrorToleranceThreshold( 0.001 );
      inverter2->SetMaxErrorToleranceThreshold( 0.1 );
      inverter2->Update();

      inverseField = inverter2->GetOutput();
      inverseField->DisconnectPipeline();
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
    if( this->m_SmoothingVariance > 0.0 )
      {
      smoothHitImage = this->SmoothImage( hitImage, this->m_SmoothingVariance );
      smoothTotalImage = this->SmoothImage( totalImage, this->m_SmoothingVariance );
      }
    else
      {
      smoothHitImage = hitImage;
      smoothTotalImage = totalImage;
      }

    ImageRegionConstIterator<RealImageType> ItSmoothHitImage( smoothHitImage,
      smoothHitImage->GetRequestedRegion() );
    ImageRegionConstIterator<RealImageType> ItSmoothTotalImage( smoothTotalImage,
      smoothTotalImage->GetRequestedRegion() );

    ItCorticalThicknessImage.GoToBegin();
    ItForwardIncrementalField.GoToBegin();
    ItSegmentationImage.GoToBegin();
    ItSmoothHitImage.GoToBegin();
    ItSmoothTotalImage.GoToBegin();
    ItVelocityField.GoToBegin();

    while( !ItSegmentationImage.IsAtEnd() )
      {
      ItVelocityField.Set( ItVelocityField.Get() + ItForwardIncrementalField.Get() );
      const typename InputImageType::PixelType grayMatterPixel =
        static_cast<typename InputImageType::PixelType>( this->m_GrayMatterLabel );
      if(  ItSegmentationImage.Get() == grayMatterPixel )
        {
        RealType thicknessValue = 0.0;
        if( ItSmoothHitImage.Get() > 0.001 )
          {
          thicknessValue = ItSmoothTotalImage.Get() / ItSmoothHitImage.Get();
          if( thicknessValue < 0.0 )
            {
            thicknessValue = 0.0;
            }
//           if( thicknessValue > this->m_ThicknessPriorEstimate )
//             {
// //             thicknessValue = this->m_ThicknessPriorEstimate;
//             RealType fraction = this->m_ThicknessPriorEstimate / thicknessValue;
//             ItVelocityField.Set( ItVelocityField.Get() * vnl_math_sqr( fraction ) );
//             }
          }

        ItCorticalThicknessImage.Set( thicknessValue );
        }

      ++ItCorticalThicknessImage;
      ++ItForwardIncrementalField;
      ++ItSmoothHitImage;
      ++ItSegmentationImage;
      ++ItSmoothTotalImage;
      ++ItVelocityField;
      }

    velocityField = this->SmoothDisplacementField( velocityField,
                                                   this->m_SmoothingVelocityFieldVariance );

    // Calculate current energy and current convergence measurement

    currentEnergy /= numberOfGrayMatterVoxels;
    priorEnergy /= priorEnergyCount;

    if( this->m_ThicknessPriorImage )
      {
      itkDebugMacro( "   PriorEnergy = " << priorEnergy );
      }
    this->m_CurrentEnergy = currentEnergy;

    convergenceMonitoring->AddEnergyValue( this->m_CurrentEnergy );
    this->m_CurrentConvergenceMeasurement = convergenceMonitoring->GetConvergenceValue();

    if( this->m_CurrentConvergenceMeasurement < this->m_ConvergenceThreshold )
      {
      isConverged = true;
      }
    reporter.CompletedStep();
    }

  // Replace the identity direction with the original direction in the outputs

  RealImagePointer warpedWhiteMatterProbabilityImage = this->WarpImage( whiteMatterProbabilityImage, inverseField );
  warpedWhiteMatterProbabilityImage->SetDirection( this->GetSegmentationImage()->GetDirection() );
  corticalThicknessImage->SetDirection( this->GetSegmentationImage()->GetDirection() );

  this->SetNthOutput( 0, corticalThicknessImage );
  this->SetNthOutput( 1, warpedWhiteMatterProbabilityImage );
}

template <class TInputImage, class TOutputImage>
typename DiReCTImageFilter<TInputImage, TOutputImage>::InputImagePointer
DiReCTImageFilter<TInputImage, TOutputImage>
::ExtractRegion( const InputImageType *segmentationImage,
                 typename DiReCTImageFilter<TInputImage, TOutputImage>::LabelType whichRegion )
{
  typedef BinaryThresholdImageFilter<InputImageType, InputImageType>
    ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( segmentationImage );
  thresholder->SetLowerThreshold( whichRegion );
  thresholder->SetUpperThreshold( whichRegion );
  thresholder->SetInsideValue( 1 );
  thresholder->SetOutsideValue( 0 );
  thresholder->Update();

  InputImagePointer thresholdRegion = thresholder->GetOutput();
  thresholdRegion->Update();
  thresholdRegion->DisconnectPipeline();

  return thresholdRegion;
}

template <class TInputImage, class TOutputImage>
typename DiReCTImageFilter<TInputImage, TOutputImage>::InputImagePointer
DiReCTImageFilter<TInputImage, TOutputImage>
::ExtractRegionalContours( const InputImageType *segmentationImage,
                           typename DiReCTImageFilter<TInputImage, TOutputImage>::LabelType whichRegion )
{
  InputImagePointer thresholdedRegion = this->ExtractRegion(
      segmentationImage, whichRegion );

  typedef BinaryContourImageFilter<InputImageType, InputImageType>
    ContourFilterType;
  typename ContourFilterType::Pointer contourFilter = ContourFilterType::New();
  contourFilter->SetInput( thresholdedRegion );
  contourFilter->SetFullyConnected( true );
  contourFilter->SetBackgroundValue( 0 );
  contourFilter->SetForegroundValue( 1 );

  InputImagePointer contours = contourFilter->GetOutput();
  contours->Update();
  contours->DisconnectPipeline();
  contours->SetRegions( segmentationImage->GetRequestedRegion() );

  return contours;
}


template <class TInputImage, class TOutputImage>
void
DiReCTImageFilter<TInputImage, TOutputImage>
::MakeThicknessImage( RealImagePointer hitImage, RealImagePointer totalImage,
  InputImagePointer segmentationImage, RealImagePointer corticalThicknessImage )
{
		RealImagePointer smoothHitImage;
		RealImagePointer smoothTotalImage;
		if( this->m_SmoothingVariance > 0.0 )
				{
				smoothHitImage = this->SmoothImage( hitImage, this->m_SmoothingVariance );
				smoothTotalImage = this->SmoothImage( totalImage, this->m_SmoothingVariance );
				}
		else
				{
				smoothHitImage = hitImage;
				smoothTotalImage = totalImage;
				}

  ImageRegionIterator<RealImageType> ItCorticalThicknessImage(
    corticalThicknessImage,
    corticalThicknessImage->GetRequestedRegion() );
  ImageRegionConstIteratorWithIndex<InputImageType> ItSegmentationImage(
    segmentationImage,
    segmentationImage->GetRequestedRegion() );

		ImageRegionConstIterator<RealImageType> ItSmoothHitImage( smoothHitImage,
				smoothHitImage->GetRequestedRegion() );
		ImageRegionConstIterator<RealImageType> ItSmoothTotalImage( smoothTotalImage,
				smoothTotalImage->GetRequestedRegion() );

		ItCorticalThicknessImage.GoToBegin();
		ItSegmentationImage.GoToBegin();
		ItSmoothHitImage.GoToBegin();
		ItSmoothTotalImage.GoToBegin();

		RealType meanThickness = 0;
		unsigned long count = 0;
		while( !ItSegmentationImage.IsAtEnd() )
				{
				const typename InputImageType::PixelType grayMatterPixel =
						static_cast<typename InputImageType::PixelType>( this->m_GrayMatterLabel );
				if(  ItSegmentationImage.Get() == grayMatterPixel )
						{
						RealType thicknessValue = 0.0;
						if( ItSmoothHitImage.Get() > 0.001 )
								{
								thicknessValue = ItSmoothTotalImage.Get() / ItSmoothHitImage.Get();
       	meanThickness += thicknessValue;
	       count++;
								if( thicknessValue < 0.0 )
										{
										thicknessValue = 0.0;
										}
								}
						ItCorticalThicknessImage.Set( thicknessValue );
						}
				++ItCorticalThicknessImage;
				++ItSmoothHitImage;
				++ItSegmentationImage;
				++ItSmoothTotalImage;
				}
}

template <class TInputImage, class TOutputImage>
typename DiReCTImageFilter<TInputImage, TOutputImage>::RealImagePointer
DiReCTImageFilter<TInputImage, TOutputImage>
::WarpImage( const RealImageType *inputImage,
             const DisplacementFieldType *displacementField )
{
  typedef WarpImageFilter<RealImageType, RealImageType, DisplacementFieldType> WarperType;
  typename WarperType::Pointer warper = WarperType::New();
  warper->SetInput( inputImage );
  warper->SetDisplacementField( displacementField );
  warper->SetEdgePaddingValue( 0 );
  warper->SetOutputSpacing( inputImage->GetSpacing() );
  warper->SetOutputOrigin( inputImage->GetOrigin() );
  warper->SetOutputDirection( inputImage->GetDirection() );

  RealImagePointer warpedImage = warper->GetOutput();
  warpedImage->Update();
  warpedImage->DisconnectPipeline();

  return warpedImage;
}

template <class TInputImage, class TOutputImage>
typename DiReCTImageFilter<TInputImage, TOutputImage>::DisplacementFieldPointer
DiReCTImageFilter<TInputImage, TOutputImage>
::SmoothDisplacementField( const DisplacementFieldType *inputField,
                           const RealType variance )
{
  typedef ImageDuplicator<DisplacementFieldType> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage( inputField );
  duplicator->Update();
  DisplacementFieldPointer outputField = duplicator->GetModifiableOutput();

  typedef VectorNeighborhoodOperatorImageFilter<DisplacementFieldType,
                                                DisplacementFieldType> SmootherType;
  typename SmootherType::Pointer smoother = SmootherType::New();

  typedef GaussianOperator<VectorValueType, ImageDimension> GaussianType;
  GaussianType gaussian;
  gaussian.SetVariance( variance );
  gaussian.SetMaximumError( 0.001 );
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    gaussian.SetDirection( d );
    gaussian.SetMaximumKernelWidth(
      outputField->GetRequestedRegion().GetSize()[d] );
    gaussian.CreateDirectional();

    smoother->SetOperator( gaussian );
    smoother->SetInput( outputField );

    outputField = smoother->GetOutput();
    outputField->Update();
    outputField->DisconnectPipeline();
    }

  // Ensure zero motion on the boundary

  RealType weight1 = 1.0;
  if( variance < 0.5 )
    {
    weight1 = 1.0 - 1.0 * ( variance / 0.5 );
    }
  RealType weight2 = 1.0 - weight1;

  typedef MultiplyByConstantImageFilter<DisplacementFieldType, RealType,
                                        DisplacementFieldType> MultiplierType;

  typename MultiplierType::Pointer multiplier1 = MultiplierType::New();
  multiplier1->SetConstant2( weight1 );
  multiplier1->SetInput1( outputField );

  typename MultiplierType::Pointer multiplier2 = MultiplierType::New();
  multiplier2->SetConstant2( weight2 );
  multiplier2->SetInput1( inputField );

  typedef AddImageFilter<DisplacementFieldType, DisplacementFieldType, DisplacementFieldType> AdderType;
  typename AdderType::Pointer adder = AdderType::New();
  adder->SetInput1( multiplier1->GetOutput() );
  adder->SetInput2( multiplier2->GetOutput() );

  outputField = adder->GetOutput();
  outputField->Update();
  outputField->DisconnectPipeline();

  VectorType zeroVector( 0.0 );

  ImageLinearIteratorWithIndex<DisplacementFieldType> It( outputField,
                                                          outputField->GetRequestedRegion() );
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    It.SetDirection( d );
    It.GoToBegin();
    while( !It.IsAtEnd() )
      {
      It.GoToBeginOfLine();
      It.Set( zeroVector );
      It.GoToEndOfLine();
      --It;
      It.Set( zeroVector );

      It.NextLine();
      }
    }

  return outputField;
}

template <class TInputImage, class TOutputImage>
typename DiReCTImageFilter<TInputImage, TOutputImage>::RealImagePointer
DiReCTImageFilter<TInputImage, TOutputImage>
::SmoothImage( const RealImageType *inputImage, const RealType variance )
{
  typedef DiscreteGaussianImageFilter<RealImageType, RealImageType> SmootherType;
  typename SmootherType::Pointer smoother = SmootherType::New();
  smoother->SetVariance( variance );
  smoother->SetUseImageSpacingOff();
  smoother->SetMaximumError( 0.01 );
  smoother->SetInput( inputImage );

  typename RealImageType::Pointer smoothImage = smoother->GetOutput();
  smoothImage->Update();
  smoothImage->DisconnectPipeline();

  return smoothImage;
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutputImage>
void
DiReCTImageFilter<TInputImage, TOutputImage>
::PrintSelf( std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  std::cout << indent << "Gray matter label = "
                   << this->m_GrayMatterLabel << std::endl;
  std::cout << indent << "White matter label = "
                   << this->m_WhiteMatterLabel << std::endl;
  std::cout << indent << "Maximum number of iterations = "
                   << this->m_MaximumNumberOfIterations << std::endl;
  std::cout << indent << "Thickness prior estimate = "
                   << this->m_ThicknessPriorEstimate << std::endl;
  std::cout << indent << "Smoothing sigma = "
                   << this->m_SmoothingVariance << std::endl;
  std::cout << indent << "Smoothing velocity field sigma = "
                   << this->m_SmoothingVelocityFieldVariance << std::endl;
  std::cout << indent << "Number of integration points = "
                   << this->m_NumberOfIntegrationPoints << std::endl;
  std::cout << indent << "Maximum number of invert displacement field iterations = "
                   << this->m_MaximumNumberOfInvertDisplacementFieldIterations << std::endl;
  std::cout << indent << "Initial gradient step = "
                   << this->m_InitialGradientStep << std::endl;
  std::cout << indent << "Current gradient step = "
                   << this->m_CurrentGradientStep << std::endl;
  std::cout << indent << "Convergence threshold = "
                   << this->m_ConvergenceThreshold << std::endl;
  std::cout << indent << "Convergence window size = "
                   << this->m_ConvergenceWindowSize << std::endl;
}
} // end namespace itk

#endif
