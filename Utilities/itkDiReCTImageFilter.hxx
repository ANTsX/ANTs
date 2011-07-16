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
#include "itkAndImageFilter.h"
#include "itkBinaryContourImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBSplineControlPointImageFunction.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkComposeDiffeomorphismsImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkGaussianOperator.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkIterationReporter.h"
#include "itkMaximumImageFilter.h"
#include "itkMultiplyByConstantImageFilter.h"
#include "itkMultiplyByConstantVectorImageFilter.h"
#include "itkOrImageFilter.h"
#include "itkPointSet.h"
#include "itkStatisticsImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkVectorNeighborhoodOperatorImageFilter.h"
#include "itkWarpImageFilter.h"

namespace itk
{
template <class TInputImage, class TOutputImage>
DiReCTImageFilter<TInputImage, TOutputImage>
::DiReCTImageFilter() :
  m_ThicknessPriorEstimate( 10.0 ),
  m_SmoothingSigma( 1.5 ),
  m_GradientStep( 0.5 ),
  m_NumberOfIntegrationPoints( 10 ),
  m_GrayMatterLabel( 2 ),
  m_WhiteMatterLabel( 3 ),
  m_MaximumNumberOfIterations( 50 ),
  m_CurrentEnergy( NumericTraits<RealType>::max() ),
  m_ConvergenceThreshold( 0.001 ),
  m_ConvergenceWindowSize( 10 )
{
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
  // Convert all input direction matrices to identities saving the original
  // directions to put back at the end of filter processing. We do this
  // because the white and gray matters reside in the same space and the
  // assumption simplifies the underlying registration code.

  typename InputImageType::DirectionType identity;
  identity.SetIdentity();
  unsigned int                                        numinputs = this->GetNumberOfInputs();
  std::vector<typename InputImageType::DirectionType> directions;
  for( unsigned int d = 0; d < numinputs; d++ )
    {
    directions.push_back(  this->GetInput( d )->GetDirection()  );
    const_cast<InputImageType *>( this->GetInput( d ) )->SetDirection( identity );
    }

  // Extract the gray and white matter segmentations and combine to form the
  // gm/wm region.  Dilate the latter region by 1 voxel.

  InputImagePointer grayMatter = this->ExtractRegion(
      this->GetSegmentationImage(), this->m_GrayMatterLabel );
  InputImagePointer whiteMatter = this->ExtractRegion(
      this->GetSegmentationImage(), this->m_WhiteMatterLabel );

  typedef AddImageFilter<InputImageType, InputImageType, InputImageType> AdderType;
  typename AdderType::Pointer adder = AdderType::New();
  adder->SetInput1( grayMatter );
  adder->SetInput2( whiteMatter );
  adder->Update();

  InputImagePointer thresholdedRegion = this->ExtractRegion(
      const_cast<const InputImageType *>( adder->GetOutput() ), 1 );

  typedef BinaryBallStructuringElement<InputPixelType, ImageDimension>
    StructuringElementType;
  typedef BinaryDilateImageFilter<InputImageType, InputImageType,
                                  StructuringElementType> DilatorType;

  StructuringElementType structuringElement;
  structuringElement.SetRadius( 1 );
  structuringElement.CreateStructuringElement();

  typename DilatorType::Pointer dilator = DilatorType::New();
  dilator->SetInput( thresholdedRegion );
  dilator->SetKernel( structuringElement );
  dilator->SetDilateValue( 1 );
  dilator->Update();

  InputImagePointer dilatedMatters = dilator->GetOutput();

  // Extract the white and gm/wm matter contours

  InputImagePointer dilatedMatterContours = this->ExtractRegionalContours(
      dilatedMatters, 1 );
  InputImagePointer whiteMatterContoursTmp = this->ExtractRegionalContours(
      this->GetSegmentationImage(), this->m_WhiteMatterLabel );

  typedef CastImageFilter<InputImageType, RealImageType> CasterType;
  typename CasterType::Pointer caster = CasterType::New();
  caster->SetInput( whiteMatterContoursTmp );
  caster->Update();
  RealImagePointer whiteMatterContours = caster->GetOutput();

  // Create mask image prior to the use of the boolean logic used in the code
  // to avoid performing

  typedef AndImageFilter<InputImageType, InputImageType, InputImageType>
    AndFilterType;
  typedef OrImageFilter<InputImageType, InputImageType, InputImageType>
    OrFilterType;

  typename OrFilterType::Pointer orFilter1 = OrFilterType::New();
  orFilter1->SetInput1( dilatedMatterContours );
  orFilter1->SetInput2( whiteMatterContoursTmp );

  typename OrFilterType::Pointer orFilter2 = OrFilterType::New();
  orFilter2->SetInput1( orFilter1->GetOutput() );
  orFilter2->SetInput2( grayMatter );

  typedef BinaryThresholdImageFilter<InputImageType, InputImageType>
    ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( this->GetSegmentationImage() );
  thresholder->SetLowerThreshold( 0 );
  thresholder->SetUpperThreshold( 0 );
  thresholder->SetInsideValue( 0 );
  thresholder->SetOutsideValue( 1 );
  thresholder->Update();

  typename AndFilterType::Pointer andFilter = AndFilterType::New();
  andFilter->SetInput1( orFilter2->GetOutput() );
  andFilter->SetInput2( thresholder->GetOutput() );
  andFilter->Update();

  InputImagePointer maskImage = andFilter->GetOutput();

  // Initialize fields and images.

  VectorType zeroVector( 0.0 );

  RealImagePointer corticalThicknessImage = RealImageType::New();
  corticalThicknessImage->CopyInformation( this->GetInput() );
  corticalThicknessImage->SetRegions( this->GetInput()->GetRequestedRegion() );
  corticalThicknessImage->Allocate();
  corticalThicknessImage->FillBuffer( 0.0 );

  VectorImagePointer forwardIncrementalField = VectorImageType::New();
  forwardIncrementalField->CopyInformation( this->GetInput() );
  forwardIncrementalField->SetRegions( this->GetInput()->GetRequestedRegion() );
  forwardIncrementalField->Allocate();

  RealImagePointer hitImage = RealImageType::New();
  hitImage->CopyInformation( this->GetInput() );
  hitImage->SetRegions( this->GetInput()->GetRequestedRegion() );
  hitImage->Allocate();

  VectorImagePointer integratedField = VectorImageType::New();
  integratedField->CopyInformation( this->GetInput() );
  integratedField->SetRegions( this->GetInput()->GetRequestedRegion() );
  integratedField->Allocate();
  integratedField->FillBuffer( zeroVector );

  VectorImagePointer inverseField = VectorImageType::New();
  inverseField->CopyInformation( this->GetInput() );
  inverseField->SetRegions( this->GetInput()->GetRequestedRegion() );
  inverseField->Allocate();

  VectorImagePointer inverseIncrementalField = VectorImageType::New();
  inverseIncrementalField->CopyInformation( this->GetInput() );
  inverseIncrementalField->SetRegions( this->GetInput()->GetRequestedRegion() );
  inverseIncrementalField->Allocate();

  RealImagePointer speedImage = RealImageType::New();
  speedImage->CopyInformation( this->GetInput() );
  speedImage->SetRegions( this->GetInput()->GetRequestedRegion() );
  speedImage->Allocate();

  RealImagePointer thicknessImage = RealImageType::New();
  thicknessImage->CopyInformation( this->GetInput() );
  thicknessImage->SetRegions( this->GetInput()->GetRequestedRegion() );
  thicknessImage->Allocate();

  RealImagePointer totalImage = RealImageType::New();
  totalImage->CopyInformation( this->GetInput() );
  totalImage->SetRegions( this->GetInput()->GetRequestedRegion() );
  totalImage->Allocate();

  VectorImagePointer velocityField = VectorImageType::New();
  velocityField->CopyInformation( this->GetInput() );
  velocityField->SetRegions( this->GetInput()->GetRequestedRegion() );
  velocityField->Allocate();
  velocityField->FillBuffer( zeroVector );

  // Instantiate most of the iterators all in one place

  ImageRegionIterator<RealImageType> ItCorticalThicknessImage(
    corticalThicknessImage,
    corticalThicknessImage->GetRequestedRegion() );
  ImageRegionConstIterator<RealImageType> ItGrayMatterProbabilityMap(
    this->GetGrayMatterProbabilityImage(),
    this->GetGrayMatterProbabilityImage()->GetRequestedRegion() );
  ImageRegionIterator<RealImageType> ItHitImage(
    hitImage,
    hitImage->GetRequestedRegion() );
  ImageRegionIterator<VectorImageType> ItForwardIncrementalField(
    forwardIncrementalField,
    forwardIncrementalField->GetRequestedRegion() );
  ImageRegionConstIterator<InputImageType> ItDilatedMatterContours(
    dilatedMatterContours,
    dilatedMatterContours->GetRequestedRegion() );
  ImageRegionIterator<VectorImageType> ItIntegratedField(
    integratedField,
    integratedField->GetRequestedRegion() );
  ImageRegionIterator<VectorImageType> ItInverseIncrementalField(
    inverseIncrementalField,
    inverseIncrementalField->GetRequestedRegion() );
  ImageRegionConstIterator<InputImageType> ItMaskImage(
    maskImage,
    maskImage->GetRequestedRegion() );
  ImageRegionConstIteratorWithIndex<InputImageType> ItSegmentationImage(
    this->GetSegmentationImage(),
    this->GetSegmentationImage()->GetRequestedRegion() );
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

  // Instantiate objects for profiling energy convergence

  typedef Vector<RealType, 1>                   ProfilePointDataType;
  typedef Image<ProfilePointDataType, 1>        CurveType;
  typedef PointSet<ProfilePointDataType, 1>     EnergyProfileType;
  typedef typename EnergyProfileType::PointType ProfilePointType;

  typename EnergyProfileType::Pointer energyProfile = EnergyProfileType::New();
  energyProfile->Initialize();

  // Instantiate the progress reporter

  IterationReporter reporter( this, 0, 1 );

  bool isConverged = false;
  this->m_CurrentConvergenceMeasurement = NumericTraits<RealType>::max();
  this->m_ElapsedIterations = 0;
  while( this->m_ElapsedIterations++ < this->m_MaximumNumberOfIterations &&
         isConverged == false )
    {
    ProfilePointDataType currentEnergy;
    currentEnergy[0] = 0.0;
    RealType numberOfGrayMatterVoxels = 0.0;

    forwardIncrementalField->FillBuffer( zeroVector );
    inverseField->FillBuffer( zeroVector );
    inverseIncrementalField->FillBuffer( zeroVector );

    hitImage->FillBuffer( 0.0 );
    totalImage->FillBuffer( 0.0 );
    thicknessImage->FillBuffer( 0.0 );

    ImageRegionIterator<VectorImageType> ItVelocityField(
      velocityField,
      velocityField->GetRequestedRegion() );

    unsigned int integrationPoint = 0;
    while( integrationPoint++ < this->m_NumberOfIntegrationPoints )
      {
      typedef ComposeDiffeomorphismsImageFilter<VectorImageType> ComposerType;
      typename ComposerType::Pointer composer = ComposerType::New();
      composer->SetDeformationField( inverseIncrementalField );
      composer->SetWarpingField( inverseField );
      composer->Update();

      inverseField = composer->GetOutput();
      inverseField->DisconnectPipeline();

      RealImagePointer warpedWhiteMatterProbabilityMap = this->WarpImage(
          this->GetWhiteMatterProbabilityImage(), inverseField );
      RealImagePointer warpedWhiteMatterContours = this->WarpImage(
          whiteMatterContours, inverseField );
      RealImagePointer warpedThicknessImage = this->WarpImage(
          thicknessImage, inverseField );

      typedef GradientRecursiveGaussianImageFilter<RealImageType, VectorImageType>
        GradientImageFilterType;
      typename GradientImageFilterType::Pointer gradientFilter =
        GradientImageFilterType::New();
      gradientFilter->SetInput( warpedWhiteMatterProbabilityMap );
      gradientFilter->SetSigma( this->m_SmoothingSigma );
      gradientFilter->Update();

      VectorImagePointer gradientImage = gradientFilter->GetOutput();

      // Instantiate the iterators all in one place

      ImageRegionIterator<VectorImageType> ItGradientImage(
        gradientImage,
        gradientImage->GetRequestedRegion() );
      ImageRegionIterator<VectorImageType> ItInverseField(
        inverseField,
        inverseField->GetRequestedRegion() );
      ImageRegionIterator<RealImageType> ItWarpedThicknessImage(
        warpedThicknessImage,
        warpedThicknessImage->GetRequestedRegion() );
      ImageRegionIterator<RealImageType> ItWarpedWhiteMatterProbabilityMap(
        warpedWhiteMatterProbabilityMap,
        warpedWhiteMatterProbabilityMap->GetRequestedRegion() );
      ImageRegionIterator<RealImageType> ItWarpedWhiteMatterContours(
        warpedWhiteMatterContours,
        warpedWhiteMatterContours->GetRequestedRegion() );

      // Generate speed image

      speedImage->FillBuffer( 0.0 );

      ItGradientImage.GoToBegin();
      ItGrayMatterProbabilityMap.GoToBegin();
      ItSegmentationImage.GoToBegin();
      ItSpeedImage.GoToBegin();
      ItWarpedWhiteMatterProbabilityMap.GoToBegin();

      while( !ItSegmentationImage.IsAtEnd() )
        {
        if( ItSegmentationImage.Get() == this->m_GrayMatterLabel )
          {
          RealType norm = ( ItGradientImage.Get() ).GetNorm();
          if( norm > 1e-3 && !vnl_math_isnan( norm ) && !vnl_math_isinf( norm ) )
            {
            ItGradientImage.Set( ItGradientImage.Get() / norm );
            }
          else
            {
            ItGradientImage.Set( zeroVector );
            }
          RealType delta = ( ItWarpedWhiteMatterProbabilityMap.Get()
                             - ItGrayMatterProbabilityMap.Get() );

          currentEnergy[0] += vnl_math_abs( delta );
          numberOfGrayMatterVoxels++;

          RealType speedValue = -1.0 * delta * ItGrayMatterProbabilityMap.Get()
            * this->m_GradientStep;
          if( vnl_math_isnan( speedValue ) || vnl_math_isinf( speedValue ) )
            {
            speedValue = 0.0;
            }
          ItSpeedImage.Set( speedValue );
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
      ItIntegratedField.GoToBegin();
      ItInverseField.GoToBegin();
      ItInverseIncrementalField.GoToBegin();
      ItMaskImage.GoToBegin();
      ItSegmentationImage.GoToBegin();
      ItSpeedImage.GoToBegin();
      ItVelocityField.GoToBegin();
      ItWhiteMatterContours.GoToBegin();

      while( !ItSegmentationImage.IsAtEnd() )
        {
        typename InputImageType::IndexType index =
          ItSegmentationImage.GetIndex();
        typename InputImageType::PixelType segmentationValue =
          ItSegmentationImage.Get();

        if( !ItMaskImage.Get() )
          {
          ItIntegratedField.Set( zeroVector );
          ItInverseField.Set( zeroVector );
          ItVelocityField.Set( zeroVector );
          }
        ItInverseIncrementalField.Set( ItVelocityField.Get() );
        ItForwardIncrementalField.Set( ItForwardIncrementalField.Get()
                                       + ItGradientImage.Get() * ItSpeedImage.Get() );

        if( segmentationValue == this->m_GrayMatterLabel ||
            segmentationValue == this->m_WhiteMatterLabel )
          {
          if( integrationPoint == 1 )
            {
            typename InputImageType::PixelType whiteMatterContoursValue =
              ItWhiteMatterContours.Get();
            hitImage->SetPixel( index, whiteMatterContoursValue );

            VectorType vector = integratedField->GetPixel( index );
            RealType   weightedNorm = vector.GetNorm() * whiteMatterContoursValue;

            thicknessImage->SetPixel( index, weightedNorm );
            totalImage->SetPixel( index, weightedNorm );
            }
          else if( segmentationValue == this->m_GrayMatterLabel )
            {
            hitImage->SetPixel( index, hitImage->GetPixel( index )
                                + warpedWhiteMatterContours->GetPixel( index ) );
            totalImage->SetPixel( index, totalImage->GetPixel( index )
                                  + warpedThicknessImage->GetPixel( index ) );
            }
          }

        ++ItForwardIncrementalField;
        ++ItGradientImage;
        ++ItIntegratedField;
        ++ItInverseField;
        ++ItInverseIncrementalField;
        ++ItMaskImage;
        ++ItSegmentationImage;
        ++ItSpeedImage;
        ++ItVelocityField;
        ++ItWhiteMatterContours;
        }

      if( integrationPoint == 1 )
        {
        integratedField->FillBuffer( zeroVector );
        }
      this->InvertDeformationField( inverseField, integratedField );
      this->InvertDeformationField( integratedField, inverseField );
      }

    // calculate the size of the solution to allow us to adjust the
    // gradient step length.

    RealType maxNorm = 0;

    typename InputImageType::SpacingType spacing = grayMatter->GetSpacing();

    ItIntegratedField.GoToBegin();
    for( ItIntegratedField.GoToBegin(); !ItIntegratedField.IsAtEnd();
         ++ItIntegratedField )
      {
      VectorType vector = ItIntegratedField.Get();
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        vector[d] = vector[d] / spacing[d];
        }
      RealType norm = vector.GetNorm();
      if( norm > maxNorm )
        {
        maxNorm = norm;
        }
      }
    std::cout << " max_norm " << maxNorm << " it " << this->m_ElapsedIterations << std::endl;
    if( this->m_ElapsedIterations == 2 )
      {
      this->m_GradientStep = this->m_GradientStep * 2.0 / maxNorm;
      velocityField->FillBuffer( zeroVector );
      }

    ItCorticalThicknessImage.GoToBegin();
    ItForwardIncrementalField.GoToBegin();
    ItHitImage.GoToBegin();
    ItSegmentationImage.GoToBegin();
    ItTotalImage.GoToBegin();
    ItVelocityField.GoToBegin();

    while( !ItSegmentationImage.IsAtEnd() )
      {
      ItVelocityField.Set( ItVelocityField.Get()
                           + ItForwardIncrementalField.Get() );
      if( ItSegmentationImage.Get() == this->m_GrayMatterLabel )
        {
        RealType thicknessValue = 0.0;
        if( ItHitImage.Get() > 0.001 )
          {
          thicknessValue = ItTotalImage.Get() / ItHitImage.Get();
          if( thicknessValue < 0.0 )
            {
            thicknessValue = 0.0;
            }
          if( thicknessValue > this->m_ThicknessPriorEstimate )
            {
            thicknessValue = this->m_ThicknessPriorEstimate;
            }
          }

        ItCorticalThicknessImage.Set( thicknessValue );
        }

      ++ItCorticalThicknessImage;
      ++ItForwardIncrementalField;
      ++ItHitImage;
      ++ItSegmentationImage;
      ++ItTotalImage;
      ++ItVelocityField;
      }

    velocityField = this->SmoothDeformationField( velocityField,
                                                  this->m_SmoothingSigma );

    // Calculate current energy and current convergence measurement

    currentEnergy[0] /= numberOfGrayMatterVoxels;
    this->m_CurrentEnergy = currentEnergy[0];

    ProfilePointType point;
    point[0] = this->m_ElapsedIterations - 1;

    energyProfile->SetPoint( this->m_ElapsedIterations - 1, point );
    energyProfile->SetPointData( this->m_ElapsedIterations - 1, currentEnergy );

    if( this->m_ElapsedIterations >= this->m_ConvergenceWindowSize )
      {
      typename CurveType::PointType    origin;
      typename CurveType::SizeType     size;
      typename CurveType::SpacingType  spacing;

      origin[0] = this->m_ElapsedIterations - this->m_ConvergenceWindowSize;
      size[0] = this->m_ConvergenceWindowSize;
      spacing[0] = 1.0;

      typedef BSplineScatteredDataPointSetToImageFilter<EnergyProfileType,
                                                        CurveType> BSplinerType;
      typename BSplinerType::Pointer bspliner = BSplinerType::New();

      typename EnergyProfileType::Pointer energyProfileWindow =
        EnergyProfileType::New();
      energyProfileWindow->Initialize();

      RealType totalEnergy = 0.0;

      unsigned int startIndex = static_cast<unsigned int>( origin[0] );
      for( unsigned int i = startIndex; i < this->m_ElapsedIterations; i++ )
        {
        ProfilePointType windowPoint;
        windowPoint[0] =
          static_cast<typename ProfilePointType::CoordRepType>( i );

        ProfilePointDataType windowEnergy;
        windowEnergy.Fill( 0.0 );
        energyProfile->GetPointData( i, &windowEnergy );

        totalEnergy += vnl_math_abs( windowEnergy[0] );
        }
      for( unsigned int i = startIndex; i < this->m_ElapsedIterations; i++ )
        {
        ProfilePointType windowPoint;
        windowPoint[0] = static_cast<typename ProfilePointType::CoordRepType>( i );

        ProfilePointDataType windowEnergy;
        windowEnergy.Fill( 0.0 );
        energyProfile->GetPointData( i, &windowEnergy );

        energyProfileWindow->SetPoint( i - startIndex, windowPoint );
        energyProfileWindow->SetPointData( i - startIndex,
                                           windowEnergy / totalEnergy );
        }

      bspliner->SetInput( energyProfileWindow );
      bspliner->SetOrigin( origin );
      bspliner->SetSpacing( spacing );
      bspliner->SetSize( size );
      bspliner->SetNumberOfLevels( 1 );
      bspliner->SetSplineOrder( 1 );
      typename BSplinerType::ArrayType ncps;
      ncps.Fill( bspliner->GetSplineOrder()[0] + 1 );
      bspliner->SetNumberOfControlPoints( ncps );
      bspliner->Update();

      typedef BSplineControlPointImageFunction<CurveType> BSplinerFunctionType;
      typename BSplinerFunctionType::Pointer bsplinerFunction =
        BSplinerFunctionType::New();
      bsplinerFunction->SetOrigin( origin );
      bsplinerFunction->SetSpacing( spacing );
      bsplinerFunction->SetSize( size );
      bsplinerFunction->SetSplineOrder( bspliner->GetSplineOrder() );
      bsplinerFunction->SetInputImage( bspliner->GetPhiLattice() );

      ProfilePointType endPoint;
      endPoint[0] = static_cast<RealType>( this->m_ElapsedIterations - 1 );
      typename BSplinerFunctionType::GradientType gradient =
        bsplinerFunction->EvaluateGradientAtParametricPoint( endPoint );
      this->m_CurrentConvergenceMeasurement = -gradient[0][0];

      if( this->m_CurrentConvergenceMeasurement < this->m_ConvergenceThreshold )
        {
        isConverged = true;
        }
      }

    reporter.CompletedStep();
    }

  this->SetNthOutput( 0, corticalThicknessImage );
  // Replace direction matrices to the inputs.
  for( unsigned int d = 0; d < this->GetNumberOfInputs(); d++ )
    {
    const_cast<InputImageType *>( this->GetInput( d ) )->
    SetDirection( directions[d] );
    }
}

template <class TInputImage, class TOutputImage>
typename DiReCTImageFilter<TInputImage, TOutputImage>::InputImagePointer
DiReCTImageFilter<TInputImage, TOutputImage>
::ExtractRegion( const InputImageType *segmentationImage,
                 unsigned int whichRegion )
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
                           unsigned int whichRegion )
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
typename DiReCTImageFilter<TInputImage, TOutputImage>::RealImagePointer
DiReCTImageFilter<TInputImage, TOutputImage>
::WarpImage( const RealImageType *inputImage,
             const VectorImageType *deformationField )
{
  typedef WarpImageFilter<RealImageType, RealImageType, VectorImageType>
    WarperType;
  typename WarperType::Pointer warper = WarperType::New();
  warper->SetInput( inputImage );
  warper->SetDeformationField( deformationField );
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
void
DiReCTImageFilter<TInputImage, TOutputImage>
::InvertDeformationField( const VectorImageType *deformationField,
                          VectorImageType *inverseField )
{
  typename VectorImageType::SpacingType spacing =
    deformationField->GetSpacing();
  VectorType spacingFactor;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    spacingFactor[d] = 1.0 / spacing[d];
    }

  RealType     maxNorm = 1.0;
  RealType     meanNorm = 1.0;
  unsigned int iteration = 0;
  while( iteration++ < 20 && maxNorm > 0.1 && meanNorm > 0.001 )
    {
    meanNorm = 0.0;
    maxNorm = 0.0;

    typedef ComposeDiffeomorphismsImageFilter<VectorImageType> ComposerType;
    typename ComposerType::Pointer composer = ComposerType::New();
    composer->SetDeformationField( deformationField );
    composer->SetWarpingField( inverseField );

    typedef MultiplyByConstantVectorImageFilter<VectorImageType, VectorType,
                                                VectorImageType> ConstantMultiplierType;
    typename ConstantMultiplierType::Pointer constantMultiplier =
      ConstantMultiplierType::New();
    constantMultiplier->SetConstantVector( spacingFactor );
    constantMultiplier->SetInput( composer->GetOutput() );

    typedef VectorMagnitudeImageFilter<VectorImageType, RealImageType>
      NormFilterType;
    typename NormFilterType::Pointer normFilter = NormFilterType::New();
    normFilter->SetInput( constantMultiplier->GetOutput() );

    typedef StatisticsImageFilter<RealImageType> StatisticsType;
    typename StatisticsType::Pointer statistics = StatisticsType::New();
    statistics->SetInput( normFilter->GetOutput() );
    statistics->Update();

    meanNorm = statistics->GetMean();
    maxNorm = statistics->GetMaximum();

    RealType epsilon = 0.5;
    if( iteration == 1 )
      {
      epsilon = 0.75;
      }
    RealType normFactor = 1.0;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      normFactor /= spacing[d];
      }

    ImageRegionIterator<VectorImageType> ItE( composer->GetOutput(),
                                              composer->GetOutput()->GetRequestedRegion() );
    ImageRegionIterator<VectorImageType> ItI( inverseField,
                                              inverseField->GetRequestedRegion() );
    for( ItI.GoToBegin(), ItE.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItE )
      {
      VectorType update = -ItE.Get();
      RealType   updateNorm = update.GetNorm();

      if( updateNorm > epsilon * maxNorm / normFactor )
        {
        update *= ( epsilon * maxNorm / ( updateNorm * normFactor ) );
        }
      ItI.Set( ItI.Get() + update * epsilon );
      }
    }
}

template <class TInputImage, class TOutputImage>
typename DiReCTImageFilter<TInputImage, TOutputImage>::VectorImagePointer
DiReCTImageFilter<TInputImage, TOutputImage>
::SmoothDeformationField( const VectorImageType *inputField,
                          const RealType variance )
{
  typedef ImageDuplicator<VectorImageType> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage( inputField );
  duplicator->Update();
  VectorImagePointer outputField = duplicator->GetOutput();

  typedef VectorNeighborhoodOperatorImageFilter<VectorImageType,
                                                VectorImageType> SmootherType;
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

  typedef MultiplyByConstantImageFilter<VectorImageType, RealType,
                                        VectorImageType> MultiplierType;

  typename MultiplierType::Pointer multiplier1 = MultiplierType::New();
  multiplier1->SetConstant2( weight1 );
  multiplier1->SetInput1( outputField );

  typename MultiplierType::Pointer multiplier2 = MultiplierType::New();
  multiplier2->SetConstant2( weight2 );
  multiplier2->SetInput1( inputField );

  typedef AddImageFilter<VectorImageType, VectorImageType, VectorImageType>
    AdderType;
  typename AdderType::Pointer adder = AdderType::New();
  adder->SetInput1( multiplier1->GetOutput() );
  adder->SetInput2( multiplier2->GetOutput() );

  outputField = adder->GetOutput();
  outputField->Update();
  outputField->DisconnectPipeline();

  VectorType zeroVector( 0.0 );

  ImageLinearIteratorWithIndex<VectorImageType> It( outputField,
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
            << this->m_SmoothingSigma << std::endl;
  std::cout << indent << "Gradient step = "
            << this->m_GradientStep << std::endl;
  std::cout << indent << "Convergence threshold = "
            << this->m_ConvergenceThreshold << std::endl;
  std::cout << indent << "Convergence window size = "
            << this->m_ConvergenceWindowSize << std::endl;
}
} // end namespace itk

#endif
