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
#ifndef __itkDiReCTImageFilter953_hxx
#define __itkDiReCTImageFilter953_hxx

#include "itkDiReCTImageFilter953.h"

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
#include "itkLabelMap.h"
#include "itkLabelImageToLabelMapFilter.h"
#include "itkMaximumImageFilter.h"
#include "itkMultiplyByConstantImageFilter.h"
#include "itkMultiplyByConstantVectorImageFilter.h"
#include "itkOrImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkVectorNeighborhoodOperatorImageFilter.h"
#include "itkWarpImageFilter.h"

#include "itkImageFileWriter.h"

namespace itk
{
template <class TInputImage, class TOutputImage>
DiReCTImageFilter953<TInputImage, TOutputImage>
::DiReCTImageFilter953() :
  m_ThicknessPriorEstimate( 10.0 ),
  m_SmoothingSigma( 1.5 ),
  m_GradientStep( 0.005 ),
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
DiReCTImageFilter953<TInputImage, TOutputImage>
::~DiReCTImageFilter953()
{
}

template <class TInputImage, class TOutputImage>
void
DiReCTImageFilter953<TInputImage, TOutputImage>
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
    directions.push_back( this->GetInput( d )->GetDirection()  );
    const_cast<InputImageType *>( this->GetInput( d ) )->SetDirection( identity );
    }

  // Extract the gray and white matter segmentations and combine to form the
  // gm/wm region.  Dilate the latter region by 1 voxel.

  typedef LabelImageToLabelMapFilter<InputImageType>     LabelImageFilterType;
  typedef typename LabelImageFilterType::OutputImageType LabelMapType;

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

  // Extract the white and gm/wm matter contours

  InputImagePointer dilatedMatterContours = this->ExtractRegionalContours(
      dilator->GetOutput(), 1 );
  InputImagePointer whiteMatterContoursTmp = this->ExtractRegionalContours(
      this->GetSegmentationImage(), this->m_WhiteMatterLabel );

  typedef CastImageFilter<InputImageType, RealImageType> CasterType;
  typename CasterType::Pointer caster = CasterType::New();
  caster->SetInput( whiteMatterContoursTmp );
  caster->Update();

  // Create mask image prior to the use of the boolean logic used in the code

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

  typename AndFilterType::Pointer andFilter = AndFilterType::New();
  andFilter->SetInput1( orFilter2->GetOutput() );
  andFilter->SetInput2( thresholder->GetOutput() );

  InputImagePointer velocityFieldMaskImage = andFilter->GetOutput();
  velocityFieldMaskImage->Update();
  velocityFieldMaskImage->DisconnectPipeline();

  typename ThresholderType::Pointer thresholderWhiteMatter
    = ThresholderType::New();
  thresholderWhiteMatter->SetInput( this->GetSegmentationImage() );
  thresholderWhiteMatter->SetLowerThreshold( this->m_WhiteMatterLabel );
  thresholderWhiteMatter->SetUpperThreshold( this->m_WhiteMatterLabel );
  thresholderWhiteMatter->SetInsideValue( 1 );
  thresholderWhiteMatter->SetOutsideValue( 0 );

  typename ThresholderType::Pointer thresholderGrayMatter
    = ThresholderType::New();
  thresholderGrayMatter->SetInput( this->GetSegmentationImage() );
  thresholderGrayMatter->SetLowerThreshold( this->m_GrayMatterLabel );
  thresholderGrayMatter->SetUpperThreshold( this->m_GrayMatterLabel );
  thresholderGrayMatter->SetInsideValue( 1 );
  thresholderGrayMatter->SetOutsideValue( 0 );

  typename OrFilterType::Pointer mattersFilter = OrFilterType::New();
  mattersFilter->SetInput1( thresholderWhiteMatter->GetOutput() );
  mattersFilter->SetInput2( thresholderGrayMatter->GetOutput() );

  InputImagePointer mattersMaskImage = mattersFilter->GetOutput();
  mattersMaskImage->Update();
  mattersMaskImage->DisconnectPipeline();

  SparseImagePointer whiteMatterContours = this->ConvertRealImageToSparseImage(
      caster->GetOutput(), mattersMaskImage );

  typedef LabelImageToLabelMapFilter<InputImageType> LabelFilterType;
  typename LabelFilterType::Pointer labelFilter = LabelFilterType::New();
  labelFilter->SetInput( andFilter->GetOutput() );
  labelFilter->SetBackgroundValue( 0 );

  typedef typename LabelFilterType::OutputImageType LabelMapType;
  typename LabelMapType::Pointer velocityMaskImage = labelFilter->GetOutput();
  velocityMaskImage->Update();
  velocityMaskImage->DisconnectPipeline();

  // Initialize fields and images.

  SparseVectorImagePointer forwardIncrementalField = SparseVectorImageType::New();
  forwardIncrementalField->Initialize();

  SparseVectorImagePointer integratedField = SparseVectorImageType::New();
  integratedField->Initialize();

  SparseVectorImagePointer inverseField = SparseVectorImageType::New();
  inverseField->Initialize();

  SparseVectorImagePointer inverseIncrementalField = SparseVectorImageType::New();
  inverseIncrementalField->Initialize();

  SparseVectorImagePointer velocityField = SparseVectorImageType::New();
  velocityField->Initialize();

  SparseImagePointer hitImage = SparseImageType::New();
  hitImage->Initialize();

  SparseImagePointer totalImage = SparseImageType::New();
  totalImage->Initialize();

  SparseImagePointer thicknessField = SparseImageType::New();
  thicknessField->Initialize();

  SparseImagePointer speedImage = SparseImageType::New();
  speedImage->Initialize();

  SparseImagePointer corticalThicknessField = SparseImageType::New();
  corticalThicknessField->Initialize();

  VectorType zeroVector( 0.0 );

  // Instantiate most of the iterators all in one place

  ImageRegionConstIteratorWithIndex<InputImageType> ItSegmentationImage(
    this->GetSegmentationImage(),
    this->GetSegmentationImage()->GetRequestedRegion() );
  ImageRegionConstIterator<RealImageType> ItGrayMatterProbabilityMap(
    this->GetGrayMatterProbabilityImage(),
    this->GetGrayMatterProbabilityImage()->GetRequestedRegion() );

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

    forwardIncrementalField->Initialize();
    inverseField->Initialize();
    inverseIncrementalField->Initialize();

    unsigned int integrationPoint = 0;
    while( integrationPoint++ < this->m_NumberOfIntegrationPoints )
      {
      itkDebugMacro( "  Integration point" << integrationPoint << " of "
                                           << this->m_NumberOfIntegrationPoints );

      VectorImagePointer inverseFieldImage =
        this->ConvertSparseVectorImageToVectorImage(
          inverseField );
      VectorImagePointer inverseIncrementalFieldImage =
        this->ConvertSparseVectorImageToVectorImage(
          inverseIncrementalField );

      typedef ComposeDiffeomorphismsImageFilter<VectorImageType> ComposerType;
      typename ComposerType::Pointer composer = ComposerType::New();
      composer->SetDisplacementField( inverseIncrementalFieldImage );
      composer->SetWarpingField( inverseFieldImage );
      composer->Update();

      RealImagePointer whiteMatterContoursImage =
        this->ConvertSparseImageToRealImage( whiteMatterContours,
                                             mattersMaskImage );
      SparseImagePointer warpedWhiteMatterContours =
        this->ConvertRealImageToSparseImage(
          this->WarpImage( whiteMatterContoursImage, composer->GetOutput() ),
          mattersMaskImage );
      RealImagePointer thicknessImage =
        this->ConvertSparseImageToRealImage( thicknessField,
                                             mattersMaskImage );
      SparseImagePointer warpedThicknessField =
        this->ConvertRealImageToSparseImage( this->WarpImage(
                                               thicknessImage, composer->GetOutput() ), mattersMaskImage );
      RealImagePointer warpedWhiteMatterProbabilityMap = this->WarpImage(
          this->GetWhiteMatterProbabilityImage(), composer->GetOutput() );

      inverseField = this->ConvertVectorImageToSparseVectorImage(
          composer->GetOutput() );

      typedef GradientRecursiveGaussianImageFilter<RealImageType, VectorImageType>
        GradientImageFilterType;
      typename GradientImageFilterType::Pointer gradientFilter =
        GradientImageFilterType::New();
      gradientFilter->SetInput( warpedWhiteMatterProbabilityMap );
      gradientFilter->SetSigma( this->m_SmoothingSigma );
      gradientFilter->Update();

      SparseVectorImagePointer gradientField =
        this->ConvertVectorImageToSparseVectorImage(
          gradientFilter->GetOutput(), mattersMaskImage );

      ImageRegionIterator<RealImageType> ItWarpedWhiteMatterProbabilityMap(
        warpedWhiteMatterProbabilityMap,
        warpedWhiteMatterProbabilityMap->GetRequestedRegion() );

      // Generate speed image

      ItGrayMatterProbabilityMap.GoToBegin();
      ItSegmentationImage.GoToBegin();
      ItWarpedWhiteMatterProbabilityMap.GoToBegin();

      unsigned long mattersCount = 0;
      while( !ItSegmentationImage.IsAtEnd() )
        {
        InputPixelType label = ItSegmentationImage.Get();
        if( label == this->m_GrayMatterLabel || label == this->m_WhiteMatterLabel )
          {
          RealType   speedValue = 0.0;
          VectorType gradientVector = zeroVector;
          if( label == this->m_GrayMatterLabel )
            {
            gradientField->GetPointData( mattersCount, &gradientVector );
            RealType norm = gradientVector.GetNorm();
            if( norm > 1e-3 && !vnl_math_isnan( norm ) && !vnl_math_isinf( norm ) )
              {
              gradientField->SetPointData( mattersCount, gradientVector / norm );
              }
            else
              {
              gradientField->SetPointData( mattersCount, zeroVector );
              }
            RealType delta = ( ItWarpedWhiteMatterProbabilityMap.Get()
                               - ItGrayMatterProbabilityMap.Get() );

            currentEnergy[0] += vnl_math_abs( delta );
            numberOfGrayMatterVoxels++;

            speedValue = -1.0 * delta * ItGrayMatterProbabilityMap.Get()
              * this->m_GradientStep;
            if( vnl_math_isnan( speedValue ) || vnl_math_isinf( speedValue ) )
              {
              speedValue = 0.0;
              }
            }
          speedImage->SetPointData( mattersCount++, speedValue );
          }
        ++ItGrayMatterProbabilityMap;
        ++ItSegmentationImage;
        ++ItWarpedWhiteMatterProbabilityMap;
        }

      // Zero outside velocity field mask

      if( integratedField->GetPoints() )
        {
        typename SparseVectorImageType::PointsContainerIterator ItP1 =
          integratedField->GetPoints()->Begin();
        typename SparseVectorImageType::PointDataContainer::Iterator ItD1 =
          integratedField->GetPointData()->Begin();
        while( ItP1 != integratedField->GetPoints()->End() )
          {
          typename SparseVectorImageType::PointType point = ItP1.Value();
          IndexType index;
          for( unsigned int d = 0; d < ImageDimension; d++ )
            {
            index[d] = static_cast<IndexValueType>( point[d] );
            }
          if( !velocityMaskImage->GetPixel( index ) )
            {
            ItD1.Value() = zeroVector;
            }
          ++ItP1;
          ++ItD1;
          }
        }

      if( inverseField->GetPoints() )
        {
        typename SparseVectorImageType::PointsContainerIterator ItP2 =
          inverseField->GetPoints()->Begin();
        typename SparseVectorImageType::PointDataContainer::Iterator ItD2 =
          inverseField->GetPointData()->Begin();
        while( ItP2 != inverseField->GetPoints()->End() )
          {
          typename SparseVectorImageType::PointType point = ItP2.Value();
          IndexType index;
          for( unsigned int d = 0; d < ImageDimension; d++ )
            {
            index[d] = static_cast<IndexValueType>( point[d] );
            }
          if( !velocityMaskImage->GetPixel( index ) )
            {
            ItD2.Value() = zeroVector;
            }
          ++ItP2;
          ++ItD2;
          }
        }

      if( velocityField->GetPoints() )
        {
        typename SparseVectorImageType::PointsContainerIterator ItP3 =
          velocityField->GetPoints()->Begin();
        typename SparseVectorImageType::PointDataContainer::Iterator ItD3 =
          velocityField->GetPointData()->Begin();
        while( ItP3 != velocityField->GetPoints()->End() )
          {
          typename SparseVectorImageType::PointType point = ItP3.Value();
          IndexType index;
          for( unsigned int d = 0; d < ImageDimension; d++ )
            {
            index[d] = static_cast<IndexValueType>( point[d] );
            }
          if( !velocityMaskImage->GetPixel( index ) )
            {
            ItD3.Value() = zeroVector;
            }
          ++ItP3;
          ++ItD3;
          }
        }

      // Copy velocity field to inverse incremental field

      inverseIncrementalField->Initialize();

      unsigned long count = 0;

      typename SparseVectorImageType::PointsContainerIterator ItP =
        velocityField->GetPoints()->Begin();
      while( ItP != velocityField->GetPoints()->End() )
        {
        inverseIncrementalField->SetPoint( count++, ItP.Value() );
        ++ItP;
        }

      count = 0;
      typename SparseVectorImageType::PointDataContainerIterator ItD =
        velocityField->GetPointData()->Begin();
      while( ItD != velocityField->GetPointData()->End() )
        {
        inverseIncrementalField->SetPointData( count++, ItD.Value() );
        ++ItD;
        }

      // Calculate objective function value

      VectorImagePointer forwardIncrementalFieldImage =
        this->ConvertSparseVectorImageToVectorImage( forwardIncrementalField );
      VectorImagePointer integratedFieldImage =
        this->ConvertSparseVectorImageToVectorImage( integratedField );

      ImageRegionIterator<VectorImageType> ItForwardIncremental(
        forwardIncrementalFieldImage,
        forwardIncrementalFieldImage->GetRequestedRegion() );
      ImageRegionIterator<VectorImageType> ItIntegrated(
        integratedFieldImage,
        integratedFieldImage->GetRequestedRegion() );

      ItForwardIncremental.GoToBegin();
      ItIntegrated.GoToBegin();
      ItSegmentationImage.GoToBegin();

      mattersCount = 0;
      while( !ItSegmentationImage.IsAtEnd() )
        {
        InputPixelType label = ItSegmentationImage.Get();
        typename InputImageType::IndexType index =
          ItSegmentationImage.GetIndex();

        RealType   speedValue = 0.0;
        VectorType gradientVector( 0.0 );
        if( label == this->m_GrayMatterLabel ||
            label == this->m_WhiteMatterLabel )
          {
          speedImage->GetPointData( mattersCount, &speedValue );
          gradientField->GetPointData( mattersCount, &gradientVector );
          }
        ItForwardIncremental.Set( ItForwardIncremental.Get() + gradientVector
                                  * speedValue );

        if( label == this->m_GrayMatterLabel ||
            label == this->m_WhiteMatterLabel )
          {
          if( integrationPoint == 1 )
            {
            RealType whiteMatterContoursValue = 0.0;
            whiteMatterContours->GetPointData(
              mattersCount, &whiteMatterContoursValue );
            hitImage->SetPointData( mattersCount, whiteMatterContoursValue );

            VectorType vector = ItIntegrated.Get();
            RealType   weightedNorm = vector.GetNorm() * whiteMatterContoursValue;

            thicknessField->SetPointData( mattersCount, weightedNorm );
            totalImage->SetPointData( mattersCount, weightedNorm );
            }
          else if( label == this->m_GrayMatterLabel )
            {
            RealType hitValue = 0.0;
            hitImage->GetPointData( mattersCount, &hitValue );

            RealType totalValue = 0.0;
            totalImage->GetPointData( mattersCount, &totalValue );

            RealType warpedWhiteMatterContoursValue = 0.0;
            warpedWhiteMatterContours->GetPointData( mattersCount,
                                                     &warpedWhiteMatterContoursValue );

            hitImage->SetPointData( mattersCount, hitValue
                                    + warpedWhiteMatterContoursValue );

            RealType warpedThicknessValue = 0.0;
            warpedThicknessField->GetPointData( mattersCount,
                                                &warpedThicknessValue );
            totalImage->SetPointData( mattersCount, totalValue
                                      + warpedThicknessValue );
            }
          mattersCount++;
          }

        ++ItSegmentationImage;
        ++ItForwardIncremental;
        ++ItIntegrated;
        }

      forwardIncrementalField = this->ConvertVectorImageToSparseVectorImage(
          forwardIncrementalFieldImage );

      if( integrationPoint == 1 )
        {
        integratedField->Initialize();
        }
      else
        {
        integratedField = this->ConvertVectorImageToSparseVectorImage(
            integratedFieldImage );
        }

      VectorImagePointer inverseFieldImage3 =
        this->ConvertSparseVectorImageToVectorImage(
          inverseField );
      VectorImagePointer integratedFieldImage2 =
        this->ConvertSparseVectorImageToVectorImage(
          integratedField );

      this->InvertDisplacementField( inverseFieldImage3, integratedFieldImage2 );
      this->InvertDisplacementField( integratedFieldImage2, inverseFieldImage3 );

      inverseField = this->ConvertVectorImageToSparseVectorImage(
          inverseFieldImage3 );
      integratedField = this->ConvertVectorImageToSparseVectorImage(
          integratedFieldImage2 );
      }

    // calculate the size of the solution to allow us to adjust the
    // gradient step length.

    RealType maxNorm = 0.0;

    typename InputImageType::SpacingType spacing =
      this->GetSegmentationImage()->GetSpacing();

    typename SparseVectorImageType::PointDataContainerIterator ItIntegratedField =
      integratedField->GetPointData()->Begin();
    while( ItIntegratedField != integratedField->GetPointData()->End() )
      {
      VectorType vector = ItIntegratedField.Value();
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        vector[d] = vector[d] / spacing[d];
        }
      RealType norm = vector.GetNorm();
      if( norm > maxNorm )
        {
        maxNorm = norm;
        }
      ++ItIntegratedField;
      }

    itkDebugMacro( "   MaxNorm = " << maxNorm );

    if( this->m_ElapsedIterations == 2 )
      {
      this->m_GradientStep = this->m_GradientStep * 1.0 / maxNorm;
      velocityField->Initialize();
      }

    VectorImagePointer forwardIncrementalFieldImage =
      this->ConvertSparseVectorImageToVectorImage( forwardIncrementalField );
    VectorImagePointer velocityFieldImage =
      this->ConvertSparseVectorImageToVectorImage( velocityField );

    ImageRegionIterator<VectorImageType> ItForwardIncremental(
      forwardIncrementalFieldImage,
      forwardIncrementalFieldImage->GetRequestedRegion() );
    ImageRegionIterator<VectorImageType> ItVelocity(
      velocityFieldImage,
      velocityFieldImage->GetRequestedRegion() );

    ItSegmentationImage.GoToBegin();
    ItForwardIncremental.GoToBegin();
    ItVelocity.GoToBegin();
    unsigned mattersCount = 0;
    while( !ItSegmentationImage.IsAtEnd() )
      {
      InputPixelType label = ItSegmentationImage.Get();

      ItVelocity.Set( ItVelocity.Get() + ItForwardIncremental.Get() );

      if( label == this->m_GrayMatterLabel )
        {
        RealType thicknessValue = 0.0;
        RealType hitValue = 0.0;
        hitImage->GetPointData( mattersCount, &hitValue );
        if( hitValue > 0.001 )
          {
          RealType totalValue = 0.0;
          totalImage->GetPointData( mattersCount, &totalValue );
          thicknessValue = totalValue / hitValue;
          if( thicknessValue < 0.0 )
            {
            thicknessValue = 0.0;
            }
          if( thicknessValue > this->m_ThicknessPriorEstimate )
            {
            thicknessValue = this->m_ThicknessPriorEstimate;
            }
          }
        corticalThicknessField->SetPointData( mattersCount, thicknessValue );
        }
      if( label == this->m_GrayMatterLabel || label == this->m_WhiteMatterLabel )
        {
        mattersCount++;
        }

      ++ItSegmentationImage;
      ++ItVelocity;
      ++ItForwardIncremental;
      }

    velocityField = this->ConvertVectorImageToSparseVectorImage(
        this->SmoothDisplacementField( velocityFieldImage,
                                       this->m_SmoothingSigma ) );

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
        windowPoint[0] =
          static_cast<typename ProfilePointType::CoordRepType>( i );

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

  RealImagePointer output = this->ConvertSparseImageToRealImage(
      corticalThicknessField, mattersMaskImage );

  this->SetNthOutput( 0, output );
  // Replace direction matrices to the inputs.
  for( unsigned int d = 0; d < this->GetNumberOfInputs(); d++ )
    {
    const_cast<InputImageType *>( this->GetInput( d ) )->
    SetDirection( directions[d] );
    }
}

template <class TInputImage, class TOutputImage>
typename DiReCTImageFilter953<TInputImage, TOutputImage>::InputImagePointer
DiReCTImageFilter953<TInputImage, TOutputImage>
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

  InputImagePointer thresholdRegion = thresholder->GetOutput();
  thresholdRegion->Update();
  thresholdRegion->DisconnectPipeline();

  return thresholdRegion;
}

template <class TInputImage, class TOutputImage>
typename DiReCTImageFilter953<TInputImage, TOutputImage>::InputImagePointer
DiReCTImageFilter953<TInputImage, TOutputImage>
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
typename DiReCTImageFilter953<TInputImage, TOutputImage>::RealImagePointer
DiReCTImageFilter953<TInputImage, TOutputImage>
::WarpImage( const RealImageType *inputImage,
             const VectorImageType *deformationField )
{
  typedef WarpImageFilter<RealImageType, RealImageType, VectorImageType>
    WarperType;
  typename WarperType::Pointer warper = WarperType::New();
  warper->SetInput( inputImage );
  warper->SetDisplacementField( deformationField );
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
DiReCTImageFilter953<TInputImage, TOutputImage>
::InvertDisplacementField( const VectorImageType *deformationField,
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
    composer->SetDisplacementField( deformationField );
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
typename DiReCTImageFilter953<TInputImage, TOutputImage>::VectorImagePointer
DiReCTImageFilter953<TInputImage, TOutputImage>
::SmoothDisplacementField( const VectorImageType *inputField,
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

template <class TInputImage, class TOutputImage>
typename DiReCTImageFilter953<TInputImage, TOutputImage>::SparseVectorImagePointer
DiReCTImageFilter953<TInputImage, TOutputImage>
::ConvertVectorImageToSparseVectorImage( const VectorImageType *inputImage,
                                         const InputImageType *maskImage = NULL )
{
  SparseVectorImagePointer sparseVectorImage = SparseVectorImageType::New();

  sparseVectorImage->Initialize();

  unsigned long count = 0;

  if( maskImage )
    {
    ImageRegionConstIterator<VectorImageType> ItI( inputImage,
                                                   inputImage->GetRequestedRegion() );
    ImageRegionConstIterator<InputImageType> ItM( maskImage,
                                                  maskImage->GetRequestedRegion() );
    for( ItI.GoToBegin(), ItM.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItM )
      {
      if( ItM.Get() != 0 )
        {
        sparseVectorImage->SetPointData( count++, ItI.Get() );
        }
      }
    }
  else
    {
    ImageRegionConstIteratorWithIndex<VectorImageType> ItI( inputImage,
                                                            inputImage->GetRequestedRegion() );
    for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
      {
      if( ( ItI.Get() ).GetSquaredNorm() > 1e-10 )
        {
        IndexType       index = ItI.GetIndex();
        SparseIndexType point;
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          point[d] =
            static_cast<typename SparseIndexType::CoordRepType>( index[d] );
          }
        sparseVectorImage->SetPoint( count, point );
        sparseVectorImage->SetPointData( count++, ItI.Get() );
        }
      }
    }

  return sparseVectorImage;
}

template <class TInputImage, class TOutputImage>
typename DiReCTImageFilter953<TInputImage, TOutputImage>::VectorImagePointer
DiReCTImageFilter953<TInputImage, TOutputImage>
::ConvertSparseVectorImageToVectorImage( const SparseVectorImageType *inputImage,
                                         const InputImageType *maskImage = NULL )
{
  VectorType zeroVector( 0.0 );

  VectorImagePointer vectorImage = VectorImageType::New();

  vectorImage->CopyInformation( this->GetSegmentationImage() );
  vectorImage->SetRegions( this->GetSegmentationImage()->GetRequestedRegion() );
  vectorImage->Allocate();
  vectorImage->FillBuffer( zeroVector );

  if( maskImage )
    {
    unsigned long count = 0;

    ImageRegionIterator<VectorImageType> ItI( vectorImage,
                                              vectorImage->GetRequestedRegion() );
    ImageRegionConstIterator<InputImageType> ItM( maskImage,
                                                  maskImage->GetRequestedRegion() );
    for( ItI.GoToBegin(), ItM.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItM )
      {
      if( ItM.Get() != 0 )
        {
        VectorType value( 0.0 );
        inputImage->GetPointData( count++, &value );
        ItI.Set( value );
        }
      }
    }
  else if( inputImage->GetPoints() != NULL &&
           inputImage->GetPointData() != NULL &&
           ( inputImage->GetPoints()->Size() == inputImage->GetPointData()->Size() ) )
    {
    typename SparseVectorImageType::PointsContainerConstIterator ItP =
      inputImage->GetPoints()->Begin();
    typename SparseVectorImageType::PointDataContainerIterator ItD =
      inputImage->GetPointData()->Begin();
    while( ItP != inputImage->GetPoints()->End() )
      {
      typename SparseVectorImageType::PointType point = ItP.Value();
      IndexType index;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        index[d] = static_cast<IndexValueType>( point[d] );
        }
      vectorImage->SetPixel( index, ItD.Value() );

      ++ItP;
      ++ItD;
      }
    }
  return vectorImage;
}

template <class TInputImage, class TOutputImage>
typename DiReCTImageFilter953<TInputImage, TOutputImage>::SparseImagePointer
DiReCTImageFilter953<TInputImage, TOutputImage>
::ConvertRealImageToSparseImage( const RealImageType *inputImage,
                                 const InputImageType *maskImage = NULL )
{
  SparseImagePointer sparseImage = SparseImageType::New();

  sparseImage->Initialize();

  unsigned long count = 0;

  if( maskImage )
    {
    ImageRegionConstIterator<RealImageType> ItI( inputImage,
                                                 inputImage->GetRequestedRegion() );
    ImageRegionConstIterator<InputImageType> ItM( maskImage,
                                                  maskImage->GetRequestedRegion() );
    for( ItI.GoToBegin(), ItM.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItM )
      {
      if( ItM.Get() != 0 )
        {
        sparseImage->SetPointData( count++, ItI.Get() );
        }
      }
    }
  else
    {
    ImageRegionConstIteratorWithIndex<RealImageType> ItI( inputImage,
                                                          inputImage->GetRequestedRegion() );
    for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
      {
      if( ItI.Get() > 1e-10 )
        {
        IndexType       index = ItI.GetIndex();
        SparseIndexType point;
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          point[d] =
            static_cast<typename SparseIndexType::CoordRepType>( index[d] );
          }
        sparseImage->SetPoint( count, point );
        sparseImage->SetPointData( count++, ItI.Get() );
        }
      }
    }

  return sparseImage;
}

template <class TInputImage, class TOutputImage>
typename DiReCTImageFilter953<TInputImage, TOutputImage>::RealImagePointer
DiReCTImageFilter953<TInputImage, TOutputImage>
::ConvertSparseImageToRealImage( const SparseImageType *inputImage,
                                 const InputImageType *maskImage = NULL )
{
  RealImagePointer realImage = RealImageType::New();

  realImage->CopyInformation( this->GetSegmentationImage() );
  realImage->SetRegions( this->GetSegmentationImage()->GetRequestedRegion() );
  realImage->Allocate();
  realImage->FillBuffer( 0.0 );

  if( maskImage )
    {
    unsigned long count = 0;

    ImageRegionIterator<RealImageType> ItI( realImage,
                                            realImage->GetRequestedRegion() );
    ImageRegionConstIterator<InputImageType> ItM( maskImage,
                                                  maskImage->GetRequestedRegion() );
    for( ItI.GoToBegin(), ItM.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItM )
      {
      if( ItM.Get() != 0 )
        {
        RealType value = 0.0;
        inputImage->GetPointData( count++, &value );
        ItI.Set( value );
        }
      }
    }
  else if( inputImage->GetPoints() != NULL &&
           inputImage->GetPointData() != NULL &&
           ( inputImage->GetPoints()->Size() == inputImage->GetPointData()->Size() ) )
    {
    typename SparseImageType::PointsContainerConstIterator ItP =
      inputImage->GetPoints()->Begin();
    typename SparseImageType::PointDataContainerIterator ItD =
      inputImage->GetPointData()->Begin();
    while( ItP != inputImage->GetPoints()->End() )
      {
      typename SparseImageType::PointType point = ItP.Value();
      IndexType index;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        index[d] = static_cast<IndexValueType>( point[d] );
        }
      realImage->SetPixel( index, ItD.Value() );

      ++ItP;
      ++ItD;
      }
    }

  return realImage;
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutputImage>
void
DiReCTImageFilter953<TInputImage, TOutputImage>
::PrintSelf( std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Gray matter label = "
     << this->m_GrayMatterLabel << std::endl;
  os << indent << "White matter label = "
     << this->m_WhiteMatterLabel << std::endl;
  os << indent << "Maximum number of iterations = "
     << this->m_MaximumNumberOfIterations << std::endl;
  os << indent << "Thickness prior estimate = "
     << this->m_ThicknessPriorEstimate << std::endl;
  os << indent << "Smoothing sigma = "
     << this->m_SmoothingSigma << std::endl;
  os << indent << "Gradient step = "
     << this->m_GradientStep << std::endl;
  os << indent << "Convergence threshold = "
     << this->m_ConvergenceThreshold << std::endl;
  os << indent << "Convergence window size = "
     << this->m_ConvergenceWindowSize << std::endl;
}
} // end namespace itk

#endif
