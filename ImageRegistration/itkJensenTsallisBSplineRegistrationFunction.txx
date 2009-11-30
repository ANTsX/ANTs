/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkJensenTsallisBSplineRegistrationFunction.txx,v $
  Language:  C++
  Date:      $Date: 2008/11/15 23:46:06 $
  Version:   $Revision: 1.18 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkJensenTsallisBSplineRegistrationFunction_txx_
#define _itkJensenTsallisBSplineRegistrationFunction_txx_

#include "itkJensenTsallisBSplineRegistrationFunction.h"

#include "itkImageLinearConstIteratorWithIndex.h"
// #include "itkBSplineControlPointImageFilter.h"

namespace itk
{
template <class TFixedImage, class TFixedPointSet,
          class TMovingImage, class TMovingPointSet,
          class TDeformationField>
JensenTsallisBSplineRegistrationFunction<TFixedImage,
                                         TFixedPointSet, TMovingImage, TMovingPointSet, TDeformationField>
::JensenTsallisBSplineRegistrationFunction()
{
  this->m_UseRegularizationTerm = false;
  this->m_UseInputAsSamples = true;
  this->m_UseAnisotropicCovariances = false;

  this->m_NumberOfFixedSamples = 100;
  this->m_FixedPointSetSigma = 1.0;
  this->m_FixedKernelSigma = 0.0;
  this->m_FixedEvaluationKNeighborhood = 50;

  this->m_NumberOfMovingSamples = 100;
  this->m_MovingPointSetSigma = 1.0;
  this->m_MovingKernelSigma = 0.0;
  this->m_MovingEvaluationKNeighborhood = 50;

  unsigned int covarianceKNeighborhood = static_cast<unsigned int>( vcl_pow( 3.0,
                                                                             static_cast<RealType>( ImageDimension ) ) )
    - 1;

  this->m_FixedCovarianceKNeighborhood = covarianceKNeighborhood;
  this->m_MovingCovarianceKNeighborhood = covarianceKNeighborhood;

  this->m_Alpha = 2.0;

//  this->m_FixedControlPointLattice = NULL;
//  this->m_MovingControlPointLattice = NULL;

  this->m_DerivativeFixedField = NULL;
  this->m_DerivativeMovingField = NULL;
  this->m_IsPointSetMetric = true;

  this->m_SplineOrder = 3;
  this->m_NumberOfLevels = 1;
  this->m_MeshResolution.Fill( 1 );
}

template <class TFixedImage, class TFixedPointSet,
          class TMovingImage, class TMovingPointSet,
          class TDeformationField>
void
JensenTsallisBSplineRegistrationFunction<TFixedImage,
                                         TFixedPointSet, TMovingImage, TMovingPointSet, TDeformationField>
::InitializeIteration( void )
{
  if( this->m_FixedKernelSigma == 0 )
    {
    double maxFixedSpacing = 0.0;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      if( this->GetFixedImage()->GetSpacing()[d] )
        {
        maxFixedSpacing = this->GetFixedImage()->GetSpacing()[d];
        }
      }
    this->m_FixedKernelSigma = 2.0 * maxFixedSpacing;
    }

  if( this->m_MovingKernelSigma == 0 )
    {
    double maxMovingSpacing = 0.0;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      if( this->GetMovingImage()->GetSpacing()[d] )
        {
        maxMovingSpacing = this->GetMovingImage()->GetSpacing()[d];
        }
      }
    this->m_MovingKernelSigma = 2.0 * maxMovingSpacing;
    }

  typename PointSetMetricType::Pointer pointSetMetric
    = PointSetMetricType::New();
  pointSetMetric->SetFixedPointSet( this->m_FixedPointSet );
  pointSetMetric->SetMovingPointSet( this->m_MovingPointSet );

  pointSetMetric->SetUseRegularizationTerm( this->m_UseRegularizationTerm );
  pointSetMetric->SetUseInputAsSamples( this->m_UseInputAsSamples );
  pointSetMetric->SetUseAnisotropicCovariances( this->m_UseAnisotropicCovariances );
  pointSetMetric->SetNumberOfFixedSamples( this->m_NumberOfFixedSamples );
  pointSetMetric->SetFixedPointSetSigma( this->m_FixedPointSetSigma );
  pointSetMetric->SetFixedKernelSigma( this->m_FixedKernelSigma );
  pointSetMetric->SetFixedCovarianceKNeighborhood( this->m_FixedCovarianceKNeighborhood );
  pointSetMetric->SetFixedEvaluationKNeighborhood( this->m_FixedEvaluationKNeighborhood );
  pointSetMetric->SetNumberOfMovingSamples( this->m_NumberOfMovingSamples );
  pointSetMetric->SetMovingPointSetSigma( this->m_MovingPointSetSigma );
  pointSetMetric->SetMovingKernelSigma( this->m_MovingKernelSigma );
  pointSetMetric->SetMovingCovarianceKNeighborhood( this->m_MovingCovarianceKNeighborhood );
  pointSetMetric->SetMovingEvaluationKNeighborhood( this->m_MovingEvaluationKNeighborhood );
  pointSetMetric->SetAlpha( this->m_Alpha );

  pointSetMetric->Initialize();

  typename PointSetMetricType::DefaultTransformType::ParametersType parameters;
  parameters.Fill( 0.0 );

  /**
   * Calculate with respect to the moving point set
   */

  pointSetMetric->SetUseWithRespectToTheMovingPointSet( true );
  typename PointSetMetricType::DerivativeType movingGradient;
  typename PointSetMetricType::MeasureType movingMeasure;

  pointSetMetric->GetValueAndDerivative( parameters, movingMeasure,
                                         movingGradient );
  this->m_Energy += movingMeasure[0];

  typename BSplinePointSetType::Pointer movingGradientPoints =
    BSplinePointSetType::New();
  movingGradientPoints->Initialize();

  typename BSplineWeightsType::Pointer movingWeights =
    BSplineWeightsType::New();
  movingWeights->Initialize();

  unsigned long count = 0;
  for( unsigned int n = 0; n < pointSetMetric->GetNumberOfValues(); n++ )
    {
    typename MovingPointSetType::PointType point;
    this->m_MovingPointSet->GetPoint( n, &point );

    typename BSplinePointSetType::PointType bsplinePoint;
    VectorType gradient;
    bsplinePoint.CastFrom( point );

    bool isInside = true;
    for( unsigned int d = 0; d < PointDimension; d++ )
      {
      if( bsplinePoint[d] <= this->GetMovingImage()->GetOrigin()[d] ||
          bsplinePoint[d] >= this->GetMovingImage()->GetOrigin()[d]
          + ( this->GetMovingImage()->GetLargestPossibleRegion().GetSize()[d] - 1 )
          * this->GetMovingImage()->GetSpacing()[d] )
        {
        isInside = false;
        break;
        }
      }

    if( isInside )
      {
      for( unsigned int d = 0; d < PointDimension; d++ )
        {
        gradient[d] = movingGradient(n, d);
        }
      movingGradientPoints->SetPoint( count, bsplinePoint );
      movingGradientPoints->SetPointData( count, gradient );
      movingWeights->InsertElement( count, 1.0 );
      count++;
      }
    }

  VectorType zeroVector;
  zeroVector.Fill( 0.0 );

  ImageLinearConstIteratorWithIndex<MovingImageType> ItM(
    this->GetMovingImage(),
    this->GetMovingImage()->GetLargestPossibleRegion() );
  for( unsigned int d = 0; d < PointDimension; d++ )
    {
    ItM.SetDirection( d );
    ItM.GoToBegin();
    while( !ItM.IsAtEnd() )
      {
      typename MovingImageType::PointType point;
      typename BSplinePointSetType::PointType bsplinePoint;

      ItM.GoToBeginOfLine();
      this->GetMovingImage()->TransformIndexToPhysicalPoint(
        ItM.GetIndex(), point );
      bsplinePoint.CastFrom( point );

      movingGradientPoints->SetPoint( count, bsplinePoint );
      movingGradientPoints->SetPointData( count, zeroVector );
      count++;

      ItM.GoToEndOfLine();
      --ItM;
      this->GetMovingImage()->TransformIndexToPhysicalPoint(
        ItM.GetIndex(), point );
      bsplinePoint.CastFrom( point );

      movingGradientPoints->SetPoint( count, bsplinePoint );
      movingGradientPoints->SetPointData( count, zeroVector );
      movingWeights->InsertElement( count, 1000.0 );
      count++;

      ItM.NextLine();
      }
    }

  ArrayType numberOfMovingControlPoints;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    numberOfMovingControlPoints[d] = this->m_SplineOrder
      + static_cast<unsigned int>( vcl_floor( 0.5 + static_cast<RealType>(
                                                this->GetMovingImage()->GetLargestPossibleRegion().GetSize()[d] )
                                              / static_cast<RealType>( this->m_MeshResolution[d] ) ) );
    }

  typename BSplineFilterType::Pointer movingBSpliner
    = BSplineFilterType::New();
  movingBSpliner->SetInput( movingGradientPoints );
  movingBSpliner->SetPointWeights( movingWeights.GetPointer() );
  movingBSpliner->SetOrigin( this->GetMovingImage()->GetOrigin() );
  movingBSpliner->SetSpacing( this->GetMovingImage()->GetSpacing() );
  movingBSpliner->SetSize(
    this->GetMovingImage()->GetLargestPossibleRegion().GetSize() );
  movingBSpliner->SetNumberOfLevels( this->m_NumberOfLevels );
  movingBSpliner->SetSplineOrder( this->m_SplineOrder );
  movingBSpliner->SetNumberOfControlPoints( numberOfMovingControlPoints );
  movingBSpliner->SetGenerateOutputImage( true );
  movingBSpliner->Update();

  this->m_DerivativeMovingField = movingBSpliner->GetOutput();
//  this->m_MovingControlPointLattice = movingBSpliner->GetPhiLattice();

  /**
   * Calculate with respect to the fixed point set
   */

  pointSetMetric->SetUseWithRespectToTheMovingPointSet( false );
  typename PointSetMetricType::DerivativeType fixedGradient;
  typename PointSetMetricType::MeasureType fixedMeasure;
  pointSetMetric->GetValueAndDerivative( parameters, fixedMeasure,
                                         fixedGradient );
  this->m_Energy += fixedMeasure[0];

  typename BSplinePointSetType::Pointer fixedGradientPoints =
    BSplinePointSetType::New();
  fixedGradientPoints->Initialize();

  typename BSplineWeightsType::Pointer fixedWeights =
    BSplineWeightsType::New();
  fixedWeights->Initialize();

  count = 0;
  for( unsigned int n = 0; n < pointSetMetric->GetNumberOfValues(); n++ )
    {
    typename FixedPointSetType::PointType point;
    this->m_FixedPointSet->GetPoint( n, &point );

    typename BSplinePointSetType::PointType bsplinePoint;
    VectorType gradient;
    bsplinePoint.CastFrom( point );

    bool isInside = true;
    for( unsigned int d = 0; d < PointDimension; d++ )
      {
      if( bsplinePoint[d] <= this->GetFixedImage()->GetOrigin()[d] ||
          bsplinePoint[d] >= this->GetFixedImage()->GetOrigin()[d]
          + ( this->GetFixedImage()->GetLargestPossibleRegion().GetSize()[d] - 1 )
          * this->GetFixedImage()->GetSpacing()[d] )
        {
        isInside = false;
        break;
        }
      }

    if( isInside )
      {
      for( unsigned int d = 0; d < PointDimension; d++ )
        {
        gradient[d] = fixedGradient(n, d);
        }

      fixedGradientPoints->SetPoint( count, bsplinePoint );
      fixedGradientPoints->SetPointData( count, gradient );
      fixedWeights->InsertElement( count, 1.0 );
      count++;
      }
    }

  ImageLinearConstIteratorWithIndex<FixedImageType> ItF(
    this->GetFixedImage(),
    this->GetFixedImage()->GetLargestPossibleRegion() );
  for( unsigned int d = 0; d < PointDimension; d++ )
    {
    ItF.SetDirection( d );
    ItF.GoToBegin();
    while( !ItF.IsAtEnd() )
      {
      typename FixedImageType::PointType point;
      typename BSplinePointSetType::PointType bsplinePoint;

      ItF.GoToBeginOfLine();
      this->GetFixedImage()->TransformIndexToPhysicalPoint(
        ItF.GetIndex(), point );
      bsplinePoint.CastFrom( point );

      fixedGradientPoints->SetPoint( count, bsplinePoint );
      fixedGradientPoints->SetPointData( count, zeroVector );
      fixedWeights->InsertElement( count, 1.0 );
      count++;

      ItF.GoToEndOfLine();
      --ItF;
      this->GetFixedImage()->TransformIndexToPhysicalPoint(
        ItF.GetIndex(), point );
      bsplinePoint.CastFrom( point );

      fixedGradientPoints->SetPoint( count, bsplinePoint );
      fixedGradientPoints->SetPointData( count, zeroVector );
      fixedWeights->InsertElement( count, 1000.0 );
      count++;

      ItF.NextLine();
      }
    }

  ArrayType numberOfFixedControlPoints;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    numberOfFixedControlPoints[d] = this->m_SplineOrder
      + static_cast<unsigned int>( vcl_floor( 0.5 + static_cast<RealType>(
                                                this->GetFixedImage()->GetLargestPossibleRegion().GetSize()[d] )
                                              / static_cast<RealType>( this->m_MeshResolution[d] ) ) );
    }

  typename BSplineFilterType::Pointer fixedBSpliner
    = BSplineFilterType::New();
  fixedBSpliner->SetInput( fixedGradientPoints );
  fixedBSpliner->SetPointWeights( fixedWeights.GetPointer() );
  fixedBSpliner->SetOrigin( this->GetFixedImage()->GetOrigin() );
  fixedBSpliner->SetSpacing( this->GetFixedImage()->GetSpacing() );
  fixedBSpliner->SetSize(
    this->GetFixedImage()->GetLargestPossibleRegion().GetSize() );
  fixedBSpliner->SetNumberOfLevels( this->m_NumberOfLevels );
  fixedBSpliner->SetSplineOrder( this->m_SplineOrder );
  fixedBSpliner->SetNumberOfControlPoints( numberOfFixedControlPoints );
  fixedBSpliner->SetGenerateOutputImage( true );
  fixedBSpliner->Update();

  this->m_DerivativeFixedField = fixedBSpliner->GetOutput();
//  this->m_FixedControlPointLattice = fixedBSpliner->GetPhiLattice();
}

template <class TFixedImage, class TFixedPointSet,
          class TMovingImage, class TMovingPointSet,
          class TDeformationField>
typename JensenTsallisBSplineRegistrationFunction<TFixedImage,
                                                  TFixedPointSet, TMovingImage, TMovingPointSet,
                                                  TDeformationField>::VectorType
JensenTsallisBSplineRegistrationFunction<TFixedImage,
                                         TFixedPointSet, TMovingImage, TMovingPointSet, TDeformationField>
::ComputeUpdate( const NeighborhoodType & neighborhood,
                 void *globalData, const FloatOffsetType & offset )
{
  if( this->m_DerivativeFixedField )
    {
    return -this->m_DerivativeFixedField->GetPixel( neighborhood.GetIndex() );
    }
  else
    {
    itkExceptionMacro( "Initialize() has not been called." );
    }

/*
  typedef BSplineControlPointImageFilter<ControlPointLatticeType,
    DeformationFieldType> BSplineControlPointImageFilterType;

  typename BSplineControlPointImageFilterType::Pointer movingBSpliner
    = BSplineControlPointImageFilterType::New();
  movingBSpliner->SetInput( this->m_MovingControlPointLattice );
  movingBSpliner->SetOrigin( this->GetMovingImage()->GetOrigin() );
  movingBSpliner->SetSpacing( this->GetMovingImage()->GetSpacing() );
  movingBSpliner->SetSize(
    this->GetMovingImage()->GetLargestPossibleRegion().GetSize() );
  movingBSpliner->SetSplineOrder( this->m_SplineOrder );

  VectorType gradient;
  movingBSpliner->EvaluateAtIndex( neighborhood.GetIndex(),
    gradient );

  return gradient;
*/
}

template <class TFixedImage, class TFixedPointSet,
          class TMovingImage, class TMovingPointSet,
          class TDeformationField>
typename JensenTsallisBSplineRegistrationFunction<TFixedImage,
                                                  TFixedPointSet, TMovingImage, TMovingPointSet,
                                                  TDeformationField>::VectorType
JensenTsallisBSplineRegistrationFunction<TFixedImage,
                                         TFixedPointSet, TMovingImage, TMovingPointSet, TDeformationField>
::ComputeUpdateInv( const NeighborhoodType & neighborhood,
                    void *globalData, const FloatOffsetType & offset )
{
  if( this->m_DerivativeMovingField )
    {
    return -this->m_DerivativeMovingField->GetPixel( neighborhood.GetIndex() );
    }
  else
    {
    itkExceptionMacro( "Initialize() has not been called." );
    }

/*
  typedef BSplineControlPointImageFilter<ControlPointLatticeType,
    DeformationFieldType> BSplineControlPointImageFilterType;

  typename BSplineControlPointImageFilterType::Pointer fixedBSpliner
    = BSplineControlPointImageFilterType::New();
  fixedBSpliner->SetInput( this->m_FixedControlPointLattice );
  fixedBSpliner->SetOrigin( this->GetFixedImage()->GetOrigin() );
  fixedBSpliner->SetSpacing( this->GetFixedImage()->GetSpacing() );
  fixedBSpliner->SetSize(
    this->GetFixedImage()->GetLargestPossibleRegion().GetSize() );
  fixedBSpliner->SetSplineOrder( this->m_SplineOrder );

  VectorType gradient;
  fixedBSpliner->EvaluateAtIndex( neighborhood.GetIndex(),
    gradient );

  return gradient;
*/
}

template <class TFixedImage, class TFixedPointSet,
          class TMovingImage, class TMovingPointSet,
          class TDeformationField>
void
JensenTsallisBSplineRegistrationFunction<TFixedImage,
                                         TFixedPointSet, TMovingImage, TMovingPointSet, TDeformationField>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Use regularization term: "
     << this->m_UseRegularizationTerm << std::endl;
  if( !this->m_UseInputAsSamples )
    {
    os << indent << "Number of fixed samples: "
       << this->m_NumberOfFixedSamples << std::endl;
    os << indent << "Number of moving samples: "
       << this->m_NumberOfMovingSamples << std::endl;
    }
  os << indent << "Alpha: "
     << this->m_Alpha << std::endl;
  os << indent << "Fixed sigma: "
     << this->m_FixedPointSetSigma << std::endl;
  os << indent << "Moving sigma: "
     << this->m_MovingPointSetSigma << std::endl;

  os << indent << "Spline order: "
     << this->m_SplineOrder << std::endl;
  os << indent << "Number of levels: "
     << this->m_NumberOfLevels << std::endl;
  os << indent << "Number of control points: "
     << this->m_MeshResolution << std::endl;
}
} // end namespace itk

#endif
