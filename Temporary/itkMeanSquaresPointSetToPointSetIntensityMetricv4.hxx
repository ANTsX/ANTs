/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkMeanSquaresPointSetToPointSetIntensityMetricv4_hxx
#define itkMeanSquaresPointSetToPointSetIntensityMetricv4_hxx

#include "itkMeanSquaresPointSetToPointSetIntensityMetricv4.h"

namespace itk
{

/** Constructor */
template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::MeanSquaresPointSetToPointSetIntensityMetricv4()
{
  this->SetUsePointSetData( true );
}

/** Destructor */
template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::~MeanSquaresPointSetToPointSetIntensityMetricv4()
{
}

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::InitializePointSets() const
{
  Superclass::InitializePointSets();

  // We need to also transform the neighborhood gradient vectors of
  // this->m_FixedPointSet and put them in this->m_FixedTransformedPointSet
  // We don't worry about doing this for the this->m_VirtualTransformedPointSet
  // as this point set consists of only the points without the corresponding
  // pixels.

  // Currently, we do not store the physical locations for the neighborhood voxels
  // in the pixel variable.  This presents a problem for transforming the gradient
  // at a neighborhood voxel location as we need the origin for a proper transformation.
  // So, for right now, we simply use the point (center of the neighborhood) as the
  // origin of each neighborhood vector.  We might need to change this.

  typedef CovariantVector<TInternalComputationValueType, Dimension> CovariantVectorType;

  typename FixedTransformType::InverseTransformBasePointer inverseTransform = this->m_FixedTransform->GetInverseTransform();

  typename FixedPointsContainer::ConstIterator It = this->m_FixedPointSet->GetPoints()->Begin();
  typename FixedPointsContainer::ConstIterator ItD = this->m_FixedPointSet->GetPointData()->Begin();

  SizeValueType numberOfVoxelsInNeighborhood = ( ItD.Value() ).GetSize() / ( 1 + Dimension );

  while( It != this->m_FixedPointSet->GetPoints()->End() )
    {
    PointType untransformedPoint = It.Value();
    PixelType neighborhoodIntensityAndGradientValues = ItD.Value();

    // transform the point into virtual space
    PointType transformPointInVirtualDomain = inverseTransform->TransformPoint( untransformedPoint );

    // Now transform each gradient vector to the space of the moving transform.
    for( SizeValueType i = 0; i < numberOfVoxelsInNeighborhood; ++i )
      {
      SizeValueType intensityIndex = i * ( ImageDimension + 1 );

      // First transform gradient vector to the space of the virtual domain
      CovariantVectorType untransformedGradient( 0.0 );
      for( SizeValueType d = 0; d < Dimension; d++ )
        {
        untransformedGradient[d] = neighborhoodIntensityAndGradientValues[intensityIndex + 1 + d];
        }

      CovariantVectorType transformedGradientInVirtualDomain =
        inverseTransform->TransformCovariantVector( untransformedGradient, untransformedPoint );

      // Second transform gradient vector to the space of the virtual domain

      CovariantVectorType transformedGradientInMovingDomain =
        this->m_MovingTransformCovariantVector( transformedGradientInVirtualDomain, transformPointInVirtualDomain );


      // Now repopulate the pixel neighborhoodIntensityAndGradientValues
      for( SizeValueType d = 0; d < Dimension; d++ )
        {
        neighborhoodIntensityAndGradientValues[intensityIndex + 1 + d] = transformedGradientInMovingDomain[d];
        }
      }

    this->m_FixedTransformedPointSet->SetPointData( It.Index(), neighborhoodIntensityAndGradientValues );

    ++It;
    ++ItD;
    }
}

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
typename MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::MeasureType
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::GetLocalNeighborhoodValue( const PointType & point, const PixelType & pixel ) const
{
  PointIdentifier pointId = this->m_MovingTransformedPointsLocator->FindClosestPoint( point );

  PointType closestPoint = this->m_MovingTransformedPointSet->GetPoint( pointId );
  PixelType closestPixel = this->m_MovingTransformedPointSet->GetPointData( pointId );

  SizeValueType numberOfVoxelsInNeighborhood = pixel.size() / ( 1 + Dimension );

  PixelType differenceIntensity( numberOfVoxelsInNeighborhood );
  for( SizeValueType i = 0; i < numberOfVoxelsInNeighborhood; ++i )
    {
    SizeValueType intensityIndex = i * ( ImageDimension + 1 );
    differenceIntensity[i] = pixel[intensityIndex] - closestPixel[intensityIndex];
    }
  const MeasureType measure = differenceIntensity.squared_magnitude();

  return measure;
}

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::GetLocalNeighborhoodValueAndDerivative( const PointType & point,
  MeasureType &measure, LocalDerivativeType & localDerivative, const PixelType & pixel ) const
{
  // Important note:  we assume that the gradients for each of the
  // neighborhood voxels that are located in the "pixel" variable
  // have been transformed according to the "fixed" transform.  This
  // is why we override the GetValue() and GetValueAndDerivative()
  // and IntializePointSets() functions.

  PointIdentifier pointId = this->m_MovingTransformedPointsLocator->FindClosestPoint( point );

  PointType closestPoint = this->m_MovingTransformedPointSet->GetPoint( pointId );
  PixelType closestPixel = this->m_MovingTransformedPointSet->GetPointData( pointId );

  SizeValueType numberOfVoxelsInNeighborhood = pixel.size() / ( 1 + Dimension );

  localDerivative.Fill( 0.0 );

  PixelType differenceIntensity( numberOfVoxelsInNeighborhood );
  for( SizeValueType i = 0; i < numberOfVoxelsInNeighborhood; ++i )
    {
    SizeValueType intensityIndex = i * ( ImageDimension + 1 );
    differenceIntensity[i] = pixel[intensityIndex] - closestPixel[intensityIndex];

    for( SizeValueType d = 0; d < Dimension; d++ )
      {
      localDerivative[d] += pixel[intensityIndex + 1 + d];
      }
    }

  measure = differenceIntensity.squared_magnitude();

  for( SizeValueType d = 0; d < Dimension; d++ )
    {
    localDerivative[d] /= static_cast<TInternalComputationValueType>( numberOfVoxelsInNeighborhood );
    }
}

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
typename MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>::MeasureType
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::GetValue() const
{
  this->InitializeForIteration();

  MeasureType value = 0.0;

  PointsConstIterator It = this->m_FixedTransformedPointSet->GetPoints()->Begin();
  // Virtual point set will be the same size as fixed point set as long as it's
  // generated from the fixed point set.
  if( this->m_VirtualTransformedPointSet->GetNumberOfPoints() != this->m_FixedTransformedPointSet->GetNumberOfPoints() )
    {
    itkExceptionMacro("Expected FixedTransformedPointSet to be the same size as VirtualTransformedPointSet.");
    }
  PointsConstIterator virtualIt = this->m_VirtualTransformedPointSet->GetPoints()->Begin();

  while( It != this->m_FixedTransformedPointSet->GetPoints()->End() )
    {
    /* Verify the virtual point is in the virtual domain.
     * If user hasn't defined a virtual space, and the active transform is not
     * a displacement field transform type, then this will always return true. */
    if( ! this->IsInsideVirtualDomain( virtualIt.Value() ) )
      {
      ++It;
      ++virtualIt;
      continue;
      }

    PixelType pixel();
    if( this->GetUsePointSetData() )
      {
      this->m_FixedPointSet->GetPointData( It.Index(), &pixel );
      }

    value += this->GetLocalNeighborhoodValue( It.Value(), pixel );
    ++virtualIt;
    ++It;
    }

  DerivativeType derivative;
  if( this->VerifyNumberOfValidPoints( value, derivative ) )
    {
    value /= static_cast<MeasureType>( this->m_NumberOfValidPoints );
    }
  this->m_Value = value;

  return value;
}

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::GetDerivative( DerivativeType & derivative ) const
{
  MeasureType value = NumericTraits<MeasureType>::ZeroValue();
  this->CalculateValueAndDerivative( value, derivative, false );
}

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::GetValueAndDerivative( MeasureType & value, DerivativeType & derivative ) const
{
  this->CalculateValueAndDerivative( value, derivative, true );
}

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::CalculateValueAndDerivative( MeasureType & value, DerivativeType & derivative, bool calculateValue ) const
{
  this->InitializeForIteration();
  derivative.SetSize( this->GetNumberOfParameters() );
  derivative.Fill( NumericTraits<DerivativeValueType>::ZeroValue() );

  value = NumericTraits<MeasureType>::ZeroValue();
  MovingTransformJacobianType  jacobian( MovingPointDimension, this->GetNumberOfLocalParameters() );
  MovingTransformJacobianType  jacobianPositional( MovingPointDimension, MovingPointDimension );

  DerivativeType localTransformDerivative( this->GetNumberOfLocalParameters() );
  localTransformDerivative.Fill( NumericTraits<DerivativeValueType>::ZeroValue() );

  // Virtual point set will be the same size as fixed point set as long as it's
  // generated from the fixed point set.
  if( this->m_VirtualTransformedPointSet->GetNumberOfPoints() != this->m_FixedTransformedPointSet->GetNumberOfPoints() )
    {
    itkExceptionMacro("Expected FixedTransformedPointSet to be the same size as VirtualTransformedPointSet.");
    }
  PointsConstIterator virtualIt = this->m_VirtualTransformedPointSet->GetPoints()->Begin();
  PointsConstIterator It = this->m_FixedTransformedPointSet->GetPoints()->Begin();
  PointsConstIterator end = this->m_FixedTransformedPointSet->GetPoints()->End();

  while( It != end )
    {
    MeasureType pointValue = NumericTraits<MeasureType>::ZeroValue();
    LocalDerivativeType pointDerivative;

    /* Verify the virtual point is in the virtual domain.
     * If user hasn't defined a virtual space, and the active transform is not
     * a displacement field transform type, then this will always return true. */
    if( ! this->IsInsideVirtualDomain( virtualIt.Value() ) )
      {
      ++It;
      ++virtualIt;
      continue;
      }

    PixelType pixel();
    if( this->GetUsePointSetData() )
      {
      this->m_FixedPointSet->GetPointData( It.Index(), &pixel );
      }

    if( calculateValue )
      {
      this->GetLocalNeighborhoodValueAndDerivative( It.Value(), pointValue, pointDerivative, pixel );
      value += pointValue;
      }
    else
      {
      pointDerivative = this->GetLocalNeighborhoodDerivative( It.Value(), pixel );
      }

    // Map into parameter space
    if( this->HasLocalSupport() )
      {
      // Reset to zero since we're not accumulating in the local-support case.
      localTransformDerivative.Fill( NumericTraits<DerivativeValueType>::ZeroValue() );
      }
    this->GetMovingTransform()->
      ComputeJacobianWithRespectToParametersCachedTemporaries(virtualIt.Value(),
                                                              jacobian,
                                                              jacobianPositional);
    for ( NumberOfParametersType par = 0; par < this->GetNumberOfLocalParameters(); par++ )
      {
      for( DimensionType d = 0; d < PointDimension; ++d )
        {
        localTransformDerivative[par] += jacobian(d, par) * pointDerivative[d];
        }
      }

    // For local-support transforms, store the per-point result
    if( this->HasLocalSupport() )
      {
      this->StorePointDerivative( virtualIt.Value(), localTransformDerivative, derivative );
      }

    ++It;
    ++virtualIt;
    }

  if( this->VerifyNumberOfValidPoints( value, derivative ) )
    {
    // For global-support transforms, average the accumulated derivative result
    if( ! this->HasLocalSupport() )
      {
      derivative = localTransformDerivative / static_cast<DerivativeValueType>(this->m_NumberOfValidPoints);
      }
    value /= static_cast<MeasureType>( this->m_NumberOfValidPoints );
    }
  this->m_Value = value;
}

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
typename LightObject::Pointer
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::InternalClone( void ) const
{
  typename Self::Pointer rval = Self::New();
  rval->SetMovingPointSet( this->m_MovingPointSet );
  rval->SetFixedPointSet( this->m_FixedPointSet );

  return rval.GetPointer();
}

/** PrintSelf method */
template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
}

} // end namespace itk

#endif
