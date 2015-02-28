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
  this->m_UsePointSetData = true;
}

/** Destructor */
template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::~MeanSquaresPointSetToPointSetIntensityMetricv4()
{
}

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
typename MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::MeasureType
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::GetLocalNeighborhoodValue( const PointType & point, const PixelType & pixel ) const
{
  PointIdentifier pointId = this->m_MovingTransformedPointsLocator->FindClosestPoint( point );

  PixelType closestPixel;
  NumericTraits<PixelType>::SetLength( closestPixel, 1 );
  if( this->m_UsePointSetData )
    {
    bool doesPointDataExist = this->m_MovingPointSet->GetPointData( pointId, &closestPixel );
    if( ! doesPointDataExist )
      {
      itkExceptionMacro( "The corresponding data for point " << point << "(pointId = " << pointId << ") does not exist." );
      }
    }

  SizeValueType numberOfVoxelsInNeighborhood = pixel.size() / ( 1 + PointDimension );
  SizeValueType centerIntensityIndex = static_cast<SizeValueType>( 0.5 * numberOfVoxelsInNeighborhood )
    * ( PointDimension + 1 );

  const MeasureType measure = vnl_math_sqr( pixel[centerIntensityIndex] - closestPixel[centerIntensityIndex] );

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

  PixelType closestPixel;
  NumericTraits<PixelType>::SetLength( closestPixel, 1 );
  if( this->m_UsePointSetData )
    {
    bool doesPointDataExist = this->m_MovingPointSet->GetPointData( pointId, &closestPixel );
    if( ! doesPointDataExist )
      {
      itkExceptionMacro( "The corresponding data for point " << point << " (pointId = " << pointId << ") does not exist." );
      }
    }

  SizeValueType numberOfVoxelsInNeighborhood = pixel.size() / ( 1 + PointDimension );
  SizeValueType centerIntensityIndex = static_cast<SizeValueType>( 0.5 * numberOfVoxelsInNeighborhood )
    * ( PointDimension + 1 );

  measure = vnl_math_sqr( pixel[centerIntensityIndex] - closestPixel[centerIntensityIndex] );

  localDerivative.Fill( 0.0 );
  for( SizeValueType d = 0; d < PointDimension; d++ )
    {
    localDerivative[d] += pixel[centerIntensityIndex + 1 + d];
    }
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
