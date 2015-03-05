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
  this->m_EuclideanDistanceSigma = std::sqrt( 5.0 );
  this->m_IntensityDistanceSigma = std::sqrt( 5.0 );
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
      itkExceptionMacro( "The corresponding data for point " << point << " (pointId = " << pointId << ") does not exist." );
      }
    }

  PointType closestPoint;
  closestPoint.Fill( 0.0 );
  closestPoint = this->m_MovingTransformedPointSet->GetPoint( pointId );

  // the probabilistic icp term
  const MeasureType euclideanDistance = point.EuclideanDistanceTo( closestPoint );
  MeasureType distanceProbability = std::exp( -0.5 * vnl_math_sqr( euclideanDistance / this->m_EuclideanDistanceSigma ) );

  SizeValueType numberOfVoxelsInNeighborhood = pixel.size() / ( 1 + PointDimension );
  SizeValueType centerIntensityIndex = static_cast<SizeValueType>( 0.5 * numberOfVoxelsInNeighborhood )
    * ( PointDimension + 1 );

  // the probabilistic intensity term
  MeasureType intensityDistance = pixel[centerIntensityIndex] - closestPixel[centerIntensityIndex];
  MeasureType intensityProbability = std::exp( -0.5 * vnl_math_sqr( intensityDistance / this->m_IntensityDistanceSigma ) );

  const MeasureType measure = ( -1.0 ) * intensityProbability * distanceProbability;

  return measure;
}

template<typename TFixedPointSet, typename TMovingPointSet, class TInternalComputationValueType>
void
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>
::GetLocalNeighborhoodValueAndDerivative( const PointType & point,
  MeasureType &measure, LocalDerivativeType & localDerivative, const PixelType & pixel ) const
{
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

  PointType closestPoint;
  closestPoint.Fill( 0.0 );
  closestPoint = this->m_MovingTransformedPointSet->GetPoint( pointId );

  // the probabilistic icp term
  const MeasureType euclideanDistance = point.EuclideanDistanceTo( closestPoint );
  MeasureType distanceProbability = std::exp( -0.5 * vnl_math_sqr( euclideanDistance / this->m_EuclideanDistanceSigma ) );

  SizeValueType numberOfVoxelsInNeighborhood = pixel.size() / ( 1 + PointDimension );
  SizeValueType centerIntensityIndex =
     static_cast<SizeValueType>( 0.5 * numberOfVoxelsInNeighborhood ) * ( PointDimension + 1 );

  // the probabilistic intensity term
  MeasureType intensityDistance = pixel[centerIntensityIndex] - closestPixel[centerIntensityIndex];
  MeasureType intensityProbability = std::exp( -0.5 * vnl_math_sqr( intensityDistance / this->m_IntensityDistanceSigma ) );

  measure = ( -1.0 ) * intensityProbability * distanceProbability;

  // total derivative is
  // d/dx( intProb * distProb ) =
  //   intProb * d/dx( distProb ) + distProb * d/dx( intProb ) =
  //   intProb * distProb * dist + distProb * intProb * intdiff =

  localDerivative = ( closestPoint - point ) * intensityProbability * distanceProbability;
  for( SizeValueType d = 0; d < PointDimension; d++ )
    {
    localDerivative[d] += intensityProbability * distanceProbability * intensityDistance *
      closestPixel[centerIntensityIndex + 1 + d];
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

  os << "Euclidean distance sigma = " << this->m_EuclideanDistanceSigma << std::endl;
  os << "intensity distance sigma = " << this->m_IntensityDistanceSigma << std::endl;
}

} // end namespace itk

#endif
