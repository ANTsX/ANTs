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


namespace itk
{

/** Constructor */
template <typename TFixedPointSet, typename TMovingPointSet, typename TInternalComputationValueType>
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>::
  MeanSquaresPointSetToPointSetIntensityMetricv4()
{
  this->m_EuclideanDistanceSigma = std::sqrt(5.0);
  this->m_IntensityDistanceSigma = std::sqrt(5.0);

  this->m_EstimateIntensityDistanceSigmaAutomatically = true;
  this->m_EstimateEuclideanDistanceSigmaAutomatically = true;

  this->m_UsePointSetData = true;
}

/** Destructor */
template <typename TFixedPointSet, typename TMovingPointSet, typename TInternalComputationValueType>
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>::
  ~MeanSquaresPointSetToPointSetIntensityMetricv4() = default;

/** Initialize the metric */
template <typename TFixedPointSet, typename TMovingPointSet, typename TInternalComputationValueType>
void
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>::
  Initialize(void)
{
  Superclass::Initialize();

  if (this->m_EstimateIntensityDistanceSigmaAutomatically)
  {
    this->EstimateIntensityDistanceSigma();
  }
  if (this->m_EstimateEuclideanDistanceSigmaAutomatically)
  {
    this->EstimateEuclideanDistanceSigma();
  }
}

template <typename TFixedPointSet, typename TMovingPointSet, typename TInternalComputationValueType>
void
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>::
  InitializePointSets() const
{
  Superclass::InitializePointSets();

  if (this->m_CalculateValueAndDerivativeInTangentSpace == true)
  {
    this->TransformMovingPointSetGradients();
    this->TransformFixedPointSetGradients();
  }
}

template <typename TFixedPointSet, typename TMovingPointSet, typename TInternalComputationValueType>
void
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>::
  TransformFixedPointSetGradients() const
{
  typename FixedTransformType::InverseTransformBasePointer inverseTransform =
    this->m_FixedTransform->GetInverseTransform();

  typename FixedPointsContainer::ConstIterator It = this->m_FixedPointSet->GetPoints()->Begin();

  while (It != this->m_FixedPointSet->GetPoints()->End())
  {
    PixelType pixel;
    NumericTraits<PixelType>::SetLength(pixel, 1);
    bool doesPointDataExist = this->m_FixedPointSet->GetPointData(It.Index(), &pixel);
    if (!doesPointDataExist)
    {
      itkExceptionMacro("The corresponding data for point " << It.Value() << " (pointId = " << It.Index()
                                                            << ") does not exist.");
    }
    SizeValueType numberOfVoxelsInNeighborhood = pixel.size() / (1 + PointDimension);

    // Here we assume that transforming the vector at the neighborhood voxel
    // is close to performing the transformation at the center voxel.
    for (SizeValueType n = 0; n < numberOfVoxelsInNeighborhood; n++)
    {
      CovariantVectorType covariantVector;

      for (unsigned int d = 0; d < PointDimension; d++)
      {
        covariantVector[d] = pixel[n * (PointDimension + 1) + d + 1];
      }

      // First, transform from fixed to virtual domain.  Then go from virtual domain
      // to moving.

      CovariantVectorType transformedCovariantVector =
        inverseTransform->TransformCovariantVector(covariantVector, It.Value());

      for (unsigned int d = 0; d < PointDimension; d++)
      {
        pixel[n * (PointDimension + 1) + d + 1] = transformedCovariantVector[d];
      }
    }

    this->m_FixedTransformedPointSet->SetPointData(It.Index(), pixel);
    ++It;
  }
}

template <typename TFixedPointSet, typename TMovingPointSet, typename TInternalComputationValueType>
void
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>::
  TransformMovingPointSetGradients() const
{
  typename MovingTransformType::InverseTransformBasePointer inverseTransform =
    this->m_MovingTransform->GetInverseTransform();

  typename MovingPointsContainer::ConstIterator It = this->m_MovingPointSet->GetPoints()->Begin();

  while (It != this->m_MovingPointSet->GetPoints()->End())
  {
    PixelType pixel;
    NumericTraits<PixelType>::SetLength(pixel, 1);
    bool doesPointDataExist = this->m_MovingPointSet->GetPointData(It.Index(), &pixel);
    if (!doesPointDataExist)
    {
      itkExceptionMacro("The corresponding data for point " << It.Value() << " (pointId = " << It.Index()
                                                            << ") does not exist.");
    }
    SizeValueType numberOfVoxelsInNeighborhood = pixel.size() / (1 + PointDimension);

    // Here we assume that transforming the vector at the neighborhood voxel
    // is close to performing the transformation at the center voxel.
    for (SizeValueType n = 0; n < numberOfVoxelsInNeighborhood; n++)
    {
      CovariantVectorType covariantVector;

      for (unsigned int d = 0; d < PointDimension; d++)
      {
        covariantVector[d] = pixel[n * (PointDimension + 1) + d + 1];
      }

      // First, transform from fixed to virtual domain.  Then go from virtual domain
      // to moving.

      CovariantVectorType transformedCovariantVector =
        inverseTransform->TransformCovariantVector(covariantVector, It.Value());

      for (unsigned int d = 0; d < PointDimension; d++)
      {
        pixel[n * (PointDimension + 1) + d + 1] = transformedCovariantVector[d];
      }
    }
    this->m_MovingTransformedPointSet->SetPointData(It.Index(), pixel);
    ++It;
  }
}

template <typename TFixedPointSet, typename TMovingPointSet, typename TInternalComputationValueType>
void
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>::
  EstimateEuclideanDistanceSigma()
{
  TInternalComputationValueType runningDistanceMean = 0.0;
  TInternalComputationValueType runningDistanceSigma = 0.0;

  if (this->m_FixedTransformedPointSet->GetNumberOfPoints() <= 1)
  {
    itkExceptionMacro("Need more than 1 point to estimate the distance sigma.");
  }

  unsigned int count = 0;

  PointsConstIterator ItF = this->m_FixedTransformedPointSet->GetPoints()->Begin();
  while (ItF != this->m_FixedTransformedPointSet->GetPoints()->End())
  {
    PointIdentifier pointId = this->m_MovingTransformedPointsLocator->FindClosestPoint(ItF.Value());

    PointType closestPoint;
    closestPoint.Fill(0.0);
    closestPoint = this->m_MovingTransformedPointSet->GetPoint(pointId);

    TInternalComputationValueType distance = closestPoint.EuclideanDistanceTo(ItF.Value());

    if (count == 0)
    {
      runningDistanceMean = distance;
      runningDistanceSigma = 0.0;
    }
    else
    {
      TInternalComputationValueType runningDistanceMeanPreviousIteration = runningDistanceMean;
      runningDistanceMean =
        runningDistanceMeanPreviousIteration +
        (distance - runningDistanceMeanPreviousIteration) / static_cast<TInternalComputationValueType>(count + 1);
      runningDistanceSigma += (distance - runningDistanceMeanPreviousIteration) * (distance - runningDistanceMean);
    }
    ++count;

    ++ItF;
  }
  this->m_EuclideanDistanceSigma = std::sqrt(runningDistanceSigma / static_cast<TInternalComputationValueType>(count));
}

template <typename TFixedPointSet, typename TMovingPointSet, typename TInternalComputationValueType>
void
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>::
  EstimateIntensityDistanceSigma()
{
  // Get the min/max intensities from the fixed point set

  PointsConstIterator ItF = this->m_FixedPointSet->GetPoints()->Begin();

  TInternalComputationValueType maxFixedIntensity = NumericTraits<TInternalComputationValueType>::NonpositiveMin();
  TInternalComputationValueType minFixedIntensity = NumericTraits<TInternalComputationValueType>::max();

  while (ItF != this->m_FixedPointSet->GetPoints()->End())
  {
    PixelType pixel;
    NumericTraits<PixelType>::SetLength(pixel, 1);
    if (this->m_UsePointSetData)
    {
      bool doesPointDataExist = this->m_FixedPointSet->GetPointData(ItF.Index(), &pixel);
      if (!doesPointDataExist)
      {
        itkExceptionMacro("The corresponding data for point " << ItF.Value() << " (pointId = " << ItF.Index()
                                                              << ") does not exist.");
      }
    }

    SizeValueType numberOfVoxelsInNeighborhood = pixel.size() / (1 + PointDimension);
    SizeValueType centerIntensityIndex =
      static_cast<SizeValueType>(0.5 * numberOfVoxelsInNeighborhood) * (PointDimension + 1);

    if (pixel[centerIntensityIndex] > maxFixedIntensity)
    {
      maxFixedIntensity = pixel[centerIntensityIndex];
    }
    else if (pixel[centerIntensityIndex] > minFixedIntensity)
    {
      minFixedIntensity = pixel[centerIntensityIndex];
    }

    ++ItF;
  }

  // Get the min/max intensities from the moving point set

  PointsConstIterator ItM = this->m_MovingPointSet->GetPoints()->Begin();

  TInternalComputationValueType maxMovingIntensity = NumericTraits<TInternalComputationValueType>::NonpositiveMin();
  TInternalComputationValueType minMovingIntensity = NumericTraits<TInternalComputationValueType>::max();

  while (ItM != this->m_MovingPointSet->GetPoints()->End())
  {
    PixelType pixel;
    NumericTraits<PixelType>::SetLength(pixel, 1);
    if (this->m_UsePointSetData)
    {
      bool doesPointDataExist = this->m_MovingPointSet->GetPointData(ItM.Index(), &pixel);
      if (!doesPointDataExist)
      {
        itkExceptionMacro("The corresponding data for point " << ItM.Value() << " (pointId = " << ItM.Index()
                                                              << ") does not exist.");
      }
    }

    SizeValueType numberOfVoxelsInNeighborhood = pixel.size() / (1 + PointDimension);
    SizeValueType centerIntensityIndex =
      static_cast<SizeValueType>(0.5 * numberOfVoxelsInNeighborhood) * (PointDimension + 1);

    if (pixel[centerIntensityIndex] > maxMovingIntensity)
    {
      maxMovingIntensity = pixel[centerIntensityIndex];
    }
    else if (pixel[centerIntensityIndex] > minMovingIntensity)
    {
      minMovingIntensity = pixel[centerIntensityIndex];
    }

    ++ItM;
  }

  // Now determine the sigma using a reasonable heuristic.

  this->m_IntensityDistanceSigma =
    std::max(maxMovingIntensity, maxFixedIntensity) - std::min(minMovingIntensity, maxMovingIntensity);
  if (Math::FloatAlmostEqual(this->m_IntensityDistanceSigma, NumericTraits<TInternalComputationValueType>::ZeroValue()))
  {
    this->m_IntensityDistanceSigma = std::max(maxMovingIntensity, maxFixedIntensity);
  }
}

template <typename TFixedPointSet, typename TMovingPointSet, typename TInternalComputationValueType>
typename MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet,
                                                        TMovingPointSet,
                                                        TInternalComputationValueType>::MeasureType
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>::
  GetLocalNeighborhoodValue(const PointType & point, const PixelType & pixel) const
{
  PointIdentifier pointId = this->m_MovingTransformedPointsLocator->FindClosestPoint(point);

  PixelType closestPixel;
  NumericTraits<PixelType>::SetLength(closestPixel, 1);
  if (this->m_UsePointSetData)
  {
    bool doesPointDataExist = false;
    if (this->m_CalculateValueAndDerivativeInTangentSpace == true)
    {
      doesPointDataExist = this->m_MovingTransformedPointSet->GetPointData(pointId, &closestPixel);
    }
    else
    {
      doesPointDataExist = this->m_MovingPointSet->GetPointData(pointId, &closestPixel);
    }
    if (!doesPointDataExist)
    {
      itkExceptionMacro("The corresponding data for point " << point << " (pointId = " << pointId
                                                            << ") does not exist.");
    }
  }

  PointType closestPoint;
  closestPoint.Fill(0.0);
  closestPoint = this->m_MovingTransformedPointSet->GetPoint(pointId);

  // the probabilistic icp term
  const MeasureType euclideanDistance = point.EuclideanDistanceTo(closestPoint);
  MeasureType       distanceProbability =
    std::exp(static_cast<MeasureType>(-0.5) * itk::Math::sqr(euclideanDistance / this->m_EuclideanDistanceSigma));

  SizeValueType numberOfVoxelsInNeighborhood = pixel.size() / (1 + PointDimension);
  SizeValueType centerIntensityIndex =
    static_cast<SizeValueType>(static_cast<MeasureType>(-0.5) * numberOfVoxelsInNeighborhood) * (PointDimension + 1);

  // the probabilistic intensity term
  MeasureType intensityDistance = pixel[centerIntensityIndex] - closestPixel[centerIntensityIndex];
  MeasureType intensityProbability =
    std::exp(static_cast<MeasureType>(-0.5) * itk::Math::sqr(intensityDistance / this->m_IntensityDistanceSigma));

  const MeasureType measure = -itk::NumericTraits<MeasureType>::OneValue() * intensityProbability * distanceProbability;

  return measure;
}

template <typename TFixedPointSet, typename TMovingPointSet, typename TInternalComputationValueType>
void
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>::
  GetLocalNeighborhoodValueAndDerivative(const PointType &     point,
                                         MeasureType &         measure,
                                         LocalDerivativeType & localDerivative,
                                         const PixelType &     pixel) const
{
  PointIdentifier pointId = this->m_MovingTransformedPointsLocator->FindClosestPoint(point);

  PixelType closestPixel;
  NumericTraits<PixelType>::SetLength(closestPixel, 1);
  if (this->m_UsePointSetData)
  {
    bool doesPointDataExist = false;
    if (this->m_CalculateValueAndDerivativeInTangentSpace == true)
    {
      doesPointDataExist = this->m_MovingTransformedPointSet->GetPointData(pointId, &closestPixel);
    }
    else
    {
      doesPointDataExist = this->m_MovingPointSet->GetPointData(pointId, &closestPixel);
    }
    if (!doesPointDataExist)
    {
      itkExceptionMacro("The corresponding data for point " << point << " (pointId = " << pointId
                                                            << ") does not exist.");
    }
  }

  PointType closestPoint;
  closestPoint.Fill(0.0);
  closestPoint = this->m_MovingTransformedPointSet->GetPoint(pointId);

  // the probabilistic icp term
  const MeasureType euclideanDistance = point.EuclideanDistanceTo(closestPoint);
  MeasureType       distanceProbability =
    std::exp(static_cast<MeasureType>(-0.5) * itk::Math::sqr(euclideanDistance / this->m_EuclideanDistanceSigma));

  SizeValueType numberOfVoxelsInNeighborhood = pixel.size() / (1 + PointDimension);
  SizeValueType centerIntensityIndex =
    static_cast<SizeValueType>(0.5 * numberOfVoxelsInNeighborhood) * (PointDimension + 1);

  // the probabilistic intensity term
  MeasureType intensityDistance = pixel[centerIntensityIndex] - closestPixel[centerIntensityIndex];
  MeasureType intensityProbability =
    std::exp(static_cast<MeasureType>(-0.5) * itk::Math::sqr(intensityDistance / this->m_IntensityDistanceSigma));

  measure = -itk::NumericTraits<MeasureType>::OneValue() * intensityProbability * distanceProbability;

  // total derivative is
  // d/dx( intProb * distProb ) =
  //   intProb * d/dx( distProb ) + distProb * d/dx( intProb ) =
  //   intProb * distProb * dist + distProb * intProb * intdiff =

  localDerivative = (closestPoint - point) * intensityProbability * distanceProbability;
  for (SizeValueType d = 0; d < PointDimension; d++)
  {
    localDerivative[d] +=
      intensityProbability * distanceProbability * intensityDistance * closestPixel[centerIntensityIndex + 1 + d];
  }
}

template <typename TFixedPointSet, typename TMovingPointSet, typename TInternalComputationValueType>
typename LightObject::Pointer
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>::
  InternalClone(void) const
{
  typename Self::Pointer rval = Self::New();
  rval->SetMovingPointSet(this->m_MovingPointSet);
  rval->SetFixedPointSet(this->m_FixedPointSet);

  return rval.GetPointer();
}

/** PrintSelf method */
template <typename TFixedPointSet, typename TMovingPointSet, typename TInternalComputationValueType>
void
MeanSquaresPointSetToPointSetIntensityMetricv4<TFixedPointSet, TMovingPointSet, TInternalComputationValueType>::
  PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << "Euclidean distance sigma = " << this->m_EuclideanDistanceSigma << std::endl;
  os << "intensity distance sigma = " << this->m_IntensityDistanceSigma << std::endl;
}

} // end namespace itk

#endif
