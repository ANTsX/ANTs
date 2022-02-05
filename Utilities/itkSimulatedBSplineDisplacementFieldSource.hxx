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
#ifndef itkSimulatedBSplineDisplacementFieldSource_hxx
#define itkSimulatedBSplineDisplacementFieldSource_hxx


#include "itkContinuousIndex.h"

namespace itk
{
template <typename TOutputImage>
SimulatedBSplineDisplacementFieldSource<TOutputImage>::SimulatedBSplineDisplacementFieldSource()
  : m_SplineOrder(3)
{
  this->m_NumberOfFittingLevels.Fill(1);
  this->m_NumberOfControlPoints.Fill(4);
  this->m_DisplacementNoiseStandardDeviation.Fill(1.0);
}

template <typename TOutputImage>
void
SimulatedBSplineDisplacementFieldSource<TOutputImage>::GenerateData()
{
  typename BSplineFilterType::Pointer bsplineFilter = BSplineFilterType::New();

  bsplineFilter->SetEstimateInverse(false);
  bsplineFilter->SetEnforceStationaryBoundary(this->GetEnforceStationaryBoundary());
  bsplineFilter->SetSplineOrder(this->m_SplineOrder);
  bsplineFilter->SetNumberOfFittingLevels(this->m_NumberOfFittingLevels);
  bsplineFilter->SetNumberOfControlPoints(this->m_NumberOfControlPoints);
  bsplineFilter->SetBSplineDomain(
    this->GetOutputOrigin(), this->GetOutputSpacing(), this->GetOutputSize(), this->GetOutputDirection());

  SizeType outputSize = this->GetOutputSize();

  typename PointSetType::Pointer randomPointSet = PointSetType::New();
  randomPointSet->Initialize();

  for (SizeValueType n = 0; n < this->GetNumberOfRandomPoints(); n++)
  {
    VectorType                                randomVector;
    ContinuousIndex<RealType, ImageDimension> randomIndex;
    for (SizeValueType d = 0; d < ImageDimension; d++)
    {
      randomIndex[d] = this->GetRandomizer()->GetUniformVariate(NumericTraits<double>::ZeroValue(),
                                                                static_cast<double>(outputSize[d] - 1));
      randomVector[d] = this->GetRandomizer()->GetNormalVariate(
        NumericTraits<double>::ZeroValue(), std::pow(this->m_DisplacementNoiseStandardDeviation[d], 2));
    }
    typename OutputImageType::PointType imagePoint;
    this->GetOutput()->TransformContinuousIndexToPhysicalPoint(randomIndex, imagePoint);

    PointType physicalPoint;
    physicalPoint.CastFrom(imagePoint);

    randomPointSet->SetPoint(n, physicalPoint);
    randomPointSet->SetPointData(n, randomVector);
  }
  bsplineFilter->SetPointSet(randomPointSet);
  bsplineFilter->Update();

  this->ProcessObject::SetNthOutput(0, bsplineFilter->GetOutput());
}

template <typename TOutputImage>
void
SimulatedBSplineDisplacementFieldSource<TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Spline order: " << this->m_SplineOrder << std::endl;
  os << indent << "Number of fitting levels: " << this->m_NumberOfFittingLevels << std::endl;
  os << indent << "Number of control points: " << this->m_NumberOfControlPoints << std::endl;
  os << indent << "Displacement noise standard deviation: " << this->m_DisplacementNoiseStandardDeviation << std::endl;
}


} // end namespace itk

#endif
