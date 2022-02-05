/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
  https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt
  for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __antsJointHistogramParzenWindowsListSampleFunction_hxx
#define __antsJointHistogramParzenWindowsListSampleFunction_hxx


#include "itkArray.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkContinuousIndex.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkDivideByConstantImageFilter.h"
#include "itkStatisticsImageFilter.h"

namespace itk
{
namespace ants
{
namespace Statistics
{
template <typename TListSample, typename TOutput, typename TCoordRep>
JointHistogramParzenWindowsListSampleFunction<TListSample, TOutput, TCoordRep>::
  JointHistogramParzenWindowsListSampleFunction()
{
  this->m_NumberOfJointHistogramBins = 32;
  this->m_Sigma = 1.0;
  this->m_UseNNforJointHistIncrements = true;
}

template <typename TListSample, typename TOutput, typename TCoordRep>
JointHistogramParzenWindowsListSampleFunction<TListSample, TOutput, TCoordRep>::
  ~JointHistogramParzenWindowsListSampleFunction()
{}

template <typename TListSample, typename TOutput, typename TCoordRep>
void
JointHistogramParzenWindowsListSampleFunction<TListSample, TOutput, TCoordRep>::IncrementJointHistogram(
  RealType     eigenvalue1,
  RealType     eigenvalue2,
  unsigned int which_hist)
{
  RealType newWeight = 1.0;

  // now define two joint histograms, one for shape, one for orientation.
  // first, the shape histogram --- 0,0 origin and spacing of 1
  if (this->m_JointHistogramImages.size() == which_hist)
  {
    typename JointHistogramImageType::SpacingType spacing;
    spacing.Fill(1);
    typename JointHistogramImageType::PointType origin;
    origin.Fill(0);
    typename JointHistogramImageType::SizeType size;
    size.Fill(this->m_NumberOfJointHistogramBins);
    typename JointHistogramImageType::DirectionType direction;
    direction.SetIdentity();
    typename JointHistogramImageType::Pointer curJHI =
      AllocImage<JointHistogramImageType>(size, spacing, origin, regions, 0);

    this->m_JointHistogramImages.push_back(curJHI);
  }

  typename JointHistogramImageType::PointType shapePoint;
  if (eigenvalue1 > 1)
  {
    eigenvalue1 = 1;
  }
  if (eigenvalue2 > 1)
  {
    eigenvalue2 = 1;
  }
  if (eigenvalue1 < 0)
  {
    eigenvalue1 = 0;
  }
  if (eigenvalue2 < 0)
  {
    eigenvalue2 = 0;
  }
  shapePoint[0] = eigenvalue1 * (this->m_NumberOfJointHistogramBins - 1);
  shapePoint[1] = eigenvalue2 * (this->m_NumberOfJointHistogramBins - 1);

  ContinuousIndex<double, 2> shapeCidx;
  this->m_JointHistogramImages[which_hist]->TransformPhysicalPointToContinuousIndex(shapePoint, shapeCidx);

  typename JointHistogramImageType::IndexType shapeIdx;

  /** Nearest neighbor increment to JH */
  if (this->m_UseNNforJointHistIncrements)
  {
    shapeIdx[0] = std::floor(shapeCidx[0] + 0.5);
    shapeIdx[1] = std::floor(shapeCidx[1] + 0.5);
    if (this->m_JointHistogramImages[which_hist]->GetLargestPossibleRegion().IsInside(shapeIdx))
    {
      RealType oldWeight = this->m_JointHistogramImages[which_hist]->GetPixel(shapeIdx);
      this->m_JointHistogramImages[which_hist]->SetPixel(shapeIdx, 1 + oldWeight);
    }
  }
  else
  {
    /** linear addition */
    shapeIdx[0] = static_cast<typename JointHistogramImageType::IndexType::IndexValueType>(std::floor(shapeCidx[0]));
    shapeIdx[1] = static_cast<typename JointHistogramImageType::IndexType::IndexValueType>(std::floor(shapeCidx[1]));
    RealType dist1 = sqrt((shapeCidx[0] - shapeIdx[0]) * (shapeCidx[0] - shapeIdx[0]) +
                          (shapeCidx[1] - shapeIdx[1]) * (shapeCidx[1] - shapeIdx[1]));
    shapeIdx[0]++;
    RealType dist2 = sqrt((shapeCidx[0] - shapeIdx[0]) * (shapeCidx[0] - shapeIdx[0]) +
                          (shapeCidx[1] - shapeIdx[1]) * (shapeCidx[1] - shapeIdx[1]));
    shapeIdx[1]++;
    RealType dist3 = sqrt((shapeCidx[0] - shapeIdx[0]) * (shapeCidx[0] - shapeIdx[0]) +
                          (shapeCidx[1] - shapeIdx[1]) * (shapeCidx[1] - shapeIdx[1]));
    shapeIdx[0]--;
    RealType dist4 = sqrt((shapeCidx[0] - shapeIdx[0]) * (shapeCidx[0] - shapeIdx[0]) +
                          (shapeCidx[1] - shapeIdx[1]) * (shapeCidx[1] - shapeIdx[1]));
    RealType distsum = dist1 + dist2 + dist3 + dist4;
    dist1 /= distsum;
    dist2 /= distsum;
    dist3 /= distsum;
    dist4 /= distsum;

    shapeIdx[0] = static_cast<typename JointHistogramImageType::IndexType::IndexValueType>(std::floor(shapeCidx[0]));
    shapeIdx[1] = static_cast<typename JointHistogramImageType::IndexType::IndexValueType>(std::floor(shapeCidx[1]));
    if (this->m_JointHistogramImages[which_hist]->GetLargestPossibleRegion().IsInside(shapeIdx))
    {
      RealType oldWeight = this->m_JointHistogramImages[which_hist]->GetPixel(shapeIdx);
      this->m_JointHistogramImages[which_hist]->SetPixel(shapeIdx, (1.0 - dist1) * newWeight + oldWeight);
    }
    shapeIdx[0]++;
    if (this->m_JointHistogramImages[which_hist]->GetLargestPossibleRegion().IsInside(shapeIdx))
    {
      RealType oldWeight = this->m_JointHistogramImages[which_hist]->GetPixel(shapeIdx);
      this->m_JointHistogramImages[which_hist]->SetPixel(shapeIdx, (1.0 - dist2) * newWeight + oldWeight);
    }
    shapeIdx[1]++;
    if (this->m_JointHistogramImages[which_hist]->GetLargestPossibleRegion().IsInside(shapeIdx))
    {
      RealType oldWeight = this->m_JointHistogramImages[which_hist]->GetPixel(shapeIdx);
      this->m_JointHistogramImages[which_hist]->SetPixel(shapeIdx, (1.0 - dist3) * newWeight + oldWeight);
    }
    shapeIdx[0]--;
    if (this->m_JointHistogramImages[which_hist]->GetLargestPossibleRegion().IsInside(shapeIdx))
    {
      RealType oldWeight = this->m_JointHistogramImages[which_hist]->GetPixel(shapeIdx);
      this->m_JointHistogramImages[which_hist]->SetPixel(shapeIdx, (1.0 - dist4) * newWeight + oldWeight);
    }
  }
  return;
}

template <typename TListSample, typename TOutput, typename TCoordRep>
void
JointHistogramParzenWindowsListSampleFunction<TListSample, TOutput, TCoordRep>::SetInputListSample(
  const InputListSampleType * ptr)
{
  this->m_ListSample = ptr;
  this->m_JointHistogramImages.clear();
  if (!this->m_ListSample)
  {
    return;
  }

  if (this->m_ListSample->Size() <= 1)
  {
    itkWarningMacro("The input list sample has <= 1 element."
                    << "Function evaluations will be equal to 0.");
    return;
  }
  typename InputListSampleType::ConstIterator It = this->m_ListSample->Begin();
  InputMeasurementVectorType                  inputMeasurement = It.GetMeasurementVector();
  unsigned int                                Dimension = inputMeasurement.Size();
  if ((Dimension % 2) != 0)
  {
    itkWarningMacro("The input list should contain 2*N images where N > 0.");
    return;
  }
  /**
   * Find the min/max values to define the histogram domain
   */
  Array<RealType> minValues(Dimension);
  minValues.Fill(NumericTraits<RealType>::max());
  Array<RealType> maxValues(Dimension);
  maxValues.Fill(NumericTraits<RealType>::NonpositiveMin());

  It = this->m_ListSample->Begin();
  while (It != this->m_ListSample->End())
  {
    InputMeasurementVectorType inputMeasurement = It.GetMeasurementVector();
    for (unsigned int d = 0; d < Dimension; d++)
    {
      if (inputMeasurement[d] < minValues[d])
      {
        minValues[d] = inputMeasurement[d];
      }
      if (inputMeasurement[d] > maxValues[d])
      {
        maxValues[d] = inputMeasurement[d];
      }
    }
    ++It;
  }

  It = this->m_ListSample->Begin();
  while (It != this->m_ListSample->End())
  {
    InputMeasurementVectorType inputMeasurement = It.GetMeasurementVector();
    /** joint-hist model for the eigenvalues */
    unsigned int jhcount = 0;
    for (unsigned int d = 0; d < Dimension; d = d + 2)
    {
      RealType value1 = (inputMeasurement[d] - minValues[d]) / (maxValues[d] - minValues[d]);
      RealType value2 = (inputMeasurement[d + 1] - minValues[d + 1]) / (maxValues[d + 1] - minValues[d + 1]);
      this->IncrementJointHistogram(value1, value2, jhcount);
      jhcount++;
    }
    ++It;
  }
  for (unsigned int d = 0; d < this->m_JointHistogramImages.size(); d++)
  {
    typedef DiscreteGaussianImageFilter<JointHistogramImageType, JointHistogramImageType> GaussianFilterType;
    typename GaussianFilterType::Pointer gaussian = GaussianFilterType::New();
    gaussian->SetInput(this->m_JointHistogramImages[d]);
    gaussian->SetVariance(this->m_Sigma * this->m_Sigma);
    gaussian->SetMaximumError(0.01);
    gaussian->SetUseImageSpacing(false);
    gaussian->Update();

    typedef StatisticsImageFilter<JointHistogramImageType> StatsFilterType;
    typename StatsFilterType::Pointer                      stats = StatsFilterType::New();
    stats->SetInput(gaussian->GetOutput());
    stats->Update();

    typedef DivideByConstantImageFilter<JointHistogramImageType, RealType, JointHistogramImageType> DividerType;
    typename DividerType::Pointer divider = DividerType::New();
    divider->SetInput(gaussian->GetOutput());
    divider->SetConstant(stats->GetSum());
    divider->Update();
    this->m_JointHistogramImages[d] = divider->GetOutput();
  }
}

template <typename TListSample, typename TOutput, typename TCoordRep>
TOutput
JointHistogramParzenWindowsListSampleFunction<TListSample, TOutput, TCoordRep>::Evaluate(
  const InputMeasurementVectorType & measurement) const
{
  try
  {
    typedef BSplineInterpolateImageFunction<JointHistogramImageType> InterpolatorType;

    RealType probability = 1.0;
    for (unsigned int d = 0; d < this->m_JointHistogramImages.size(); d++)
    {
      typename JointHistogramImageType::PointType point;
      point[0] = measurement[d];

      typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
      interpolator->SetSplineOrder(3);
      interpolator->SetInputImage(this->m_JointHistogramImages[d]);
      if (interpolator->IsInsideBuffer(point))
      {
        probability *= interpolator->Evaluate(point);
      }
      else
      {
        return 0;
      }
    }
    return probability;
  }
  catch (...)
  {
    return 0;
  }
}

/**
 * Standard "PrintSelf" method
 */
template <typename TListSample, typename TOutput, typename TCoordRep>
void
JointHistogramParzenWindowsListSampleFunction<TListSample, TOutput, TCoordRep>::PrintSelf(std::ostream & os,
                                                                                          Indent         indent) const
{
  os << indent << "Sigma: " << this->m_Sigma << std::endl;
  os << indent << "Number of histogram bins: " << this->m_NumberOfJointHistogramBins << std::endl;
}
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#endif
