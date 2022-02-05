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
#ifndef __antsBoxPlotQuantileListSampleFilter_hxx
#define __antsBoxPlotQuantileListSampleFilter_hxx


#include "itkDenseFrequencyContainer2.h"
#include "itkHistogram.h"
#include "itkSampleToHistogramFilter.h"

namespace itk
{
namespace ants
{
namespace Statistics
{
template <typename TScalarListSample>
BoxPlotQuantileListSampleFilter<TScalarListSample>::BoxPlotQuantileListSampleFilter()
{
  this->AllocateOutput();
  this->GetOutput()->SetMeasurementVectorSize(1);

  this->m_OutlierHandling = Winsorize;
  this->m_WhiskerScalingFactor = 1.5;
  this->m_LowerPercentile = 0.25;
  this->m_UpperPercentile = 0.75;
}

template <typename TScalarListSample>
BoxPlotQuantileListSampleFilter<TScalarListSample>::~BoxPlotQuantileListSampleFilter() = default;

template <typename TScalarListSample>
void
BoxPlotQuantileListSampleFilter<TScalarListSample>::GenerateData()
{
  if (this->GetInput()->GetMeasurementVectorSize() != 1)
  {
    itkExceptionMacro("The input sample must be univariate.");
  }
  if (this->m_LowerPercentile >= this->m_UpperPercentile)
  {
    itkExceptionMacro("Lower percentile must be less than upper percentile.");
  }

  const unsigned int scalarMeasurementVectorSize = this->GetInput()->GetMeasurementVectorSize();
  this->GetOutput()->SetMeasurementVectorSize(scalarMeasurementVectorSize);

  /**
   * Initialize the histogram in preparation
   */

  typedef itk::Statistics::Histogram<RealType, itk::Statistics::DenseFrequencyContainer2> HistogramType;
  typedef itk::Statistics::SampleToHistogramFilter<ScalarListSampleType, HistogramType>   SampleFilterType;
  typename SampleFilterType::HistogramSizeType                                            histogramSize(1);
  histogramSize.Fill(200);

  typename SampleFilterType::Pointer sampleFilter = SampleFilterType::New();
  sampleFilter->SetInput(this->GetInput());
  sampleFilter->SetHistogramSize(histogramSize);
  sampleFilter->Update();

  RealType lowerQuantile = sampleFilter->GetOutput()->Quantile(0, this->m_LowerPercentile);
  RealType upperQuantile = sampleFilter->GetOutput()->Quantile(0, this->m_UpperPercentile);

  RealType upperBound = upperQuantile + this->m_WhiskerScalingFactor * (upperQuantile - lowerQuantile);
  RealType lowerBound = lowerQuantile - this->m_WhiskerScalingFactor * (upperQuantile - lowerQuantile);

  typename ScalarListSampleType::ConstIterator It = this->GetInput()->Begin();
  It = this->GetInput()->Begin();
  while (It != this->GetInput()->End())
  {
    MeasurementVectorType                                inputMeasurement = It.GetMeasurementVector();
    typename ScalarListSampleType::MeasurementVectorType outputMeasurement;
    outputMeasurement.SetSize(scalarMeasurementVectorSize);
    if (static_cast<RealType>(inputMeasurement[0]) < lowerBound ||
        static_cast<RealType>(inputMeasurement[0]) > upperBound)
    {
      this->m_OutlierInstanceIdentifiers.push_back(It.GetInstanceIdentifier());
      if (this->m_OutlierHandling == None)
      {
        outputMeasurement[0] = inputMeasurement[0];
        this->GetOutput()->PushBack(outputMeasurement);
      }
      // else trim from the output
    }
    else
    {
      outputMeasurement[0] = inputMeasurement[0];
      this->GetOutput()->PushBack(outputMeasurement);
    }
    ++It;
  }

  if (this->m_OutlierHandling == Winsorize)
  {
    /** Retabulate the histogram with the outliers removed */
    typename SampleFilterType::Pointer sampleFilter2 = SampleFilterType::New();
    sampleFilter2->SetInput(this->GetOutput());
    sampleFilter2->SetHistogramSize(histogramSize);
    sampleFilter2->Update();

    RealType lowerQuantile2 = sampleFilter2->GetOutput()->Quantile(0, this->m_LowerPercentile);
    RealType upperQuantile2 = sampleFilter2->GetOutput()->Quantile(0, this->m_UpperPercentile);

    RealType upperBound2 = upperQuantile2 + this->m_WhiskerScalingFactor * (upperQuantile2 - lowerQuantile2);
    RealType lowerBound2 = lowerQuantile2 - this->m_WhiskerScalingFactor * (upperQuantile2 - lowerQuantile2);

    this->GetOutput()->Clear();

    It = this->GetInput()->Begin();
    while (It != this->GetInput()->End())
    {
      MeasurementVectorType                                inputMeasurement = It.GetMeasurementVector();
      typename ScalarListSampleType::MeasurementVectorType outputMeasurement;
      outputMeasurement.SetSize(scalarMeasurementVectorSize);
      outputMeasurement[0] = inputMeasurement[0];
      if (static_cast<RealType>(inputMeasurement[0]) < lowerBound)
      {
        outputMeasurement[0] = lowerBound2;
      }
      else if (static_cast<RealType>(inputMeasurement[0]) > upperBound)
      {
        outputMeasurement[0] = upperBound2;
      }
      this->GetOutput()->PushBack(outputMeasurement);
      ++It;
    }
  }
}

template <typename TScalarListSample>
void
BoxPlotQuantileListSampleFilter<TScalarListSample>::PrintSelf(std::ostream & os, Indent indent) const
{
  os << indent << "Percentile Bounds: [" << this->m_LowerPercentile << ", " << this->m_UpperPercentile << "]"
     << std::endl;
  os << indent << "Whisker scaling factor: " << this->m_WhiskerScalingFactor << std::endl;
  os << indent << "Outlier handling: ";
  if (this->m_OutlierHandling == None)
  {
    os << "None" << std::endl;
  }
  if (this->m_OutlierHandling == Trim)
  {
    os << "Trim" << std::endl;
  }
  if (this->m_OutlierHandling == Winsorize)
  {
    os << "Winsorize" << std::endl;
  }
}
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#endif
