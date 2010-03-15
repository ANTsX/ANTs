/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBoxPlotQuantileListSampleFilter.txx,v $
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBoxPlotQuantileListSampleFilter_txx
#define __itkBoxPlotQuantileListSampleFilter_txx

#include "itkBoxPlotQuantileListSampleFilter.h"

#include "itkHistogram.h"
#include "itkListSampleToHistogramFilter.h"

namespace itk
{
namespace Statistics
{
template <class TScalarListSample>
BoxPlotQuantileListSampleFilter<TScalarListSample>
::BoxPlotQuantileListSampleFilter()
{
  this->AllocateOutput();
  this->GetOutput()->SetMeasurementVectorSize( 1 );

  this->m_OutlierHandling = Winsorize;
  this->m_WhiskerScalingFactor = 1.5;
  this->m_LowerPercentile = 0.25;
  this->m_UpperPercentile = 0.75;
}

template <class TScalarListSample>
BoxPlotQuantileListSampleFilter<TScalarListSample>
::~BoxPlotQuantileListSampleFilter()
{
}

template <class TScalarListSample>
void
BoxPlotQuantileListSampleFilter<TScalarListSample>
::GenerateData()
{
  if( this->GetInput()->GetMeasurementVectorSize() != 1 )
    {
    itkExceptionMacro( "The input sample must be univariate." );
    }
  if( this->m_LowerPercentile >= this->m_UpperPercentile )
    {
    itkExceptionMacro( "Lower percentile must be less than upper percentile." );
    }

  const unsigned int scalarMeasurementVectorSize =
    this->GetOutput()->GetMeasurementVectorSize();
  this->GetOutput()->SetMeasurementVectorSize( scalarMeasurementVectorSize );

  /**
   * Initialize the histogram in preparation for the
   */

  typedef Histogram<RealType, 1> HistogramType;
  typename HistogramType::Pointer histogram = HistogramType::New();

  typename HistogramType::MeasurementVectorType minimumValue;
  minimumValue[0] = NumericTraits<RealType>::max();
  typename HistogramType::MeasurementVectorType maximumValue;
  maximumValue[0] = NumericTraits<RealType>::NonpositiveMin();
  typename HistogramType::SizeType histogramSize;
  histogramSize.Fill( 200 );

  typename ScalarListSampleType::ConstIterator It = this->GetInput()->Begin();
  while( It != this->GetInput()->End() )
    {
    MeasurementVectorType inputMeasurement = It.GetMeasurementVector();
    if( inputMeasurement[0] < minimumValue[0] )
      {
      minimumValue[0] = inputMeasurement[0];
      }
    if( inputMeasurement[0] > maximumValue[0] )
      {
      maximumValue[0] = inputMeasurement[0];
      }
    ++It;
    }

  histogram->Initialize( histogramSize, minimumValue, maximumValue );

  typedef ListSampleToHistogramFilter<ScalarListSampleType, HistogramType>
    SampleFilterType;
  typename SampleFilterType::Pointer sampleFilter = SampleFilterType::New();
  sampleFilter->SetHistogram( histogram );
  sampleFilter->SetListSample( this->GetInput() );
  sampleFilter->Update();

  RealType lowerQuantile = histogram->Quantile( 0, this->m_LowerPercentile );
  RealType upperQuantile = histogram->Quantile( 0, this->m_UpperPercentile );

  RealType upperBound = upperQuantile
    + this->m_WhiskerScalingFactor * ( upperQuantile - lowerQuantile );
  RealType lowerBound = lowerQuantile
    - this->m_WhiskerScalingFactor * ( upperQuantile - lowerQuantile );

  It = this->GetInput()->Begin();
  while( It != this->GetInput()->End() )
    {
    MeasurementVectorType inputMeasurement = It.GetMeasurementVector();
    typename ScalarListSampleType::MeasurementVectorType outputMeasurement;
    outputMeasurement.SetSize( scalarMeasurementVectorSize );
    if( inputMeasurement[0] < lowerBound || inputMeasurement[0] > upperBound )
      {
      this->m_OutlierInstanceIdentifiers.push_back( It.GetInstanceIdentifier() );
      if( this->m_OutlierHandling == None )
        {
        outputMeasurement[0] = inputMeasurement[0];
        this->GetOutput()->PushBack( outputMeasurement );
        }
      // else trim from the output
      }
    else
      {
      outputMeasurement[0] = inputMeasurement[0];
      this->GetOutput()->PushBack( outputMeasurement );
      }
    ++It;
    }

  if( this->m_OutlierHandling == Winsorize )
    {
    /** Retabulate the histogram with the outliers removed */
    minimumValue[0] = NumericTraits<RealType>::max();
    maximumValue[0] = NumericTraits<RealType>::NonpositiveMin();

    typename ScalarListSampleType::Iterator ItO = this->GetOutput()->Begin();
    while( ItO != this->GetOutput()->End() )
      {
      MeasurementVectorType outputMeasurement = ItO.GetMeasurementVector();
      if( outputMeasurement[0] < minimumValue[0] )
        {
        minimumValue[0] = outputMeasurement[0];
        }
      if( outputMeasurement[0] > maximumValue[0] )
        {
        maximumValue[0] = outputMeasurement[0];
        }
      ++ItO;
      }

    typename HistogramType::Pointer histogram2 = HistogramType::New();
    histogram2->Initialize( histogramSize, minimumValue, maximumValue );

    typename SampleFilterType::Pointer sampleFilter2 = SampleFilterType::New();
    sampleFilter2->SetHistogram( histogram2 );
    sampleFilter2->SetListSample( this->GetOutput() );
    sampleFilter2->Update();

    RealType lowerQuantile2 = histogram2->Quantile( 0, this->m_LowerPercentile );
    RealType upperQuantile2 = histogram2->Quantile( 0, this->m_UpperPercentile );

    RealType upperBound2 = upperQuantile2
      + this->m_WhiskerScalingFactor * ( upperQuantile2 - lowerQuantile2 );
    RealType lowerBound2 = lowerQuantile2
      - this->m_WhiskerScalingFactor * ( upperQuantile2 - lowerQuantile2 );

    this->GetOutput()->Clear();

    It = this->GetInput()->Begin();
    while( It != this->GetInput()->End() )
      {
      MeasurementVectorType inputMeasurement = It.GetMeasurementVector();
      typename ScalarListSampleType::MeasurementVectorType outputMeasurement;
      outputMeasurement.SetSize( scalarMeasurementVectorSize );
      outputMeasurement[0] = inputMeasurement[0];
      if( inputMeasurement[0] < lowerBound )
        {
        outputMeasurement[0] = lowerBound2;
        }
      else if( inputMeasurement[0] > upperBound )
        {
        outputMeasurement[0] = upperBound2;
        }
      this->GetOutput()->PushBack( outputMeasurement );
      ++It;
      }
    }
}

template <class TScalarListSample>
void
BoxPlotQuantileListSampleFilter<TScalarListSample>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  os << indent << "Percentile Bounds: ["
     << this->m_LowerPercentile << ", " << this->m_UpperPercentile << "]"
     << std::endl;
  os << indent << "Whisker scaling factor: "
     << this->m_WhiskerScalingFactor << std::endl;
  os << indent << "Outlier handling: ";
  if( this->m_OutlierHandling == None )
    {
    os << "None" << std::endl;
    }
  if( this->m_OutlierHandling == Trim )
    {
    os << "Trim" << std::endl;
    }
  if( this->m_OutlierHandling == Winsorize )
    {
    os << "Winsorize" << std::endl;
    }
}
} // end of namespace Statistics
} // end of namespace itk

#endif
