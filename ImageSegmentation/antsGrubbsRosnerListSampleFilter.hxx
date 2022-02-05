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
#ifndef __antsGrubbsRosnerListSampleFilter_hxx
#define __antsGrubbsRosnerListSampleFilter_hxx


#include "itkTDistribution.h"

namespace itk
{
namespace ants
{
namespace Statistics
{
template <typename TScalarListSample>
GrubbsRosnerListSampleFilter<TScalarListSample>::GrubbsRosnerListSampleFilter()
{
  this->AllocateOutput();
  this->GetOutput()->SetMeasurementVectorSize(1);

  this->m_OutlierHandling = Winsorize;
  this->m_WinsorizingLevel = 0.10;
  this->m_SignificanceLevel = 0.05;
}

template <typename TScalarListSample>
GrubbsRosnerListSampleFilter<TScalarListSample>::~GrubbsRosnerListSampleFilter() = default;

template <typename TScalarListSample>
void
GrubbsRosnerListSampleFilter<TScalarListSample>::GenerateData()
{
  if (this->GetInput()->GetMeasurementVectorSize() != 1)
  {
    itkExceptionMacro("The input sample must be univariate.");
  }

  const unsigned int scalarMeasurementVectorSize = this->GetOutput()->GetMeasurementVectorSize();

  this->GetOutput()->SetMeasurementVectorSize(scalarMeasurementVectorSize);

  /**
   * A common hueristic is that Grubbs-Rosner outlier removal does not work for
   * sample sizes less than or equal to 6.
   */
  if (this->GetInput()->Size() <= 6)
  {
    typename ScalarListSampleType::ConstIterator It = this->GetInput()->Begin();
    while (It != this->GetInput()->End())
    {
      MeasurementVectorType inputMeasurement = It.GetMeasurementVector();
      MeasurementVectorType outputMeasurement;
      outputMeasurement.SetSize(scalarMeasurementVectorSize);
      for (unsigned int d = 0; d < scalarMeasurementVectorSize; d++)
      {
        outputMeasurement[d] = inputMeasurement[d];
      }
      this->GetOutput()->PushBack(outputMeasurement);
      ++It;
    }

    return;
  }

  /**
   * Otherwise, iterate through the input list, removing t
   */

  RealType mean = 0.0;
  RealType variance = 0.0;
  RealType count = 0.0;

  typename ScalarListSampleType::ConstIterator It = this->GetInput()->Begin();
  while (It != this->GetInput()->End())
  {
    MeasurementVectorType inputMeasurement = It.GetMeasurementVector();

    count += NumericTraits<RealType>::OneValue();
    variance += (count - NumericTraits<RealType>::OneValue()) *
                itk::Math::sqr(static_cast<RealType>(inputMeasurement[0]) - mean) / count;
    mean = mean + (static_cast<RealType>(inputMeasurement[0]) - mean) / count;
    ++It;
  }

  variance /= (count - NumericTraits<RealType>::OneValue());

  bool outlierFound = true;
  this->m_OutlierInstanceIdentifiers.clear();
  while (outlierFound == true && (this->GetInput()->Size() - this->m_OutlierInstanceIdentifiers.size() > 6))
  {
    outlierFound = false;
    InstanceIdentifierType id = this->FindMaximumNonOutlierDeviationValue(mean, variance);
    if (this->GetInput()->GetFrequency(id) > 0)
    {
      MeasurementVectorType measurement = this->GetInput()->GetMeasurementVector(id);
      outlierFound = this->IsMeasurementAnOutlier(
        measurement[0], mean, variance, this->GetInput()->Size() - this->m_OutlierInstanceIdentifiers.size());
      if (outlierFound)
      {
        /** Retabulate the variance and mean by removing the previous estimate */
        RealType count2 = this->GetInput()->Size() - this->m_OutlierInstanceIdentifiers.size();
        mean = (mean * count2 - static_cast<RealType>(measurement[0])) / (count2 - NumericTraits<RealType>::OneValue());
        variance = (count2 - 1.0) * variance - (count2 - NumericTraits<RealType>::OneValue()) *
                                                 itk::Math::sqr(static_cast<RealType>(measurement[0]) - mean) / count2;
        variance /= (count2 - static_cast<RealType>(2.0));
        this->m_OutlierInstanceIdentifiers.push_back(id);
      }
    }
  }

  RealType lowerWinsorBound = 0.0;
  RealType upperWinsorBound = 0.0;
  if (this->m_OutlierHandling == Winsorize)
  {
    typename itk::Statistics::TDistribution::Pointer tdistribution = itk::Statistics::TDistribution::New();
    RealType                                         t = tdistribution->EvaluateInverseCDF(
      1.0 - 0.5 * this->m_WinsorizingLevel, this->GetInput()->Size() - this->m_OutlierInstanceIdentifiers.size());

    lowerWinsorBound = mean - t * std::sqrt(variance);
    upperWinsorBound = mean + t * std::sqrt(variance);
  }

  It = this->GetInput()->Begin();
  while (It != this->GetInput()->End())
  {
    MeasurementVectorType inputMeasurement = It.GetMeasurementVector();
    MeasurementVectorType outputMeasurement;
    outputMeasurement.SetSize(scalarMeasurementVectorSize);

    if (this->m_OutlierHandling == None ||
        std::find(this->m_OutlierInstanceIdentifiers.begin(),
                  this->m_OutlierInstanceIdentifiers.end(),
                  It.GetInstanceIdentifier()) == this->m_OutlierInstanceIdentifiers.end())
    {
      outputMeasurement[0] = inputMeasurement[0];
      this->GetOutput()->PushBack(outputMeasurement);
    }
    else if (this->m_OutlierHandling == Winsorize)
    {
      if (static_cast<RealType>(inputMeasurement[0]) < lowerWinsorBound)
      {
        outputMeasurement[0] = lowerWinsorBound;
      }
      else
      {
        outputMeasurement[0] = upperWinsorBound;
      }
      this->GetOutput()->PushBack(outputMeasurement);
    }
    ++It;
  }
}

template <typename TScalarListSample>
typename GrubbsRosnerListSampleFilter<TScalarListSample>::InstanceIdentifierType
GrubbsRosnerListSampleFilter<TScalarListSample>::FindMaximumNonOutlierDeviationValue(RealType mean,
                                                                                     RealType itkNotUsed(variance))
{
  RealType               maximumDeviation = 0.0;
  InstanceIdentifierType maximumID = NumericTraits<InstanceIdentifierType>::max();

  typename ScalarListSampleType::ConstIterator It = this->GetInput()->Begin();
  while (It != this->GetInput()->End())
  {
    MeasurementVectorType  inputMeasurement = It.GetMeasurementVector();
    InstanceIdentifierType inputID = It.GetInstanceIdentifier();

    if (std::find(this->m_OutlierInstanceIdentifiers.begin(), this->m_OutlierInstanceIdentifiers.end(), inputID) ==
        this->m_OutlierInstanceIdentifiers.end())
    {
      if (Math::abs(static_cast<RealType>(inputMeasurement[0]) - mean) > maximumDeviation)
      {
        maximumDeviation = Math::abs(static_cast<RealType>(inputMeasurement[0]) - mean);
        maximumID = inputID;
      }
    }
    ++It;
  }

  return maximumID;
}

template <typename TScalarListSample>
bool
GrubbsRosnerListSampleFilter<TScalarListSample>::IsMeasurementAnOutlier(RealType      x,
                                                                        RealType      mean,
                                                                        RealType      variance,
                                                                        unsigned long N)
{
  /**
   * The Grubb critical two-sided value is defined to be
   * (N-1)/sqrt(N)*sqrt( t*t / (N-2+t*t) ) where t is at the
   * (alpha / (2N)) signficance level with N-2 degrees of freedom.
   */

  RealType sig = this->m_SignificanceLevel / (2.0 * static_cast<RealType>(N));

  typename itk::Statistics::TDistribution::Pointer tdistribution = itk::Statistics::TDistribution::New();

  RealType t = tdistribution->EvaluateInverseCDF(1.0 - sig, N - 2);

  RealType nu = static_cast<RealType>(N - 1);
  RealType g = nu / std::sqrt(nu + 1.0) * std::sqrt(t * t / (nu - 1 + t * t));

  return g < (itk::Math::abs(x - mean) / std::sqrt(variance));
}

template <typename TScalarListSample>
void
GrubbsRosnerListSampleFilter<TScalarListSample>::PrintSelf(std::ostream & os, Indent indent) const
{
  os << indent << "Significance level: " << this->m_SignificanceLevel << std::endl;
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
    os << "Winsorize";
    os << " (level = " << this->m_WinsorizingLevel << ")" << std::endl;
  }
  if (this->m_OutlierInstanceIdentifiers.size() > 0)
  {
    os << indent << "Outlier Identifiers: " << std::endl;
    for (unsigned int d = 0; d < this->m_OutlierInstanceIdentifiers.size(); d++)
    {
      os << indent << "   " << this->m_OutlierInstanceIdentifiers[d] << std::endl;
    }
  }
  else
  {
    os << indent << "There are no outliers." << std::endl;
  }
}
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#endif
