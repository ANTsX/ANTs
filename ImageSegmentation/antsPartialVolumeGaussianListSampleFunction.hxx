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
#ifndef __antsPartialVolumeGaussianListSampleFunction_hxx
#define __antsPartialVolumeGaussianListSampleFunction_hxx


#include "itkMeanSampleFilter.h"
#include "itkWeightedMeanSampleFilter.h"

namespace itk
{
namespace ants
{
namespace Statistics
{
template <typename TListSample, typename TOutput, typename TCoordRep>
PartialVolumeGaussianListSampleFunction<TListSample, TOutput, TCoordRep>::PartialVolumeGaussianListSampleFunction()
{
  this->m_Gaussian = GaussianType::New();

  this->m_IsCalculated[0] = false;
  this->m_IsCalculated[1] = false;
}

template <typename TListSample, typename TOutput, typename TCoordRep>
PartialVolumeGaussianListSampleFunction<TListSample, TOutput, TCoordRep>::~PartialVolumeGaussianListSampleFunction() =
  default;

template <typename TListSample, typename TOutput, typename TCoordRep>
void
PartialVolumeGaussianListSampleFunction<TListSample, TOutput, TCoordRep>::SetIndexedInputListSample(
  const unsigned int          d,
  const InputListSampleType * ptr)
{
  Superclass::SetIndexedInputListSample(d, ptr);

  if (d > 1)
  {
    itkExceptionMacro("This class only requires two input list samples.");
  }

  if (!this->GetInputListSample(d))
  {
    return;
  }
  else
  {
    this->CalculateGaussianParametersFromListSample(
      this->GetInputListSample(d), this->GetListSampleWeights(d), this->m_Mean[d]);

    this->m_IsCalculated[d] = true;
  }
  if (this->m_IsCalculated[0] && this->m_IsCalculated[1])
  {
    this->CalculateGaussianParameters();
  }
}

template <typename TListSample, typename TOutput, typename TCoordRep>
void
PartialVolumeGaussianListSampleFunction<TListSample, TOutput, TCoordRep>::CalculateGaussianParametersFromListSample(
  const InputListSampleType *       listSample,
  const ListSampleWeightArrayType * weights,
  MeanType &                        mean)
{
  if (!listSample)
  {
    return;
  }

  if (listSample->Size() > 1)
  {
    if (weights->Size() == listSample->Size())
    {
      typedef typename itk::Statistics::WeightedMeanSampleFilter<InputListSampleType> MeanCalculatorType;
      typename MeanCalculatorType::Pointer meanCalculator = MeanCalculatorType::New();

      meanCalculator->SetWeights(*weights);
      meanCalculator->SetInput(listSample);
      meanCalculator->Update();

      NumericTraits<MeanType>::SetLength(mean, listSample->GetMeasurementVectorSize());
      for (unsigned int d = 0; d < listSample->GetMeasurementVectorSize(); d++)
      {
        mean[d] = meanCalculator->GetMean()[d];
      }
    }
    else
    {
      typedef itk::Statistics::MeanSampleFilter<InputListSampleType> MeanCalculatorType;
      typename MeanCalculatorType::Pointer                           meanCalculator = MeanCalculatorType::New();
      meanCalculator->SetInput(listSample);
      meanCalculator->Update();

      NumericTraits<MeanType>::SetLength(mean, listSample->GetMeasurementVectorSize());
      for (unsigned int d = 0; d < listSample->GetMeasurementVectorSize(); d++)
      {
        mean[d] = meanCalculator->GetMean()[d];
      }
    }
  }
  else
  {
    itkWarningMacro("The input list sample has <= 1 element.");
  }
}

template <typename TListSample, typename TOutput, typename TCoordRep>
void
PartialVolumeGaussianListSampleFunction<TListSample, TOutput, TCoordRep>::CalculateGaussianParameters()
{
  if (this->m_Mean[0].Size() != this->m_Mean[1].Size())
  {
    itkExceptionMacro("Mean sizes are unequal.");
  }

  MeanType mean;
  NumericTraits<MeanType>::SetLength(mean, this->m_Mean[0].Size());

  CovarianceType covariance;
  covariance.SetSize(mean.Size(), mean.Size());
  covariance.SetIdentity();
  for (unsigned int d = 0; d < mean.Size(); d++)
  {
    mean[d] = 0.5 * (this->m_Mean[0][d] + this->m_Mean[1][d]);
    covariance(d, d) = 1.0 / 12.0 * itk::Math::sqr(this->m_Mean[0][d]) +
                       -1.0 / 6.0 * this->m_Mean[0][d] * this->m_Mean[1][d] +
                       1.0 / 12.0 * itk::Math::sqr(this->m_Mean[1][d]);
  }

  this->m_Gaussian->SetMean(mean);
  this->m_Gaussian->SetCovariance(covariance);
}

template <typename TListSample, typename TOutput, typename TCoordRep>
TOutput
PartialVolumeGaussianListSampleFunction<TListSample, TOutput, TCoordRep>::Evaluate(
  const InputMeasurementVectorType & measurement) const
{
  if (this->m_IsCalculated[0] && this->m_IsCalculated[1])
  {
    try
    {
      return this->m_Gaussian->Evaluate(measurement);
    }
    catch (...)
    {
      return 0.0;
    }
  }
  else
  {
    return 0.0;
  }
}

/**
 * Standard "PrintSelf" method
 */
template <typename TListSample, typename TOutput, typename TCoordRep>
void
PartialVolumeGaussianListSampleFunction<TListSample, TOutput, TCoordRep>::PrintSelf(std::ostream & os,
                                                                                    Indent         indent) const
{
  os << indent << "mean = " << this->m_Gaussian->GetMean() << ", ";

  CovarianceType covariance = this->m_Gaussian->GetCovariance();
  os << "covariance = [";
  for (unsigned int r = 0; r < covariance.Rows(); r++)
  {
    for (unsigned int c = 0; c < covariance.Cols() - 1; c++)
    {
      os << covariance(r, c) << ", ";
    }
    if (r == covariance.Rows() - 1)
    {
      os << covariance(r, covariance.Cols() - 1) << "]" << std::endl;
    }
    else
    {
      os << covariance(r, covariance.Cols() - 1) << "; ";
    }
  }
}
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#endif
