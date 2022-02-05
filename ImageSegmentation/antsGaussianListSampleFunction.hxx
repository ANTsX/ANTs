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
#ifndef __antsGaussianListSampleFunction_hxx
#define __antsGaussianListSampleFunction_hxx


#include "itkCovarianceSampleFilter.h"
#include "itkMeanSampleFilter.h"
#include "itkWeightedCovarianceSampleFilter.h"
#include "itkWeightedMeanSampleFilter.h"

namespace itk
{
namespace ants
{
namespace Statistics
{
template <typename TListSample, typename TOutput, typename TCoordRep>
GaussianListSampleFunction<TListSample, TOutput, TCoordRep>::GaussianListSampleFunction()
{
  this->m_Gaussian = GaussianType::New();
}

template <typename TListSample, typename TOutput, typename TCoordRep>
GaussianListSampleFunction<TListSample, TOutput, TCoordRep>::~GaussianListSampleFunction() = default;

template <typename TListSample, typename TOutput, typename TCoordRep>
void
GaussianListSampleFunction<TListSample, TOutput, TCoordRep>::SetInputListSample(const InputListSampleType * ptr)
{
  Superclass::SetInputListSample(ptr);

  if (!this->GetInputListSample())
  {
    return;
  }

  if (this->GetInputListSample()->Size() > 1)
  {
    if (this->GetListSampleWeights()->Size() == this->GetInputListSample()->Size())
    {
      typedef typename itk::Statistics::WeightedCovarianceSampleFilter<InputListSampleType> CovarianceCalculatorType;
      typename CovarianceCalculatorType::Pointer covarianceCalculator = CovarianceCalculatorType::New();

      covarianceCalculator->SetWeights(*this->GetListSampleWeights());
      covarianceCalculator->SetInput(this->GetInputListSample());
      covarianceCalculator->Update();

      typename GaussianType::MeanVectorType mean;
      NumericTraits<typename GaussianType::MeanVectorType>::SetLength(
        mean, this->GetInputListSample()->GetMeasurementVectorSize());
      for (unsigned int d = 0; d < this->GetInputListSample()->GetMeasurementVectorSize(); d++)
      {
        mean[d] = covarianceCalculator->GetMean()[d];
      }
      this->m_Gaussian->SetMean(mean);
      this->m_Gaussian->SetCovariance(covarianceCalculator->GetCovarianceMatrix());
    }
    else
    {
      typedef itk::Statistics::CovarianceSampleFilter<InputListSampleType> CovarianceCalculatorType;
      typename CovarianceCalculatorType::Pointer covarianceCalculator = CovarianceCalculatorType::New();
      covarianceCalculator->SetInput(this->GetInputListSample());
      covarianceCalculator->Update();

      typename GaussianType::MeanVectorType mean;
      NumericTraits<typename GaussianType::MeanVectorType>::SetLength(
        mean, this->GetInputListSample()->GetMeasurementVectorSize());
      for (unsigned int d = 0; d < this->GetInputListSample()->GetMeasurementVectorSize(); d++)
      {
        mean[d] = covarianceCalculator->GetMean()[d];
      }
      this->m_Gaussian->SetMean(mean);
      this->m_Gaussian->SetCovariance(covarianceCalculator->GetCovarianceMatrix());
    }

    // Check to see if the covariance matrix is nonsingular

    vnl_matrix_inverse<double> inv_cov((this->m_Gaussian->GetCovariance()).GetVnlMatrix());
    double                     det = inv_cov.determinant_magnitude();

    if (det < 0.0)
    {
      itkExceptionMacro("Determinant of the covariance < 0");
    }
    else if (det < 1.0e-6)
    {
      itkExceptionMacro("Covariance is singular (determinant = " << det << " < 1.0e-6)");
    }
  }
  else
  {
    itkWarningMacro("The input list sample has <= 1 element.  "
                    << "Function evaluations will be equal to 0.");
  }
}

template <typename TListSample, typename TOutput, typename TCoordRep>
TOutput
GaussianListSampleFunction<TListSample, TOutput, TCoordRep>::Evaluate(
  const InputMeasurementVectorType & measurement) const
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

/**
 * Standard "PrintSelf" method
 */
template <typename TListSample, typename TOutput, typename TCoordRep>
void
GaussianListSampleFunction<TListSample, TOutput, TCoordRep>::PrintSelf(std::ostream & os, Indent indent) const
{
  os << indent << "mean = " << this->m_Gaussian->GetMean() << ", ";

  typename GaussianType::CovarianceMatrixType covariance = this->m_Gaussian->GetCovariance();
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
