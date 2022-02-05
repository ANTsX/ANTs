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
#ifndef __antsLogEuclideanGaussianListSampleFunction_hxx
#define __antsLogEuclideanGaussianListSampleFunction_hxx


#include "itkDecomposeTensorFunction.h"

#include "vnl/vnl_trace.h"

namespace itk
{
namespace ants
{
namespace Statistics
{
template <typename TListSample, typename TOutput, typename TCoordRep>
LogEuclideanGaussianListSampleFunction<TListSample, TOutput, TCoordRep>::LogEuclideanGaussianListSampleFunction() =
  default;

template <typename TListSample, typename TOutput, typename TCoordRep>
LogEuclideanGaussianListSampleFunction<TListSample, TOutput, TCoordRep>::~LogEuclideanGaussianListSampleFunction() =
  default;

template <typename TListSample, typename TOutput, typename TCoordRep>
void
LogEuclideanGaussianListSampleFunction<TListSample, TOutput, TCoordRep>::SetInputListSample(
  const InputListSampleType * ptr)
{
  Superclass::SetInputListSample(ptr);

  if (!this->GetInputListSample())
  {
    return;
  }

  if (this->GetInputListSample()->Size() > 1)
  {
    RealType     L = static_cast<RealType>(this->GetInputListSample()->GetMeasurementVectorSize());
    unsigned int D = static_cast<unsigned int>(
      static_cast<RealType>(0.5) * (-NumericTraits<RealType>::OneValue() +
                                    std::sqrt(NumericTraits<RealType>::OneValue() + static_cast<RealType>(8.0) * L)));
    this->m_MeanTensor.SetSize(D, D);
    this->m_MeanTensor.Fill(0.0);

    unsigned long N = 0;
    RealType      totalWeight = 0.0;

    typename InputListSampleType::ConstIterator It = this->GetInputListSample()->Begin();
    while (It != this->GetInputListSample()->End())
    {
      InputMeasurementVectorType measurement = It.GetMeasurementVector();

      TensorType T(D, D);

      unsigned int index = 0;
      for (unsigned int i = 0; i < D; i++)
      {
        for (unsigned int j = i; j < D; j++)
        {
          T(i, j) = measurement(index++);
          T(j, i) = T(i, j);
        }
      }
      T = this->LogTensorTransform(T);

      RealType weight = 1.0;
      if (this->GetListSampleWeights()->Size() == this->GetInputListSample()->Size())
      {
        weight = (*this->GetListSampleWeights())[N++];
      }
      totalWeight += weight;
      this->m_MeanTensor += (T * weight);
      ++It;
    }

    if (totalWeight > NumericTraits<RealType>::ZeroValue())
    {
      this->m_MeanTensor /= totalWeight;
    }
    this->m_MeanTensor = this->ExpTensorTransform(this->m_MeanTensor);

    /**
     * Now calculate the dispersion (i.e. variance)
     */
    this->m_Dispersion = 0.0;

    N = 0;

    It = this->GetInputListSample()->Begin();
    while (It != this->GetInputListSample()->End())
    {
      InputMeasurementVectorType measurement = It.GetMeasurementVector();

      TensorType T(D, D);

      unsigned int index = 0;
      for (unsigned int i = 0; i < D; i++)
      {
        for (unsigned int j = i; j < D; j++)
        {
          T(i, j) = measurement(index++);
          T(j, i) = T(i, j);
        }
      }
      RealType distance = this->CalculateTensorDistance(T, this->m_MeanTensor);

      RealType weight = 1.0;
      if (this->GetListSampleWeights()->Size() == this->GetInputListSample()->Size())
      {
        weight = (*this->GetListSampleWeights())[N++];
      }

      this->m_Dispersion += (weight * itk::Math::sqr(distance));
      ++It;
    }

    this->m_Dispersion /= static_cast<RealType>(N);
  }
  else
  {
    itkWarningMacro("The input list sample has <= 1 element."
                    << "Function evaluations will be equal to 0.");
  }
}

template <typename TListSample, typename TOutput, typename TCoordRep>
typename LogEuclideanGaussianListSampleFunction<TListSample, TOutput, TCoordRep>::TensorType
LogEuclideanGaussianListSampleFunction<TListSample, TOutput, TCoordRep>::LogTensorTransform(const TensorType & T) const
{
  TensorType V;
  TensorType W;

  TensorType Tc = T;

  typedef DecomposeTensorFunction<TensorType> DecomposerType;
  typename DecomposerType::Pointer            decomposer = DecomposerType::New();
  decomposer->EvaluateSymmetricEigenDecomposition(Tc, W, V);
  for (unsigned int i = 0; i < W.Rows(); i++)
  {
    if (W(i, i) > NumericTraits<typename TensorType::ComponentType>::ZeroValue())
    {
      W(i, i) = std::log(W(i, i));
    }
    else
    {
      W(i, i) = NumericTraits<typename TensorType::ComponentType>::ZeroValue();
    }
  }
  W *= V.GetTranspose();
  TensorType logT = V * W;
  return logT;
}

template <typename TListSample, typename TOutput, typename TCoordRep>
typename LogEuclideanGaussianListSampleFunction<TListSample, TOutput, TCoordRep>::TensorType
LogEuclideanGaussianListSampleFunction<TListSample, TOutput, TCoordRep>::ExpTensorTransform(const TensorType & T) const
{
  TensorType V;
  TensorType W;

  TensorType Tc = T;

  typedef DecomposeTensorFunction<TensorType> DecomposerType;
  typename DecomposerType::Pointer            decomposer = DecomposerType::New();
  decomposer->EvaluateSymmetricEigenDecomposition(Tc, W, V);
  for (unsigned int i = 0; i < W.Rows(); i++)
  {
    W(i, i) = std::exp(W(i, i));
  }
  W *= V.GetTranspose();
  TensorType expT = V * W;
  return expT;
}

template <typename TListSample, typename TOutput, typename TCoordRep>
typename LogEuclideanGaussianListSampleFunction<TListSample, TOutput, TCoordRep>::RealType
LogEuclideanGaussianListSampleFunction<TListSample, TOutput, TCoordRep>::CalculateTensorDistance(
  const TensorType & S,
  const TensorType & T) const
{
  TensorType logS = this->LogTensorTransform(S);
  TensorType logT = this->LogTensorTransform(T);

  TensorType diff = logS - logT;
  TensorType diffSq = diff * diff;
  RealType   distance = std::sqrt(vnl_trace((diffSq).GetVnlMatrix()));

  //  RealType distance = ( ( logS - logT ).GetVnlMatrix() ).frobenius_norm();
  return distance;
}

template <typename TListSample, typename TOutput, typename TCoordRep>
TOutput
LogEuclideanGaussianListSampleFunction<TListSample, TOutput, TCoordRep>::Evaluate(
  const InputMeasurementVectorType & measurement) const
{
  unsigned int D = this->m_MeanTensor.Rows();

  TensorType T(D, D);

  unsigned int index = 0;

  for (unsigned int i = 0; i < D; i++)
  {
    for (unsigned int j = i; j < D; j++)
    {
      T(i, j) = measurement(index++);
      T(j, i) = T(i, j);
    }
  }
  RealType distance = this->CalculateTensorDistance(T, this->m_MeanTensor);
  RealType preFactor =
    NumericTraits<RealType>::OneValue() /
    (std::sqrt(static_cast<RealType>(2.0) * static_cast<RealType>(itk::Math::pi) * this->m_Dispersion));
  RealType probability =
    preFactor * std::exp(static_cast<RealType>(-0.5) * itk::Math::sqr(distance) / this->m_Dispersion);

  return probability;
}

/**
 * Standard "PrintSelf" method
 */
template <typename TListSample, typename TOutput, typename TCoordRep>
void
LogEuclideanGaussianListSampleFunction<TListSample, TOutput, TCoordRep>::PrintSelf(std::ostream & os,
                                                                                   Indent         indent) const
{
  os << indent << "Mean tensor = [";
  for (unsigned int r = 0; r < this->m_MeanTensor.Rows(); r++)
  {
    for (unsigned int c = 0; c < this->m_MeanTensor.Cols() - 1; c++)
    {
      os << this->m_MeanTensor(r, c) << ", ";
    }
    if (r == this->m_MeanTensor.Rows() - 1)
    {
      os << this->m_MeanTensor(r, this->m_MeanTensor.Cols() - 1) << "], ";
    }
    else
    {
      os << this->m_MeanTensor(r, this->m_MeanTensor.Cols() - 1) << "; ";
    }
  }
  os << "Dispersion (variance) = " << this->m_Dispersion << std::endl;
}
} // end of namespace Statistics
} // end of namespace ants
} // end of namespace itk

#endif
