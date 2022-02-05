/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkManifoldParzenWindowsPointSetFunction_hxx
#define __itkManifoldParzenWindowsPointSetFunction_hxx


#include "vnl/vnl_vector.h"
#include "itkMath.h"

namespace itk
{
template <typename TPointSet, typename TOutput, typename TCoordRep>
ManifoldParzenWindowsPointSetFunction<TPointSet, TOutput, TCoordRep>::ManifoldParzenWindowsPointSetFunction()
{
  this->m_BucketSize = 4;
  this->m_CovarianceKNeighborhood = 5;

  this->m_EvaluationKNeighborhood = 50;

  this->m_SamplePoints = nullptr;
  this->m_KdTreeGenerator = nullptr;

  this->m_RegularizationSigma = 1.0;
  this->m_KernelSigma = 1.0;

  this->m_Normalize = true;
  this->m_UseAnisotropicCovariances = true;

  this->m_Randomizer = RandomizerType::New();
  this->m_Randomizer->SetSeed();
}

template <typename TPointSet, typename TOutput, typename TCoordRep>
ManifoldParzenWindowsPointSetFunction<TPointSet, TOutput, TCoordRep>::~ManifoldParzenWindowsPointSetFunction()
{}

template <typename TPointSet, typename TOutput, typename TCoordRep>
void
ManifoldParzenWindowsPointSetFunction<TPointSet, TOutput, TCoordRep>::SetInputPointSet(const InputPointSetType * ptr)
{
  this->m_PointSet = ptr;

  /**
   * Generate KdTree and create set of gaussians from input point set
   */
  this->m_SamplePoints = SampleType::New();
  this->m_SamplePoints->SetMeasurementVectorSize(Dimension);

  std::vector<typename GaussianType::Pointer> inputGaussians;
  inputGaussians.resize(this->GetInputPointSet()->GetNumberOfPoints());
  this->m_Gaussians.resize(this->GetInputPointSet()->GetNumberOfPoints());

  MeasurementVectorType mv;

  unsigned long                count = 0;
  PointsContainerConstIterator It = this->GetInputPointSet()->GetPoints()->Begin();
  while (It != this->GetInputPointSet()->GetPoints()->End())
  {
    PointType point = It.Value();

    typename GaussianType::MeanType mean(Dimension);
    for (unsigned int d = 0; d < Dimension; d++)
    {
      mv[d] = point[d];
      mean[d] = point[d];
    }
    this->m_SamplePoints->PushBack(mv);

    inputGaussians[count] = GaussianType::New();
    inputGaussians[count]->SetGenerateRandomSamples(false);
    inputGaussians[count]->SetSigma(this->m_KernelSigma);
    inputGaussians[count]->SetMean(mean);

    count++;
    ++It;
  }

  this->m_KdTreeGenerator = TreeGeneratorType::New();
  m_KdTreeGenerator->SetSample(this->m_SamplePoints);
  m_KdTreeGenerator->SetBucketSize(this->m_BucketSize);
  m_KdTreeGenerator->Update();

  /**
   * Calculate covariance matrices
   */

  It = this->GetInputPointSet()->GetPoints()->Begin();
  while (It != this->GetInputPointSet()->GetPoints()->End())
  {
    PointType     point = It.Value();
    unsigned long index = It.Index();

    this->m_Gaussians[index] = GaussianType::New();
    this->m_Gaussians[index]->SetGenerateRandomSamples(true);
    this->m_Gaussians[index]->SetMean(inputGaussians[index]->GetMean());

    if (this->m_CovarianceKNeighborhood > 0 && this->m_UseAnisotropicCovariances)
    {
      CovarianceMatrixType Cout(Dimension, Dimension);
      Cout.Fill(0);

      MeasurementVectorType queryPoint;
      for (unsigned int d = 0; d < Dimension; d++)
      {
        queryPoint[d] = point[d];
      }

      typename TreeGeneratorType::KdTreeType ::InstanceIdentifierVectorType neighbors;
      this->m_KdTreeGenerator->GetOutput()->Search(queryPoint, this->m_CovarianceKNeighborhood, neighbors);

      RealType denominator = 0.0;
      for (unsigned int j = 0; j < this->m_CovarianceKNeighborhood; j++)
      {
        if (neighbors[j] != index && neighbors[j] < this->GetInputPointSet()->GetNumberOfPoints())
        {
          MeasurementVectorType neighbor = this->m_KdTreeGenerator->GetOutput()->GetMeasurementVector(neighbors[j]);

          RealType kernelValue = inputGaussians[index]->Evaluate(neighbor);

          denominator += kernelValue;
          if (kernelValue > 0.0)
          {
            for (unsigned int m = 0; m < Dimension; m++)
            {
              for (unsigned int n = m; n < Dimension; n++)
              {
                RealType covariance = kernelValue * (neighbor[m] - queryPoint[m]) * (neighbor[n] - queryPoint[n]);
                Cout(m, n) += covariance;
                Cout(n, m) += covariance;
              }
            }
          }
        }
      }
      if (this->m_Normalize && denominator > 0.0)
      {
        Cout /= denominator;
      }
      else
      {
        Cout /= static_cast<RealType>(this->m_CovarianceKNeighborhood);
      }
      for (unsigned int m = 0; m < Dimension; m++)
      {
        Cout(m, m) += (this->m_RegularizationSigma * this->m_RegularizationSigma);
      }

      this->m_Gaussians[index]->SetCovariance(Cout);
    }
    else
    {
      this->m_Gaussians[index]->SetSigma(this->m_RegularizationSigma);
    }
    ++It;
  }
}

template <typename TPointSet, typename TOutput, typename TCoordRep>
void
ManifoldParzenWindowsPointSetFunction<TPointSet, TOutput, TCoordRep>::GenerateKdTree()
{
  /**
   * Generate KdTree and create set of gaussians from input point set
   */
  this->m_SamplePoints = SampleType::New();
  this->m_SamplePoints->SetMeasurementVectorSize(Dimension);

  MeasurementVectorType mv;

  typename GaussianContainerType::const_iterator it;
  for (it = this->m_Gaussians.begin(); it != this->m_Gaussians.end(); ++it)
  {
    typename GaussianType::MeanType mean = (*it)->GetMean();
    for (unsigned int d = 0; d < Dimension; d++)
    {
      mv[d] = mean[d];
    }
    this->m_SamplePoints->PushBack(mv);
  }

  this->m_KdTreeGenerator = TreeGeneratorType::New();
  this->m_KdTreeGenerator->SetSample(this->m_SamplePoints);
  this->m_KdTreeGenerator->SetBucketSize(this->m_BucketSize);
  this->m_KdTreeGenerator->Update();
}

template <typename TPointSet, typename TOutput, typename TCoordRep>
TOutput
ManifoldParzenWindowsPointSetFunction<TPointSet, TOutput, TCoordRep>::Evaluate(const InputPointType & point) const
{
  if (!this->m_KdTreeGenerator)
  {
    MeasurementVectorType measurement;
    for (unsigned int d = 0; d < Dimension; d++)
    {
      measurement[d] = point[d];
    }

    OutputType                                     sum = 0.0;
    typename GaussianContainerType::const_iterator it;
    for (it = this->m_Gaussians.begin(); it != this->m_Gaussians.end(); ++it)
    {
      sum += static_cast<OutputType>((*it)->Evaluate(measurement));
    }
    return static_cast<OutputType>(sum / static_cast<OutputType>(this->m_Gaussians.size()));
  }
  else
  {
    MeasurementVectorType queryPoint;
    for (unsigned int d = 0; d < Dimension; d++)
    {
      queryPoint[d] = point[d];
    }

    unsigned int numberOfNeighbors =
      std::min(this->m_EvaluationKNeighborhood, static_cast<unsigned int>(this->m_Gaussians.size()));

    OutputType sum = 0.0;

    if (numberOfNeighbors == this->m_Gaussians.size())
    {
      for (unsigned int j = 0; j < this->m_Gaussians.size(); j++)
      {
        sum += static_cast<OutputType>(this->m_Gaussians[j]->Evaluate(queryPoint));
      }
    }
    else
    {
      typename TreeGeneratorType::KdTreeType ::InstanceIdentifierVectorType neighbors;
      this->m_KdTreeGenerator->GetOutput()->Search(queryPoint, static_cast<unsigned int>(numberOfNeighbors), neighbors);
      for (unsigned int j = 0; j < numberOfNeighbors; j++)
      {
        sum += static_cast<OutputType>(this->m_Gaussians[neighbors[j]]->Evaluate(queryPoint));
      }
    }
    return static_cast<OutputType>(sum / static_cast<OutputType>(this->m_Gaussians.size()));
  }
}

template <typename TPointSet, typename TOutput, typename TCoordRep>
typename ManifoldParzenWindowsPointSetFunction<TPointSet, TOutput, TCoordRep>::NeighborhoodIdentifierType
ManifoldParzenWindowsPointSetFunction<TPointSet, TOutput, TCoordRep>::GetNeighborhoodIdentifiers(
  MeasurementVectorType point,
  unsigned int          numberOfNeighbors)
{
  if (numberOfNeighbors > this->m_KdTreeGenerator->GetOutput()->Size())
  {
    numberOfNeighbors = this->m_KdTreeGenerator->GetOutput()->Size();
  }

  NeighborhoodIdentifierType neighbors;
  this->m_KdTreeGenerator->GetOutput()->Search(point, numberOfNeighbors, neighbors);
  return neighbors;
}

template <typename TPointSet, typename TOutput, typename TCoordRep>
typename ManifoldParzenWindowsPointSetFunction<TPointSet, TOutput, TCoordRep>::NeighborhoodIdentifierType
ManifoldParzenWindowsPointSetFunction<TPointSet, TOutput, TCoordRep>::GetNeighborhoodIdentifiers(
  InputPointType point,
  unsigned int   numberOfNeighbors)
{
  MeasurementVectorType queryPoint(Dimension);

  for (unsigned int d = 0; d < Dimension; d++)
  {
    queryPoint[d] = point[d];
  }
  return this->GetNeighborhoodIdentifiers(queryPoint, numberOfNeighbors);
}

template <typename TPointSet, typename TOutput, typename TCoordRep>
typename ManifoldParzenWindowsPointSetFunction<TPointSet, TOutput, TCoordRep>::PointType
ManifoldParzenWindowsPointSetFunction<TPointSet, TOutput, TCoordRep>::GenerateRandomSample()
{
  typename GaussianType::MeasurementVectorType gaussianSample;
  gaussianSample =
    this->m_Gaussians[this->m_Randomizer->GetIntegerVariate(this->GetInputPointSet()->GetNumberOfPoints() - 1)]
      ->GenerateRandomSample();

  PointType sample;
  for (unsigned int d = 0; d < Dimension; d++)
  {
    sample[d] = gaussianSample[d];
  }

  return sample;
}

/**
 * Standard "PrintSelf" method
 */
template <typename TPointSet, typename TOutput, typename TCoordRep>
void
ManifoldParzenWindowsPointSetFunction<TPointSet, TOutput, TCoordRep>::PrintSelf(std::ostream & os, Indent indent) const
{
  os << indent << "Covariance: " << this->m_CovarianceKNeighborhood << std::endl;
  os << indent << "Evaluation: " << this->m_EvaluationKNeighborhood << std::endl;
  os << indent << "Regularization sigma: " << this->m_RegularizationSigma << std::endl;
  os << indent << "Kernel sigma: " << this->m_KernelSigma << std::endl;
  os << indent << "Bucket size: " << this->m_BucketSize << std::endl;
  os << indent << "Normalize: " << this->m_Normalize << std::endl;
  os << indent << "Use anisotropic covariances: " << this->m_UseAnisotropicCovariances << std::endl;
}
} // end namespace itk

#endif
