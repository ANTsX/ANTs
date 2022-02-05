/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkJensenHavrdaCharvatTsallisPointSetMetric_hxx
#define __itkJensenHavrdaCharvatTsallisPointSetMetric_hxx


namespace itk
{
template <typename TPointSet>
JensenHavrdaCharvatTsallisPointSetMetric<TPointSet>::JensenHavrdaCharvatTsallisPointSetMetric()
{
  this->m_UseRegularizationTerm = false;
  this->m_UseInputAsSamples = true;
  this->m_UseAnisotropicCovariances = false;

  this->m_NumberOfFixedSamples = 100;
  this->m_FixedPointSetSigma = 1.0;
  this->m_FixedKernelSigma = 10.0;
  this->m_FixedCovarianceKNeighborhood = 5;
  this->m_FixedEvaluationKNeighborhood = 50;

  this->m_NumberOfMovingSamples = 100;
  this->m_MovingPointSetSigma = 1.0;
  this->m_MovingKernelSigma = 10.0;
  this->m_MovingCovarianceKNeighborhood = 5;
  this->m_MovingEvaluationKNeighborhood = 50;

  this->m_Alpha = 2.0;
  this->m_UseWithRespectToTheMovingPointSet = true;

  typename DefaultTransformType::Pointer transform = DefaultTransformType::New();
  transform->SetIdentity();

  Superclass::SetTransform(transform);
}

/** Initialize the metric */
template <typename TPointSet>
void
JensenHavrdaCharvatTsallisPointSetMetric<TPointSet>::Initialize(void) throw(ExceptionObject)
{
  Superclass::Initialize();

  /**
   * Initialize the fixed points
   */
  this->m_FixedDensityFunction = DensityFunctionType::New();
  this->m_FixedDensityFunction->SetBucketSize(4);
  this->m_FixedDensityFunction->SetKernelSigma(this->m_FixedKernelSigma);
  this->m_FixedDensityFunction->SetRegularizationSigma(this->m_FixedPointSetSigma);
  this->m_FixedDensityFunction->SetNormalize(true);
  this->m_FixedDensityFunction->SetUseAnisotropicCovariances(this->m_UseAnisotropicCovariances);
  this->m_FixedDensityFunction->SetCovarianceKNeighborhood(this->m_FixedCovarianceKNeighborhood);
  this->m_FixedDensityFunction->SetEvaluationKNeighborhood(this->m_FixedEvaluationKNeighborhood);
  this->m_FixedDensityFunction->SetInputPointSet(this->m_FixedPointSet);

  if (!this->m_UseInputAsSamples)
  {
    this->m_FixedSamplePoints = PointSetType::New();
    this->m_FixedSamplePoints->Initialize();
    for (unsigned long i = 0; i < this->m_NumberOfFixedSamples; i++)
    {
      this->m_FixedSamplePoints->SetPoint(i, this->m_FixedDensityFunction->GenerateRandomSample());
    }
  }

  /**
   * Initialize the moving points
   */
  this->m_MovingDensityFunction = DensityFunctionType::New();
  this->m_MovingDensityFunction->SetBucketSize(4);
  this->m_MovingDensityFunction->SetKernelSigma(this->m_MovingKernelSigma);
  this->m_MovingDensityFunction->SetRegularizationSigma(this->m_MovingPointSetSigma);
  this->m_MovingDensityFunction->SetNormalize(true);
  this->m_MovingDensityFunction->SetUseAnisotropicCovariances(this->m_UseAnisotropicCovariances);
  this->m_MovingDensityFunction->SetCovarianceKNeighborhood(this->m_MovingCovarianceKNeighborhood);
  this->m_MovingDensityFunction->SetEvaluationKNeighborhood(this->m_MovingEvaluationKNeighborhood);
  this->m_MovingDensityFunction->SetInputPointSet(this->m_MovingPointSet);

  if (!this->m_UseInputAsSamples)
  {
    this->m_MovingSamplePoints = PointSetType::New();
    this->m_MovingSamplePoints->Initialize();
    for (unsigned long i = 0; i < this->m_NumberOfMovingSamples; i++)
    {
      this->m_MovingSamplePoints->SetPoint(i, this->m_MovingDensityFunction->GenerateRandomSample());
    }
  }
}

/** Return the number of values, i.e the number of points in the moving set */
template <typename TPointSet>
unsigned int
JensenHavrdaCharvatTsallisPointSetMetric<TPointSet>::GetNumberOfValues() const
{
  if (this->m_UseWithRespectToTheMovingPointSet)
  {
    if (this->m_MovingPointSet)
    {
      return this->m_MovingPointSet->GetPoints()->Size();
    }
  }
  else
  {
    if (this->m_FixedPointSet)
    {
      return this->m_FixedPointSet->GetPoints()->Size();
    }
  }

  return 0;
}

/** Get the match Measure */
template <typename TPointSet>
typename JensenHavrdaCharvatTsallisPointSetMetric<TPointSet>::MeasureType
JensenHavrdaCharvatTsallisPointSetMetric<TPointSet>::GetValue(const TransformParametersType & parameters) const
{
  /**
   * Only identity transform is valid
   */
  //  this->SetTransformParameters( parameters );

  PointSetPointer points[2];
  PointSetPointer samples[2];

  typename DensityFunctionType::Pointer densityFunctions[2];

  if (this->m_UseWithRespectToTheMovingPointSet)
  {
    points[0] = const_cast<PointSetType *>(static_cast<const PointSetType *>(this->m_FixedPointSet));
    points[1] = const_cast<PointSetType *>(static_cast<const PointSetType *>(this->m_MovingPointSet));

    if (this->m_UseInputAsSamples)
    {
      samples[0] = points[0];
      samples[1] = points[1];
    }
    else
    {
      samples[0] = this->m_FixedSamplePoints;
      samples[1] = this->m_MovingSamplePoints;
    }
    densityFunctions[0] = this->m_FixedDensityFunction;
    densityFunctions[1] = this->m_MovingDensityFunction;
  }
  else
  {
    points[1] = const_cast<PointSetType *>(static_cast<const PointSetType *>(this->m_FixedPointSet));
    points[0] = const_cast<PointSetType *>(static_cast<const PointSetType *>(this->m_MovingPointSet));

    if (this->m_UseInputAsSamples)
    {
      samples[1] = points[1];
      samples[0] = points[0];
    }
    else
    {
      samples[1] = this->m_FixedSamplePoints;
      samples[0] = this->m_MovingSamplePoints;
    }
    densityFunctions[1] = this->m_FixedDensityFunction;
    densityFunctions[0] = this->m_MovingDensityFunction;
  }

  RealType totalNumberOfPoints =
    static_cast<RealType>(points[0]->GetNumberOfPoints()) + static_cast<RealType>(points[1]->GetNumberOfPoints());
  RealType totalNumberOfSamples =
    static_cast<RealType>(samples[0]->GetNumberOfPoints()) + static_cast<RealType>(samples[1]->GetNumberOfPoints());

  MeasureType measure;
  measure.SetSize(1);
  measure.Fill(0);

  RealType energyTerm1 = 0.0;
  RealType energyTerm2 = 0.0;

  /**
   * first term
   */
  RealType prefactor = -1.0 / totalNumberOfSamples;
  if (this->m_Alpha != 1.0)
  {
    prefactor /= (this->m_Alpha - 1.0);
  }
  typename PointSetType::PointsContainerConstIterator It = samples[0]->GetPoints()->Begin();
  while (It != samples[0]->GetPoints()->End())
  {
    PointType samplePoint = It.Value();

    RealType probabilityStar =
      //         densityFunctions[0]->Evaluate( samplePoint ) *
      //         static_cast<RealType>( points[0]->GetNumberOfPoints() ) +
      densityFunctions[1]->Evaluate(samplePoint) * static_cast<RealType>(points[1]->GetNumberOfPoints());
    probabilityStar /= totalNumberOfPoints;

    if (probabilityStar == 0)
    {
      ++It;
      continue;
    }

    if (this->m_Alpha == 1.0)
    {
      energyTerm1 += std::log(probabilityStar);
    }
    else
    {
      energyTerm1 += std::pow(probabilityStar, static_cast<RealType>(this->m_Alpha - 1.0));
    }
    ++It;
  }

  if (this->m_Alpha != 1.0)
  {
    energyTerm1 -= 1.0;
  }
  energyTerm1 *= prefactor;

  /**
   * second term, i.e. regularization term
   */
  if (this->m_UseRegularizationTerm)
  {
    RealType prefactor2 = -static_cast<RealType>(points[1]->GetNumberOfPoints()) /
                          (totalNumberOfPoints * static_cast<RealType>(samples[1]->GetNumberOfPoints()));
    if (this->m_Alpha != 1.0)
    {
      prefactor2 /= (this->m_Alpha - 1.0);
    }
    typename PointSetType::PointsContainerConstIterator It = samples[1]->GetPoints()->Begin();
    while (It != samples[1]->GetPoints()->End())
    {
      PointType samplePoint = It.Value();

      RealType probability = densityFunctions[1]->Evaluate(samplePoint);

      if (probability == 0)
      {
        ++It;
        continue;
      }

      if (this->m_Alpha == 1.0)
      {
        energyTerm2 += (prefactor2 * std::log(probability));
      }
      else
      {
        energyTerm2 += (prefactor2 * std::pow(probability, static_cast<RealType>(this->m_Alpha - 1.0)));
      }
      ++It;
    }

    if (this->m_Alpha != 1.0)
    {
      energyTerm2 -= 1.0;
    }
    energyTerm2 *= prefactor2;
  }

  measure[0] = energyTerm1 - energyTerm2;

  return measure;
}

/** Get the Derivative Measure */
template <typename TPointSet>
void
JensenHavrdaCharvatTsallisPointSetMetric<TPointSet>::GetDerivative(const TransformParametersType & parameters,
                                                                   DerivativeType &                derivative) const
{
  /**
   * Only identity transform is valid
   */
  //  this->SetTransformParameters( parameters );

  PointSetPointer points[2];
  PointSetPointer samples[2];

  typename DensityFunctionType::Pointer densityFunctions[2];

  unsigned int kNeighborhood;

  if (this->m_UseWithRespectToTheMovingPointSet)
  {
    points[0] = const_cast<PointSetType *>(static_cast<const PointSetType *>(this->m_FixedPointSet));
    points[1] = const_cast<PointSetType *>(static_cast<const PointSetType *>(this->m_MovingPointSet));

    if (this->m_UseInputAsSamples)
    {
      samples[0] = points[0];
      samples[1] = points[1];
    }
    else
    {
      samples[0] = this->m_FixedSamplePoints;
      samples[1] = this->m_MovingSamplePoints;
    }
    densityFunctions[0] = this->m_FixedDensityFunction;
    densityFunctions[1] = this->m_MovingDensityFunction;

    kNeighborhood = this->m_MovingEvaluationKNeighborhood;
  }
  else
  {
    points[1] = const_cast<PointSetType *>(static_cast<const PointSetType *>(this->m_FixedPointSet));
    points[0] = const_cast<PointSetType *>(static_cast<const PointSetType *>(this->m_MovingPointSet));

    if (this->m_UseInputAsSamples)
    {
      samples[1] = points[1];
      samples[0] = points[0];
    }
    else
    {
      samples[1] = this->m_FixedSamplePoints;
      samples[0] = this->m_MovingSamplePoints;
    }
    densityFunctions[1] = this->m_FixedDensityFunction;
    densityFunctions[0] = this->m_MovingDensityFunction;

    kNeighborhood = this->m_FixedEvaluationKNeighborhood;
  }
  RealType totalNumberOfPoints =
    static_cast<RealType>(points[0]->GetNumberOfPoints()) + static_cast<RealType>(points[1]->GetNumberOfPoints());

  RealType totalNumberOfSamples =
    static_cast<RealType>(samples[0]->GetNumberOfPoints()) + static_cast<RealType>(samples[1]->GetNumberOfPoints());

  derivative.SetSize(points[1]->GetPoints()->Size(), PointDimension);
  derivative.Fill(0);

  /**
   * first term
   */

  RealType prefactor = 1.0 / (totalNumberOfSamples * totalNumberOfPoints);

  typename PointSetType::PointsContainerConstIterator It = samples[0]->GetPoints()->Begin();
  while (It != samples[0]->GetPoints()->End())
  {
    PointType fixedSamplePoint = It.Value();

    RealType probabilityStar =
      //         densityFunctions[0]->Evaluate( fixedSamplePoint ) *
      //         static_cast<RealType>( points[0]->GetNumberOfPoints() ) +
      densityFunctions[1]->Evaluate(fixedSamplePoint) * static_cast<RealType>(points[1]->GetNumberOfPoints());
    probabilityStar /= totalNumberOfPoints;

    if (probabilityStar == 0)
    {
      ++It;
      continue;
    }

    RealType probabilityStarFactor = std::pow(probabilityStar, static_cast<RealType>(2.0 - this->m_Alpha));

    typename GaussianType::MeasurementVectorType sampleMeasurement;
    for (unsigned int d = 0; d < PointDimension; d++)
    {
      sampleMeasurement[d] = fixedSamplePoint[d];
    }

    typename DensityFunctionType::NeighborhoodIdentifierType neighbors =
      densityFunctions[1]->GetNeighborhoodIdentifiers(sampleMeasurement, kNeighborhood);
    for (unsigned int n = 0; n < neighbors.size(); n++)
    {
      RealType gaussian = densityFunctions[1]->GetGaussian(neighbors[n])->Evaluate(sampleMeasurement);

      if (gaussian == 0)
      {
        continue;
      }

      typename GaussianType::MeanType mean = densityFunctions[1]->GetGaussian(neighbors[n])->GetMean();
      for (unsigned int d = 0; d < PointDimension; d++)
      {
        mean[d] -= fixedSamplePoint[d];
      }

      if (this->m_UseAnisotropicCovariances)
      {
        typename GaussianType::MatrixType Ci = densityFunctions[1]->GetGaussian(neighbors[n])->GetInverseCovariance();
        mean = Ci * mean;
      }
      else
      {
        mean /= itk::Math::sqr(densityFunctions[1]->GetGaussian(neighbors[n])->GetSigma());
      }

      mean *= (prefactor * gaussian / probabilityStarFactor);
      for (unsigned int d = 0; d < PointDimension; d++)
      {
        derivative(neighbors[n], d) += mean[d];
      }
    }
    ++It;
  }

  /**
   * second term, i.e. regularization term
   */
  if (this->m_UseRegularizationTerm)
  {
    RealType prefactor2 = -1.0 / (static_cast<RealType>(samples[1]->GetNumberOfPoints()) * totalNumberOfPoints);

    typename PointSetType::PointsContainerConstIterator It = samples[1]->GetPoints()->Begin();
    while (It != samples[1]->GetPoints()->End())
    {
      PointType movingSamplePoint = It.Value();

      RealType probability = densityFunctions[1]->Evaluate(movingSamplePoint);

      if (probability == 0)
      {
        ++It;
        continue;
      }

      RealType probabilityFactor = std::pow(probability, static_cast<RealType>(2.0 - this->m_Alpha));
      probabilityFactor *= (samples[1]->GetNumberOfPoints() / totalNumberOfSamples);

      typename GaussianType::MeasurementVectorType sampleMeasurement;
      for (unsigned int d = 0; d < PointDimension; d++)
      {
        sampleMeasurement[d] = movingSamplePoint[d];
      }

      typename DensityFunctionType::NeighborhoodIdentifierType neighbors =
        densityFunctions[1]->GetNeighborhoodIdentifiers(sampleMeasurement, kNeighborhood);
      for (unsigned int i = 0; i < neighbors.size(); i++)
      {
        RealType gaussian = densityFunctions[1]->GetGaussian(neighbors[i])->Evaluate(sampleMeasurement);
        if (gaussian == 0)
        {
          continue;
        }

        typename GaussianType::MeanType mean = densityFunctions[1]->GetGaussian(neighbors[i])->GetMean();
        for (unsigned int d = 0; d < PointDimension; d++)
        {
          mean[d] -= movingSamplePoint[d];
        }

        if (this->m_UseAnisotropicCovariances)
        {
          typename GaussianType::MatrixType Ci = densityFunctions[1]->GetGaussian(neighbors[i])->GetInverseCovariance();
          mean = Ci * mean;
        }
        else
        {
          mean /= itk::Math::sqr(densityFunctions[1]->GetGaussian(neighbors[i])->GetSigma());
        }

        mean *= (prefactor2 * gaussian / probabilityFactor);
        for (unsigned int d = 0; d < PointDimension; d++)
        {
          derivative(neighbors[i], d) += mean[d];
        }
      }
      ++It;
    }
  }
}

/** Get both the match Measure and theDerivative Measure  */
template <typename TPointSet>
void
JensenHavrdaCharvatTsallisPointSetMetric<TPointSet>::GetValueAndDerivative(const TransformParametersType & parameters,
                                                                           MeasureType &                   value,
                                                                           DerivativeType & derivative) const
{
  /**
   * Only identity transform is valid
   */
  //  this->SetTransformParameters( parameters );

  PointSetPointer points[2];
  PointSetPointer samples[2];

  typename DensityFunctionType::Pointer densityFunctions[2];

  unsigned int kNeighborhood;

  if (this->m_UseWithRespectToTheMovingPointSet)
  {
    points[0] = const_cast<PointSetType *>(static_cast<const PointSetType *>(this->m_FixedPointSet));
    points[1] = const_cast<PointSetType *>(static_cast<const PointSetType *>(this->m_MovingPointSet));

    if (this->m_UseInputAsSamples)
    {
      samples[0] = points[0];
      samples[1] = points[1];
    }
    else
    {
      samples[0] = this->m_FixedSamplePoints;
      samples[1] = this->m_MovingSamplePoints;
    }
    densityFunctions[0] = this->m_FixedDensityFunction;
    densityFunctions[1] = this->m_MovingDensityFunction;

    kNeighborhood = this->m_MovingEvaluationKNeighborhood;
  }
  else
  {
    points[1] = const_cast<PointSetType *>(static_cast<const PointSetType *>(this->m_FixedPointSet));
    points[0] = const_cast<PointSetType *>(static_cast<const PointSetType *>(this->m_MovingPointSet));

    if (this->m_UseInputAsSamples)
    {
      samples[1] = points[1];
      samples[0] = points[0];
    }
    else
    {
      samples[1] = this->m_FixedSamplePoints;
      samples[0] = this->m_MovingSamplePoints;
    }
    densityFunctions[1] = this->m_FixedDensityFunction;
    densityFunctions[0] = this->m_MovingDensityFunction;

    kNeighborhood = this->m_FixedEvaluationKNeighborhood;
  }
  RealType totalNumberOfPoints =
    static_cast<RealType>(points[0]->GetNumberOfPoints()) + static_cast<RealType>(points[1]->GetNumberOfPoints());

  RealType totalNumberOfSamples =
    static_cast<RealType>(samples[0]->GetNumberOfPoints()) + static_cast<RealType>(samples[1]->GetNumberOfPoints());

  derivative.SetSize(points[1]->GetPoints()->Size(), PointDimension);
  derivative.Fill(0);

  value.SetSize(1);
  value.Fill(0);

  /**
   * first term
   */
  RealType energyTerm1 = 0.0;
  RealType energyTerm2 = 0.0;

  RealType prefactor[2];
  prefactor[0] = -1.0 / totalNumberOfSamples;
  if (this->m_Alpha != 1.0)
  {
    prefactor[0] /= (this->m_Alpha - 1.0);
  }
  prefactor[1] = 1.0 / (totalNumberOfSamples * totalNumberOfPoints);

  typename PointSetType::PointsContainerConstIterator It = samples[0]->GetPoints()->Begin();
  while (It != samples[0]->GetPoints()->End())
  {
    PointType fixedSamplePoint = It.Value();

    RealType probabilityStar =
      //       densityFunctions[0]->Evaluate( fixedSamplePoint ) *
      //       static_cast<RealType>( points[0]->GetNumberOfPoints() ) +
      densityFunctions[1]->Evaluate(fixedSamplePoint) * static_cast<RealType>(points[1]->GetNumberOfPoints());

    probabilityStar /= totalNumberOfPoints;

    if (probabilityStar == 0)
    {
      ++It;
      continue;
    }

    if (this->m_Alpha == 1.0)
    {
      energyTerm1 += (prefactor[0] * std::log(probabilityStar));
    }
    else
    {
      energyTerm1 += (prefactor[0] * std::pow(probabilityStar, static_cast<RealType>(this->m_Alpha - 1.0)));
    }

    RealType probabilityStarFactor = std::pow(probabilityStar, static_cast<RealType>(2.0 - this->m_Alpha));

    typename GaussianType::MeasurementVectorType sampleMeasurement;
    for (unsigned int d = 0; d < PointDimension; d++)
    {
      sampleMeasurement[d] = fixedSamplePoint[d];
    }

    typename DensityFunctionType::NeighborhoodIdentifierType neighbors =
      densityFunctions[1]->GetNeighborhoodIdentifiers(sampleMeasurement, kNeighborhood);
    for (unsigned int n = 0; n < neighbors.size(); n++)
    {
      RealType gaussian = densityFunctions[1]->GetGaussian(neighbors[n])->Evaluate(sampleMeasurement);

      if (gaussian == 0)
      {
        continue;
      }

      typename GaussianType::MeanType mean = densityFunctions[1]->GetGaussian(neighbors[n])->GetMean();
      for (unsigned int d = 0; d < PointDimension; d++)
      {
        mean[d] -= fixedSamplePoint[d];
      }

      if (this->m_UseAnisotropicCovariances)
      {
        typename GaussianType::MatrixType Ci = densityFunctions[1]->GetGaussian(neighbors[n])->GetInverseCovariance();
        mean = Ci * mean;
      }
      else
      {
        mean /= itk::Math::sqr(densityFunctions[1]->GetGaussian(neighbors[n])->GetSigma());
      }

      mean *= (prefactor[1] * gaussian / probabilityStarFactor);
      for (unsigned int d = 0; d < PointDimension; d++)
      {
        derivative(neighbors[n], d) += mean[d];
      }
    }
    ++It;
  }

  if (this->m_Alpha != 1.0)
  {
    energyTerm1 -= 1.0;
  }
  energyTerm1 *= prefactor[0];

  /**
   * second term, i.e. regularization term
   */
  if (this->m_UseRegularizationTerm)
  {
    RealType prefactor2[2];
    prefactor2[0] = -static_cast<RealType>(points[1]->GetNumberOfPoints()) /
                    (totalNumberOfPoints * static_cast<RealType>(samples[1]->GetNumberOfPoints()));
    prefactor2[1] = -1.0 / (static_cast<RealType>(samples[1]->GetNumberOfPoints()) * totalNumberOfPoints);
    if (this->m_Alpha != 1.0)
    {
      prefactor2[0] /= (this->m_Alpha - 1.0);
    }

    typename PointSetType::PointsContainerConstIterator It = samples[1]->GetPoints()->Begin();
    while (It != samples[1]->GetPoints()->End())
    {
      PointType movingSamplePoint = It.Value();

      RealType probability = densityFunctions[1]->Evaluate(movingSamplePoint);

      if (probability == 0)
      {
        ++It;
        continue;
      }

      if (this->m_Alpha == 1.0)
      {
        energyTerm2 += (prefactor2[0] * std::log(probability));
      }
      else
      {
        energyTerm2 += (prefactor2[0] * std::pow(probability, static_cast<RealType>(this->m_Alpha - 1.0)));
      }

      RealType probabilityFactor = std::pow(probability, static_cast<RealType>(2.0 - this->m_Alpha));
      probabilityFactor *= (samples[1]->GetNumberOfPoints() / totalNumberOfSamples);

      typename GaussianType::MeasurementVectorType sampleMeasurement;
      for (unsigned int d = 0; d < PointDimension; d++)
      {
        sampleMeasurement[d] = movingSamplePoint[d];
      }

      typename DensityFunctionType::NeighborhoodIdentifierType neighbors =
        densityFunctions[1]->GetNeighborhoodIdentifiers(sampleMeasurement, kNeighborhood);
      for (unsigned int i = 0; i < neighbors.size(); i++)
      {
        RealType gaussian = densityFunctions[1]->GetGaussian(neighbors[i])->Evaluate(sampleMeasurement);
        if (gaussian == 0)
        {
          continue;
        }

        typename GaussianType::MeanType mean = densityFunctions[1]->GetGaussian(neighbors[i])->GetMean();
        for (unsigned int d = 0; d < PointDimension; d++)
        {
          mean[d] -= movingSamplePoint[d];
        }

        if (this->m_UseAnisotropicCovariances)
        {
          typename GaussianType::MatrixType Ci = densityFunctions[1]->GetGaussian(neighbors[i])->GetInverseCovariance();
          mean = Ci * mean;
        }
        else
        {
          mean /= itk::Math::sqr(densityFunctions[1]->GetGaussian(neighbors[i])->GetSigma());
        }

        mean *= (prefactor2[1] * gaussian / probabilityFactor);
        for (unsigned int d = 0; d < PointDimension; d++)
        {
          derivative(neighbors[i], d) += mean[d];
        }
      }
      ++It;
    }

    if (this->m_Alpha != 1.0)
    {
      energyTerm2 -= 1.0;
    }
    energyTerm2 *= prefactor2[0];
  }

  value[0] = energyTerm1 - energyTerm2;
}

template <typename TPointSet>
void
JensenHavrdaCharvatTsallisPointSetMetric<TPointSet>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Use with respect to the moving point set: " << this->m_UseWithRespectToTheMovingPointSet
     << std::endl;
  os << indent << "Use regularization term: " << this->m_UseRegularizationTerm << std::endl;
  os << indent << "Alpha: " << this->m_Alpha << std::endl;

  os << indent << "Fixed sigma: " << this->m_FixedPointSetSigma << std::endl;
  os << indent << "Moving sigma: " << this->m_MovingPointSetSigma << std::endl;

  if (!this->m_UseInputAsSamples)
  {
    os << indent << "Number of fixed samples: " << this->m_NumberOfFixedSamples << std::endl;
    os << indent << "Number of moving samples: " << this->m_NumberOfMovingSamples << std::endl;
  }
  else
  {
    os << indent << "Use input points as samples." << std::endl;
  }

  if (this->m_UseAnisotropicCovariances)
  {
    os << indent << "Fixed kernel sigma: " << this->m_FixedKernelSigma << std::endl;
    os << indent << "Moving kernel sigma: " << this->m_MovingKernelSigma << std::endl;
    os << indent << "Fixed covariance k-neighborhood: " << this->m_FixedCovarianceKNeighborhood << std::endl;
    os << indent << "Moving covariance k-neighborhood: " << this->m_MovingCovarianceKNeighborhood << std::endl;
  }
  else
  {
    os << indent << "Isotropic covariances are used." << std::endl;
  }
}
} // end namespace itk

#endif
