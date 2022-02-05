/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkSyNDemonsRegistrationFunction_hxx_
#define _itkSyNDemonsRegistrationFunction_hxx_

#include "itkSyNDemonsRegistrationFunction.h"
#include "itkMacro.h"
#include "itkMath.h"

namespace itk
{
/*
 * Default constructor
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
SyNDemonsRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::SyNDemonsRegistrationFunction()
{
  RadiusType   r;
  unsigned int j;

  for (j = 0; j < ImageDimension; j++)
  {
    r[j] = 0;
  }
  this->SetRadius(r);

  m_TimeStep = 1.0;
  m_DenominatorThreshold = 1e-9;
  m_IntensityDifferenceThreshold = 0.001;
  this->SetMovingImage(nullptr);
  this->SetFixedImage(nullptr);
  m_FixedImageSpacing.Fill(1.0);
  m_FixedImageOrigin.Fill(0.0);
  m_Normalizer = 1.0;
  m_FixedImageGradientCalculator = GradientCalculatorType::New();

  typename DefaultInterpolatorType::Pointer interp = DefaultInterpolatorType::New();

  m_MovingImageInterpolator = static_cast<InterpolatorType *>(interp.GetPointer());

  m_Metric = NumericTraits<double>::max();
  m_SumOfSquaredDifference = 0.0;
  m_NumberOfPixelsProcessed = 0L;
  m_RMSChange = NumericTraits<double>::max();
  m_SumOfSquaredChange = 0.0;

  m_MovingImageGradientCalculator = MovingImageGradientCalculatorType::New();
  m_UseMovingImageGradient = false;
  m_UseSSD = false;
}

/*
 * Standard "PrintSelf" method.
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
void
SyNDemonsRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::PrintSelf(std::ostream & os,
                                                                                        Indent         indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "MovingImageIterpolator: ";
  os << m_MovingImageInterpolator.GetPointer() << std::endl;
  os << indent << "FixedImageGradientCalculator: ";
  os << m_FixedImageGradientCalculator.GetPointer() << std::endl;
  os << indent << "DenominatorThreshold: ";
  os << m_DenominatorThreshold << std::endl;
  os << indent << "IntensityDifferenceThreshold: ";
  os << m_IntensityDifferenceThreshold << std::endl;

  os << indent << "UseMovingImageGradient: ";
  os << m_UseMovingImageGradient << std::endl;

  os << indent << "Metric: ";
  os << m_Metric << std::endl;
  os << indent << "SumOfSquaredDifference: ";
  os << m_SumOfSquaredDifference << std::endl;
  os << indent << "NumberOfPixelsProcessed: ";
  os << m_NumberOfPixelsProcessed << std::endl;
  os << indent << "RMSChange: ";
  os << m_RMSChange << std::endl;
  os << indent << "SumOfSquaredChange: ";
  os << m_SumOfSquaredChange << std::endl;
}

/**
 *
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
void
SyNDemonsRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::SetIntensityDifferenceThreshold(
  double threshold)
{
  m_IntensityDifferenceThreshold = threshold;
}

/**
 *
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
double
SyNDemonsRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::GetIntensityDifferenceThreshold() const
{
  return m_IntensityDifferenceThreshold;
}

/*
 * Set the function state values before each iteration
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
void
SyNDemonsRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::InitializeIteration()
{
  //  std::cout << " INIT ITER " << std::endl;
  if (!this->GetMovingImage() || !this->GetFixedImage() || !m_MovingImageInterpolator)
  {
    itkExceptionMacro(<< "MovingImage, FixedImage and/or Interpolator not set");
  }
  // cache fixed image information
  m_FixedImageSpacing = this->GetFixedImage()->GetSpacing();
  m_FixedImageOrigin = this->GetFixedImage()->GetOrigin();

  // compute the normalizer
  m_Normalizer = 0.0;
  for (unsigned int k = 0; k < ImageDimension; k++)
  {
    m_Normalizer += m_FixedImageSpacing[k] * m_FixedImageSpacing[k];
  }
  m_Normalizer /= static_cast<double>(ImageDimension);

  this->m_Energy = 0;

  // setup gradient calculator
  m_FixedImageGradientCalculator->SetInputImage(this->GetFixedImage());
  m_MovingImageGradientCalculator->SetInputImage(this->GetMovingImage());

  // setup moving image interpolator
  m_MovingImageInterpolator->SetInputImage(this->GetMovingImage());

  // initialize metric computation variables
  m_SumOfSquaredDifference = 0.0;
  m_NumberOfPixelsProcessed = 0L;
  m_SumOfSquaredChange = 0.0;
}

/*
 * Compute update at a specify neighbourhood
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
typename SyNDemonsRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::PixelType
SyNDemonsRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::ComputeUpdate(
  const NeighborhoodType & it,
  void *                   gd,
  const FloatOffsetType &  itkNotUsed(offset))
{
  PixelType    update;
  unsigned int j;

  IndexType index = it.GetIndex();

  // Get fixed image related information
  double              fixedValue;
  CovariantVectorType gradient;
  CovariantVectorType mgradient;
  double              gradientSquaredMagnitude = 0;

  // Note: no need to check the index is within
  // fixed image buffer. This is done by the external filter.
  fixedValue = (double)this->GetFixedImage()->GetPixel(index);
  double movingValue = (double)this->GetMovingImage()->GetPixel(index);

  //  if (fixedValue > 0)std::cout << " fxv  " << fixedValue << " movingValue " << movingValue << std::endl;

  gradient = m_FixedImageGradientCalculator->EvaluateAtIndex(index);

  mgradient = m_MovingImageGradientCalculator->EvaluateAtIndex(index);
  for (j = 0; j < ImageDimension; j++)
  {
    if (this->m_UseMovingImageGradient)
    {
      gradient[j] = gradient[j] + mgradient[j];
    }
    gradientSquaredMagnitude += itk::Math::sqr(gradient[j]);
  }

  /**
   * Compute Update.
   * In the original equation the denominator is defined as (g-f)^2 + grad_mag^2.
   * However there is a mismatch in units between the two terms.
   * The units for the second term is intensity^2/mm^2 while the
   * units for the first term is intensity^2. This mismatch is particularly
   * problematic when the fixed image does not have unit spacing.
   * In this implemenation, we normalize the first term by a factor K,
   * such that denominator = (g-f)^2/K + grad_mag^2
   * where K = mean square spacing to compensate for the mismatch in units.
   */
  double speedValue = fixedValue - movingValue;
  if (fabs(speedValue) < static_cast<double>(this->m_RobustnessParameter))
  {
    speedValue = 0;
  }

  // update the metric
  GlobalDataStruct * globalData = reinterpret_cast<GlobalDataStruct *>(gd);
  if (globalData)
  {
    globalData->m_SumOfSquaredDifference += itk::Math::sqr(speedValue);
    globalData->m_NumberOfPixelsProcessed += 1;
  }

  double denominator = itk::Math::sqr(speedValue) / m_Normalizer + gradientSquaredMagnitude;
  this->m_Energy += speedValue * speedValue;
  if (m_UseSSD)
  {
    denominator = 1;
  }
  if (itk::Math::abs(speedValue) < m_IntensityDifferenceThreshold || denominator < m_DenominatorThreshold)
  {
    for (j = 0; j < ImageDimension; j++)
    {
      update[j] = 0.0;
    }
    return update;
  }
  for (j = 0; j < ImageDimension; j++)
  {
    update[j] = speedValue * gradient[j] / denominator;
    if (globalData)
    {
      globalData->m_SumOfSquaredChange += static_cast<double>(itk::Math::sqr(update[j]));
    }
  }

  return update;
}

/*
 * Compute update at a specify neighbourhood
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
typename SyNDemonsRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::PixelType
SyNDemonsRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::ComputeUpdateInv(
  const NeighborhoodType & it,
  void *                   gd,
  const FloatOffsetType &  itkNotUsed(offset))
{
  PixelType    update;
  unsigned int j;

  IndexType index = it.GetIndex();

  // Get fixed image related information
  double              fixedValue;
  CovariantVectorType gradient;
  double              gradientSquaredMagnitude = 0;

  // Note: no need to check the index is within
  // fixed image buffer. This is done by the external filter.
  fixedValue = (double)this->GetFixedImage()->GetPixel(index);
  double movingValue = (double)this->GetMovingImage()->GetPixel(index);

  //  if (fixedValue > 0)std::cout << " fxv  " << fixedValue << " movingValue " << movingValue << std::endl;

  //    gradient = m_FixedImageGradientCalculator->EvaluateAtIndex( index );

  gradient = m_MovingImageGradientCalculator->EvaluateAtIndex(index);
  for (j = 0; j < ImageDimension; j++)
  {
    gradientSquaredMagnitude += itk::Math::sqr(gradient[j]);
  }

  /**
   * Compute Update.
   * In the original equation the denominator is defined as (g-f)^2 + grad_mag^2.
   * However there is a mismatch in units between the two terms.
   * The units for the second term is intensity^2/mm^2 while the
   * units for the first term is intensity^2. This mismatch is particularly
   * problematic when the fixed image does not have unit spacing.
   * In this implemenation, we normalize the first term by a factor K,
   * such that denominator = (g-f)^2/K + grad_mag^2
   * where K = mean square spacing to compensate for the mismatch in units.
   */
  double speedValue = movingValue - fixedValue;
  if (std::fabs(speedValue) < static_cast<double>(this->m_RobustnessParameter))
  {
    speedValue = 0;
  }

  // update the metric
  GlobalDataStruct * globalData = reinterpret_cast<GlobalDataStruct *>(gd);
  if (globalData)
  {
    globalData->m_SumOfSquaredDifference += itk::Math::sqr(speedValue);
    globalData->m_NumberOfPixelsProcessed += 1;
  }

  double denominator = itk::Math::sqr(speedValue) / m_Normalizer + gradientSquaredMagnitude;

  if (itk::Math::abs(speedValue) < m_IntensityDifferenceThreshold || denominator < m_DenominatorThreshold)
  {
    for (j = 0; j < ImageDimension; j++)
    {
      update[j] = 0.0;
    }
    return update;
  }
  for (j = 0; j < ImageDimension; j++)
  {
    update[j] = speedValue * gradient[j] / denominator;
    if (globalData)
    {
      globalData->m_SumOfSquaredChange += static_cast<double>(itk::Math::sqr(update[j]));
    }
  }

  return update;
}

/*
 * Update the metric and release the per-thread-global data.
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
void
SyNDemonsRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::ReleaseGlobalDataPointer(void * gd) const
{
  GlobalDataStruct * globalData = reinterpret_cast<GlobalDataStruct *>(gd);

  m_SumOfSquaredDifference += globalData->m_SumOfSquaredDifference;
  m_NumberOfPixelsProcessed += globalData->m_NumberOfPixelsProcessed;
  m_SumOfSquaredChange += globalData->m_SumOfSquaredChange;
  if (m_NumberOfPixelsProcessed)
  {
    m_Metric = m_SumOfSquaredDifference / static_cast<double>(m_NumberOfPixelsProcessed);
    m_RMSChange = std::sqrt(m_SumOfSquaredChange / static_cast<double>(m_NumberOfPixelsProcessed));
  }

  delete globalData;
}
} // end namespace itk

#endif
