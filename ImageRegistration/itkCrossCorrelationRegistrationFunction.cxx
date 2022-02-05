/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkCrossCorrelationRegistrationFunction_hxx_
#define _itkCrossCorrelationRegistrationFunction_hxx_

#include "itkCrossCorrelationRegistrationFunction.h"
#include "itkMacro.h"
#include "itkMath.h"
#include "itkImageFileWriter.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkMeanImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkImageFileWriter.h"

#include <deque>

namespace itk
{
/*
 * Default constructor
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
CrossCorrelationRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::
  CrossCorrelationRegistrationFunction()
{
  m_AvgMag = 0;
  m_Iteration = 0;
  RadiusType   r;
  unsigned int j;
  for (j = 0; j < ImageDimension; j++)
  {
    r[j] = 2;
  }
  this->SetRadius(r);
  this->m_Energy = 0.0;
  m_TimeStep = 1.0;
  m_DenominatorThreshold = 1e-9;
  m_IntensityDifferenceThreshold = 0.001;
  Superclass::m_MovingImage = nullptr;
  m_MetricGradientImage = nullptr;
  Superclass::m_FixedImage = nullptr;
  m_FixedImageSpacing.Fill(1.0);
  m_FixedImageOrigin.Fill(0.0);
  m_FixedImageGradientCalculator = GradientCalculatorType::New();
  binaryimage = nullptr;
  m_FullyRobust = false;
  m_MovingImageGradientCalculator = GradientCalculatorType::New();

  typename DefaultInterpolatorType::Pointer interp = DefaultInterpolatorType::New();

  m_MovingImageInterpolator = static_cast<InterpolatorType *>(interp.GetPointer());
  for (auto & finitediffimage : finitediffimages)
  {
    finitediffimage = nullptr;
  }

  m_NumberOfHistogramBins = 32;

  m_FixedImageMask = nullptr;
  m_MovingImageMask = nullptr;
}

/*
 * Standard "PrintSelf" method.
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
void
CrossCorrelationRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::PrintSelf(std::ostream & os,
                                                                                               Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  /*
    os << indent << "MovingImageIterpolator: ";
    os << m_MovingImageInterpolator.GetPointer() << std::endl;
    os << indent << "FixedImageGradientCalculator: ";
    os << m_FixedImageGradientCalculator.GetPointer() << std::endl;
    os << indent << "DenominatorThreshold: ";
    os << m_DenominatorThreshold << std::endl;
    os << indent << "IntensityDifferenceThreshold: ";
    os << m_IntensityDifferenceThreshold << std::endl;
  */
}

/*
 * Set the function state values before each iteration
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
void
CrossCorrelationRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::InitializeIteration()
{
  if (!Superclass::m_MovingImage || !Superclass::m_FixedImage || !m_MovingImageInterpolator)
  {
    itkExceptionMacro(<< "MovingImage, FixedImage and/or Interpolator not set");
    throw ExceptionObject(__FILE__, __LINE__);
  }

  // cache fixed image information
  m_FixedImageSpacing = Superclass::m_FixedImage->GetSpacing();
  m_FixedImageOrigin = Superclass::m_FixedImage->GetOrigin();

  // setup gradient calculator
  m_FixedImageGradientCalculator->SetInputImage(Superclass::m_FixedImage);
  m_MovingImageGradientCalculator->SetInputImage(Superclass::m_MovingImage);

  // setup moving image interpolator
  m_MovingImageInterpolator->SetInputImage(Superclass::m_MovingImage);

  m_MetricTotal = 0.0;
  this->m_Energy = 0.0;

  // compute the normalizer
  m_Normalizer = 0.0;
  for (unsigned int k = 0; k < ImageDimension; k++)
  {
    m_Normalizer += static_cast<float>(itk::Math::sqr(m_FixedImageSpacing[k]));
  }
  m_Normalizer /= static_cast<float>(ImageDimension);

  bool makeimg = false;
  if (m_Iteration == 0)
  {
    makeimg = true;
  }
  else if (!finitediffimages[0])
  {
    makeimg = true;
  }
  else
  {
    for (unsigned int dd = 0; dd < ImageDimension; dd++)
    {
      if (finitediffimages[0]->GetLargestPossibleRegion().GetSize()[dd] !=
          this->GetFixedImage()->GetLargestPossibleRegion().GetSize()[dd])
      {
        makeimg = true;
      }
    }
  }

  if (makeimg)
  {
    finitediffimages[0] = this->MakeImage();
    finitediffimages[1] = this->MakeImage();
    finitediffimages[2] = this->MakeImage();
    finitediffimages[3] = this->MakeImage();
    finitediffimages[4] = this->MakeImage();
  }

  // float sig=15.;

  RadiusType r;
  for (unsigned int j = 0; j < ImageDimension; j++)
  {
    r[j] = this->GetRadius()[j];
  }

  typedef itk::ImageRegionIteratorWithIndex<TFixedImage> Iterator;
  Iterator tIter(this->GetFixedImage(), this->GetFixedImage()->GetLargestPossibleRegion());

  // compute local means
  //  typedef itk::ImageRegionIteratorWithIndex<MetricImageType> Iterator;

  //
  // The following change was made to speed up the correlation calculation.
  //

  typedef std::deque<float> SumQueueType;

  SumQueueType Qsuma2;
  SumQueueType Qsumb2;
  SumQueueType Qsuma;
  SumQueueType Qsumb;
  SumQueueType Qsumab;
  SumQueueType Qcount;

  ImageLinearConstIteratorWithIndex<MetricImageType> outIter(this->finitediffimages[0],
                                                             this->finitediffimages[0]->GetLargestPossibleRegion());
  outIter.SetDirection(0);
  outIter.GoToBegin();
  while (!outIter.IsAtEnd())
  {
    // Push the zeros onto the stack that are outsized the image boundary at
    // the beginning of the line.
    Qsuma2 = SumQueueType(r[0], 0.0);
    Qsumb2 = SumQueueType(r[0], 0.0);
    Qsuma = SumQueueType(r[0], 0.0);
    Qsumb = SumQueueType(r[0], 0.0);
    Qsumab = SumQueueType(r[0], 0.0);
    Qcount = SumQueueType(r[0], 0.0);

    NeighborhoodIterator<MetricImageType> hoodIt(
      this->GetRadius(), this->finitediffimages[0], this->finitediffimages[0]->GetLargestPossibleRegion());
    IndexType oindex = outIter.GetIndex();
    hoodIt.SetLocation(oindex);
    unsigned int hoodlen = hoodIt.Size();
    // Now add the rest of the values from each hyperplane
    for (unsigned int i = r[0]; i < (2 * r[0] + 1); i++)
    {
      float suma2 = 0.0;
      float sumb2 = 0.0;
      float suma = 0.0;
      float sumb = 0.0;
      float sumab = 0.0;
      float count = 0.0;
      for (unsigned int indct = i; indct < hoodlen; indct += (2 * r[0] + 1))
      {
        bool isInBounds = true;
        hoodIt.GetPixel(indct, isInBounds);
        IndexType index = hoodIt.GetIndex(indct);

        if (!isInBounds || (this->m_FixedImageMask && this->m_FixedImageMask->GetPixel(index) <
                                                        static_cast<typename FixedImageType::PixelType>(0.25)))
        {
          continue;
        }

        float a = this->GetFixedImage()->GetPixel(index);
        float b = this->GetMovingImage()->GetPixel(index);

        suma2 += a * a;
        sumb2 += b * b;
        suma += a;
        sumb += b;
        sumab += a * b;
        count += itk::NumericTraits<float>::OneValue();
      }

      Qsuma2.push_back(suma2);
      Qsumb2.push_back(sumb2);
      Qsuma.push_back(suma);
      Qsumb.push_back(sumb);
      Qsumab.push_back(sumab);
      Qcount.push_back(count);
    }

    while (!outIter.IsAtEndOfLine())
    {
      // Test to see if there are any voxels we need to handle in the current
      // window.

      float suma2 = 0.0;
      float sumb2 = 0.0;
      float suma = 0.0;
      float sumb = 0.0;
      float sumab = 0.0;
      float count = 0.0;

      typename SumQueueType::iterator itcount = Qcount.begin();
      while (itcount != Qcount.end())
      {
        count += *itcount;
        ++itcount;
      }

      // If there are values, we need to calculate the different quantities
      if (count > 0)
      {
        typename SumQueueType::iterator ita2 = Qsuma2.begin();
        typename SumQueueType::iterator itb2 = Qsumb2.begin();
        typename SumQueueType::iterator ita = Qsuma.begin();
        typename SumQueueType::iterator itb = Qsumb.begin();
        typename SumQueueType::iterator itab = Qsumab.begin();

        while (ita2 != Qsuma2.end())
        {
          suma2 += *ita2;
          sumb2 += *itb2;
          suma += *ita;
          sumb += *itb;
          sumab += *itab;

          ++ita2;
          ++itb2;
          ++ita;
          ++itb;
          ++itab;
        }

        float fixedMean = suma / count;
        float movingMean = sumb / count;

        float sff = suma2 - fixedMean * suma - fixedMean * suma + count * fixedMean * fixedMean;
        float smm = sumb2 - movingMean * sumb - movingMean * sumb + count * movingMean * movingMean;
        float sfm = sumab - movingMean * suma - fixedMean * sumb + count * movingMean * fixedMean;

        IndexType _oindex = outIter.GetIndex();

        float val = this->GetFixedImage()->GetPixel(_oindex) - fixedMean;
        this->finitediffimages[0]->SetPixel(_oindex, val);
        val = this->GetMovingImage()->GetPixel(_oindex) - movingMean;
        this->finitediffimages[1]->SetPixel(_oindex, val);
        this->finitediffimages[2]->SetPixel(_oindex, sfm); // A
        this->finitediffimages[3]->SetPixel(_oindex, sff); // B
        this->finitediffimages[4]->SetPixel(_oindex, smm); // C
      }

      // Increment the iterator and check to see if we're at the end of the
      // line.  If so, go to the next line.  Otherwise, add the
      // the values for the next hyperplane.
      ++outIter;

      if (!outIter.IsAtEndOfLine())
      {
        hoodIt.SetLocation(outIter.GetIndex());

        suma2 = 0.0;
        sumb2 = 0.0;
        suma = 0.0;
        sumb = 0.0;
        sumab = 0.0;
        count = 0.0;
        for (unsigned int indct = 2 * r[0]; indct < hoodlen; indct += (2 * r[0] + 1))
        {
          bool isInBounds = true;
          hoodIt.GetPixel(indct, isInBounds);
          IndexType index = hoodIt.GetIndex(indct);

          if (!isInBounds || (this->m_FixedImageMask && this->m_FixedImageMask->GetPixel(index) <
                                                          static_cast<typename FixedImageType::PixelType>(0.25)))
          {
            continue;
          }

          float a = this->GetFixedImage()->GetPixel(index);
          float b = this->GetMovingImage()->GetPixel(index);

          suma2 += a * a;
          sumb2 += b * b;
          suma += a;
          sumb += b;
          sumab += a * b;
          count += itk::NumericTraits<float>::OneValue();
        }

        Qsuma2.push_back(suma2);
        Qsumb2.push_back(sumb2);
        Qsuma.push_back(suma);
        Qsumb.push_back(sumb);
        Qsumab.push_back(sumab);
        Qcount.push_back(count);
      }

      Qsuma2.pop_front();
      Qsumb2.pop_front();
      Qsuma.pop_front();
      Qsumb.pop_front();
      Qsumab.pop_front();
      Qcount.pop_front();
    }

    outIter.NextLine();
  }

  // m_FixedImageGradientCalculator->SetInputImage(finitediffimages[0]);

  m_MaxMag = 0.0;
  m_MinMag = 9.e9;
  m_AvgMag = 0.0;
  m_Iteration++;
}

/*
 * Compute the ncc metric everywhere
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
typename TDisplacementField::PixelType
CrossCorrelationRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::ComputeMetricAtPairB(
  IndexType oindex,
  typename TDisplacementField::PixelType /* vec */)
{
  typename TDisplacementField::PixelType deriv;
  deriv.Fill(0.0);
  double sff = 0.0;
  double smm = 0.0;
  double sfm = 0.0;
  //  double fixedValue;
  //  double movingValue;
  CovariantVectorType gradI;
  if (this->m_FixedImageMask)
  {
    if (this->m_FixedImageMask->GetPixel(oindex) < static_cast<typename FixedImageType::PixelType>(0.25))
    {
      return deriv;
    }
  }

  sfm = finitediffimages[2]->GetPixel(oindex);
  sff = finitediffimages[3]->GetPixel(oindex);
  smm = finitediffimages[4]->GetPixel(oindex);
  if (itk::Math::FloatAlmostEqual(sff, itk::NumericTraits<double>::ZeroValue()) ||
      itk::Math::FloatAlmostEqual(smm, itk::NumericTraits<double>::ZeroValue()))
  {
    return deriv;
  }

  this->localCrossCorrelation = 0;
  if (sff * smm > 1.e-5)
  {
    this->localCrossCorrelation = sfm * sfm / (sff * smm);
  }
  IndexType index = oindex; // hoodIt.GetIndex(indct);
  gradI = m_FixedImageGradientCalculator->EvaluateAtIndex(index);
  //    gradJ = m_MovingImageGradientCalculator->EvaluateAtIndex( index );

  float Ji = finitediffimages[1]->GetPixel(index);
  float Ii = finitediffimages[0]->GetPixel(index);

  m_TEMP = 2.0 * sfm / (sff * smm) * (static_cast<double>(Ji) - sfm / sff * static_cast<double>(Ii));
  for (unsigned int qq = 0; qq < ImageDimension; qq++)
  {
    deriv[qq] -= static_cast<float>(2.0 * sfm / (sff * smm) *
                                    (static_cast<double>(Ji) - sfm / sff * static_cast<double>(Ii)) * gradI[qq]);
    //        derivinv[qq]-=2.0*sfm/(sff*smm)*( Ii - sfm/smm*Ji )*gradJ[qq];
  }

  //  if ( this->localCrossCorrelation*(-1.0) < this->m_RobustnessParameter) deriv.Fill(0);
  //  if ( this->localCrossCorrelation*(-1.0) < this->m_RobustnessParameter) {
  //  std::cout << " localC " << this->localCrossCorrelation << std::endl; }
  if (this->localCrossCorrelation < 1)
  {
    this->m_Energy -= this->localCrossCorrelation;
  }
  return deriv; // localCrossCorrelation;
}

/*
 * Compute the ncc metric everywhere
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
typename TDisplacementField::PixelType
CrossCorrelationRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::ComputeMetricAtPairC(
  IndexType oindex,
  typename TDisplacementField::PixelType /* vec */)
{
  typename TDisplacementField::PixelType deriv;
  deriv.Fill(0.0);
  double              sff = 0.0;
  double              smm = 0.0;
  double              sfm = 0.0;
  CovariantVectorType gradJ;
  if (this->m_FixedImageMask)
  {
    if (this->m_FixedImageMask->GetPixel(oindex) < static_cast<typename TFixedImage::PixelType>(0.25))
    {
      return deriv;
    }
  }

  sfm = finitediffimages[2]->GetPixel(oindex);
  sff = finitediffimages[3]->GetPixel(oindex);
  smm = finitediffimages[4]->GetPixel(oindex);

  if (itk::Math::FloatAlmostEqual(sff, itk::NumericTraits<double>::ZeroValue()) ||
      itk::Math::FloatAlmostEqual(smm, itk::NumericTraits<double>::ZeroValue()))
  {
    return deriv;
  }

  IndexType index = oindex; // hoodIt.GetIndex(indct);
  if (itk::Math::FloatAlmostEqual(sff, itk::NumericTraits<double>::ZeroValue()))
  {
    sff = itk::NumericTraits<double>::OneValue();
  }
  if (itk::Math::FloatAlmostEqual(smm, itk::NumericTraits<double>::ZeroValue()))
  {
    smm = itk::NumericTraits<double>::OneValue();
  }

  // /gradI = m_FixedImageGradientCalculator->EvaluateAtIndex( index );
  gradJ = m_MovingImageGradientCalculator->EvaluateAtIndex(index);

  float Ji = finitediffimages[1]->GetPixel(index);
  float Ii = finitediffimages[0]->GetPixel(index);
  for (unsigned int qq = 0; qq < ImageDimension; qq++)
  {
    deriv[qq] -= static_cast<float>(2.0 * sfm / (sff * smm) *
                                    (static_cast<double>(Ii) - sfm / smm * static_cast<double>(Ji)) * gradJ[qq]);
    // deriv[qq]   -=2.0*sfm/(sff*smm)*( Ji - sfm/sff*Ii )*gradI[qq];
  }

  if (!itk::Math::FloatAlmostEqual(sff * smm, itk::NumericTraits<double>::ZeroValue()))
  {
    this->localCrossCorrelation = sfm * sfm / (sff * smm);
  }
  else
  {
    this->localCrossCorrelation = 1.0;
  }
  //  if ( this->localCrossCorrelation*(-1.0) < this->m_RobustnessParameter) deriv.Fill(0);

  return deriv; // localCrossCorrelation;
}
} // end namespace itk

#endif
