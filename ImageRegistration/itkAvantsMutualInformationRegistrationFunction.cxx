/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkAvantsMutualInformationRegistrationFunction_hxx
#define _itkAvantsMutualInformationRegistrationFunction_hxx

#include "antsAllocImage.h"
#include "itkAvantsMutualInformationRegistrationFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkCovariantVector.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageIterator.h"
#include "itkMath.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::
  AvantsMutualInformationRegistrationFunction()
{
  this->Superclass::m_NormalizeGradient = true;
  this->m_NumberOfSpatialSamples = 5000;
  this->m_NumberOfHistogramBins = 50;

  //  this->SetComputeGradient(false); // don't use the default gradient for no

  // Initialize PDFs to NULL
  m_JointPDF = nullptr;

  m_OpticalFlow = false;
  typename TransformType::Pointer transformer = TransformType::New();
  this->SetTransform(transformer);

  typename BSplineInterpolatorType::Pointer interpolator = BSplineInterpolatorType::New();
  this->SetInterpolator(interpolator);

  m_FixedImageMask = nullptr;
  m_MovingImageMask = nullptr;

  // Initialize memory
  m_MovingImageNormalizedMin = 0.0;
  m_FixedImageNormalizedMin = 0.0;
  m_MovingImageTrueMin = 0.0;
  m_MovingImageTrueMax = 0.0;
  m_FixedImageBinSize = 0.0;
  m_MovingImageBinSize = 0.0;
  m_CubicBSplineDerivativeKernel = nullptr;
  m_BSplineInterpolator = nullptr;
  m_DerivativeCalculator = nullptr;
  m_NumberOfParameters = ImageDimension;

  m_FixedImageGradientCalculator = GradientCalculatorType::New();
  m_MovingImageGradientCalculator = GradientCalculatorType::New();
  this->m_Padding = 2;

  typename DefaultInterpolatorType::Pointer interp = DefaultInterpolatorType::New();
  typename DefaultInterpolatorType::Pointer interp2 = DefaultInterpolatorType::New();

  m_MovingImageInterpolator = static_cast<InterpolatorType *>(interp.GetPointer());
  m_FixedImageInterpolator = static_cast<InterpolatorType *>(interp2.GetPointer());
  m_Interpolator = static_cast<InterpolatorType *>(interp.GetPointer());

  this->m_RobustnessParameter = -1.e19;
}

/**
 * Print out internal information about this class
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::PrintSelf(
  std::ostream & os,
  Indent         indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "NumberOfSpatialSamples: ";
  os << m_NumberOfSpatialSamples << std::endl;
  os << indent << "NumberOfHistogramBins: ";
  os << m_NumberOfHistogramBins << std::endl;

  // Debugging information
  os << indent << "NumberOfParameters: ";
  os << m_NumberOfParameters << std::endl;
  os << indent << "FixedImageNormalizedMin: ";
  os << m_FixedImageNormalizedMin << std::endl;
  os << indent << "MovingImageNormalizedMin: ";
  os << m_MovingImageNormalizedMin << std::endl;
  os << indent << "MovingImageTrueMin: ";
  os << m_MovingImageTrueMin << std::endl;
  os << indent << "MovingImageTrueMax: ";
  os << m_MovingImageTrueMax << std::endl;
  os << indent << "FixedImageBinSize: ";
  os << m_FixedImageBinSize << std::endl;
  os << indent << "MovingImageBinSize: ";
  os << m_MovingImageBinSize << std::endl;
  os << indent << "InterpolatorIsBSpline: ";
  os << m_InterpolatorIsBSpline << std::endl;
}

/**
 * Initialize
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::InitializeIteration()
{
  m_CubicBSplineKernel = CubicBSplineFunctionType::New();
  m_CubicBSplineDerivativeKernel = CubicBSplineDerivativeFunctionType::New();
  this->m_Energy = 0;
  pdfinterpolator = pdfintType::New();
  dpdfinterpolator = dpdfintType::New();
  pdfinterpolator2 = pdfintType2::New();
  pdfinterpolator3 = pdfintType2::New();
  m_DerivativeCalculator = DerivativeFunctionType::New();

  //  this->ComputeMetricImage();

  /*
  bool makenewimage=false;
  if (!this->m_MetricImage ) makenewimage=true;
  else if (imagesize[0] != this->m_MetricImage->GetLargestPossibleRegion().GetSize()[0])
    makenewimage = true;
  else this->m_MetricImage->FillBuffer(0);
  if (makenewimage)
  {
    this->m_MetricImage = TFixedImage::New();
    this->m_MetricImage->SetLargestPossibleRegion(img->GetLargestPossibleRegion()  );
    this->m_MetricImage->SetBufferedRegion(img->GetLargestPossibleRegion());
    this->m_MetricImage->SetSpacing(img->GetSpacing());
    this->m_MetricImage->SetOrigin(img->GetOrigin());
    this->m_MetricImage->AllocateInitialized();
  }
  */
  m_FixedImageGradientCalculator->SetInputImage(this->m_FixedImage);
  m_MovingImageGradientCalculator->SetInputImage(this->m_MovingImage);
  m_FixedImageInterpolator->SetInputImage(this->m_FixedImage);
  m_Interpolator->SetInputImage(this->m_MovingImage);

  m_FixedImageSpacing = this->m_FixedImage->GetSpacing();
  m_FixedImageOrigin = this->m_FixedImage->GetOrigin();
  m_Normalizer = 0.0;
  m_NumberOfSpatialSamples = 1;
  for (unsigned int k = 0; k < ImageDimension; k++)
  {
    m_Normalizer += static_cast<float>(itk::Math::sqr(m_FixedImageSpacing[k]));
    m_NumberOfSpatialSamples *= this->m_FixedImage->GetLargestPossibleRegion().GetSize()[k];
  }
  m_Normalizer /= static_cast<float>(ImageDimension);

  /**
   * Compute binsize for the histograms.
   *
   * The binsize for the image intensities needs to be adjusted so that
   * we can avoid dealing with boundary conditions using the cubic
   * spline as the Parzen window.  We do this by increasing the size
   * of the bins so that the joint histogram becomes "padded" at the
   * borders. Because we are changing the binsize,
   * we also need to shift the minimum by the padded amount in order to
   * avoid minimum values filling in our padded region.
   *
   * Note that there can still be non-zero bin values in the padded region,
   * it's just that these bins will never be a central bin for the Parzen
   * window.
   *
  double fixedImageMax = 1.0;
  double fixedImageMin = 0.0;
  double movingImageMax = 1.0;
  double movingImageMin = 0.0;
  m_MovingImageTrueMin = movingImageMin;
  m_MovingImageTrueMax = movingImageMax;

   */

  double movingImageMin = NumericTraits<double>::max();
  double movingImageMax = NumericTraits<double>::NonpositiveMin();
  double fixedImageMin = NumericTraits<double>::max();
  double fixedImageMax = NumericTraits<double>::NonpositiveMin();

  typedef ImageRegionConstIterator<MovingImageType> MovingIteratorType;
  MovingIteratorType movingImageIterator(this->m_MovingImage, this->m_MovingImage->GetBufferedRegion());
  for (movingImageIterator.GoToBegin(); !movingImageIterator.IsAtEnd(); ++movingImageIterator)
  {
    bool takesample = true;

    if (this->m_FixedImageMask)
    {
      if (this->m_FixedImageMask->GetPixel(movingImageIterator.GetIndex()) <
          static_cast<typename FixedImageType::PixelType>(1.e-6))
      {
        takesample = false;
      }
    }

    if (takesample)
    {
      double sample = static_cast<double>(movingImageIterator.Get());
      double fsample = static_cast<double>(this->m_FixedImage->GetPixel(movingImageIterator.GetIndex()));

      if (sample < movingImageMin)
      {
        movingImageMin = sample;
      }

      if (sample > movingImageMax)
      {
        movingImageMax = sample;
      }

      if (fsample < fixedImageMin)
      {
        fixedImageMin = fsample;
      }

      if (fsample > fixedImageMax)
      {
        fixedImageMax = fsample;
      }
    }
  }
  this->m_MovingImageTrueMax = movingImageMax;
  this->m_FixedImageTrueMax = fixedImageMax;
  this->m_MovingImageTrueMin = movingImageMin;
  this->m_FixedImageTrueMin = fixedImageMin;

  // Instantiate a region, index, size
  JointPDFRegionType jointPDFRegion;
  JointPDFIndexType  jointPDFIndex;
  JointPDFSizeType   jointPDFSize;

  // the jointPDF is of size NumberOfBins x NumberOfBins
  jointPDFSize.Fill(m_NumberOfHistogramBins);
  jointPDFIndex.Fill(0);
  jointPDFRegion.SetIndex(jointPDFIndex);
  jointPDFRegion.SetSize(jointPDFSize);
  this->m_JointPDF = AllocImage<JointPDFType>(jointPDFRegion);

  // By setting these values, the joint histogram physical locations will correspond to intensity values.
  JointPDFSpacingType spacing;
  spacing[0] = 1 / (double)(this->m_NumberOfHistogramBins - (double)this->m_Padding * 2 - 1);
  spacing[1] = spacing[0];
  this->m_JointPDF->SetSpacing(spacing);
  this->m_JointPDFSpacing = this->m_JointPDF->GetSpacing();
  JointPDFPointType origin;
  origin[0] = this->m_JointPDFSpacing[0] * (double)this->m_Padding * (-1.0);
  origin[1] = origin[0];
  this->m_JointPDF->SetOrigin(origin);

  // Instantiate a region, index, size
  typedef typename MarginalPDFType::RegionType MarginalPDFRegionType;
  typedef typename MarginalPDFType::SizeType   MarginalPDFSizeType;
  MarginalPDFRegionType                        marginalPDFRegion;
  MarginalPDFIndexType                         marginalPDFIndex;
  MarginalPDFSizeType                          marginalPDFSize;

  // the marginalPDF is of size NumberOfBins x NumberOfBins
  marginalPDFSize.Fill(m_NumberOfHistogramBins);
  marginalPDFIndex.Fill(0);
  marginalPDFRegion.SetIndex(marginalPDFIndex);
  marginalPDFRegion.SetSize(marginalPDFSize);
  // do the same thing for the marginal pdfs
  this->m_FixedImageMarginalPDF = AllocImage<MarginalPDFType>(marginalPDFRegion);
  this->m_MovingImageMarginalPDF = AllocImage<MarginalPDFType>(marginalPDFRegion);

  // By setting these values, the marginal histogram physical locations will correspond to intensity values.
  typename MarginalPDFType::PointType fixedorigin;
  typename MarginalPDFType::PointType movingorigin;
  fixedorigin[0] = origin[0];
  movingorigin[0] = origin[1];
  this->m_FixedImageMarginalPDF->SetOrigin(fixedorigin);
  this->m_MovingImageMarginalPDF->SetOrigin(movingorigin);
  typename MarginalPDFType::SpacingType mspacing;
  mspacing[0] = spacing[0];
  this->m_FixedImageMarginalPDF->SetSpacing(mspacing);
  mspacing[0] = spacing[1];
  this->m_MovingImageMarginalPDF->SetSpacing(mspacing);

  // For the derivatives of the joint PDF define a region starting from {0,0,0}
  // with size {m_NumberOfParameters,m_NumberOfHistogramBins,
  // m_NumberOfHistogramBins}. The dimension represents transform parameters,
  // fixed image parzen window index and moving image parzen window index,
  // respectively.

  m_NormalizeMetric = 1.0;
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    m_NormalizeMetric *= this->m_FixedImage->GetLargestPossibleRegion().GetSize()[i];
  }

  this->GetProbabilities();
  this->ComputeMutualInformation();

  pdfinterpolator->SetInputImage(m_JointPDF);
  pdfinterpolator2->SetInputImage(m_FixedImageMarginalPDF);
  pdfinterpolator3->SetInputImage(m_MovingImageMarginalPDF);
  /*  pdfinterpolator->SetSplineOrder(3);
  pdfinterpolator2->SetSplineOrder(3);
  pdfinterpolator3->SetSplineOrder(3);
  */
}

/**
 * Get the both Value and Derivative Measure
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::GetProbabilities()
{
  typedef ImageRegionConstIteratorWithIndex<FixedImageType> RandomIterator;
  RandomIterator randIter(this->m_FixedImage, this->m_FixedImage->GetLargestPossibleRegion());
  this->m_FixedImageMarginalPDF->FillBuffer(0);
  this->m_MovingImageMarginalPDF->FillBuffer(0);

  // Reset the joint pdfs to zero
  m_JointPDF->FillBuffer(0.0);

  // unsigned long  nSamples = 0;
  RandomIterator iter(this->m_FixedImage, this->m_FixedImage->GetLargestPossibleRegion());
  for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
  {
    bool takesample = true;
    if (this->m_FixedImageMask)
    {
      if (this->m_FixedImageMask->GetPixel(iter.GetIndex()) < static_cast<typename FixedImageType::PixelType>(1.e-6))
      {
        takesample = false;
      }
    }

    if (takesample)
    {
      // Get sampled index
      FixedImageIndexType index = iter.GetIndex();
      /*    bool inimage=true;
          for (unsigned int dd=0; dd<ImageDimension; dd++)
            {
              if ( index[dd] < 1 ||    index[dd] >= static_cast<typename IndexType::IndexValueType>(imagesize[dd]-1) )
                inimage=false;
            }
      */
      double movingImageValue = this->GetMovingParzenTerm(this->m_MovingImage->GetPixel(index));
      double fixedImageValue = this->GetFixedParzenTerm(this->m_FixedImage->GetPixel(index));

      /** add the paired intensity points to the joint histogram */
      JointPDFPointType jointPDFpoint;
      this->ComputeJointPDFPoint(fixedImageValue, movingImageValue, jointPDFpoint);
      JointPDFIndexType jointPDFIndex;
      jointPDFIndex.Fill(0);
      jointPDFIndex = this->m_JointPDF->TransformPhysicalPointToIndex(jointPDFpoint);
      this->m_JointPDF->SetPixel(jointPDFIndex, this->m_JointPDF->GetPixel(jointPDFIndex) + 1);

      // ++nSamples;
    }
  }

  /**
   * Normalize the PDFs, compute moving image marginal PDF
   *
   */
  typedef ImageRegionIterator<JointPDFType> JointPDFIteratorType;
  JointPDFIteratorType                      jointPDFIterator(m_JointPDF, m_JointPDF->GetBufferedRegion());

  // Compute joint PDF normalization factor (to ensure joint PDF sum adds to 1.0)
  double jointPDFSum = 0.0;
  jointPDFIterator.GoToBegin();
  while (!jointPDFIterator.IsAtEnd())
  {
    float temp = jointPDFIterator.Get();
    jointPDFSum += static_cast<double>(temp);
    ++jointPDFIterator;
  }

  // of derivatives
  if (itk::Math::FloatAlmostEqual(jointPDFSum, itk::NumericTraits<double>::ZeroValue()))
  {
    itkExceptionMacro("Joint PDF summed to zero");
  }

  // Normalize the PDF bins
  jointPDFIterator.GoToEnd();
  while (!jointPDFIterator.IsAtBegin())
  {
    --jointPDFIterator;
    jointPDFIterator.Value() /= static_cast<PDFValueType>(jointPDFSum);
  }

  constexpr bool smoothjh = true;
  if (smoothjh)
  {
    typedef DiscreteGaussianImageFilter<JointPDFType, JointPDFType> dgtype;
    typename dgtype::Pointer                                        dg = dgtype::New();
    dg->SetInput(this->m_JointPDF);
    dg->SetVariance(1.5);
    dg->SetUseImageSpacing(false);
    dg->SetMaximumError(.01f);
    dg->Update();
    this->m_JointPDF = dg->GetOutput();
  }

  // Compute moving image marginal PDF by summing over fixed image bins.
  typedef ImageLinearIteratorWithIndex<JointPDFType> JointPDFLinearIterator;
  JointPDFLinearIterator                             linearIter(m_JointPDF, m_JointPDF->GetBufferedRegion());
  linearIter.SetDirection(0);
  linearIter.GoToBegin();
  unsigned int fixedIndex = 0;
  while (!linearIter.IsAtEnd())
  {
    double sum = 0.0;
    while (!linearIter.IsAtEndOfLine())
    {
      sum += static_cast<double>(linearIter.Get());
      ++linearIter;
    }

    MarginalPDFIndexType mind;
    mind[0] = fixedIndex;
    m_FixedImageMarginalPDF->SetPixel(mind, static_cast<PDFValueType>(sum));
    linearIter.NextLine();
    ++fixedIndex;
  }

  linearIter.SetDirection(1);
  linearIter.GoToBegin();
  unsigned int movingIndex = 0;
  while (!linearIter.IsAtEnd())
  {
    double sum = 0.0;
    while (!linearIter.IsAtEndOfLine())
    {
      sum += static_cast<double>(linearIter.Get());
      ++linearIter;
    }

    MarginalPDFIndexType mind;
    mind[0] = movingIndex;
    m_MovingImageMarginalPDF->SetPixel(mind, static_cast<PDFValueType>(sum));
    linearIter.NextLine();
    ++movingIndex;
  }
}

/**
 * Get the both Value and Derivative Measure
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
double
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::GetValueAndDerivative(
  IndexType oindex,
  MeasureType & /* valuei */,
  DerivativeType & /* derivative1 */,
  DerivativeType & /* derivative2 */)
{
  double         value = 0;
  DerivativeType zero(ImageDimension);

  zero.Fill(0);

  double movingImageValue = this->GetMovingParzenTerm(this->m_MovingImage->GetPixel(oindex));
  double fixedImageValue = this->GetFixedParzenTerm(this->m_FixedImage->GetPixel(oindex));

  JointPDFPointType pdfind;
  this->ComputeJointPDFPoint(fixedImageValue, movingImageValue, pdfind);
  const double jointPDFValue = pdfinterpolator->Evaluate(pdfind);
  const double dJPDF = this->ComputeJointPDFDerivative(pdfind, 0, 0);

  typename pdfintType2::ContinuousIndexType mind;
  mind[0] = pdfind[0];
  const double fixedImagePDFValue = pdfinterpolator2->Evaluate(mind);
  const double dFmPDF = this->ComputeFixedImageMarginalPDFDerivative(mind, 0);

  const double eps = 1.e-16;
  if (jointPDFValue > eps && (fixedImagePDFValue) > 0)
  {
    const double pRatio = std::log(jointPDFValue) - std::log(fixedImagePDFValue);
    const double term1 = dJPDF * pRatio;
    const double term2 = std::log((double)2) * dFmPDF * jointPDFValue / fixedImagePDFValue;
    value = (term2 - term1);
  } // end if-block to check non-zero bin contribution
  else
  {
    value = 0;
  }
  return value;
}

template <typename TFixedImage, typename TMovingImage, typename TDisplacementField>
double
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDisplacementField>::GetValueAndDerivativeInv(
  IndexType oindex,
  MeasureType & /* valuei */,
  DerivativeType & /* derivative1 */,
  DerivativeType & /* derivative2 */)
{
  double         value = 0;
  DerivativeType zero(ImageDimension);

  zero.Fill(0);

  double movingImageValue = this->GetMovingParzenTerm(this->m_MovingImage->GetPixel(oindex));
  double fixedImageValue = this->GetFixedParzenTerm(this->m_FixedImage->GetPixel(oindex));

  JointPDFPointType pdfind;
  this->ComputeJointPDFPoint(fixedImageValue, movingImageValue, pdfind);
  const double jointPDFValue = pdfinterpolator->Evaluate(pdfind);
  const double dJPDF = this->ComputeJointPDFDerivative(pdfind, 0, 1);

  typename pdfintType2::ContinuousIndexType mind;
  mind[0] = pdfind[1];
  const double movingImagePDFValue = pdfinterpolator3->EvaluateAtContinuousIndex(mind);
  const double dMmPDF = this->ComputeMovingImageMarginalPDFDerivative(mind, 0);

  const double eps = 1.e-16;
  if (jointPDFValue > eps && (movingImagePDFValue) > 0)
  {
    const double pRatio = std::log(jointPDFValue) - std::log(movingImagePDFValue);
    const double term1 = dJPDF * pRatio;
    const double term2 = std::log((double)2) * dMmPDF * jointPDFValue / movingImagePDFValue;
    value = (term2 - term1);
  } // end if-block to check non-zero bin contribution
  else
  {
    value = 0;
  }

  return value;
}
} // end namespace itk

#endif
