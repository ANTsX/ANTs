/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkN3MRIBiasFieldCorrectionImageFilter_hxx
#define __itkN3MRIBiasFieldCorrectionImageFilter_hxx


#include "itkDivideImageFilter.h"
#include "itkExpImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkIterationReporter.h"
#include "itkLBFGSBOptimizer.h"
#include "itkLogImageFilter.h"
#include "itkSubtractImageFilter.h"

#include "vnl/algo/vnl_fft_1d.h"
#include "vnl/vnl_complex_traits.h"

namespace itk
{
/**
 * N3BiasFieldScaleCostFunction class definitions
 */
template <typename TInputImage, typename TBiasFieldImage, typename TMaskImage, typename TConfidenceImage>
N3BiasFieldScaleCostFunction<TInputImage, TBiasFieldImage, TMaskImage, TConfidenceImage>::N3BiasFieldScaleCostFunction()
{
  this->m_InputImage = nullptr;
  this->m_BiasFieldImage = nullptr;
  this->m_MaskImage = nullptr;
  this->m_ConfidenceImage = nullptr;
}

template <typename TInputImage, typename TBiasFieldImage, typename TMaskImage, typename TConfidenceImage>
N3BiasFieldScaleCostFunction<TInputImage, TBiasFieldImage, TMaskImage, TConfidenceImage>::
  ~N3BiasFieldScaleCostFunction() = default;

template <typename TInputImage, typename TBiasFieldImage, typename TMaskImage, typename TConfidenceImage>
typename N3BiasFieldScaleCostFunction<TInputImage, TBiasFieldImage, TMaskImage, TConfidenceImage>::MeasureType
N3BiasFieldScaleCostFunction<TInputImage, TBiasFieldImage, TMaskImage, TConfidenceImage>::GetValue(
  const ParametersType & parameters) const
{
  ImageRegionConstIterator<TInputImage>     ItI(this->m_InputImage, this->m_InputImage->GetRequestedRegion());
  ImageRegionConstIterator<TBiasFieldImage> ItB(this->m_BiasFieldImage, this->m_BiasFieldImage->GetRequestedRegion());

  MeasureType mu = 0.0;
  MeasureType N = 0.0;
  for (ItI.GoToBegin(), ItB.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItB)
  {
    if ((!this->m_MaskImage ||
         this->m_MaskImage->GetPixel(ItI.GetIndex()) != NumericTraits<typename TMaskImage::PixelType>::ZeroValue()) &&
        (!this->m_ConfidenceImage || this->m_ConfidenceImage->GetPixel(ItI.GetIndex()) >
                                       NumericTraits<typename TConfidenceImage::PixelType>::ZeroValue()))
    {
      mu += static_cast<MeasureType>(ItI.Get()) /
            (parameters[0] * (static_cast<MeasureType>(ItB.Get()) - NumericTraits<MeasureType>::OneValue()) +
             NumericTraits<MeasureType>::OneValue());
      N += NumericTraits<MeasureType>::OneValue();
    }
  }
  mu /= N;

  MeasureType value = NumericTraits<MeasureType>::ZeroValue();
  for (ItI.GoToBegin(), ItB.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItB)
  {
    if ((!this->m_MaskImage ||
         this->m_MaskImage->GetPixel(ItI.GetIndex()) != NumericTraits<typename TMaskImage::PixelType>::ZeroValue()) &&
        (!this->m_ConfidenceImage || this->m_ConfidenceImage->GetPixel(ItI.GetIndex()) >
                                       NumericTraits<typename TConfidenceImage::PixelType>::ZeroValue()))
    {
      value += Math::sqr((static_cast<MeasureType>(ItI.Get()) /
                          (static_cast<MeasureType>(parameters[0]) *
                             (static_cast<MeasureType>(ItB.Get()) - NumericTraits<MeasureType>::OneValue()) +
                           NumericTraits<MeasureType>::OneValue())) /
                           mu -
                         NumericTraits<MeasureType>::OneValue());
    }
  }
  value /= (N - NumericTraits<MeasureType>::OneValue());

  return value;
}

template <typename TInputImage, typename TBiasFieldImage, typename TMaskImage, typename TConfidenceImage>
void
N3BiasFieldScaleCostFunction<TInputImage, TBiasFieldImage, TMaskImage, TConfidenceImage>::GetDerivative(
  const ParametersType & parameters,
  DerivativeType &       derivative) const
{
  ImageRegionConstIterator<TInputImage>     ItI(this->m_InputImage, this->m_InputImage->GetRequestedRegion());
  ImageRegionConstIterator<TBiasFieldImage> ItB(this->m_BiasFieldImage, this->m_BiasFieldImage->GetRequestedRegion());

  MeasureType mu = 0.0;
  MeasureType dmu = 0.0;
  MeasureType N = 0.0;
  for (ItI.GoToBegin(), ItB.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItB)
  {
    if ((!this->m_MaskImage ||
         this->m_MaskImage->GetPixel(ItI.GetIndex()) != NumericTraits<typename TMaskImage::PixelType>::ZeroValue()) &&
        (!this->m_ConfidenceImage || this->m_ConfidenceImage->GetPixel(ItI.GetIndex()) >
                                       NumericTraits<typename TConfidenceImage::PixelType>::ZeroValue()))
    {
      MeasureType d = parameters[0] * (static_cast<MeasureType>(ItB.Get()) - 1.0) + 1.0;
      mu += (static_cast<MeasureType>(ItI.Get()) / d);
      dmu += -static_cast<MeasureType>(ItI.Get()) * (static_cast<MeasureType>(ItB.Get()) - 1.0) / d;
      N += 1.0;
    }
  }
  mu /= N;
  dmu /= N;

  MeasureType value = 0.0;
  for (ItI.GoToBegin(), ItB.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItB)
  {
    if ((!this->m_MaskImage ||
         this->m_MaskImage->GetPixel(ItI.GetIndex()) != NumericTraits<typename TMaskImage::PixelType>::ZeroValue()) &&
        (!this->m_ConfidenceImage || this->m_ConfidenceImage->GetPixel(ItI.GetIndex()) >
                                       NumericTraits<typename TConfidenceImage::PixelType>::ZeroValue()))
    {
      MeasureType d = static_cast<MeasureType>(parameters[0]) *
                        (static_cast<MeasureType>(ItB.Get()) - NumericTraits<MeasureType>::OneValue()) +
                      NumericTraits<MeasureType>::OneValue();
      MeasureType t = static_cast<MeasureType>(ItI.Get()) / d;
      MeasureType dt = -t * (static_cast<MeasureType>(ItB.Get()) - NumericTraits<MeasureType>::OneValue());
      value += ((t / mu - NumericTraits<MeasureType>::OneValue()) * (dt / mu - dmu * t / (Math::sqr(mu))));
    }
  }
  derivative.SetSize(1);
  derivative(0) = static_cast<MeasureType>(2.0) * value / (N - NumericTraits<MeasureType>::OneValue());
}

template <typename TInputImage, typename TBiasFieldImage, typename TMaskImage, typename TConfidenceImage>
unsigned int
N3BiasFieldScaleCostFunction<TInputImage, TBiasFieldImage, TMaskImage, TConfidenceImage>::GetNumberOfParameters() const
{
  return NumericTraits<unsigned int>::OneValue();
}

/**
 * N3MRIBiasFieldCorrectionImageFilter class definitions
 */

template <typename TInputImage, typename TMaskImage, typename TOutputImage>
N3MRIBiasFieldCorrectionImageFilter<TInputImage, TMaskImage, TOutputImage>::N3MRIBiasFieldCorrectionImageFilter()
{
  this->SetNumberOfRequiredInputs(1);

  this->m_NumberOfHistogramBins = 200;
  this->m_WienerFilterNoise = 0.01;
  this->m_BiasFieldFullWidthAtHalfMaximum = 0.15;

  this->m_MaximumNumberOfIterations = 50;
  this->m_ConvergenceThreshold = 0.001;

  this->m_SplineOrder = 3;
  this->m_NumberOfFittingLevels.Fill(4);
  this->m_NumberOfControlPoints.Fill(4);

  this->m_UseOptimalBiasFieldScaling = true;
  this->m_BiasFieldScaling = 1.0;
}

template <typename TInputImage, typename TMaskImage, typename TOutputImage>
void
N3MRIBiasFieldCorrectionImageFilter<TInputImage, TMaskImage, TOutputImage>::GenerateData()
{
  /**
   * Calculate the log of the input image.
   */
  typename RealImageType::Pointer logInputImage = RealImageType::New();

  typedef ExpImageFilter<RealImageType, RealImageType>  ExpImageFilterType;
  typedef LogImageFilter<InputImageType, RealImageType> LogFilterType;

  typename LogFilterType::Pointer logFilter1 = LogFilterType::New();
  logFilter1->SetInput(this->GetInput());
  logFilter1->Update();
  logInputImage = logFilter1->GetOutput();

  /**
   * Remove possible nans/infs from the log input image.
   */
  ImageRegionIteratorWithIndex<RealImageType> It(logInputImage, logInputImage->GetRequestedRegion());
  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    if ((!this->GetMaskImage() ||
         this->GetMaskImage()->GetPixel(It.GetIndex()) != NumericTraits<MaskPixelType>::ZeroValue()) &&
        (!this->GetConfidenceImage() ||
         this->GetConfidenceImage()->GetPixel(It.GetIndex()) > NumericTraits<RealType>::ZeroValue()))
    {
      if (std::isnan(It.Get()) || std::isinf(It.Get()) || It.Get() < NumericTraits<RealType>::ZeroValue())
      {
        It.Set(NumericTraits<RealType>::ZeroValue());
      }
    }
  }

  /**
   * Provide an initial log bias field of zeros
   */
  typename RealImageType::Pointer logBiasField = AllocImage<RealImageType>(this->GetInput()->GetRequestedRegion(), 0.0);
  logBiasField->SetOrigin(this->GetInput()->GetOrigin());
  logBiasField->SetSpacing(this->GetInput()->GetSpacing());
  logBiasField->SetDirection(this->GetInput()->GetDirection());

  /**
   * Iterate until convergence or iterative exhaustion.
   */
  IterationReporter reporter(this, 0, 1);

  bool isConverged = false;
  this->m_ElapsedIterations = 0;
  while (!isConverged && this->m_ElapsedIterations++ < this->m_MaximumNumberOfIterations)
  {
    typedef SubtractImageFilter<RealImageType, RealImageType, RealImageType> SubtracterType;

    typename SubtracterType::Pointer subtracter1 = SubtracterType::New();
    subtracter1->SetInput1(logInputImage);
    subtracter1->SetInput2(logBiasField);
    subtracter1->Update();

    typename RealImageType::Pointer logSharpenedImage = this->SharpenImage(subtracter1->GetOutput());

    typename SubtracterType::Pointer subtracter2 = SubtracterType::New();
    subtracter2->SetInput1(logInputImage);
    subtracter2->SetInput2(logSharpenedImage);
    subtracter2->Update();

    typename RealImageType::Pointer newLogBiasField = this->SmoothField(subtracter2->GetOutput());

    this->m_CurrentConvergenceMeasurement = this->CalculateConvergenceMeasurement(logBiasField, newLogBiasField);
    isConverged = (this->m_CurrentConvergenceMeasurement < this->m_ConvergenceThreshold);

    logBiasField = newLogBiasField;

    reporter.CompletedStep();
  }

  typedef ExpImageFilter<RealImageType, RealImageType> ExpImageFilterType;
  typename ExpImageFilterType::Pointer                 expFilter = ExpImageFilterType::New();
  expFilter->SetInput(logBiasField);
  expFilter->Update();

  /**
   * Calculate the optimal scaling
   */
  if (this->m_UseOptimalBiasFieldScaling)
  {
    this->m_BiasFieldScaling = this->CalculateOptimalBiasFieldScaling(expFilter->GetOutput());

    ImageRegionIterator<RealImageType> ItE(expFilter->GetOutput(), expFilter->GetOutput()->GetRequestedRegion());
    for (ItE.GoToBegin(); !ItE.IsAtEnd(); ++ItE)
    {
      ItE.Set(this->m_BiasFieldScaling * (ItE.Get() - NumericTraits<RealType>::OneValue()) +
              NumericTraits<RealType>::OneValue());
    }
  }

  /**
   * Divide the input image by the bias field to get the final image.
   */
  typedef DivideImageFilter<InputImageType, RealImageType, OutputImageType> DividerType;
  typename DividerType::Pointer                                             divider = DividerType::New();
  divider->SetInput1(this->GetInput());
  divider->SetInput2(expFilter->GetOutput());
  divider->Update();

  this->SetNthOutput(0, divider->GetOutput());
}

template <typename TInputImage, typename TMaskImage, typename TOutputImage>
typename N3MRIBiasFieldCorrectionImageFilter<TInputImage, TMaskImage, TOutputImage>::RealImageType::Pointer
N3MRIBiasFieldCorrectionImageFilter<TInputImage, TMaskImage, TOutputImage>::SharpenImage(
  typename N3MRIBiasFieldCorrectionImageFilter<TInputImage, TMaskImage, TOutputImage>::RealImageType::Pointer
    unsharpenedImage)
{
  /**
   * Build the histogram for the uncorrected image.  Store copy
   * in a vnl_vector to utilize vnl FFT routines.  Note that variables
   * in real space are denoted by a single uppercase letter whereas their
   * frequency counterparts are indicated by a trailing lowercase 'f'.
   */
  RealType binMaximum = NumericTraits<RealType>::NonpositiveMin();
  RealType binMinimum = NumericTraits<RealType>::max();

  ImageRegionIterator<RealImageType> ItU(unsharpenedImage, unsharpenedImage->GetLargestPossibleRegion());
  for (ItU.GoToBegin(); !ItU.IsAtEnd(); ++ItU)
  {
    if ((!this->GetMaskImage() ||
         this->GetMaskImage()->GetPixel(ItU.GetIndex()) != NumericTraits<MaskPixelType>::ZeroValue()) &&
        (!this->GetConfidenceImage() ||
         this->GetConfidenceImage()->GetPixel(ItU.GetIndex()) > NumericTraits<RealType>::ZeroValue()))
    {
      RealType pixel = ItU.Get();
      if (pixel > binMaximum)
      {
        binMaximum = pixel;
      }
      else if (pixel < binMinimum)
      {
        binMinimum = pixel;
      }
    }
  }
  RealType histogramSlope = (binMaximum - binMinimum) / static_cast<RealType>(this->m_NumberOfHistogramBins - 1);

  /**
   * Create the intensity profile (within the masked region, if applicable)
   * using a triangular parzen windowing scheme.
   */
  vnl_vector<RealType> H(this->m_NumberOfHistogramBins, NumericTraits<RealType>::ZeroValue());
  for (ItU.GoToBegin(); !ItU.IsAtEnd(); ++ItU)
  {
    if ((!this->GetMaskImage() ||
         this->GetMaskImage()->GetPixel(ItU.GetIndex()) != NumericTraits<MaskPixelType>::ZeroValue()) &&
        (!this->GetConfidenceImage() ||
         this->GetConfidenceImage()->GetPixel(ItU.GetIndex()) > NumericTraits<RealType>::ZeroValue()))
    {
      RealType pixel = ItU.Get();

      float        cidx = (static_cast<RealType>(pixel) - binMinimum) / histogramSlope;
      unsigned int idx = itk::Math::floor(cidx);
      RealType     offset = cidx - static_cast<RealType>(idx);

      if (Math::FloatAlmostEqual(offset, NumericTraits<RealType>::ZeroValue()))
      {
        H[idx] += NumericTraits<RealType>::OneValue();
      }
      else if (idx < this->m_NumberOfHistogramBins - 1)
      {
        H[idx] += NumericTraits<RealType>::OneValue() - offset;
        H[idx + 1] += offset;
      }
    }
  }

  /**
   * Determine information about the intensity histogram and zero-pad
   * histogram to a power of 2.
   */

  RealType exponent =
    std::ceil(std::log(static_cast<RealType>(this->m_NumberOfHistogramBins)) / static_cast<RealType>(std::log(2.0))) +
    1;
  unsigned int paddedHistogramSize = static_cast<unsigned int>(std::pow(static_cast<RealType>(2.0), exponent) +
                                                               static_cast<RealType>(static_cast<RealType>(0.5)));
  unsigned int histogramOffset =
    static_cast<unsigned int>(static_cast<RealType>(0.5) * (paddedHistogramSize - this->m_NumberOfHistogramBins));

  using FFTComputationType = double;
  using FFTComplexType = std::complex<FFTComputationType>;

  vnl_vector<FFTComplexType> V(paddedHistogramSize, FFTComplexType(0.0, 0.0));
  for (unsigned int n = 0; n < this->m_NumberOfHistogramBins; n++)
  {
    V[n + histogramOffset] = H[n];
  }

  /**
   * Instantiate the 1-d vnl fft routine
   */
  vnl_fft_1d<FFTComputationType> fft(paddedHistogramSize);

  vnl_vector<FFTComplexType> Vf(V);
  fft.fwd_transform(Vf);

  /**
   * Create the Gaussian filter.
   */
  RealType scaledFWHM = this->m_BiasFieldFullWidthAtHalfMaximum / histogramSlope;
  RealType expFactor = static_cast<RealType>(4.0) * static_cast<RealType>(std::log(2.0)) / Math::sqr(scaledFWHM);
  RealType scaleFactor = static_cast<RealType>(2.0) *
                         std::sqrt(static_cast<RealType>(std::log(2.0)) / static_cast<RealType>(Math::pi)) / scaledFWHM;

  vnl_vector<FFTComplexType> F(paddedHistogramSize, FFTComplexType(0.0, 0.0));
  F[0] = FFTComplexType(scaleFactor, 0.0);
  unsigned int halfSize = static_cast<unsigned int>(0.5 * paddedHistogramSize);
  for (unsigned int n = 1; n <= halfSize; n++)
  {
    F[n] = F[paddedHistogramSize - n] =
      FFTComplexType(scaleFactor * std::exp(-Math::sqr(static_cast<RealType>(n)) * expFactor), 0.0);
  }
  if (paddedHistogramSize % 2 == 0)
  {
    F[halfSize] = FFTComplexType(static_cast<FFTComputationType>(scaleFactor) *
                                   std::exp(0.25 * -Math::sqr(static_cast<FFTComputationType>(paddedHistogramSize)) *
                                            static_cast<FFTComputationType>(expFactor)),
                                 NumericTraits<FFTComputationType>::ZeroValue());
  }

  vnl_vector<FFTComplexType> Ff(F);
  fft.fwd_transform(Ff);

  /**
   * Create the Wiener deconvolution filter.
   */
  vnl_vector<FFTComplexType> Gf(paddedHistogramSize);

  const auto wienerNoiseValue = static_cast<FFTComputationType>(this->m_WienerFilterNoise);

  for (unsigned int n = 0; n < paddedHistogramSize; n++)
  {
    FFTComplexType c = vnl_complex_traits<FFTComplexType>::conjugate(Ff[n]);
    Gf[n] = c / (c * Ff[n] + wienerNoiseValue);
  }

  vnl_vector<FFTComplexType> Uf(paddedHistogramSize);
  for (unsigned int n = 0; n < paddedHistogramSize; n++)
  {
    Uf[n] = Vf[n] * Gf[n].real();
  }

  vnl_vector<FFTComplexType> U(Uf);
  fft.bwd_transform(U);
  for (unsigned int n = 0; n < paddedHistogramSize; n++)
  {
    U[n] = FFTComplexType(std::max(U[n].real(), static_cast<FFTComputationType>(0.0)), 0.0);
  }

  /**
   * Compute mapping E(u|v)
   */
  vnl_vector<FFTComplexType> numerator(paddedHistogramSize);
  for (unsigned int n = 0; n < paddedHistogramSize; n++)
  {
    numerator[n] =
      FFTComplexType((static_cast<FFTComputationType>(binMinimum) +
                      (static_cast<FFTComputationType>(n) - static_cast<FFTComputationType>(histogramOffset)) *
                        static_cast<FFTComputationType>(histogramSlope)) *
                       U[n].real(),
                     0.0);
  }
  fft.fwd_transform(numerator);
  for (unsigned int n = 0; n < paddedHistogramSize; n++)
  {
    numerator[n] *= Ff[n];
  }
  fft.bwd_transform(numerator);

  vnl_vector<FFTComplexType> denominator(U);
  fft.fwd_transform(denominator);
  for (unsigned int n = 0; n < paddedHistogramSize; n++)
  {
    denominator[n] *= Ff[n];
  }
  fft.bwd_transform(denominator);

  vnl_vector<RealType> E(paddedHistogramSize);
  for (unsigned int n = 0; n < paddedHistogramSize; n++)
  {
    E[n] = numerator[n].real() / denominator[n].real();
    if (std::isinf(E[n]) || std::isnan(E[n]))
    {
      E[n] = 0.0;
    }
  }

  /**
   * Remove the zero-padding from the mapping
   */
  E = E.extract(this->m_NumberOfHistogramBins, histogramOffset);

  /**
   * Sharpen the image with the new mapping, E(u|v)
   */
  typename RealImageType::Pointer sharpenedImage = AllocImage<RealImageType>(unsharpenedImage, 0.0);

  ImageRegionIterator<RealImageType> ItC(sharpenedImage, sharpenedImage->GetLargestPossibleRegion());
  for (ItU.GoToBegin(), ItC.GoToBegin(); !ItU.IsAtEnd(); ++ItU, ++ItC)
  {
    if ((!this->GetMaskImage() ||
         this->GetMaskImage()->GetPixel(ItU.GetIndex()) != NumericTraits<MaskPixelType>::ZeroValue()) &&
        (!this->GetConfidenceImage() ||
         this->GetConfidenceImage()->GetPixel(ItU.GetIndex()) > NumericTraits<RealType>::ZeroValue()))
    {
      float        cidx = (ItU.Get() - binMinimum) / histogramSlope;
      unsigned int idx = Math::floor(cidx);

      RealType correctedPixel = 0;
      if (idx < E.size() - 1)
      {
        correctedPixel = E[idx] + (E[idx + 1] - E[idx]) * (cidx - static_cast<RealType>(idx));
      }
      else
      {
        correctedPixel = E[E.size() - 1];
      }

      ItC.Set(correctedPixel);
    }
  }

  return sharpenedImage;
}

template <typename TInputImage, typename TMaskImage, typename TOutputImage>
typename N3MRIBiasFieldCorrectionImageFilter<TInputImage, TMaskImage, TOutputImage>::RealImageType::Pointer
N3MRIBiasFieldCorrectionImageFilter<TInputImage, TMaskImage, TOutputImage>::SmoothField(
  typename N3MRIBiasFieldCorrectionImageFilter<TInputImage, TMaskImage, TOutputImage>::RealImageType::Pointer
    fieldEstimate)
{
  /**
   * Get original direction and change to identity temporarily for the
   * b-spline fitting.
   */
  typename RealImageType::DirectionType direction = fieldEstimate->GetDirection();
  typename RealImageType::DirectionType identity;
  identity.SetIdentity();
  fieldEstimate->SetDirection(identity);

  typename PointSetType::Pointer fieldPoints = PointSetType::New();
  fieldPoints->Initialize();

  typename BSplineFilterType::WeightsContainerType::Pointer weights = BSplineFilterType::WeightsContainerType::New();
  weights->Initialize();

  ImageRegionConstIteratorWithIndex<RealImageType> It(fieldEstimate, fieldEstimate->GetRequestedRegion());
  unsigned int                                     N = 0;
  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    if ((!this->GetMaskImage() ||
         this->GetMaskImage()->GetPixel(It.GetIndex()) != NumericTraits<MaskPixelType>::ZeroValue()) &&
        (!this->GetConfidenceImage() ||
         this->GetConfidenceImage()->GetPixel(It.GetIndex()) > NumericTraits<RealType>::ZeroValue()))
    {
      typename PointSetType::PointType point;
      fieldEstimate->TransformIndexToPhysicalPoint(It.GetIndex(), point);

      ScalarType scalar;
      scalar[0] = It.Get();

      fieldPoints->SetPointData(N, scalar);
      fieldPoints->SetPoint(N, point);
      if (this->GetConfidenceImage())
      {
        weights->InsertElement(N, this->GetConfidenceImage()->GetPixel(It.GetIndex()));
      }
      else
      {
        weights->InsertElement(N, 1.0);
      }
      N++;
    }
  }
  fieldEstimate->SetDirection(direction);

  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
  bspliner->SetOrigin(fieldEstimate->GetOrigin());
  bspliner->SetSpacing(fieldEstimate->GetSpacing());
  bspliner->SetSize(fieldEstimate->GetLargestPossibleRegion().GetSize());
  bspliner->SetDirection(fieldEstimate->GetDirection());
  bspliner->SetGenerateOutputImage(true);
  bspliner->SetNumberOfLevels(this->m_NumberOfFittingLevels);
  bspliner->SetSplineOrder(this->m_SplineOrder);
  bspliner->SetNumberOfControlPoints(this->m_NumberOfControlPoints);
  bspliner->SetInput(fieldPoints);
  bspliner->SetPointWeights(weights);
  bspliner->Update();

  /**
   * Save the bias field control points in case the user wants to
   * reconstruct the bias field.
   */
  this->m_LogBiasFieldControlPointLattice = bspliner->GetPhiLattice();

  typename RealImageType::Pointer smoothField = AllocImage<RealImageType>(fieldEstimate->GetLargestPossibleRegion());
  smoothField->SetOrigin(fieldEstimate->GetOrigin());
  smoothField->SetSpacing(fieldEstimate->GetSpacing());
  smoothField->SetDirection(direction);

  ImageRegionIterator<ScalarImageType> ItB(bspliner->GetOutput(), bspliner->GetOutput()->GetLargestPossibleRegion());
  ImageRegionIterator<RealImageType>   ItF(smoothField, smoothField->GetLargestPossibleRegion());
  for (ItB.GoToBegin(), ItF.GoToBegin(); !ItB.IsAtEnd(); ++ItB, ++ItF)
  {
    ItF.Set(ItB.Get()[0]);
  }

  return smoothField;
}

template <typename TInputImage, typename TMaskImage, typename TOutputImage>
typename N3MRIBiasFieldCorrectionImageFilter<TInputImage, TMaskImage, TOutputImage>::RealType
N3MRIBiasFieldCorrectionImageFilter<TInputImage, TMaskImage, TOutputImage>::CalculateConvergenceMeasurement(
  typename RealImageType::Pointer fieldEstimate1,
  typename RealImageType::Pointer fieldEstimate2)
{
  typedef SubtractImageFilter<RealImageType, RealImageType, RealImageType> SubtracterType;
  typename SubtracterType::Pointer                                         subtracter = SubtracterType::New();
  subtracter->SetInput1(fieldEstimate1);
  subtracter->SetInput2(fieldEstimate2);
  subtracter->Update();

  /**
   * Calculate statistics over the mask region
   */
  RealType mu = 0.0;
  RealType sigma = 0.0;
  RealType N = 0.0;

  ImageRegionConstIteratorWithIndex<RealImageType> It(subtracter->GetOutput(),
                                                      subtracter->GetOutput()->GetLargestPossibleRegion());
  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    if (!this->GetMaskImage() ||
        this->GetMaskImage()->GetPixel(It.GetIndex()) != NumericTraits<MaskPixelType>::ZeroValue())
    {
      RealType pixel = It.Get();
      N += NumericTraits<RealType>::OneValue();

      if (N > NumericTraits<RealType>::OneValue())
      {
        sigma = sigma + Math::sqr(pixel - mu) * (N - NumericTraits<RealType>::OneValue()) / N;
      }
      mu = mu * (NumericTraits<RealType>::OneValue() - NumericTraits<RealType>::OneValue() / N) + pixel / N;
    }
  }
  sigma = std::sqrt(sigma / (N - NumericTraits<RealType>::OneValue()));

  /**
   * Although Sled's paper proposes convergence determination via
   * the coefficient of variation, the actual mnc implementation
   * utilizes the standard deviation as the convergence measurement.
   */
  return sigma;
}

template <typename TInputImage, typename TMaskImage, typename TOutputImage>
typename N3MRIBiasFieldCorrectionImageFilter<TInputImage, TMaskImage, TOutputImage>::RealType
N3MRIBiasFieldCorrectionImageFilter<TInputImage, TMaskImage, TOutputImage>::CalculateOptimalBiasFieldScaling(
  typename RealImageType::Pointer biasField)
{
  /**
   * This section is not described in Sled's paper but rather stems from
   * our own experience with N4ITK and the resulting innovation.  For an initial
   * B-spline mesh of large resolution and a low number of fitting levels,
   * although the shape of of the bias field appears correect, the scale is too
   * small. This section finds an optimal scaling by minimizing the coefficient
   * of variation over masked region.
   */

  typedef N3BiasFieldScaleCostFunction<InputImageType, RealImageType, MaskImageType, RealImageType>
                                          ScaleCostFunctionType;
  typename ScaleCostFunctionType::Pointer scaleCostFunction = ScaleCostFunctionType::New();

  scaleCostFunction->SetInputImage(const_cast<InputImageType *>(this->GetInput()));
  scaleCostFunction->SetBiasFieldImage(biasField);
  scaleCostFunction->SetMaskImage(const_cast<MaskImageType *>(this->GetMaskImage()));
  scaleCostFunction->SetConfidenceImage(const_cast<RealImageType *>(this->GetConfidenceImage()));

  typename LBFGSBOptimizer::BoundSelectionType boundSelection;
  boundSelection.SetSize(1);
  boundSelection.Fill(1); // only set a lower bound on the scale factor
  typename LBFGSBOptimizer::BoundValueType lowerBound;
  lowerBound.SetSize(1);
  lowerBound.Fill(1.0);
  typename LBFGSBOptimizer::BoundValueType upperBound;
  upperBound.SetSize(1);
  upperBound.Fill(NumericTraits<typename LBFGSBOptimizer::BoundValueType::ValueType>::max());
  typename LBFGSBOptimizer::ParametersType initialParameters;
  initialParameters.SetSize(1);
  initialParameters.Fill(1.0);

  typename LBFGSBOptimizer::Pointer optimizer = LBFGSBOptimizer::New();
  optimizer->SetMinimize(true);
  optimizer->SetCostFunction(scaleCostFunction);
  optimizer->SetInitialPosition(initialParameters);
  optimizer->SetCostFunctionConvergenceFactor(1e1);
  optimizer->SetLowerBound(lowerBound);
  optimizer->SetUpperBound(upperBound);
  optimizer->SetBoundSelection(boundSelection);
  optimizer->SetProjectedGradientTolerance(1e-10);

  optimizer->StartOptimization();

  typename LBFGSBOptimizer::ParametersType finalParameters = optimizer->GetCurrentPosition();

  return finalParameters[0];
}

template <typename TInputImage, typename TMaskImage, typename TOutputImage>
void
N3MRIBiasFieldCorrectionImageFilter<TInputImage, TMaskImage, TOutputImage>::PrintSelf(std::ostream & os,
                                                                                      Indent         indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Number of histogram bins: " << this->m_NumberOfHistogramBins << std::endl;
  os << indent << "Wiener filter noise: " << this->m_WienerFilterNoise << std::endl;
  os << indent << "Bias field FWHM: " << this->m_BiasFieldFullWidthAtHalfMaximum << std::endl;
  os << indent << "Maximum number of iterations: " << this->m_MaximumNumberOfIterations << std::endl;
  os << indent << "Convergence threshold: " << this->m_ConvergenceThreshold << std::endl;
  os << indent << "Spline order: " << this->m_SplineOrder << std::endl;
  os << indent << "Number of fitting levels: " << this->m_NumberOfFittingLevels << std::endl;
  os << indent << "Number of control points: " << this->m_NumberOfControlPoints << std::endl;
  if (this->m_UseOptimalBiasFieldScaling)
  {
    os << indent << "Optimal bias field scaling: " << this->m_BiasFieldScaling << std::endl;
  }
}
} // end namespace itk
#endif
