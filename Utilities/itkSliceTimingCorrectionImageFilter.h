/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkSliceTimingCorrectionImageFilter_h
#define __itkSliceTimingCorrectionImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkExtrapolateImageFunction.h"
#include "itkInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"

namespace itk
{
/** \class SliceTimingCorrectionImageFilter
 * \brief Finds difference signal from alternating signal
 *
 * This filter is templated over the input image type and the output image
 * type. Each signal is interpolated over the entire range of the
 * subtraction dimension. The output image is the difference between
 * the two intepolated signals.
 *
 * \ingroup GeometricTransform
 * \ingroup MultiThreaded
 * \ingroup Streamed
 *
 * \author Jeffrey Duda
 *
 */
template <typename TInputImage, typename TOutputImage>
class SliceTimingCorrectionImageFilter final : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef SliceTimingCorrectionImageFilter              Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(SliceTimingCorrectionImageFilter);

  /** Compiler can't inherit typedef? */
  typedef typename Superclass::InputImageType  InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;
  typedef typename InputImageType::Pointer     InputImagePointer;
  typedef typename OutputImageType::Pointer    OutputImagePointer;
  typedef typename InputImageType::RegionType  InputImageRegionType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef LinearInterpolateImageFunction<InputImageType> DefaultInterpolatorType;
  typedef typename DefaultInterpolatorType::Pointer      DefaultInterpolatorPointerType;

  typedef InterpolateImageFunction<InputImageType, double> InterpolatorType;
  typedef typename InterpolatorType::Pointer               InterpolatorPointerType;

  // typedef typename InputImageType::SpacingType::ValueType     TimingType;
  typedef double TimingType;

  /** Compiler can't inherit ImageDimension enumeration? */
  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int OutputImageDimension = TOutputImage::ImageDimension;

  itkGetMacro(TimeDimension, unsigned int);
  itkSetMacro(TimeDimension, unsigned int);

  itkSetMacro(SliceDimension, unsigned int);
  itkGetMacro(SliceDimension, unsigned int);

  itkGetMacro(IndexPadding, unsigned int);
  itkSetMacro(IndexPadding, unsigned int);

  itkSetMacro(SliceTiming, TimingType);
  itkGetMacro(SliceTiming, TimingType);

  itkSetMacro(ExtrapolateEdges, bool);
  itkGetMacro(ExtrapolateEdges, bool);

  /** Set the interpolator function.  The default is
   * LinearInterpolateImageFunction<InputImageType,
   * TInterpolatorPrecisionType>. WindowedSincInterpolateImageFunction
   * is often suggested in the literature. Some
   * other options are NearestNeighborInterpolateImageFunction
   * (useful for binary masks and other images with a small number of
   * possible pixel values), and BSplineInterpolateImageFunction
   * (which provides a higher order of interpolation).  */
  itkSetObjectMacro(Interpolator, InterpolatorType);

  /** Get a pointer to the interpolator function. */
  itkGetConstObjectMacro(Interpolator, InterpolatorType);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputConvertibleToOutputCheck,
                  (Concept::Convertible<typename TInputImage::PixelType, typename TOutputImage::PixelType>));
  itkConceptMacro(DimensionCheck, (Concept::SameDimension<InputImageDimension, OutputImageDimension>));

  /** End concept checking */
#endif
protected:
  SliceTimingCorrectionImageFilter();
  ~SliceTimingCorrectionImageFilter() override {}

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /** Override VeriyInputInformation() to add the additional check
   * that all inputs have the same number of components.
   *
   * \sa ProcessObject::VerifyInputInformation
   */
  void
  VerifyInputInformation() const override;

  /** Overrides GenerateOutputInformation() in order to produce
   * an image which has a different information than the first input.
   * \sa ProcessObject::GenerateOutputInformaton() */
  void
  GenerateOutputInformation() override;

  /** Overrides GenerateInputRequestedRegion() in order to inform
   * the pipeline execution model of different input requested regions
   * than the output requested region.
   * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
  // virtual void GenerateInputRequestedRegion();

  /** This method is used to set the state of the filter before
   * multi-threading. */
  void
  BeforeThreadedGenerateData() override;

  /** SliceTimingCorrectionImageFilter can be implemented as a multithreaded filter.
   * \sa ImageSource::ThreadedGenerateData(),
   *     ImageSource::GenerateData() */
  void
  ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId) override;

private:
  SliceTimingCorrectionImageFilter(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  /** IndexValueType is used to switch among the inputs and
   * is used as the index value of the new dimension */
  typedef unsigned int IndexValueType;

  unsigned int m_TimeDimension;

  unsigned int m_SliceDimension;

  unsigned int m_IndexPadding;

  TimingType m_SliceTiming;

  InterpolatorPointerType m_Interpolator; // Image function for
                                          // interpolation

  bool m_ExtrapolateEdges;

  // ExtrapolatorPointerType m_Extrapolator;      // Image function for
  // extrapolation
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkSliceTimingCorrectionImageFilter.hxx"
#endif

#endif
