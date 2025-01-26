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
#ifndef __itkAlternatingValueSimpleSubtractionImageFilter_h
#define __itkAlternatingValueSimpleSubtractionImageFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{
/** \class AlternatingValueSimpleSubtractionImageFilter
 * \brief Finds difference signal from alternating signal
 *
 * This filter is templated over the input image type and the output image
 * type. The output image is the difference between a two alternating
 * signals that alternate over the subtraction dimension.
 *
 * \ingroup GeometricTransform
 * \ingroup MultiThreaded
 * \ingroup Streamed
 *
 * \author Jeffrey Duda
 *
 */
template <typename TInputImage, typename TOutputImage>
class AlternatingValueSimpleSubtractionImageFilter final : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef AlternatingValueSimpleSubtractionImageFilter  Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(AlternatingValueSimpleSubtractionImageFilter);

  /** Compiler can't inherit typedef? */
  typedef typename Superclass::InputImageType  InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;
  typedef typename InputImageType::Pointer     InputImagePointer;
  typedef typename OutputImageType::Pointer    OutputImagePointer;
  typedef typename InputImageType::RegionType  InputImageRegionType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  /** Compiler can't inherit ImageDimension enumeration? */
  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int OutputImageDimension = TOutputImage::ImageDimension;

  /** Set/Get dimension to subtract over */
  itkGetMacro(SubtractionDimension, unsigned int);
  itkSetMacro(SubtractionDimension, unsigned int);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputConvertibleToOutputCheck,
                  (Concept::Convertible<typename TInputImage::PixelType, typename TOutputImage::PixelType>));
  itkConceptMacro(DimensionCheck, (Concept::SameDimension<InputImageDimension, OutputImageDimension>));
  /** End concept checking */
#endif
protected:
  AlternatingValueSimpleSubtractionImageFilter();
  ~AlternatingValueSimpleSubtractionImageFilter() override = default;

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

  /** AlternatingValueSimpleSubtractionImageFilter can be implemented as a multithreaded filter.
   * \sa ImageSource::ThreadedGenerateData(),
   *     ImageSource::GenerateData() */
  void
  ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId) override;

private:
  AlternatingValueSimpleSubtractionImageFilter(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  /** IndexValueType is used to switch among the inputs and
   * is used as the index value of the new dimension */
  typedef unsigned int IndexValueType;

  unsigned int m_SubtractionDimension;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkAlternatingValueSimpleSubtractionImageFilter.hxx"
#endif

#endif
