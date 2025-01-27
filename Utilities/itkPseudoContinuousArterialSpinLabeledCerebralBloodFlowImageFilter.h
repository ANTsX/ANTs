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
#ifndef __itkPseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter_h
#define __itkPseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{
/** \class PseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter
 * \brief Join N-D images into an (N+1)-D image
 *
 * This filter is templated over the input image type and the output image
 * type. The pixel type of them must be the same and the input dimension
 * must be less than the output dimension.
 * When the input images are N-dimensinal, they are joined in order and
 * the size of the N+1'th dimension of the output is same as the number of
 * the inputs. The spacing and the origin (where the first input is placed)
 * for the N+1'th dimension is specified in this filter. The output image
 * informations for the first N dimensions are taken from the first input.
 * Note that all the inputs should have the same information.
 *
 * \ingroup GeometricTransform
 * \ingroup MultiThreaded
 * \ingroup Streamed
 *
 * \author Hideaki Hiraki
 *
 * Contributed in the users list
 * http://public.kitware.com/pipermail/insight-users/2004-February/006542.html
 *
 * \ingroup ITKImageCompose
 */
template <typename TInputImage, typename TReferenceImage, typename TOutputImage>
class PseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter
  : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef PseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage>                   Superclass;
  typedef SmartPointer<Self>                                              Pointer;
  typedef SmartPointer<const Self>                                        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(PseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter);

  /** Compiler can't inherit typedef? */
  typedef TInputImage                             InputImageType;
  typedef typename Superclass::OutputImageType    OutputImageType;
  typedef TReferenceImage                         ReferenceImageType;
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename OutputImageType::Pointer       OutputImagePointer;
  typedef typename ReferenceImageType::Pointer    ReferenceImagePointer;
  typedef typename InputImageType::RegionType     InputImageRegionType;
  typedef typename OutputImageType::RegionType    OutputImageRegionType;
  typedef typename ReferenceImageType::RegionType ReferenceImageRegionType;

  /** Compiler can't inherit ImageDimension enumeration? */
  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int OutputImageDimension = TOutputImage::ImageDimension;
  static constexpr unsigned int ReferenceImageDimension = TReferenceImage::ImageDimension;

  itkGetMacro(TI1, float);
  itkSetMacro(TI1, float);

  itkGetMacro(TI2, float);
  itkSetMacro(TI2, float);

  itkGetMacro(T1blood, float);
  itkSetMacro(T1blood, float);

  itkGetMacro(Alpha, float);
  itkSetMacro(Alpha, float);

  itkGetMacro(Lambda, float);
  itkSetMacro(Lambda, float);

  itkGetMacro(SliceDelay, float);
  itkSetMacro(SliceDelay, float);

  void
  SetDifferenceImage(const InputImageType * image);

  void
  SetReferenceImage(const ReferenceImageType * image);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputConvertibleToOutputCheck,
                  (Concept::Convertible<typename TInputImage::PixelType, typename TOutputImage::PixelType>));
  /** End concept checking */
#endif
protected:
  PseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter();
  ~PseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter() {}

  typename TInputImage::ConstPointer
  GetDifferenceImage();
  typename TReferenceImage::ConstPointer
  GetReferenceImage();

  void
  PrintSelf(std::ostream & os, Indent indent) const;

  /** PseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter can be implemented as a multithreaded filter.
   * \sa ImageSource::ThreadedGenerateData(),
   *     ImageSource::GenerateData() */
  virtual void
  ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId);

private:
  PseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  /** IndexValueType is used to switch among the inputs and
   * is used as the index value of the new dimension */
  typedef unsigned int IndexValueType;

  float m_TI1;
  float m_TI2;
  float m_T1blood;
  float m_Lambda;
  float m_Alpha;
  float m_SliceDelay;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkPseudoContinuousArterialSpinLabeledCerebralBloodFlowImageFilter.hxx"
#endif

#endif
