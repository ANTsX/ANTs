/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDeterminantTensorImageFilter_h
#define __itkDeterminantTensorImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkMatrix.h"
#include "itkVector.h"

namespace itk
{
/** \class DeterminantTensorImageFilter
 *
 */
template <typename TInputImage,
          typename TRealType = float,
          typename TOutputImage = Image<TRealType, TInputImage::ImageDimension>>
class ITK_EXPORT DeterminantTensorImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef DeterminantTensorImageFilter                  Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkOverrideGetNameOfClassMacro(DeterminantTensorImageFilter);

  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  typedef typename TOutputImage::PixelType OutputPixelType;
  typedef typename TInputImage::PixelType  InputPixelType;

  /** Image typedef support */
  typedef TInputImage                          InputImageType;
  typedef TOutputImage                         OutputImageType;
  typedef typename InputImageType::Pointer     InputImagePointer;
  typedef typename OutputImageType::Pointer    OutputImagePointer;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  /** The dimensionality of the input and output images. */
  static constexpr unsigned int ImageDimension = TOutputImage::ImageDimension;

  /** Define the data type and the vector of data type used in calculations. */
  typedef TRealType RealType;

protected:
  DeterminantTensorImageFilter();
  ~DeterminantTensorImageFilter() override = default;

  /** DeterminantTensorImageFilter can be implemented as a
   * multithreaded filter (we're only using vnl_det(), which is trivially
   * thread safe).  Therefore, this implementation provides a
   * ThreadedGenerateData() routine which is called for each
   * processing thread. The output image data is allocated
   * automatically by the superclass prior to calling
   * ThreadedGenerateData().  ThreadedGenerateData can only write to
   * the portion of the output image specified by the parameter
   * "outputRegionForThread"
   *
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData() */
  void
  ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId) override;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  DeterminantTensorImageFilter(const Self &) = delete;
  void
  operator=(const Self &) = delete;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkDeterminantTensorImageFilter.hxx"
#endif

#endif
