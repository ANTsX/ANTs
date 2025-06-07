/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkOptimalSharpeningImageFilter_h
#define __itkOptimalSharpeningImageFilter_h

#include "itkNumericTraits.h"
#include "itkImageToImageFilter.h"
#include "itkImage.h"

namespace itk
{
/**
 * \class OptimalSharpeningImageFilter
 *
 * This filter sharpens an image using a Laplacian. OptimalSharpening
 * highlights regions of rapid intensity change and therefore
 * highlights or enhances the edges.  The result is an image that
 * appears more in focus.
 *
 * \par The OptimalSharpening at each pixel location is computed by
 * convolution with the itk::LaplacianOperator.
 *
 * \par Inputs and Outputs
 * The input to this filter is a scalar-valued itk::Image of arbitrary
 * dimension. The output is a scalar-valued itk::Image.
 * \sa Image
 * \sa Neighborhood
 * \sa NeighborhoodOperator
 * \sa NeighborhoodIterator
 * \sa LaplacianOperator
 *
 * \ingroup ImageFeatureExtraction */
template <typename TInputImage, typename TOutputImage>
class OptimalSharpeningImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard "Self" & Superclass typedef.   */
  typedef OptimalSharpeningImageFilter                  Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;

  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  typedef typename TOutputImage::PixelType                  OutputPixelType;
  typedef typename TOutputImage::InternalPixelType          OutputInternalPixelType;
  typedef typename NumericTraits<OutputPixelType>::RealType RealType;
  typedef typename TInputImage::PixelType                   InputPixelType;
  typedef typename TInputImage::InternalPixelType           InputInternalPixelType;
  static constexpr unsigned int ImageDimension = TOutputImage::ImageDimension;

  /** Image typedef support. */
  typedef TInputImage                      InputImageType;
  typedef TOutputImage                     OutputImageType;
  typedef typename InputImageType::Pointer InputImagePointer;

  /** Smart pointer typedef support.   */
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods)  */
  itkOverrideGetNameOfClassMacro(OptimalSharpeningImageFilter);

  /** Method for creation through the object factory.  */
  itkNewMacro(Self);

  /** OptimalSharpeningImageFilter needs a larger input requested
   * region than the output requested region (larger in the direction
   * of the derivative).  As such, OptimalSharpeningImageFilter
   * needs to provide an implementation for
   * GenerateInputRequestedRegion() in order to inform the pipeline
   * execution model.
   *
   * \sa ImageToImageFilter::GenerateInputRequestedRegion()  */
  virtual void
  GenerateInputRequestedRegion();

  void
  SetSValue(float s)
  {
    this->m_SValue = s;
  }

  /** Use the image spacing information in calculations. Use this option if you
   *  want derivatives in physical space. Default is UseImageSpacingOn. */
  void
  UseImageSpacingOn()
  {
    this->SetUseImageSpacing(true);
  }

  /** Ignore the image spacing. Use this option if you want derivatives in
      isotropic pixel space.  Default is UseImageSpacingOn. */
  void
  UseImageSpacingOff()
  {
    this->SetUseImageSpacing(false);
  }

  /** Set/Get whether or not the filter will use the spacing of the input
      image in its calculations */
  itkSetMacro(UseImageSpacing, bool);
  itkGetMacro(UseImageSpacing, bool);

protected:
  OptimalSharpeningImageFilter()
  {
    m_UseImageSpacing = true;
    m_SValue = 0.5;
  }

  virtual ~OptimalSharpeningImageFilter() = default;

  /** Standard pipeline method. While this class does not implement a
   * ThreadedGenerateData(), its GenerateData() delegates all
   * calculations to an NeighborhoodOperatorImageFilter.  Since the
   * NeighborhoodOperatorImageFilter is multithreaded, this filter is
   * multithreaded by default.   */
  void
  GenerateData();

  void
  PrintSelf(std::ostream &, Indent) const;

private:
  OptimalSharpeningImageFilter(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  bool  m_UseImageSpacing;
  float m_SValue;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkOptimalSharpeningImageFilter.hxx"
#endif

#endif
