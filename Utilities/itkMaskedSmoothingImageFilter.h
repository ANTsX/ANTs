/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMaskedSmoothingImageFilter_h
#define __itkMaskedSmoothingImageFilter_h

#include "itkImageToImageFilter.h"

#include <vnl/vnl_sparse_matrix.h>


namespace itk
{
/** \class MaskedSmoothImageFilter
 * \brief Use a mask to spatially smooth a scalar image or displacement
 * field.  This also includes longitudinal images.
 *
 * \author Nick Tustison, Brian Avants
 *
 */

template <typename TInputImage,
          typename TMaskImage = Image<unsigned char, TInputImage::ImageDimension>,
          class TOutputImage = TInputImage>
class MaskedSmoothingImageFilter final : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  using Self = MaskedSmoothingImageFilter;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Extract dimension from input and output image. */
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;

  /** Convenient typedefs for simplifying declarations. */
  using InputImageType = TInputImage;
  using InputImagePointer = typename InputImageType::Pointer;
  using InputPixelType = typename InputImageType::PixelType;
  using RegionType = typename InputImageType::RegionType;
  using IndexType = typename InputImageType::IndexType;

  using MaskImageType = TMaskImage;
  using MaskPixelType = typename MaskImageType::PixelType;

  using OutputImageType = TOutputImage;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using OutputPixelType = typename OutputImageType::PixelType;

  using RealType = float;
  using RealImageType = Image<RealType, ImageDimension>;
  using RealImagePointer = typename RealImageType::Pointer;
  using SparseMatrixType = vnl_sparse_matrix<RealType>;

  /**
   * Set mask image function.  If a binary mask image is specified, only
   * those input image voxels inside the mask image values are used.
   */
  void
  SetMaskImage(const MaskImageType * mask)
  {
    this->SetNthInput(1, const_cast<MaskImageType *>(mask));
  }
  void
  SetInput2(const MaskImageType * mask)
  {
    this->SetMaskImage(mask);
  }

  /**
   * Get mask image function.
   */
  const MaskImageType *
  GetMaskImage() const
  {
    return static_cast<const MaskImageType *>(this->ProcessObject::GetInput(1));
  }

  /**
   * Set/Get the sparse image neighborhood radius.  Default = 2.
   */
  itkSetMacro(SparseImageNeighborhoodRadius, unsigned int);
  itkGetConstMacro(SparseImageNeighborhoodRadius, unsigned int);

  /**
   * Set/Get the variance for spatial regularization.
   */
  itkSetMacro(SmoothingVariance, RealType);
  itkGetConstMacro(SmoothingVariance, RealType);

  /**
   * Set/Get the variance for time regularization.
   */
  itkSetMacro(TimeSmoothingVariance, RealType);
  itkGetConstMacro(TimeSmoothingVariance, RealType);

  /**
   * Set/Get the time point values.  Default = no special value.
   */
  void
  SetTimePoints(std::vector<RealType> timePoints)
  {
    this->m_TimePoints = timePoints;
    this->Modified();
  }
  itkGetConstMacro(TimePoints, std::vector<RealType>);

protected:
  MaskedSmoothingImageFilter();
  ~MaskedSmoothingImageFilter() override;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  GenerateData() override;

private:
  unsigned int     m_SparseImageNeighborhoodRadius;
  SparseMatrixType m_SparseMatrix;
  RealImagePointer m_SparseMatrixIndexImage;
  RealType         m_SmoothingVariance;

  std::vector<RealType> m_TimePoints;
  RealType              m_TimeSmoothingVariance;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkMaskedSmoothingImageFilter.hxx"
#endif

#endif
