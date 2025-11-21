/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPreservationOfVectorDirectionImageFilter_h
#define __itkPreservationOfVectorDirectionImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkMatrix.h"
#include "itkNumericTraits.h"
#include "itkVector.h"
#include "itkCompositeTransform.h"

namespace itk
{
/** PreservationOfPrincipalDirectionImageFilter
 *
 * Reorients Vector images, preserving their magnitude.
 *
 *
 *
 */
template <typename TVectorImage>
class PreservationOfVectorDirectionImageFilter final
  : public ImageToImageFilter<TVectorImage, TVectorImage>
{
public:
  using Self = PreservationOfVectorDirectionImageFilter;
  using Superclass = ImageToImageFilter<TVectorImage, TVectorImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(PreservationOfVectorDirectionImageFilter);

  using InputImageType = TVectorImage;
  using OutputImageType = TVectorImage;

  using InputPixelType = typename InputImageType::PixelType;
  using RealType = typename InputPixelType::ValueType;

  static constexpr unsigned int ImageDimension = TVectorImage::ImageDimension;

  // use double for vnl matrices instead of RealType, for compatiblity with direction
  // matrices, which are always double
  using MatrixType = vnl_matrix_fixed<double, ImageDimension, ImageDimension>;
  using VectorMatrixType = vnl_matrix_fixed<double, ImageDimension, 1>;

  using RealTypeImageType = Image<RealType, ImageDimension>;

  using CompositeTransformType = CompositeTransform<RealType, ImageDimension>;
  using CompositeTransformPointer = typename CompositeTransformType::Pointer;

  using JacobianMatrixType = typename CompositeTransformType::InverseJacobianPositionType;

  using InputImagePointer = typename InputImageType::ConstPointer;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using OutputPixelType = typename OutputImageType::PixelType;

  using InputImageRegionType = typename InputImageType::RegionType;
  using OutputImageRegionType = typename OutputImageType::RegionType;

  using InputSizeType = typename InputImageType::SizeType;
  using OutputSizeType = typename OutputImageType::SizeType;
  using InputIndexType = typename InputImageType::IndexType;
  using OutputIndexType = typename OutputImageType::IndexType;

  itkSetMacro(CompositeTransform, CompositeTransformPointer);
  itkGetMacro(CompositeTransform, CompositeTransformPointer);

  // if input vectors are in physical space, we leave the output in physical space also
  itkSetMacro(InputVectorsInPhysicalSpace, bool);
  itkGetMacro(InputVectorsInPhysicalSpace, bool);

protected:
  PreservationOfVectorDirectionImageFilter();
  ~PreservationOfVectorDirectionImageFilter() override = default;

  void PrintSelf(std::ostream & os, Indent indent) const override;
  void GenerateData() override;

private:
  PreservationOfVectorDirectionImageFilter(const Self &) = delete;
  void operator=(const Self &) = delete;

  CompositeTransformPointer m_CompositeTransform;

  bool m_InputVectorsInPhysicalSpace{ false };

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkPreservationOfVectorDirectionImageFilter.cxx"
#endif

#endif
