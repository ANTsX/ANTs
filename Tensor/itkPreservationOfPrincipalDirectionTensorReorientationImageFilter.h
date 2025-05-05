/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPreservationOfPrincipalDirectionTensorReorientationImageFilter_h
#define __itkPreservationOfPrincipalDirectionTensorReorientationImageFilter_h

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
 * Reorients tensor images to preserve the principal directions of diffusion.
 *
 * This filter rebases tensors into physical space using the image direction
 * transform, and then applies the reorientation computed by the composite
 * transform.
 *
 */
template <typename TTensorImage>
class PreservationOfPrincipalDirectionTensorReorientationImageFilter final
  : public ImageToImageFilter<TTensorImage, TTensorImage>
{
public:
  using Self = PreservationOfPrincipalDirectionTensorReorientationImageFilter;
  using Superclass = ImageToImageFilter<TTensorImage, TTensorImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(PreservationOfPrincipalDirectionTensorReorientationImageFilter);

  using InputImageType = TTensorImage;
  using OutputImageType = TTensorImage;

  using InputPixelType = typename InputImageType::PixelType;
  using RealType = typename InputPixelType::ValueType;

  using MatrixType = Matrix<RealType, 3, 3>;
  using VariableMatrixType = VariableSizeMatrix<RealType>;

  static constexpr unsigned int ImageDimension = TTensorImage::ImageDimension;

  using RealTypeImageType = Image<RealType, ImageDimension>;

  using CompositeTransformType = CompositeTransform<RealType, ImageDimension>;
  using CompositeTransformPointer = typename CompositeTransformType::Pointer;

  using VnlMatrixType = vnl_matrix<RealType>;
  using VnlVectorType = vnl_vector<RealType>;

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

protected:
  PreservationOfPrincipalDirectionTensorReorientationImageFilter();
  ~PreservationOfPrincipalDirectionTensorReorientationImageFilter() override = default;

  void PrintSelf(std::ostream & os, Indent indent) const override;
  void GenerateData() override;

private:
  PreservationOfPrincipalDirectionTensorReorientationImageFilter(const Self &) = delete;
  void operator=(const Self &) = delete;

  CompositeTransformPointer m_CompositeTransform;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkPreservationOfPrincipalDirectionTensorReorientationImageFilter.cxx"
#endif

#endif
