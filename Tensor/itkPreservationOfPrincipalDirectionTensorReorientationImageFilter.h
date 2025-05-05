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
#include "itkDisplacementFieldTransform.h"

namespace itk
{
/** \class PreservationOfPrincipalDirectionImageFilter
 * \brief Applies an averaging filter to an image
 *
 * Computes an image where a given pixel is the mean value of the
 * the pixels in a neighborhood about the corresponding input pixel.
 *
 * A mean filter is one of the family of linear filters.
 *
 * \sa Image
 * \sa Neighborhood
 * \sa NeighborhoodOperator
 * \sa NeighborhoodIterator
 *
 * \ingroup IntensityImageFilters
 */
template <typename TTensorImage, typename TVectorImage>
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
  using DisplacementFieldType = TVectorImage;

  using DisplacementFieldPointer = typename DisplacementFieldType::Pointer;
  using VectorType = typename DisplacementFieldType::PixelType;
  using RealType = typename VectorType::RealValueType;

  using DisplacementFieldTransformType = DisplacementFieldTransform<double, 3>;
  using DisplacementFieldTransformPointer = typename DisplacementFieldTransformType::Pointer;

  using MatrixType = Matrix<RealType, 3, 3>;
  using VariableMatrixType = VariableSizeMatrix<RealType>;

  static constexpr unsigned int ImageDimension = TTensorImage::ImageDimension;

  using RealTypeImageType = Image<RealType, ImageDimension>;
  using AffineTransformType = MatrixOffsetTransformBase<RealType, ImageDimension, ImageDimension>;
  using AffineTransformPointer = typename AffineTransformType::Pointer;
  using InverseTransformType = typename AffineTransformType::InverseTransformBaseType;
  using InverseTransformPointer = typename InverseTransformType::Pointer;

  using VnlMatrixType = vnl_matrix<RealType>;
  using VnlVectorType = vnl_vector<RealType>;

  using InputImagePointer = typename InputImageType::ConstPointer;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using InputPixelType = typename InputImageType::PixelType;
  using OutputPixelType = typename OutputImageType::PixelType;
  using InputRealType = typename InputPixelType::ValueType;

  using InputImageRegionType = typename InputImageType::RegionType;
  using OutputImageRegionType = typename OutputImageType::RegionType;

  using InputSizeType = typename InputImageType::SizeType;
  using OutputSizeType = typename OutputImageType::SizeType;
  using InputIndexType = typename InputImageType::IndexType;
  using OutputIndexType = typename OutputImageType::IndexType;

  itkSetMacro(DisplacementField, DisplacementFieldPointer);
  itkGetMacro(DisplacementField, DisplacementFieldPointer);

  void SetAffineTransform(AffineTransformPointer aff)
  {
    this->m_AffineTransform = aff;
    this->m_UseAffine = true;
  }

  itkGetMacro(AffineTransform, AffineTransformPointer);

protected:
  PreservationOfPrincipalDirectionTensorReorientationImageFilter();
  ~PreservationOfPrincipalDirectionTensorReorientationImageFilter() override = default;

  void PrintSelf(std::ostream & os, Indent indent) const override;
  void GenerateData() override;

private:
  PreservationOfPrincipalDirectionTensorReorientationImageFilter(const Self &) = delete;
  void operator=(const Self &) = delete;

  DisplacementFieldPointer m_DisplacementField;
  DisplacementFieldTransformPointer m_DisplacementTransform;
  AffineTransformPointer m_AffineTransform;
  bool m_UseAffine = false;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkPreservationOfPrincipalDirectionTensorReorientationImageFilter.cxx"
#endif

#endif
