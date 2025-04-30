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
  typedef TTensorImage InputImageType;
  typedef TTensorImage OutputImageType;
  typedef TVectorImage DisplacementFieldType;

  typedef typename DisplacementFieldType::Pointer   DisplacementFieldPointer;
  typedef typename DisplacementFieldType::PixelType VectorType;
  typedef typename VectorType::RealValueType        RealType;

  typedef itk::DisplacementFieldTransform<double, 3> DisplacementFieldTransformType;

  typedef typename DisplacementFieldTransformType::Pointer DisplacementFieldTransformPointer;

  typedef Matrix<RealType, 3, 3> MatrixType;
  // typedef Vector<RealType, 3> VectorType;
  typedef VariableSizeMatrix<RealType> VariableMatrixType;

  static constexpr unsigned int ImageDimension = TTensorImage::ImageDimension;

  typedef itk::Image<RealType, ImageDimension> RealTypeImageType;

  typedef itk::MatrixOffsetTransformBase<RealType, ImageDimension, ImageDimension> AffineTransformType;

  typedef typename AffineTransformType::Pointer AffineTransformPointer;

  typedef typename AffineTransformType::InverseTransformBaseType InverseTransformType;

  typedef typename InverseTransformType::Pointer InverseTransformPointer;

  typedef vnl_matrix<RealType> VnlMatrixType;
  typedef vnl_vector<RealType> VnlVectorType;
  typedef vnl_vector<RealType> vvec;

  /** Standard class typedefs. */
  typedef PreservationOfPrincipalDirectionTensorReorientationImageFilter Self;
  typedef ImageToImageFilter<InputImageType, OutputImageType>            Superclass;
  typedef SmartPointer<Self>                                             Pointer;
  typedef SmartPointer<const Self>                                       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  itkSetMacro(DisplacementField, DisplacementFieldPointer);
  itkGetMacro(DisplacementField, DisplacementFieldPointer);

  void
  SetAffineTransform(AffineTransformPointer aff)
  {
    this->m_AffineTransform = aff;
    this->m_UseAffine = true;
  }

  itkGetMacro(AffineTransform, AffineTransformPointer);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(PreservationOfPrincipalDirectionTensorReorientationImageFilter);

  /** Image typedef support. */
  typedef typename InputImageType::ConstPointer InputImagePointer;
  typedef typename OutputImageType::Pointer     OutputImagePointer;
  typedef typename InputImageType::PixelType    InputPixelType;
  typedef typename OutputImageType::PixelType   OutputPixelType;
  typedef typename InputPixelType::ValueType    InputRealType;

  typedef typename InputImageType::RegionType  InputImageRegionType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef typename InputImageType::SizeType   InputSizeType;
  typedef typename OutputImageType::SizeType  OutputSizeType;
  typedef typename InputImageType::IndexType  InputIndexType;
  typedef typename OutputImageType::IndexType OutputIndexType;

protected:
  PreservationOfPrincipalDirectionTensorReorientationImageFilter();
  ~PreservationOfPrincipalDirectionTensorReorientationImageFilter() override = default;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /** PreservationOfPrincipalDirectionTensorReorientationImageFilter can be implemented as a multithreaded filter.
   * Therefore, this implementation provides a ThreadedGenerateData()
   * routine which is called for each processing thread. The output
   * image data is allocated automatically by the superclass prior to
   * calling ThreadedGenerateData().  ThreadedGenerateData can only
   * write to the portion of the output image specified by the
   * parameter "outputRegionForThread"
   *
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData() */
  void
  GenerateData() override;


private:
  PreservationOfPrincipalDirectionTensorReorientationImageFilter(const Self &) = delete;
  void
  operator=(const Self &) = delete;

  DisplacementFieldPointer m_DisplacementField;

  DisplacementFieldTransformPointer m_DisplacementTransform;

  AffineTransformPointer m_DirectionTransform;

  AffineTransformPointer m_AffineTransform;

  bool m_UseAffine;

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkPreservationOfPrincipalDirectionTensorReorientationImageFilter.cxx"
#endif

#endif
