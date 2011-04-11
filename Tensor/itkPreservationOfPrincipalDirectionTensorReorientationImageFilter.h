/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkPreservationOfPrincipalDirectionTensorReorientationImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009/01/27 19:56:28 $
  Version:   $Revision: 1.1 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPreservationOfPrincipalDirectionTensorReorientationImageFilter_h
#define __itkPreservationOfPrincipalDirectionTensorReorientationImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkVectorInterpolateImageFunction.h"
#include "itkImage.h"
#include "itkMatrix.h"
#include "itkNumericTraits.h"
#include "itkVector.h"
#include "itkSymmetricSecondRankTensor.h"

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
class ITK_EXPORT PreservationOfPrincipalDirectionTensorReorientationImageFilter :
  public         ImageToImageFilter<TTensorImage, TTensorImage>
{
public:

  typedef TTensorImage InputImageType;
  typedef TTensorImage OutputImageType;
  typedef TVectorImage DeformationFieldType;

  typedef typename DeformationFieldType::Pointer   DeformationFieldPointer;
  typedef typename DeformationFieldType::PixelType VectorType;

  typedef Matrix<float, 3, 3> MatrixType;
  // typedef Vector<float, 3> VectorType;
  typedef VariableSizeMatrix<float> VariableMatrixType;

  itkStaticConstMacro(ImageDimension, unsigned int,
                      TTensorImage::ImageDimension);

  typedef itk::Image<float, ImageDimension> FloatImageType;

  typedef itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension> AffineTransformType;

  typedef typename AffineTransformType::Pointer AffineTransformPointer;

  typedef typename AffineTransformType::InputVectorType TransformInputVectorType;

  typedef typename AffineTransformType::OutputVectorType TransformOutputVectorType;

  typedef typename AffineTransformType::InverseTransformBaseType InverseTransformType;

  typedef typename InverseTransformType::Pointer InverseTransformPointer;

  //  typedef Vector<float, 6> TensorType;
  // typedef itk::SymmetricSecondRankTensor< float, 3 >  TensorType;
  // typedef Image<TensorType, ImageDimension> TensorImageType;
  // typedef typename TensorImageType::Pointer TensorImagePointer;
  typedef typename InputImageType::PixelType InputImagePixelType;
  typedef InputImagePixelType                TensorType;

  typedef vnl_matrix<float> VnlMatrixType;
  typedef vnl_vector<float> VnlVectorType;
  typedef vnl_vector<float> vvec;

  /** Standard class typedefs. */
  typedef PreservationOfPrincipalDirectionTensorReorientationImageFilter Self;
  typedef ImageToImageFilter<InputImageType, OutputImageType>            Superclass;
  typedef SmartPointer<Self>                                             Pointer;
  typedef SmartPointer<const Self>                                       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  itkSetMacro(DeformationField, DeformationFieldPointer);
  itkGetMacro(DeformationField, DeformationFieldPointer);

  void SetAffineTransform( AffineTransformPointer aff )
  {
    this->m_AffineTransform = aff;
    this->m_UseAffine = true;
  }

  itkGetMacro( AffineTransform, AffineTransformPointer );

  itkSetMacro( UseImageDirection, bool );
  itkGetMacro( UseImageDirection, bool );

  /** Run-time type information (and related methods). */
  itkTypeMacro(PreservationOfPrincipalDirectionTensorReorientationImageFilter, ImageToImageFilter);

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
  virtual ~PreservationOfPrincipalDirectionTensorReorientationImageFilter()
  {
  }

  void PrintSelf(std::ostream& os, Indent indent) const;

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
  void GenerateData( void );

private:
  PreservationOfPrincipalDirectionTensorReorientationImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);                                                 // purposely not implemented

  AffineTransformPointer GetLocalDeformation(DeformationFieldPointer, typename DeformationFieldType::IndexType );

  TensorType ApplyReorientation(InverseTransformPointer, TensorType );

  void DirectionCorrectTransform( AffineTransformPointer, AffineTransformPointer );

  DeformationFieldPointer m_DeformationField;

  AffineTransformPointer m_DirectionTransform;

  AffineTransformPointer m_AffineTransform;

  InverseTransformPointer m_InverseAffineTransform;

  bool m_UseAffine;

  bool m_UseImageDirection;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPreservationOfPrincipalDirectionTensorReorientationImageFilter.cxx"
#endif

#endif
