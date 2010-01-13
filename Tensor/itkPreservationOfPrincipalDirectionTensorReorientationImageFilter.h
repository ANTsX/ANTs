/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPreservationOfPrincipalDirectionTensorReorientationImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009/01/27 19:56:28 $
  Version:   $Revision: 1.1 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

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

  typedef typename DeformationFieldType::Pointer DeformationFieldPointer;

  typedef Matrix<double, 3, 3> MatrixType;
  typedef Vector<double, 3>    VectorType;

  typedef Vector<double, 6> FakeVectorType;
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TTensorImage::ImageDimension);

  typedef itk::Image<float, ImageDimension> FloatImageType;

  //  typedef Vector<float, 6> TensorType;
  typedef itk::SymmetricSecondRankTensor<float, 3> TensorType;
  typedef Image<TensorType, ImageDimension>        TensorImageType;
  typedef typename TensorImageType::Pointer        TensorImagePointer;
  typedef vnl_matrix<float>                        VnlMatrixType;
  typedef vnl_vector<float>                        vvec;

  /** Standard class typedefs. */
  typedef PreservationOfPrincipalDirectionTensorReorientationImageFilter Self;
  typedef ImageToImageFilter<InputImageType, OutputImageType>            Superclass;
  typedef SmartPointer<Self>                                             Pointer;
  typedef SmartPointer<const Self>                                       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  itkSetMacro(DeformationField, DeformationFieldPointer);
  itkGetMacro(DeformationField, DeformationFieldPointer);

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

  InputPixelType PPD(typename InputImageType::IndexType ind,  VnlMatrixType jac);

  vnl_matrix<float>  FindRotation( vnl_vector<float>,  vnl_vector<float> );

  TensorType ReorientTensors( TensorType fixedTens,  TensorType movingTens,  typename TTensorImage::IndexType index);

  // ::ReorientTensors( Vector<float, 6>  fixedTens, Vector<float, 6> movingTens,  typename TTensorImage::IndexType
  // index)

  TensorImagePointer m_ReferenceImage;
  typename FloatImageType::Pointer m_ThetaImage;
  typename FloatImageType::Pointer m_PsiImage;
  typename FloatImageType::Pointer m_AngleEnergyImage;
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

  typename DeformationFieldType::PixelType TransformVector(DeformationFieldType* field,
                                                           typename FloatImageType::IndexType index )
  {
    typename DeformationFieldType::PixelType vec = field->GetPixel(index);
    typename DeformationFieldType::PixelType newvec;
    newvec.Fill(0);
    for( unsigned int row = 0; row < ImageDimension; row++ )
      {
      for( unsigned int col = 0; col < ImageDimension; col++ )
        {
        newvec[row] += vec[col] * field->GetDirection()[row][col];
        }
      }

    return newvec;
  }

private:
  PreservationOfPrincipalDirectionTensorReorientationImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);                                                 // purposely not implemented

  DeformationFieldPointer m_DeformationField;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPreservationOfPrincipalDirectionTensorReorientationImageFilter.cxx"
#endif

#endif
