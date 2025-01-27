/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDeformationFieldGradientTensorImageFilter_h
#define __itkDeformationFieldGradientTensorImageFilter_h

#include "itkConstNeighborhoodIterator.h"
#include "itkImageToImageFilter.h"
#include "itkMatrix.h"
#include "itkVector.h"

namespace itk
{
/** \class DeformationFieldGradientTensorImageFilter
 *
 */
template <typename TInputImage,
          typename TRealType = float,
          typename TOutputImage =
            Image<itk::Matrix<TRealType, TInputImage::ImageDimension, TInputImage::PixelType::Dimension>,
                  TInputImage::ImageDimension>>
class DeformationFieldGradientTensorImageFilter final : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef DeformationFieldGradientTensorImageFilter     Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkOverrideGetNameOfClassMacro(DeformationFieldGradientTensorImageFilter);

  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  typedef typename TOutputImage::PixelType OutputPixelType;
  typedef typename TInputImage::PixelType  InputPixelType;

  /** Image typedef support */
  typedef TInputImage                       InputImageType;
  typedef TOutputImage                      OutputImageType;
  typedef typename InputImageType::Pointer  InputImagePointer;
  typedef typename OutputImageType::Pointer OutputImagePointer;

  /** The dimensionality of the input and output images. */
  static constexpr unsigned int ImageDimension = TOutputImage::ImageDimension;

  /** Length of the vector pixel type of the input image. */
  static constexpr unsigned int VectorDimension = InputPixelType::Dimension;

  /** Define the data type and the vector of data type used in calculations. */
  typedef TRealType                                         RealType;
  typedef Image<RealType, ImageDimension>                   RealImageType;
  typedef Vector<RealType, VectorDimension>                 RealVectorType;
  typedef Image<RealVectorType, ImageDimension>             RealVectorImageType;
  typedef Matrix<RealType, ImageDimension, VectorDimension> RealMatrixType;


  /** Type of the iterator that will be used to move through the image.  Also
      the type which will be passed to the evaluate function */
  typedef ConstNeighborhoodIterator<RealVectorImageType>     ConstNeighborhoodIteratorType;
  typedef typename ConstNeighborhoodIteratorType::RadiusType RadiusType;

  /** Superclass typedefs. */
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

  /** DeformationFieldGradientTensorImageFilter needs a larger input requested
   * region than the output requested region (larger by the kernel
   * size to calculate derivatives).  As such,
   * DeformationFieldGradientTensorImageFilter needs to provide an
   * implementation for GenerateInputRequestedRegion() in order to inform the
   * pipeline execution model.
   *
   * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
  void
  GenerateInputRequestedRegion() override;

  itkSetClampMacro(Order, unsigned int, 1, 2);
  itkGetConstReferenceMacro(Order, unsigned int);

  itkSetMacro(UseCenteredDifference, bool);
  itkGetConstReferenceMacro(UseCenteredDifference, bool);
  itkBooleanMacro(UseCenteredDifference);

  itkSetMacro(UseImageSpacing, bool);
  itkGetConstReferenceMacro(UseImageSpacing, bool);
  itkBooleanMacro(UseImageSpacing);

  itkSetMacro(CalculateJacobian, bool);
  itkGetConstReferenceMacro(CalculateJacobian, bool);
  itkBooleanMacro(CalculateJacobian);

  /** Get access to the input image casted as real pixel values */
  itkGetConstObjectMacro(RealValuedInputImage, RealVectorImageType);

protected:
  DeformationFieldGradientTensorImageFilter();
  ~DeformationFieldGradientTensorImageFilter() override = default;

  /** Do any necessary casting/copying of the input data.  Input pixel types
     whose value types are not real number types must be cast to real number
     types.*/
  void
  BeforeThreadedGenerateData() override;

  /** DeformationFieldGradientTensorImageFilter can be implemented as a
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
  bool                                       m_UseImageSpacing;
  bool                                       m_UseCenteredDifference;
  bool                                       m_CalculateJacobian;
  unsigned int                               m_Order;
  RadiusType                                 m_NeighborhoodRadius;
  Vector<RealType, ImageDimension>           m_DerivativeWeights;
  typename RealVectorImageType::ConstPointer m_RealValuedInputImage;

  RealMatrixType
  EvaluateAtNeighborhood(const ConstNeighborhoodIteratorType &) const;

  DeformationFieldGradientTensorImageFilter(const Self &) = delete;
  void
  operator=(const Self &) = delete;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkDeformationFieldGradientTensorImageFilter.hxx"
#endif

#endif
