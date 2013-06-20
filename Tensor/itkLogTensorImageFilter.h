/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkLogTensorImageFilter.h,v $
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
#ifndef __itkLogTensorImageFilter_h
#define __itkLogTensorImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkNumericTraits.h"
#include "TensorFunctions.h"

namespace itk
{
/** \class LogTensorImageFilter
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
template <class TInputImage, class TOutputImage>
class LogTensorImageFilter :
  public         ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Extract dimension from input and output image. */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage  InputImageType;
  typedef TOutputImage OutputImageType;

  /** Standard class typedefs. */
  typedef LogTensorImageFilter                                Self;
  typedef ImageToImageFilter<InputImageType, OutputImageType> Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(LogTensorImageFilter, ImageToImageFilter);

  /** Image typedef support. */
  typedef typename InputImageType::ConstPointer InputImagePointer;
  typedef typename OutputImageType::Pointer     OutputImagePointer;
  typedef typename InputImageType::PixelType    InputPixelType;
  typedef typename OutputImageType::PixelType   OutputPixelType;
  typedef typename InputPixelType::ValueType    InputRealType;
  typedef typename OutputPixelType::ValueType   OutputRealType;

  // typedef typename DiffusionTensor3D<float> TensorType;

  // typedef typename TensorType::EigenValuesArrayType EigenValuesArrayType;
  // typedef typename TensorType::EigenVectorsMatrixType EigenVectorsMatrixType;

  typedef typename InputImageType::RegionType  InputImageRegionType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef typename InputImageType::SizeType   InputSizeType;
  typedef typename OutputImageType::SizeType  OutputSizeType;
  typedef typename InputImageType::IndexType  InputIndexType;
  typedef typename OutputImageType::IndexType OutputIndexType;

  /** Set the radius of the neighborhood used to compute the mean. */
  // itkSetMacro(Radius, InputSizeType);

  /** Get the radius of the neighborhood used to compute the mean */
  // itkGetConstReferenceMacro(Radius, InputSizeType);
protected:
  LogTensorImageFilter();
  virtual ~LogTensorImageFilter()
  {
  }

  void PrintSelf(std::ostream& os, Indent indent) const;

  /** LogTensorImageFilter can be implemented as a multithreaded filter.
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
  LogTensorImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);       // purposely not implemented
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLogTensorImageFilter.hxx"
#endif

#endif
