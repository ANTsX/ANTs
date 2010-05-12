/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkGridImageSource.h,v $
  Language:  C++
  Date:      $Date: 2008/11/15 23:46:06 $
  Version:   $Revision: 1.16 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGridImageSource_h
#define __itkGridImageSource_h

#include "itkImageSource.h"
#include "itkFixedArray.h"
#include "itkKernelFunction.h"
#include "itkVectorContainer.h"

#include "vnl/vnl_vector.h"

namespace itk
{
/** \class GridImageSource
 * \brief Generate an n-dimensional image of a grid.
 *
 * GridImageSource generates an image of a grid.
 *
 * The output image may be of any dimension.
 *
 * \ingroup DataSources
 */
template <typename TOutputImage>
class ITK_EXPORT GridImageSource : public ImageSource<TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef GridImageSource           Self;
  typedef ImageSource<TOutputImage> Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( GridImageSource, ImageSource );

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  typedef double RealType;

  /** Dimensionality of the output image */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TOutputImage::ImageDimension );

  /** Typedef for the output image types. */
  typedef TOutputImage                       ImageType;
  typedef typename TOutputImage::RegionType  ImageRegionType;
  typedef typename TOutputImage::PixelType   PixelType;
  typedef typename TOutputImage::SpacingType SpacingType;
  typedef typename TOutputImage::PointType   OriginType;
  typedef typename TOutputImage::SizeType    SizeType;

  /** Other convenient types. */
  typedef FixedArray<RealType,
                     itkGetStaticConstMacro( ImageDimension )>              ArrayType;
  typedef FixedArray<bool,
                     itkGetStaticConstMacro( ImageDimension )>              BoolArrayType;
  typedef vnl_vector<RealType>                           PixelArrayType;
  typedef VectorContainer<unsigned long, PixelArrayType> PixelArrayContainerType;

  /** Gets and sets for the output image. */
  itkSetMacro( Spacing, SpacingType );
  itkGetConstMacro( Spacing, SpacingType );
  itkSetMacro( Origin, OriginType );
  itkGetConstMacro( Origin, OriginType );
  itkSetMacro( Size, SizeType );
  itkGetConstMacro( Size, SizeType );

  /** Gets and sets for grid parameters */
  itkSetObjectMacro( KernelFunction, KernelFunction );
  itkGetObjectMacro( KernelFunction, KernelFunction );
  itkSetMacro( Sigma, ArrayType );
  itkGetConstMacro( Sigma, ArrayType );
  itkSetMacro( GridSpacing, ArrayType );
  itkGetConstMacro( GridSpacing, ArrayType );
  itkSetMacro( GridOffset, ArrayType );
  itkGetConstMacro( GridOffset, ArrayType );
  itkSetMacro( WhichDimensions, BoolArrayType );
  itkGetConstMacro( WhichDimensions, BoolArrayType );
  itkSetMacro( Scale, RealType );
  itkGetConstMacro( Scale, RealType );
protected:
  GridImageSource();
  ~GridImageSource();
  void PrintSelf(std::ostream& os, Indent indent) const;

  virtual void ThreadedGenerateData(const ImageRegionType &
                                    outputRegionForThread, int threadId );

  virtual void BeforeThreadedGenerateData();

  virtual void GenerateOutputInformation();

private:
  GridImageSource(const GridImageSource &); // purposely not implemented
  void operator=(const GridImageSource &);  // purposely not implemented

  /** Parameters for the output image. */

  SizeType    m_Size;       // size
  SpacingType m_Spacing;    // spacing
  OriginType  m_Origin;     // origin

  /** Parameters for the grid. */

  /** Internal variable to speed up the calculation of pixel values */
  typename PixelArrayContainerType::Pointer m_PixelArrays;

  /** The kernel function used to create the grid */
  typename KernelFunction::Pointer m_KernelFunction;

  /** The standard deviation of the gaussians
      or width of the box functions. */
  ArrayType m_Sigma;

  /** The grid spacing of the peaks. */
  ArrayType m_GridSpacing;

  /** The grid spacing of the peaks. */
  ArrayType m_GridOffset;

  /** Which dimensions which are gridded. */
  BoolArrayType m_WhichDimensions;

  /** A scale factor multiplied by the true value of the grid. */
  RealType m_Scale;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGridImageSource.txx"
#endif

#endif
