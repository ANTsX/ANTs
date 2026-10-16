/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/ANTsX/ANTs/blob/main/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGeneralToBSplineDisplacementFieldFilter_h
#define __itkGeneralToBSplineDisplacementFieldFilter_h

#include "ITKGetterSetterMacroShims.h"
#include "itkImageToImageFilter.h"

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkPointSet.h"

namespace itk
{
template <typename TInputImage, typename TOutputImage = TInputImage>
class GeneralToBSplineDisplacementFieldFilter final : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef GeneralToBSplineDisplacementFieldFilter       Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;

  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Extract dimension from input and output image. */
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;

  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage  InputImageType;
  typedef TOutputImage OutputImageType;

  typedef typename InputImageType::PixelType  InputPixelType;
  typedef typename OutputImageType::PixelType OutputPixelType;

  typedef typename InputPixelType::ValueType InputPixelComponentType;

  typedef InputPixelType                  VectorType;
  typedef InputPixelComponentType         RealType;
  typedef Image<RealType, ImageDimension> RealImageType;

  typedef PointSet<InputPixelType, Self::ImageDimension>         PointSetType;
  typedef BSplineScatteredDataPointSetToImageFilter<PointSetType, OutputImageType> BSplineFilterType;
  typedef typename BSplineFilterType::ArrayType                                    ArrayType;

  //  itkSetMacro( ConfidenceImage, RealImageType );
  //  itkGetConstMacro( ConfidenceImage, RealImageType );

  itkNonVirtualSetMacro(NumberOfControlPoints, ArrayType);
  itkNonVirtualGetConstMacro(NumberOfControlPoints, ArrayType);

  itkNonVirtualSetMacro(NumberOfLevels, unsigned int);
  itkNonVirtualGetConstMacro(NumberOfLevels, unsigned int);

  itkNonVirtualSetMacro(SplineOrder, unsigned int);
  itkNonVirtualGetConstMacro(SplineOrder, unsigned int);

  itkNonVirtualSetMacro(IgnorePixelValue, InputPixelType);
  itkNonVirtualGetConstMacro(IgnorePixelValue, InputPixelType);

protected:
  GeneralToBSplineDisplacementFieldFilter();
  ~GeneralToBSplineDisplacementFieldFilter() override;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  GenerateData() override;

private:
  //   typename RealImageType::Pointer              m_ConfidenceImage;

  InputPixelType m_IgnorePixelValue;
  unsigned int   m_NumberOfLevels;
  unsigned int   m_SplineOrder;
  ArrayType      m_NumberOfControlPoints;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkGeneralToBSplineDisplacementFieldFilter.hxx"
#endif

#endif
