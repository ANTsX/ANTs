/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGeneralToBSplineDeformationFieldFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/11/15 23:46:06 $
  Version:   $Revision: 1.16 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGeneralToBSplineDeformationFieldFilter_h
#define __itkGeneralToBSplineDeformationFieldFilter_h

#include "itkImageToImageFilter.h"

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkPointSet.h"

namespace itk
{
template <class TInputImage, class TOutputImage = TInputImage>
class GeneralToBSplineDeformationFieldFilter :
  public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef GeneralToBSplineDeformationFieldFilter Self;
  typedef ImageToImageFilter<
    TInputImage,
    TOutputImage>           Superclass;

  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Extract dimension from input and output image. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TInputImage::ImageDimension );

  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage  InputImageType;
  typedef TOutputImage OutputImageType;

  typedef typename InputImageType::PixelType  InputPixelType;
  typedef typename OutputImageType::PixelType OutputPixelType;

  typedef typename InputPixelType::ValueType InputPixelComponentType;

  typedef InputPixelType                  VectorType;
  typedef InputPixelComponentType         RealType;
  typedef Image<RealType, ImageDimension> RealImageType;

  typedef PointSet<InputPixelType,
                   itkGetStaticConstMacro( ImageDimension )>   PointSetType;
  typedef BSplineScatteredDataPointSetToImageFilter
  <PointSetType, OutputImageType>             BSplineFilterType;
  typedef typename BSplineFilterType::ArrayType ArrayType;

//  itkSetMacro( ConfidenceImage, RealImageType );
//  itkGetConstMacro( ConfidenceImage, RealImageType );

  itkSetMacro( NumberOfControlPoints, ArrayType );
  itkGetConstMacro( NumberOfControlPoints, ArrayType );

  itkSetMacro( NumberOfLevels, unsigned int );
  itkGetConstMacro( NumberOfLevels, unsigned int );

  itkSetMacro( SplineOrder, unsigned int );
  itkGetConstMacro( SplineOrder, unsigned int );

  itkSetMacro( IgnorePixelValue, InputPixelType );
  itkGetConstMacro( IgnorePixelValue, InputPixelType );

  itkSetMacro( UseFFDRegularization, bool );
  itkGetConstMacro( UseFFDRegularization, bool );
  itkBooleanMacro( UseFFDRegularization );
protected:

  GeneralToBSplineDeformationFieldFilter();
  virtual ~GeneralToBSplineDeformationFieldFilter();

  void PrintSelf( std::ostream& os, Indent indent ) const;

  void GenerateData();

private:
//   typename RealImageType::Pointer              m_ConfidenceImage;

  InputPixelType m_IgnorePixelValue;
  unsigned int   m_NumberOfLevels;
  unsigned int   m_SplineOrder;
  ArrayType      m_NumberOfControlPoints;

  bool m_UseFFDRegularization;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGeneralToBSplineDeformationFieldFilter.txx"
#endif

#endif
