/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVectorNormImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-04-01 14:36:18 $
  Version:   $Revision: 1.25 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVectorNormImageFilter_h
#define __itkVectorNormImageFilter_h

#include "itkUnaryFunctorImageFilter.h"

namespace itk
{
/** \class VectorNormImageFilter
 * \brief Computes the Norm(x) or SquaredNorm(x) pixel-wise
 *
 *
 * \ingroup IntensityImageFilters  Multithreaded
 *
 */

namespace Function
{
template <class TInput, class TOutput>
class VectorNorm
{
public:
  VectorNorm()
  {
  };
  ~VectorNorm()
  {
  };
  bool operator!=( const VectorNorm & ) const
  {
    return false;
  }

  bool operator==( const VectorNorm & other ) const
  {
    return !(*this != other);
  }

#if defined(_MSC_VER) && (_MSC_VER == 1300)
#pragma optimize("g",off)
#endif
  inline TOutput operator()( const TInput & A ) const
  {
    return (TOutput)( (double)A.GetNorm() );
  }
};
}
template <class TInputImage, class TOutputImage>
class ITK_EXPORT VectorNormImageFilter :
  public
  UnaryFunctorImageFilter<TInputImage, TOutputImage,
                          Function::VectorNorm<
                            typename TInputImage::PixelType,
                            typename TOutputImage::PixelType>   >
{
public:
  /** Standard class typedefs. */
  typedef VectorNormImageFilter Self;
  typedef UnaryFunctorImageFilter<
      TInputImage, TOutputImage,
      Function::VectorNorm<typename TInputImage::PixelType,
                           typename TOutputImage::PixelType> >  Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(VectorNormImageFilter,
               UnaryFunctorImageFilter);
protected:
  VectorNormImageFilter()
  {
  }

  virtual ~VectorNormImageFilter()
  {
  }

private:
  VectorNormImageFilter( const Self & ); // purposely not implemented
  void operator=( const Self & );        // purposely not implemented
};
#if defined(_MSC_VER) && (_MSC_VER == 1300)
#pragma optimize("",on)
#endif
} // end namespace itk

#endif
