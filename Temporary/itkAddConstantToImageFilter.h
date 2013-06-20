/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkAddConstantToImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-02-24 19:03:14 $
  Version:   $Revision: 1.6 $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkAddConstantToImageFilter_h
#define __itkAddConstantToImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "itkNumericTraits.h"

namespace itk
{
/** \class AddConstantToImageFilter
 *
 * \brief Add a constant to all input pixels.
 *
 * This filter is templated over the input image type
 * and the output image type.
 *
 * \author Tom Vercauteren, INRIA & Mauna Kea Technologies
 *
 * Based on filters from the Insight Journal paper:
 * http://hdl.handle.net/1926/510
 *
 * \ingroup IntensityImageFilters  Multithreaded
 * \sa UnaryFunctorImageFilter
 */
namespace Functor
{
template <class TInput, class TConstant, class TOutput>
class AddConstantTo
{
public:
  AddConstantTo() : m_Constant(NumericTraits<TConstant>::One)
  {
  };
  ~AddConstantTo()
  {
  };
  bool operator!=( const AddConstantTo & other ) const
  {
    return !(*this == other);
  }

  bool operator==( const AddConstantTo & other ) const
  {
    return other.m_Constant == m_Constant;
  }

  inline TOutput operator()( const TInput & A ) const
  {
    // Because the user has to specify the constant we don't
    // check if the cte is not 0;
    return static_cast<TOutput>( A + m_Constant );
  }

  void SetConstant(TConstant ct)
  {
    this->m_Constant = ct;
  }

  const TConstant & GetConstant() const
  {
    return m_Constant;
  }

  TConstant m_Constant;
};
}

template <class TInputImage, class TConstant, class TOutputImage>
class AddConstantToImageFilter :
  public
  UnaryFunctorImageFilter<TInputImage, TOutputImage,
                          Functor::AddConstantTo<
                            typename TInputImage::PixelType, TConstant,
                            typename TOutputImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef AddConstantToImageFilter Self;
  typedef UnaryFunctorImageFilter<
      TInputImage, TOutputImage,
      Functor::AddConstantTo<
        typename TInputImage::PixelType, TConstant,
        typename TOutputImage::PixelType>   >             Superclass;

  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(AddConstantToImageFilter, UnaryFunctorImageFilter);

  /** Set the constant that will be used to multiply all the image
    * pixels */
  void SetConstant(TConstant ct)
  {
    if( ct != this->GetFunctor().GetConstant() )
      {
      this->GetFunctor().SetConstant(ct);
      this->Modified();
      }
  }

  const TConstant & GetConstant() const
  {
    return this->GetFunctor().GetConstant();
  }

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputConvertibleToOutputCheck,
                  (Concept::Convertible<typename TInputImage::PixelType,
                                        typename TOutputImage::PixelType> ) );
  itkConceptMacro(Input1Input2OutputAddOperatorCheck,
                  (Concept::AdditiveOperators<typename TInputImage::PixelType,
                                              TConstant,
                                              typename TOutputImage::PixelType> ) );
  /** End concept checking */
#endif
protected:
  AddConstantToImageFilter()
  {
  };
  virtual ~AddConstantToImageFilter()
  {
  };

  void PrintSelf(std::ostream & os, Indent indent) const
  {
    Superclass::PrintSelf(os, indent);

    os << indent << "Constant: "
       << static_cast<typename NumericTraits<TConstant>::PrintType>(this->GetConstant() )
       << std::endl;
  }

private:
  AddConstantToImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);           // purposely not implemented
};
} // end namespace itk

#endif
