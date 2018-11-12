/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkDivideByConstantImageFilter_h
#define __itkDivideByConstantImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "itkNumericTraits.h"

namespace itk
{
/** \class DivideByConstantImageFilter
 *
 * \brief Divide input pixels by a constant.
 *
 * This filter is templated over the input image type
 * and the output image type.
 *
 * \author Tom Vercauteren, INRIA & Mauna Kea Technologies
 *
 * This implementation was taken from the Insight Journal paper:
 * http://hdl.handle.net/1926/510
 *
 * \ingroup IntensityImageFilters  Multithreaded
 * \sa UnaryFunctorImageFilter
 */
namespace Functor
{
template <typename TInput, typename TConstant, typename TOutput>
class DivideByConstant
{
public:
  DivideByConstant() : m_Constant(NumericTraits<TConstant>::OneValue())
  {
  };
  ~DivideByConstant() = default;
  bool operator!=( const DivideByConstant & other ) const
  {
    return !(*this == other);
  }

  bool operator==( const DivideByConstant & other ) const
  {
    return other.m_Constant == m_Constant;
  }

  inline TOutput operator()( const TInput & A ) const
  {
    // Because the user has to specify the constant we don't
    // check if the constant is not too small (i.e. almost equal to zero);
    return static_cast<TOutput>( A / m_Constant );
  }

  void SetConstant(TConstant ct)
  {
    if( ct == NumericTraits<TConstant>::ZeroValue() )
      {
      itkGenericExceptionMacro(
        << "The constant value used as denominator should not be set to zero");
      }

    this->m_Constant = ct;
  }

  const TConstant & GetConstant() const
  {
    return m_Constant;
  }

  TConstant m_Constant;
};
}

template <typename TInputImage, typename TConstant, typename TOutputImage>
class DivideByConstantImageFilter :
  public
  UnaryFunctorImageFilter<TInputImage, TOutputImage,
                          Functor::DivideByConstant<
                            typename TInputImage::PixelType, TConstant,
                            typename TOutputImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef DivideByConstantImageFilter Self;
  typedef UnaryFunctorImageFilter<
      TInputImage, TOutputImage, Functor::DivideByConstant<
        typename TInputImage::PixelType, TConstant,
        typename TOutputImage::PixelType>   >                  Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DivideByConstantImageFilter, UnaryFunctorImageFilter);

  /** Set the constant value that will be used for dividing all the image
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
  // The following concept check doesn't seem to work with vector immages
  // itkConceptMacro(Input1Input2OutputDivisionOperatorsCheck,
  //                (Concept::DivisionOperators<typename TInputImage::PixelType,
  //                 TConstant,
  //                 typename TOutputImage::PixelType>));
  /** End concept checking */
#endif
protected:
  DivideByConstantImageFilter() = default;
  ~DivideByConstantImageFilter() override = default;

  void PrintSelf(std::ostream & os, Indent indent) const override
  {
    Superclass::PrintSelf(os, indent);

    os << indent << "Constant: "
       << static_cast<typename NumericTraits<TConstant>::PrintType>(this->GetConstant() )
       << std::endl;
  }

private:
  DivideByConstantImageFilter(const Self &) = delete;
  void operator=(const Self &) = delete;
};
} // end namespace itk

#endif
