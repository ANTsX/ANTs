/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkMultiplyByConstantVectorImageFilter_h
#define __itkMultiplyByConstantVectorImageFilter_h

#include "itkUnaryFunctorImageFilter.h"

namespace itk
{
/** \class MultiplyByConstantVectorImageFilter
 *
 * \brief Multiply input pixels by a constant.
 *
 * This filter is templated over the input image type
 * and the output image type.
 *
 * \ingroup IntensityImageFilters  Multithreaded
 * \sa UnaryFunctorImageFilter
 */
namespace Functor
{
template <typename TInput, typename TConstantVector, typename TOutput>
class MultiplyByConstantVector
{
public:
  MultiplyByConstantVector() {}

  ~MultiplyByConstantVector() {}

  bool
  operator!=(const MultiplyByConstantVector & other) const
  {
    return !(*this == other);
  }

  bool
  operator==(const MultiplyByConstantVector & other) const
  {
    return other.m_ConstantVector == m_ConstantVector;
  }

  inline TOutput
  operator()(const TInput & A) const
  {
    // Because the user has to specify the constant we don't
    // check if the cte is not 0;

    TConstantVector value;

    for (unsigned int i = 0; i < m_ConstantVector.GetVectorDimension(); i++)
    {
      value[i] = A[i] * m_ConstantVector[i];
    }
    return value;
  }

  void
  SetConstantVector(TConstantVector ct)
  {
    this->m_ConstantVector = ct;
  }

  const TConstantVector &
  GetConstantVector() const
  {
    return m_ConstantVector;
  }

  TConstantVector m_ConstantVector;
};
} // namespace Functor

template <typename TInputImage, typename TConstantVector, typename TOutputImage>
class MultiplyByConstantVectorImageFilter
  : public UnaryFunctorImageFilter<TInputImage,
                                   TOutputImage,
                                   Functor::MultiplyByConstantVector<typename TInputImage::PixelType,
                                                                     TConstantVector,
                                                                     typename TOutputImage::PixelType>>
{
public:
  /** Standard class typedefs. */
  typedef MultiplyByConstantVectorImageFilter Self;
  typedef UnaryFunctorImageFilter<TInputImage,
                                  TOutputImage,
                                  Functor::MultiplyByConstantVector<typename TInputImage::PixelType,
                                                                    TConstantVector,
                                                                    typename TOutputImage::PixelType>>
    Superclass;

  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(MultiplyByConstantVectorImageFilter);

  /** Set the constant that will be used to multiply all the image pixels */
  void
  SetConstantVector(TConstantVector ct)
  {
    if (ct != this->GetFunctor().GetConstantVector())
    {
      this->GetFunctor().SetConstantVector(ct);
      this->Modified();
    }
  }

  const TConstantVector &
  GetConstantVector() const
  {
    return this->GetFunctor().GetConstantVector();
  }

protected:
  MultiplyByConstantVectorImageFilter() {}

  virtual ~MultiplyByConstantVectorImageFilter() {}

  void
  PrintSelf(std::ostream & os, Indent indent) const
  {
    Superclass::PrintSelf(os, indent);

    os << indent << "Constant: " << this->GetConstantVector() << std::endl;
  }

private:
  MultiplyByConstantVectorImageFilter(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented
};
} // end namespace itk

#endif
