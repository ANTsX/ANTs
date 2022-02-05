/*=========================================================================

  Program:   Advanced Normalization Tools

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 https://github.com/stnava/ANTs/blob/master/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkDecomposeTensorFunction2_h
#define itkDecomposeTensorFunction2_h

#include "itkVariableSizeMatrix.h"

namespace itk
{
/** \class DecomposeTensorFunction2
 *
 */
template <typename TInput, typename TRealType = float, typename TOutput = itk::VariableSizeMatrix<TRealType>>
class DecomposeTensorFunction2
{
public:
  /** Standard class typedefs. */
  typedef DecomposeTensorFunction2 Self;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  typedef TInput  InputMatrixType;
  typedef TOutput OutputMatrixType;

  /** Define the data type and the vector of data type used in calculations. */
  typedef TRealType RealType;

  // Wrappers for vnl routines
  void
  EvaluateEigenDecomposition(InputMatrixType &, OutputMatrixType &, OutputMatrixType &);

  void
  EvaluateSymmetricEigenDecomposition(InputMatrixType &, OutputMatrixType &, OutputMatrixType &);

  void
  EvaluateQRDecomposition(InputMatrixType &, OutputMatrixType &, OutputMatrixType &);

  void
  EvaluateSVDDecomposition(InputMatrixType &, OutputMatrixType &, OutputMatrixType &, OutputMatrixType &);

  void
  EvaluateSVDEconomyDecomposition(InputMatrixType &, OutputMatrixType &, OutputMatrixType &);

  void
  EvaluateLeftPolarDecomposition(InputMatrixType &, OutputMatrixType &, OutputMatrixType &);

  void
  EvaluateRightPolarDecomposition(InputMatrixType &, OutputMatrixType &, OutputMatrixType &);

  void
  EvaluateCholeskyDecomposition(InputMatrixType &, OutputMatrixType &);

  RealType
  EvaluateDeterminant(InputMatrixType &);

  DecomposeTensorFunction2();
  virtual ~DecomposeTensorFunction2() = default;

protected:
  void
  PrintSelf(std::ostream & os, Indent indent) const;

private:
  DecomposeTensorFunction2(const Self &) = delete;
  void
  operator=(const Self &) = delete;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkDecomposeTensorFunction2.hxx"
#endif

#endif
